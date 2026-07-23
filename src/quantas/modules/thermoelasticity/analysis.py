# -*- coding: utf-8 -*-

"""Pressure-temperature and depth-path reconstruction for thermoelastic QSA."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import RegularGridInterpolator

from quantas.core.physics.elasticity import (
    elastic_component_definitions,
    evaluate_stability_field,
    stiffness_component_linear_coefficients,
)
from quantas.core.physics.units import (
    convert_pressure,
    convert_temperature,
    convert_volume,
)
from quantas.modules.thermoelasticity.fitting import evaluate_component_predictions
from quantas.modules.thermoelasticity.models import (
    ElasticComponentFit,
    ElasticVolumeSeries,
    ReferenceEOSFit,
    ThermoelasticDepthProfile,
    ThermoelasticOptions,
    ThermoelasticProfileResult,
)


FloatArray: TypeAlias = NDArray[np.float64]


def qha_grid_in_standard_units(
    temperature: ArrayLike,
    pressure: ArrayLike,
    volume: ArrayLike,
    *,
    temperature_unit: str,
    pressure_unit: str,
    volume_unit: str,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Convert a QHA grid to K, GPa and angstrom-cubed cell volume."""
    t = np.asarray(
        convert_temperature(temperature, temperature_unit, "K"),
        dtype=np.float64,
    )
    p = np.asarray(
        convert_pressure(pressure, pressure_unit, "GPa"),
        dtype=np.float64,
    )
    v = np.asarray(convert_volume(volume, volume_unit, "A"), dtype=np.float64)
    if t.ndim != 1 or p.ndim != 1 or v.shape != (t.size, p.size):
        raise ValueError(
            "QHA temperature, pressure, and volume arrays are inconsistent"
        )
    if np.any(np.diff(t) <= 0.0) or np.any(np.diff(p) <= 0.0):
        raise ValueError(
            "QHA temperature and pressure grids must be strictly increasing"
        )
    return t, p, v


def normalized_cell_mass(series: ElasticVolumeSeries) -> tuple[float, float]:
    """Return normalized-cell mass and its relative spread from SOEC points.

    Returns
    -------
    tuple of float
        Median mass in kg and maximum relative deviation across sampled points.
    """
    masses = series.densities * series.volumes * 1.0e-30
    median = float(np.median(masses))
    spread = float(np.max(np.abs(masses - median)) / median)
    return median, spread


def density_from_volume(volume: ArrayLike, mass_kg: float) -> FloatArray:
    """Calculate density in kg m^-3 from cell mass and angstrom-cubed volume."""
    values = np.asarray(volume, dtype=np.float64)
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("volumes must be finite and positive")
    if not np.isfinite(mass_kg) or mass_kg <= 0.0:
        raise ValueError("mass_kg must be finite and positive")
    return np.asarray(mass_kg / (values * 1.0e-30), dtype=np.float64)


def reconstruct_stiffness_grid(
    *,
    symmetry: str,
    labels: tuple[str, ...],
    values: FloatArray,
    covariance: FloatArray,
) -> tuple[FloatArray, FloatArray]:
    """Reconstruct full Voigt matrices and one-sigma uncertainties.

    The covariance is used when propagating derived relations such as
    ``C66=(C11-C12)/2`` so shared EOS and volume uncertainty correlations are
    retained.
    """
    coefficients = stiffness_component_linear_coefficients(symmetry, labels)
    matrix = np.einsum("ijc,...c->...ij", coefficients, values)
    variance = np.einsum(
        "ijc,...cd,ijd->...ij",
        coefficients,
        covariance,
        coefficients,
    )
    sigma = np.sqrt(np.clip(variance, 0.0, None))
    return np.asarray(matrix, dtype=np.float64), np.asarray(sigma, dtype=np.float64)


def independent_component_labels(symmetry: str) -> tuple[str, ...]:
    """Return the canonical independent component order for a symmetry."""
    return tuple(
        definition.label
        for definition in elastic_component_definitions(symmetry)
        if not definition.derived
    )


def interpolate_qha_volume_along_profile(
    temperature_grid: ArrayLike,
    pressure_grid: ArrayLike,
    volume_grid: ArrayLike,
    profile: ThermoelasticDepthProfile,
) -> tuple[FloatArray, NDArray[np.bool_]]:
    """Interpolate QHA equilibrium volume along a geological profile.

    Linear extrapolation is returned outside the rectangular QHA grid together
    with a mask. The caller applies the selected failure/warning policy.
    """
    temperature = np.asarray(temperature_grid, dtype=np.float64)
    pressure = np.asarray(pressure_grid, dtype=np.float64)
    volume = np.asarray(volume_grid, dtype=np.float64)
    if volume.shape != (temperature.size, pressure.size):
        raise ValueError("volume_grid must have shape (nT, nP)")
    points = np.column_stack((profile.temperature, profile.pressure))
    outside = (
        (profile.temperature < temperature[0])
        | (profile.temperature > temperature[-1])
        | (profile.pressure < pressure[0])
        | (profile.pressure > pressure[-1])
    )
    interpolator = RegularGridInterpolator(
        (temperature, pressure),
        volume,
        method="linear",
        bounds_error=False,
        fill_value=None,
    )
    values = np.asarray(interpolator(points), dtype=np.float64)
    return values, np.asarray(outside, dtype=np.bool_)


def interpolate_qha_uncertainty_along_profile(
    temperature_grid: ArrayLike,
    pressure_grid: ArrayLike,
    sigma_grid: ArrayLike | None,
    profile: ThermoelasticDepthProfile,
) -> FloatArray | None:
    """Interpolate an optional QHA volume-uncertainty grid."""
    if sigma_grid is None:
        return None
    temperature = np.asarray(temperature_grid, dtype=np.float64)
    pressure = np.asarray(pressure_grid, dtype=np.float64)
    sigma = np.asarray(sigma_grid, dtype=np.float64)
    if sigma.shape != (temperature.size, pressure.size):
        return None
    interpolator = RegularGridInterpolator(
        (temperature, pressure),
        sigma,
        method="linear",
        bounds_error=False,
        fill_value=None,
    )
    return np.asarray(
        interpolator(np.column_stack((profile.temperature, profile.pressure))),
        dtype=np.float64,
    )


def evaluate_depth_profile(
    profile: ThermoelasticDepthProfile,
    *,
    temperature_grid: FloatArray,
    pressure_grid: FloatArray,
    volume_grid: FloatArray,
    sigma_volume_grid: FloatArray | None,
    series: ElasticVolumeSeries,
    reference_eos: ReferenceEOSFit,
    component_fits: Mapping[str, ElasticComponentFit],
    labels: tuple[str, ...],
    options: ThermoelasticOptions,
    mass_kg: float,
) -> ThermoelasticProfileResult:
    """Evaluate fitted elastic tensors along one depth-pressure-temperature path."""
    volume, qha_extrapolated = interpolate_qha_volume_along_profile(
        temperature_grid,
        pressure_grid,
        volume_grid,
        profile,
    )
    sigma_volume = interpolate_qha_uncertainty_along_profile(
        temperature_grid,
        pressure_grid,
        sigma_volume_grid,
        profile,
    )
    values, sigma_values, covariance = evaluate_component_predictions(
        component_fits,
        labels,
        volume,
        reference_eos,
        options,
        sigma_volume=sigma_volume,
    )
    stiffness, sigma_stiffness = reconstruct_stiffness_grid(
        symmetry=series.elastic_symmetry,
        labels=labels,
        values=values,
        covariance=covariance,
    )
    stability = evaluate_stability_field(
        stiffness,
        tolerance=options.stability_tolerance,
    )
    lower, upper = series.volume_bounds
    elastic_extrapolated = np.asarray(
        (volume < lower) | (volume > upper) | ~np.isfinite(volume),
        dtype=np.bool_,
    )
    return ThermoelasticProfileResult(
        name=profile.name,
        depth=profile.depth,
        pressure=profile.pressure,
        temperature=profile.temperature,
        volume=volume,
        density=density_from_volume(volume, mass_kg),
        independent_stiffness=values,
        sigma_independent_stiffness=sigma_values,
        independent_stiffness_covariance=covariance,
        stiffness_isothermal=stiffness,
        sigma_stiffness_isothermal=sigma_stiffness,
        qha_extrapolation_mask=qha_extrapolated,
        elastic_extrapolation_mask=elastic_extrapolated,
        stability=stability,
        metadata={
            **profile.metadata,
            "qha_interpolation": "bilinear_regular_grid",
            "qha_extrapolated_points": int(np.count_nonzero(qha_extrapolated)),
            "elastic_extrapolated_points": int(np.count_nonzero(elastic_extrapolated)),
            "mechanically_stable_points": int(np.count_nonzero(stability.stable_mask)),
            "mechanically_unstable_points": int(
                np.count_nonzero(stability.unstable_mask)
            ),
            "mechanical_stability_indeterminate_points": int(
                np.count_nonzero(stability.indeterminate_mask)
            ),
        },
    )


__all__ = [
    "density_from_volume",
    "evaluate_depth_profile",
    "independent_component_labels",
    "interpolate_qha_uncertainty_along_profile",
    "interpolate_qha_volume_along_profile",
    "normalized_cell_mass",
    "qha_grid_in_standard_units",
    "reconstruct_stiffness_grid",
]
