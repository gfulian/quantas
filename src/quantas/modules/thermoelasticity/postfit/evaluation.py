# -*- coding: utf-8 -*-

"""Post-fit pressure-temperature reconstruction from thermoelastic archives.

This module evaluates the already fitted quasi-static elastic model without
repeating either the static EOS fit or the independent-component fits.  The
QHA equilibrium-volume surface archived by the fitting run supplies
``V(P, T)`` through rectilinear interpolation or controlled extrapolation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from quantas.core.numerics import (
    validated_axis,
)
from quantas.core.physics.elasticity import (
    AdiabaticStiffnessFieldResult,
    StabilityFieldResult,
    evaluate_stability_field,
)
from quantas.modules.thermoelasticity.adiabatic import (
    calculate_adiabatic_field,
)
from quantas.modules.thermoelasticity.analysis import (
    density_from_volume,
    reconstruct_stiffness_grid,
)
from quantas.modules.thermoelasticity.fitting import evaluate_component_predictions
from quantas.modules.thermoelasticity.models import (
    ThermoelasticDepthProfile,
    ThermoelasticOptions,
    ThermoelasticProfileResult,
    ThermoelasticResult,
)


from .interpolation import (
    _interpolate_adiabatic_inputs_grid,
    _interpolate_adiabatic_inputs_points,
    interpolate_archived_grid,
    interpolate_archived_points,
)
from .policies import (
    _elastic_volume_bounds,
    _enforce_extrapolation_policy,
    _enforce_stability_policy,
    _normalized_cell_mass,
    _required_metadata_text,
    _stability_summary,
)


@dataclass(slots=True)
class _EvaluatedTensorFields:
    """Internal common tensor fields for one requested state collection."""

    volume: np.ndarray
    sigma_volume: np.ndarray | None
    qha_extrapolated: np.ndarray
    elastic_extrapolated: np.ndarray
    values: np.ndarray
    sigma_values: np.ndarray
    covariance: np.ndarray
    stiffness: np.ndarray
    sigma_stiffness: np.ndarray
    stability: StabilityFieldResult
    adiabatic: AdiabaticStiffnessFieldResult | None


def evaluate_thermoelastic_grid(
    source: ThermoelasticResult,
    temperature: ArrayLike,
    pressure: ArrayLike,
    *,
    options: ThermoelasticOptions,
) -> ThermoelasticResult:
    """Evaluate fitted tensors on a requested rectangular P--T grid."""
    target_temperature = validated_axis(temperature, "temperature")
    target_pressure = validated_axis(pressure, "pressure")
    volume, qha_extrapolated = interpolate_archived_grid(
        source.temperature,
        source.pressure,
        source.equilibrium_volume,
        target_temperature,
        target_pressure,
    )
    sigma_volume = _interpolate_optional_grid(
        source,
        source.sigma_equilibrium_volume,
        target_temperature,
        target_pressure,
    )
    adiabatic_inputs = _interpolate_adiabatic_inputs_grid(
        source,
        target_temperature,
        target_pressure,
    )
    fields = _evaluate_tensor_fields(
        source,
        volume=volume,
        sigma_volume=sigma_volume,
        qha_extrapolated=qha_extrapolated,
        temperature=np.broadcast_to(target_temperature[:, np.newaxis], volume.shape),
        adiabatic_inputs=adiabatic_inputs,
        options=options,
        context="requested pressure-temperature grid",
    )
    mass_kg = _normalized_cell_mass(source)
    return ThermoelasticResult(
        jobname=source.jobname,
        reference_eos=source.reference_eos,
        component_fits=source.component_fits,
        independent_labels=source.independent_labels,
        temperature=target_temperature,
        pressure=target_pressure,
        equilibrium_volume=fields.volume,
        density=density_from_volume(fields.volume, mass_kg),
        independent_stiffness=fields.values,
        sigma_independent_stiffness=fields.sigma_values,
        independent_stiffness_covariance=fields.covariance,
        stiffness_isothermal=fields.stiffness,
        sigma_stiffness_isothermal=fields.sigma_stiffness,
        extrapolation_mask=fields.elastic_extrapolated,
        sigma_equilibrium_volume=fields.sigma_volume,
        qha_extrapolation_mask=fields.qha_extrapolated,
        profiles={},
        isochoric_heat_capacity_cell=(
            None if adiabatic_inputs is None else adiabatic_inputs.heat_capacity
        ),
        sigma_isochoric_heat_capacity_cell=(
            None if adiabatic_inputs is None else adiabatic_inputs.sigma_heat_capacity
        ),
        thermal_expansion_tensor=(
            None
            if adiabatic_inputs is None
            else adiabatic_inputs.thermal_expansion_tensor
        ),
        sigma_thermal_expansion_tensor=(
            None
            if adiabatic_inputs is None
            else adiabatic_inputs.sigma_thermal_expansion_tensor
        ),
        stiffness_adiabatic=(
            None if fields.adiabatic is None else fields.adiabatic.stiffness
        ),
        sigma_stiffness_adiabatic=(
            None if fields.adiabatic is None else fields.adiabatic.sigma_stiffness
        ),
        adiabatic_correction=(
            None if fields.adiabatic is None else fields.adiabatic.correction
        ),
        adiabatic_thermal_stress=(
            None if fields.adiabatic is None else fields.adiabatic.thermal_stress
        ),
        adiabatic_valid_mask=(
            None if fields.adiabatic is None else fields.adiabatic.valid_mask
        ),
        stability=fields.stability,
        completed=True,
        metadata=_grid_analysis_metadata(
            source,
            fields=fields,
        ),
    )


def evaluate_thermoelastic_profile(
    source: ThermoelasticResult,
    profile: ThermoelasticDepthProfile,
    *,
    options: ThermoelasticOptions,
) -> ThermoelasticProfileResult:
    """Evaluate one depth-dependent pressure-temperature profile from fits."""
    volume, qha_extrapolated = interpolate_archived_points(
        source.temperature,
        source.pressure,
        source.equilibrium_volume,
        profile.temperature,
        profile.pressure,
    )
    sigma_volume = _interpolate_optional_points(
        source,
        source.sigma_equilibrium_volume,
        profile.temperature,
        profile.pressure,
    )
    adiabatic_inputs = _interpolate_adiabatic_inputs_points(
        source,
        profile.temperature,
        profile.pressure,
    )
    context = f"depth profile '{profile.name}'"
    fields = _evaluate_tensor_fields(
        source,
        volume=volume,
        sigma_volume=sigma_volume,
        qha_extrapolated=qha_extrapolated,
        temperature=profile.temperature,
        adiabatic_inputs=adiabatic_inputs,
        options=options,
        context=context,
    )
    return ThermoelasticProfileResult(
        name=profile.name,
        depth=profile.depth,
        pressure=profile.pressure,
        temperature=profile.temperature,
        volume=fields.volume,
        density=density_from_volume(fields.volume, _normalized_cell_mass(source)),
        independent_stiffness=fields.values,
        sigma_independent_stiffness=fields.sigma_values,
        independent_stiffness_covariance=fields.covariance,
        stiffness_isothermal=fields.stiffness,
        sigma_stiffness_isothermal=fields.sigma_stiffness,
        qha_extrapolation_mask=fields.qha_extrapolated,
        elastic_extrapolation_mask=fields.elastic_extrapolated,
        stiffness_adiabatic=(
            None if fields.adiabatic is None else fields.adiabatic.stiffness
        ),
        sigma_stiffness_adiabatic=(
            None if fields.adiabatic is None else fields.adiabatic.sigma_stiffness
        ),
        adiabatic_correction=(
            None if fields.adiabatic is None else fields.adiabatic.correction
        ),
        adiabatic_thermal_stress=(
            None if fields.adiabatic is None else fields.adiabatic.thermal_stress
        ),
        adiabatic_valid_mask=(
            None if fields.adiabatic is None else fields.adiabatic.valid_mask
        ),
        stability=fields.stability,
        metadata={
            **profile.metadata,
            "qha_interpolation": "piecewise-linear rectilinear",
            "qha_extrapolated_points": int(np.count_nonzero(fields.qha_extrapolated)),
            "elastic_extrapolated_points": int(
                np.count_nonzero(fields.elastic_extrapolated)
            ),
            "mechanical_stability": _stability_summary(fields.stability),
            "adiabatic_conversion_available": bool(fields.adiabatic is not None),
            "adiabatic_valid_points": (
                0
                if fields.adiabatic is None
                else int(np.count_nonzero(fields.adiabatic.valid_mask))
            ),
        },
    )


def _evaluate_tensor_fields(
    source: ThermoelasticResult,
    *,
    volume: np.ndarray,
    sigma_volume: np.ndarray | None,
    qha_extrapolated: np.ndarray,
    temperature: np.ndarray,
    adiabatic_inputs,
    options: ThermoelasticOptions,
    context: str,
) -> _EvaluatedTensorFields:
    """Evaluate common isothermal, adiabatic, and stability fields."""
    labels = source.independent_labels
    values, sigma_values, covariance = evaluate_component_predictions(
        source.component_fits,
        labels,
        volume,
        source.reference_eos,
        options,
        sigma_volume=sigma_volume,
    )
    stiffness, sigma_stiffness = reconstruct_stiffness_grid(
        symmetry=_required_metadata_text(source, "elastic_symmetry"),
        labels=labels,
        values=values,
        covariance=covariance,
    )
    stability = evaluate_stability_field(
        stiffness, tolerance=options.stability_tolerance
    )
    adiabatic = calculate_adiabatic_field(
        stiffness_isothermal=stiffness,
        sigma_stiffness_isothermal=sigma_stiffness,
        temperature_k=temperature,
        volume_a3=volume,
        sigma_volume_a3=sigma_volume,
        inputs=adiabatic_inputs,
        mode=options.adiabatic_mode,
        propagate_uncertainty=options.propagate_adiabatic_uncertainty,
    )
    lower, upper = _elastic_volume_bounds(source)
    elastic_extrapolated = np.asarray(
        (~np.isfinite(volume)) | (volume < lower) | (volume > upper),
        dtype=np.bool_,
    )
    _enforce_extrapolation_policy(
        qha_extrapolated | elastic_extrapolated,
        options.extrapolation_policy,
        context,
    )
    _enforce_stability_policy(
        stability.unstable_mask | stability.indeterminate_mask,
        options.stability_policy,
        context,
    )
    return _EvaluatedTensorFields(
        volume=volume,
        sigma_volume=sigma_volume,
        qha_extrapolated=qha_extrapolated,
        elastic_extrapolated=elastic_extrapolated,
        values=values,
        sigma_values=sigma_values,
        covariance=covariance,
        stiffness=stiffness,
        sigma_stiffness=sigma_stiffness,
        stability=stability,
        adiabatic=adiabatic,
    )


def _interpolate_optional_grid(
    source: ThermoelasticResult,
    field: np.ndarray | None,
    temperature: np.ndarray,
    pressure: np.ndarray,
) -> np.ndarray | None:
    """Interpolate one optional archived field on a rectangular grid."""
    if field is None:
        return None
    value, _ = interpolate_archived_grid(
        source.temperature, source.pressure, field, temperature, pressure
    )
    return np.maximum(value, 0.0)


def _interpolate_optional_points(
    source: ThermoelasticResult,
    field: np.ndarray | None,
    temperature: np.ndarray,
    pressure: np.ndarray,
) -> np.ndarray | None:
    """Interpolate one optional archived field at paired states."""
    if field is None:
        return None
    value, _ = interpolate_archived_points(
        source.temperature, source.pressure, field, temperature, pressure
    )
    return np.maximum(value, 0.0)


def _grid_analysis_metadata(
    source: ThermoelasticResult,
    *,
    fields: _EvaluatedTensorFields,
) -> dict[str, object]:
    """Build deterministic provenance for one rectangular grid analysis."""
    return {
        **source.metadata,
        "analysis_kind": "pressure-temperature-grid",
        "source_grid_temperature_min_K": float(source.temperature[0]),
        "source_grid_temperature_max_K": float(source.temperature[-1]),
        "source_grid_pressure_min_GPa": float(source.pressure[0]),
        "source_grid_pressure_max_GPa": float(source.pressure[-1]),
        "qha_extrapolated_points": int(np.count_nonzero(fields.qha_extrapolated)),
        "elastic_extrapolated_points": int(
            np.count_nonzero(fields.elastic_extrapolated)
        ),
        "qha_volume_interpolation": "piecewise-linear rectilinear",
        "qha_grid_reconstructed": True,
        "mechanical_stability": _stability_summary(fields.stability),
        "adiabatic_conversion_available": bool(fields.adiabatic is not None),
        "adiabatic_valid_points": (
            0
            if fields.adiabatic is None
            else int(np.count_nonzero(fields.adiabatic.valid_mask))
        ),
    }


__all__ = ["evaluate_thermoelastic_grid", "evaluate_thermoelastic_profile"]
