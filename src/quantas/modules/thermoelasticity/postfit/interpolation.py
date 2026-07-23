# -*- coding: utf-8 -*-

"""Post-fit pressure-temperature reconstruction from thermoelastic archives.

This module evaluates the already fitted quasi-static elastic model without
repeating either the static EOS fit or the independent-component fits.  The
QHA equilibrium-volume surface archived by the fitting run supplies
``V(P, T)`` through rectilinear interpolation or controlled extrapolation.
"""

from __future__ import annotations

from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from quantas.core.numerics import (
    RectilinearFieldInterpolator,
)
from quantas.modules.thermoelasticity.adiabatic import (
    StandardAdiabaticInputs,
)
from quantas.modules.thermoelasticity.models import (
    ThermoelasticResult,
)


FloatArray: TypeAlias = NDArray[np.float64]


def _interpolate_adiabatic_inputs_grid(
    source: ThermoelasticResult,
    target_temperature: FloatArray,
    target_pressure: FloatArray,
) -> StandardAdiabaticInputs | None:
    """Interpolate archived SI adiabatic inputs to a rectangular grid."""
    if (
        source.isochoric_heat_capacity_cell is None
        or source.thermal_expansion_tensor is None
    ):
        return None
    cv, _ = interpolate_archived_grid(
        source.temperature,
        source.pressure,
        source.isochoric_heat_capacity_cell,
        target_temperature,
        target_pressure,
    )
    alpha = np.empty((target_temperature.size, target_pressure.size, 3, 3))
    for i in range(3):
        for j in range(3):
            alpha[..., i, j], _ = interpolate_archived_grid(
                source.temperature,
                source.pressure,
                source.thermal_expansion_tensor[..., i, j],
                target_temperature,
                target_pressure,
            )
    sigma_cv = None
    if source.sigma_isochoric_heat_capacity_cell is not None:
        sigma_cv, _ = interpolate_archived_grid(
            source.temperature,
            source.pressure,
            source.sigma_isochoric_heat_capacity_cell,
            target_temperature,
            target_pressure,
        )
        sigma_cv = np.maximum(sigma_cv, 0.0)
    sigma_alpha = None
    if source.sigma_thermal_expansion_tensor is not None:
        sigma_alpha = np.empty_like(alpha)
        for i in range(3):
            for j in range(3):
                sigma_alpha[..., i, j], _ = interpolate_archived_grid(
                    source.temperature,
                    source.pressure,
                    source.sigma_thermal_expansion_tensor[..., i, j],
                    target_temperature,
                    target_pressure,
                )
        np.maximum(sigma_alpha, 0.0, out=sigma_alpha)
    return StandardAdiabaticInputs(
        heat_capacity=cv,
        thermal_expansion_tensor=alpha,
        sigma_heat_capacity=sigma_cv,
        sigma_thermal_expansion_tensor=sigma_alpha,
        metadata={"interpolation": "piecewise-linear rectilinear"},
    )


def _interpolate_adiabatic_inputs_points(
    source: ThermoelasticResult,
    target_temperature: FloatArray,
    target_pressure: FloatArray,
) -> StandardAdiabaticInputs | None:
    """Interpolate archived SI adiabatic inputs to paired states."""
    if (
        source.isochoric_heat_capacity_cell is None
        or source.thermal_expansion_tensor is None
    ):
        return None
    cv, _ = interpolate_archived_points(
        source.temperature,
        source.pressure,
        source.isochoric_heat_capacity_cell,
        target_temperature,
        target_pressure,
    )
    alpha = np.empty((target_temperature.size, 3, 3))
    for i in range(3):
        for j in range(3):
            alpha[:, i, j], _ = interpolate_archived_points(
                source.temperature,
                source.pressure,
                source.thermal_expansion_tensor[..., i, j],
                target_temperature,
                target_pressure,
            )
    sigma_cv = None
    if source.sigma_isochoric_heat_capacity_cell is not None:
        sigma_cv, _ = interpolate_archived_points(
            source.temperature,
            source.pressure,
            source.sigma_isochoric_heat_capacity_cell,
            target_temperature,
            target_pressure,
        )
        sigma_cv = np.maximum(sigma_cv, 0.0)
    sigma_alpha = None
    if source.sigma_thermal_expansion_tensor is not None:
        sigma_alpha = np.empty_like(alpha)
        for i in range(3):
            for j in range(3):
                sigma_alpha[:, i, j], _ = interpolate_archived_points(
                    source.temperature,
                    source.pressure,
                    source.sigma_thermal_expansion_tensor[..., i, j],
                    target_temperature,
                    target_pressure,
                )
        np.maximum(sigma_alpha, 0.0, out=sigma_alpha)
    return StandardAdiabaticInputs(
        heat_capacity=cv,
        thermal_expansion_tensor=alpha,
        sigma_heat_capacity=sigma_cv,
        sigma_thermal_expansion_tensor=sigma_alpha,
        metadata={"interpolation": "piecewise-linear rectilinear"},
    )


def interpolate_archived_grid(
    source_temperature: ArrayLike,
    source_pressure: ArrayLike,
    source_values: ArrayLike,
    target_temperature: ArrayLike,
    target_pressure: ArrayLike,
) -> tuple[FloatArray, NDArray[np.bool_]]:
    """Interpolate a rectilinear field through the shared numerical engine."""
    return RectilinearFieldInterpolator(
        source_temperature, source_pressure, source_values
    ).evaluate_grid(target_temperature, target_pressure)


def interpolate_archived_points(
    source_temperature: ArrayLike,
    source_pressure: ArrayLike,
    source_values: ArrayLike,
    target_temperature: ArrayLike,
    target_pressure: ArrayLike,
) -> tuple[FloatArray, NDArray[np.bool_]]:
    """Interpolate paired states through the shared numerical engine."""
    return RectilinearFieldInterpolator(
        source_temperature, source_pressure, source_values
    ).evaluate_points(target_temperature, target_pressure)


__all__ = ["interpolate_archived_grid", "interpolate_archived_points"]
