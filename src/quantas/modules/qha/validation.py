# -*- coding: utf-8 -*-

"""Validation utilities for quasi-harmonic workflow results.

The functions in this module evaluate thermodynamic consistency, numerical
completeness and agreement between alternative QHA workflows. They operate on
structured result objects and do not depend on command-line or graphical
interfaces.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable

import numpy as np
from scipy.constants import gas_constant

from quantas.core.physics.units import convert_energy_per_temperature
from quantas.modules.qha.models import QHAInput, QHAResult


@dataclass(slots=True)
class PropertyDifference:
    """Numerical difference between two QHA property arrays.

    Parameters
    ----------
    property_name : str
        QHA result attribute name.
    maximum_absolute : float
        Maximum absolute difference.
    maximum_relative : float
        Maximum symmetric relative difference.
    root_mean_square : float
        Root-mean-square absolute difference.
    compared_points : int
        Number of finite points included in the comparison.
    maximum_absolute_temperature, maximum_absolute_pressure : float
        Pressure-temperature location of the maximum absolute difference.
    maximum_relative_temperature, maximum_relative_pressure : float
        Pressure-temperature location of the maximum relative difference.
    """

    property_name: str
    maximum_absolute: float
    maximum_relative: float
    root_mean_square: float
    compared_points: int
    maximum_absolute_temperature: float
    maximum_absolute_pressure: float
    maximum_relative_temperature: float
    maximum_relative_pressure: float


@dataclass(slots=True)
class QHAValidationSummary:
    """Summary of numerical and thermodynamic validation checks.

    Parameters
    ----------
    scheme : str
        Thermodynamic workflow identifier.
    minimization : str
        Volume-minimization method.
    completed : bool
        Whether the workflow completed.
    total_points : int
        Total number of pressure-temperature states.
    valid_points : int
        Number of valid states.
    finite_properties : bool
        Whether all checked properties are finite on valid states.
    volume_decreases_with_pressure : bool
        Whether volume decreases monotonically with pressure at fixed
        temperature.
    positive_bulk_moduli : bool
        Whether ``KT`` and ``KS`` are positive on valid states.
    zero_kelvin_consistency : bool | None
        Whether ``KS = KT`` and ``Cp-Cv = 0`` at zero kelvin, when zero kelvin
        is present.
    cp_not_below_cv : bool | None
        Whether ``Cp >= Cv`` within numerical tolerance.
    dulong_petit_ratio : float | None
        High-temperature ``Cv`` divided by the Dulong-Petit limit at the
        lowest pressure.
    volumes_below_sampled_range : int
        Number of equilibrium volumes below the sampled interval.
    volumes_above_sampled_range : int
        Number of equilibrium volumes above the sampled interval.
    extrema : dict
        Minimum and maximum values of selected properties.
    warnings : list of str
        Validation warnings.
    """

    scheme: str
    minimization: str
    completed: bool
    total_points: int
    valid_points: int
    finite_properties: bool
    volume_decreases_with_pressure: bool
    positive_bulk_moduli: bool
    zero_kelvin_consistency: bool | None
    cp_not_below_cv: bool | None
    dulong_petit_ratio: float | None
    volumes_below_sampled_range: int
    volumes_above_sampled_range: int
    extrema: dict[str, tuple[float, float]] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        """Return a serializable representation of the summary.

        Returns
        -------
        dict
            Validation fields and warnings.
        """
        return {
            "scheme": self.scheme,
            "minimization": self.minimization,
            "completed": self.completed,
            "total_points": self.total_points,
            "valid_points": self.valid_points,
            "finite_properties": self.finite_properties,
            "volume_decreases_with_pressure": self.volume_decreases_with_pressure,
            "positive_bulk_moduli": self.positive_bulk_moduli,
            "zero_kelvin_consistency": self.zero_kelvin_consistency,
            "cp_not_below_cv": self.cp_not_below_cv,
            "dulong_petit_ratio": self.dulong_petit_ratio,
            "volumes_below_sampled_range": self.volumes_below_sampled_range,
            "volumes_above_sampled_range": self.volumes_above_sampled_range,
            "extrema": dict(self.extrema),
            "warnings": list(self.warnings),
        }


DEFAULT_COMPARISON_PROPERTIES: tuple[str, ...] = (
    "equilibrium_volume",
    "isothermal_bulk_modulus",
    "bulk_modulus_derivative",
    "adiabatic_bulk_modulus",
    "thermal_expansion",
    "isochoric_heat_capacity",
    "isobaric_heat_capacity",
    "heat_capacity_difference",
    "free_energy",
)


def validate_qha_result(
    result: QHAResult,
    input_data: QHAInput,
    *,
    numerical_tolerance: float = 1.0e-10,
) -> QHAValidationSummary:
    """Validate one QHA result on its pressure-temperature grid.

    Parameters
    ----------
    result : QHAResult
        Result to validate.
    input_data : QHAInput
        Input data used to establish sampled-volume and normalization limits.
    numerical_tolerance : float, optional
        Absolute tolerance used for thermodynamic inequalities and zero-kelvin
        identities.

    Returns
    -------
    QHAValidationSummary
        Structured validation summary.

    Raises
    ------
    ValueError
        If the result does not contain pressure-temperature grids or
        equilibrium volumes.
    """
    if result.temperature is None or result.pressure is None:
        raise ValueError("QHA result does not contain temperature and pressure grids")
    if result.equilibrium_volume is None:
        raise ValueError("QHA result does not contain equilibrium volumes")

    temperature = np.asarray(result.temperature, dtype=np.float64)
    pressure = np.asarray(result.pressure, dtype=np.float64)
    volume = np.asarray(result.equilibrium_volume, dtype=np.float64)
    expected_shape = (temperature.size, pressure.size)
    if volume.shape != expected_shape:
        raise ValueError("equilibrium volume has an incompatible shape")

    valid = (
        np.asarray(result.valid_mask, dtype=bool)
        if result.valid_mask is not None
        else np.ones(expected_shape, dtype=bool)
    )
    if valid.shape != expected_shape:
        valid = np.ones(expected_shape, dtype=bool)

    checked_names = (
        "equilibrium_volume",
        "isothermal_bulk_modulus",
        "bulk_modulus_derivative",
        "adiabatic_bulk_modulus",
        "thermal_expansion",
        "isochoric_heat_capacity",
        "isobaric_heat_capacity",
        "heat_capacity_difference",
        "free_energy",
    )
    arrays: dict[str, np.ndarray] = {}
    finite_properties = True
    extrema: dict[str, tuple[float, float]] = {}
    for name in checked_names:
        value = getattr(result, name, None)
        if value is None:
            continue
        array = np.asarray(value, dtype=np.float64)
        if array.shape != expected_shape:
            continue
        arrays[name] = array
        selected = array[valid]
        if selected.size:
            finite_properties &= bool(np.all(np.isfinite(selected)))
            extrema[name] = (float(np.nanmin(selected)), float(np.nanmax(selected)))

    pressure_differences = np.diff(volume, axis=1)
    volume_decreases = bool(
        pressure.size < 2 or np.all(pressure_differences <= numerical_tolerance)
    )

    positive_bulk_moduli = True
    for name in ("isothermal_bulk_modulus", "adiabatic_bulk_modulus"):
        bulk_array = arrays.get(name)
        if bulk_array is not None:
            positive_bulk_moduli &= bool(np.all(bulk_array[valid] > 0.0))

    zero_kelvin_consistency: bool | None = None
    zero_indices = np.where(np.isclose(temperature, 0.0, atol=numerical_tolerance))[0]
    if zero_indices.size:
        index = int(zero_indices[0])
        checks: list[bool] = []
        kt = arrays.get("isothermal_bulk_modulus")
        ks = arrays.get("adiabatic_bulk_modulus")
        correction = arrays.get("heat_capacity_difference")
        if kt is not None and ks is not None:
            checks.append(
                bool(
                    np.allclose(
                        ks[index],
                        kt[index],
                        rtol=1.0e-10,
                        atol=numerical_tolerance,
                    )
                )
            )
        if correction is not None:
            checks.append(
                bool(np.all(np.abs(correction[index]) <= numerical_tolerance))
            )
        zero_kelvin_consistency = all(checks) if checks else None

    cp_not_below_cv: bool | None = None
    cp = arrays.get("isobaric_heat_capacity")
    cv = arrays.get("isochoric_heat_capacity")
    if cp is not None and cv is not None:
        cp_not_below_cv = bool(np.all(cp[valid] + numerical_tolerance >= cv[valid]))

    sampled = np.asarray(input_data.volume, dtype=np.float64)
    below = int(np.sum(volume[valid] < np.min(sampled))) if sampled.size else 0
    above = int(np.sum(volume[valid] > np.max(sampled))) if sampled.size else 0

    dulong_petit_ratio = _dulong_petit_ratio(result, input_data, valid)

    warnings: list[str] = []
    if int(np.sum(valid)) != int(valid.size):
        warnings.append("one or more pressure-temperature states are invalid")
    if not finite_properties:
        warnings.append("one or more result properties contain non-finite values")
    if not volume_decreases:
        warnings.append("equilibrium volume is not monotonic with pressure")
    if not positive_bulk_moduli:
        warnings.append("one or more bulk moduli are non-positive")
    if zero_kelvin_consistency is False:
        warnings.append("zero-kelvin thermodynamic identities are not satisfied")
    if cp_not_below_cv is False:
        warnings.append("Cp is below Cv beyond numerical tolerance")
    if below or above:
        warnings.append(
            f"{below + above} equilibrium volumes lie outside the sampled interval"
        )
    if dulong_petit_ratio is not None and not (0.90 <= dulong_petit_ratio <= 1.05):
        warnings.append(
            "the high-temperature Cv value is not close to the Dulong-Petit limit"
        )

    metadata = result.metadata if isinstance(result.metadata, dict) else {}
    return QHAValidationSummary(
        scheme=str(metadata.get("scheme", "unknown")),
        minimization=str(metadata.get("minimization", "unknown")),
        completed=bool(result.completed),
        total_points=int(valid.size),
        valid_points=int(np.sum(valid)),
        finite_properties=finite_properties,
        volume_decreases_with_pressure=volume_decreases,
        positive_bulk_moduli=positive_bulk_moduli,
        zero_kelvin_consistency=zero_kelvin_consistency,
        cp_not_below_cv=cp_not_below_cv,
        dulong_petit_ratio=dulong_petit_ratio,
        volumes_below_sampled_range=below,
        volumes_above_sampled_range=above,
        extrema=extrema,
        warnings=warnings,
    )


def compare_qha_results(
    first: QHAResult,
    second: QHAResult,
    *,
    properties: Iterable[str] = DEFAULT_COMPARISON_PROPERTIES,
) -> list[PropertyDifference]:
    """Compare corresponding property arrays from two QHA results.

    Parameters
    ----------
    first, second : QHAResult
        Results defined on the same pressure-temperature grid.
    properties : iterable of str, optional
        Result attributes to compare.

    Returns
    -------
    list of PropertyDifference
        Difference metrics for available compatible arrays.

    Raises
    ------
    ValueError
        If temperature or pressure grids differ.
    """
    _assert_matching_grids(first, second)
    differences: list[PropertyDifference] = []
    for name in properties:
        first_value = getattr(first, name, None)
        second_value = getattr(second, name, None)
        if first_value is None or second_value is None:
            continue
        first_array = np.asarray(first_value, dtype=np.float64)
        second_array = np.asarray(second_value, dtype=np.float64)
        if first_array.shape != second_array.shape:
            continue
        mask = np.isfinite(first_array) & np.isfinite(second_array)
        if not np.any(mask):
            continue
        difference_full = np.full(first_array.shape, np.nan, dtype=np.float64)
        denominator_full = np.full(first_array.shape, np.nan, dtype=np.float64)
        difference_full[mask] = np.abs(first_array[mask] - second_array[mask])
        denominator_full[mask] = np.maximum(
            np.abs(first_array[mask]),
            np.abs(second_array[mask]),
        )
        absolute = difference_full[mask]
        denominator = denominator_full[mask]
        scale = float(np.max(denominator))
        relative_floor = max(scale * 1.0e-8, float(np.finfo(np.float64).tiny))
        significant_full = mask & (denominator_full >= relative_floor)
        relative_full = np.full(first_array.shape, np.nan, dtype=np.float64)
        relative_full[significant_full] = (
            difference_full[significant_full] / denominator_full[significant_full]
        )
        absolute_index = np.unravel_index(
            int(np.nanargmax(difference_full)), difference_full.shape
        )
        if np.any(significant_full):
            relative_index = np.unravel_index(
                int(np.nanargmax(relative_full)), relative_full.shape
            )
            maximum_relative = float(relative_full[relative_index])
        else:
            relative_index = absolute_index
            maximum_relative = 0.0
        temperature = np.asarray(first.temperature, dtype=np.float64)
        pressure = np.asarray(first.pressure, dtype=np.float64)
        differences.append(
            PropertyDifference(
                property_name=name,
                maximum_absolute=float(difference_full[absolute_index]),
                maximum_relative=maximum_relative,
                root_mean_square=float(np.sqrt(np.mean(absolute**2))),
                compared_points=int(absolute.size),
                maximum_absolute_temperature=float(temperature[absolute_index[0]]),
                maximum_absolute_pressure=float(pressure[absolute_index[1]]),
                maximum_relative_temperature=float(temperature[relative_index[0]]),
                maximum_relative_pressure=float(pressure[relative_index[1]]),
            )
        )
    return differences


def _assert_matching_grids(first: QHAResult, second: QHAResult) -> None:
    """Validate that two results share pressure and temperature grids."""
    for name in ("temperature", "pressure"):
        first_value = getattr(first, name, None)
        second_value = getattr(second, name, None)
        if first_value is None or second_value is None:
            raise ValueError(f"both results must contain a {name} grid")
        if not np.array_equal(np.asarray(first_value), np.asarray(second_value)):
            raise ValueError(f"QHA {name} grids do not match")


def _dulong_petit_ratio(
    result: QHAResult,
    input_data: QHAInput,
    valid: np.ndarray,
) -> float | None:
    """Return the high-temperature Cv/Dulong-Petit ratio."""
    if result.isochoric_heat_capacity is None or result.temperature is None:
        return None
    cv = np.asarray(result.isochoric_heat_capacity, dtype=np.float64)
    temperature = np.asarray(result.temperature, dtype=np.float64)
    if cv.shape != valid.shape or not temperature.size:
        return None
    it = int(np.argmax(temperature))
    valid_pressures = np.where(valid[it])[0]
    if not valid_pressures.size:
        return None
    ip = int(valid_pressures[0])

    metadata = result.metadata if isinstance(result.metadata, dict) else {}
    units = (
        metadata.get("units", {}) if isinstance(metadata.get("units", {}), dict) else {}
    )
    source_unit = str(
        units.get(
            "heat_capacity",
            f"{units.get('energy', 'Ha')} cell^-1 K^-1",
        )
    )
    try:
        cv_j_mol_cell = float(
            convert_energy_per_temperature(
                cv[it, ip],
                source_unit,
                "J mol^-1 K^-1",
            )
        )
    except (ValueError, NotImplementedError):
        return None

    formula_units = max(int(input_data.formula_units), 1)
    if "mol" not in source_unit.lower():
        cv_j_mol_formula = cv_j_mol_cell / formula_units
    else:
        cv_j_mol_formula = cv_j_mol_cell
    atoms_per_formula = float(input_data.natoms) / formula_units
    limit = 3.0 * atoms_per_formula * gas_constant
    if limit <= 0.0:
        return None
    return cv_j_mol_formula / limit
