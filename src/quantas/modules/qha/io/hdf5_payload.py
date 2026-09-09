# -*- coding: utf-8 -*-

"""HDF5 payload readers and writers for QHA workflow results.

This module owns only the QHA-specific scientific payload stored under
``/results`` and module-specific fit diagnostics.  The generic Quantas HDF5
envelope is provided by :mod:`quantas.io.hdf5`.
"""

from __future__ import annotations

from typing import Any

import h5py
import numpy as np

from quantas.core.math.fitting import FitQuality, FitResult, FitStatus
from quantas.io.hdf5 import (
    decode_scalar,
    numeric_sort_key,
    read_group_mapping,
    read_node,
    write_array_dataset,
    write_mapping,
    write_value,
)
from quantas.models import InputData
from quantas.models.kieffer import KiefferThermodynamicContribution
from quantas.modules.qha.models import QHAFailedPoint, QHAFitRecord, QHAResult


QHA_DATASETS: dict[str, tuple[str, str]] = {
    "volume": ("V", "Sampled unit-cell volumes"),
    "static_energy": ("U0", "Static energies"),
    "temperature": ("T", "Temperature"),
    "pressure": ("P", "Pressure"),
    "equilibrium_volume": ("VT", "Equilibrium volume"),
    "zero_point_energy": ("Uzp", "Zero-point energy"),
    "thermal_energy": ("Uth", "Thermal energy"),
    "internal_energy": ("Utot", "Internal energy"),
    "entropy": ("S", "Entropy"),
    "vibrational_free_energy": ("Fvib", "Vibrational Helmholtz free energy"),
    "free_energy": ("F", "Helmholtz free energy"),
    "isochoric_heat_capacity": ("Cv", "Isochoric heat capacity"),
    "isobaric_heat_capacity": ("Cp", "Isobaric heat capacity"),
    "heat_capacity_difference": ("Cp-Cv", "Heat-capacity correction"),
    "isothermal_bulk_modulus": ("KT", "Isothermal bulk modulus"),
    "adiabatic_bulk_modulus": ("KS", "Adiabatic bulk modulus"),
    "bulk_modulus_derivative": ("Kp", "Bulk-modulus pressure derivative"),
    "thermal_expansion": ("alphaV", "Volumetric thermal expansion"),
    "thermal_expansion_mixed": (
        "alphaV_mixed",
        "Volumetric thermal expansion from mixed derivative",
    ),
    "thermal_expansion_mode": (
        "alphaV_mode",
        "Volumetric thermal expansion from mode Gruneisen values",
    ),
    "thermal_expansion_numerical": (
        "alphaV_numerical",
        "Volumetric thermal expansion from equilibrium-volume gradient",
    ),
    "thermal_expansion_source": (
        "alphaV_source",
        "Thermal-expansion source code",
    ),
    "equilibrium_lattice": (
        "lattice",
        "Equilibrium crystallographic direct-lattice matrices",
    ),
    "lattice_parameters": (
        "cell_parameters",
        "Equilibrium cell parameters a, b, c, alpha, beta, gamma",
    ),
    "lattice_parameter_derivatives": (
        "cell_parameter_dT",
        "Temperature derivatives of equilibrium cell parameters",
    ),
    "axial_thermal_expansion": (
        "alphaABC",
        "Axial thermal-expansion coefficients alpha_a, alpha_b, alpha_c",
    ),
    "thermal_expansion_tensor": (
        "alpha_tensor",
        "Cartesian thermal-expansion tensor",
    ),
    "structural_extrapolation_mask": (
        "structural_extrapolation",
        "Structural-path extrapolation mask",
    ),
    "enthalpy": ("H", "Enthalpy"),
    "gibbs_free_energy": ("G", "Gibbs free energy"),
    "gruneisen": ("gruneisen", "Thermodynamic Gruneisen parameter"),
    "mode_weighted_gruneisen": (
        "gamma_mode",
        "Heat-capacity-weighted mode Gruneisen parameter",
    ),
    "mode_gruneisen": ("gamma_qnu", "Mode-resolved Gruneisen parameters"),
}


def write_qha_payload(h5: h5py.File, result: QHAResult) -> h5py.Group:
    """Write the QHA-specific result payload.

    Parameters
    ----------
    h5 : h5py.File
        Open destination file.
    result : QHAResult
        Quasi-harmonic result payload.

    Returns
    -------
    h5py.Group
        Created ``results`` group.
    """
    group = h5.create_group("results")
    group.attrs["jobname"] = result.jobname
    group.attrs["completed"] = bool(result.completed)

    metadata_group = group.create_group("metadata")
    write_mapping(metadata_group, result.metadata)
    for attribute in QHA_DATASETS:
        value = getattr(result, attribute, None)
        if value is not None:
            write_array_dataset(group, attribute, value)
    if result.valid_mask is not None:
        write_array_dataset(group, "valid_mask", result.valid_mask)
    if result.kieffer_sampled_contribution is not None:
        _write_kieffer_contribution(group, result.kieffer_sampled_contribution)

    uncertainties = group.create_group("uncertainties")
    for key, value in result.uncertainties.items():
        write_array_dataset(uncertainties, key, value)
    return group


def write_qha_fit_diagnostics(diagnostics: h5py.Group, result: QHAResult) -> None:
    """Write QHA fit records and failed pressure-temperature points.

    Parameters
    ----------
    diagnostics : h5py.Group
        Generic diagnostics group created by the shared HDF5 envelope.
    result : QHAResult
        QHA result containing fit records and failed points.
    """
    fit_group = diagnostics.create_group("fit_records")
    for index, record in enumerate(result.diagnostics_as_dict()):
        write_value(fit_group, str(index), record)

    failed_group = diagnostics.create_group("failed_points")
    for index, point in enumerate(result.failed_points_as_dict()):
        write_value(failed_group, str(index), point)


def read_qha_payload(group: h5py.Group) -> QHAResult:
    """Read a QHA result from the modern ``results`` group.

    Parameters
    ----------
    group : h5py.Group
        HDF5 group containing QHA arrays.

    Returns
    -------
    QHAResult
        Reconstructed QHA result object.
    """
    reverse = {attr: attr for attr in QHA_DATASETS}
    arrays: dict[str, Any] = {}
    for attr in reverse:
        if attr in group:
            arrays[attr] = read_node(group[attr])
    if "valid_mask" in group:
        arrays["valid_mask"] = np.asarray(read_node(group["valid_mask"]), dtype=bool)

    metadata_group = (
        group["metadata"] if "metadata" in group else group.get("qha_metadata")
    )
    metadata = read_group_mapping(metadata_group) if metadata_group is not None else {}
    uncertainties = _read_uncertainties(group.file)
    fit_records, failed_points = _read_diagnostics(group.file)

    return QHAResult(
        jobname=str(decode_scalar(group.attrs.get("jobname", "Unknown"))),
        temperature=arrays.get("temperature"),
        pressure=arrays.get("pressure"),
        volume=arrays.get("volume"),
        static_energy=arrays.get("static_energy"),
        equilibrium_volume=arrays.get("equilibrium_volume"),
        zero_point_energy=arrays.get("zero_point_energy"),
        thermal_energy=arrays.get("thermal_energy"),
        internal_energy=arrays.get("internal_energy"),
        entropy=arrays.get("entropy"),
        vibrational_free_energy=arrays.get("vibrational_free_energy"),
        free_energy=arrays.get("free_energy"),
        isochoric_heat_capacity=arrays.get("isochoric_heat_capacity"),
        isobaric_heat_capacity=arrays.get("isobaric_heat_capacity"),
        heat_capacity_difference=arrays.get("heat_capacity_difference"),
        isothermal_bulk_modulus=arrays.get("isothermal_bulk_modulus"),
        adiabatic_bulk_modulus=arrays.get("adiabatic_bulk_modulus"),
        bulk_modulus_derivative=arrays.get("bulk_modulus_derivative"),
        thermal_expansion=arrays.get("thermal_expansion"),
        thermal_expansion_mixed=arrays.get("thermal_expansion_mixed"),
        thermal_expansion_mode=arrays.get("thermal_expansion_mode"),
        thermal_expansion_numerical=arrays.get("thermal_expansion_numerical"),
        thermal_expansion_source=arrays.get("thermal_expansion_source"),
        equilibrium_lattice=arrays.get("equilibrium_lattice"),
        lattice_parameters=arrays.get("lattice_parameters"),
        lattice_parameter_derivatives=arrays.get("lattice_parameter_derivatives"),
        axial_thermal_expansion=arrays.get("axial_thermal_expansion"),
        thermal_expansion_tensor=arrays.get("thermal_expansion_tensor"),
        structural_extrapolation_mask=arrays.get("structural_extrapolation_mask"),
        enthalpy=arrays.get("enthalpy"),
        gibbs_free_energy=arrays.get("gibbs_free_energy"),
        gruneisen=arrays.get("gruneisen"),
        mode_weighted_gruneisen=arrays.get("mode_weighted_gruneisen"),
        mode_gruneisen=arrays.get("mode_gruneisen"),
        kieffer_sampled_contribution=_read_kieffer_contribution(group),
        uncertainties=uncertainties,
        fit_records=fit_records,
        failed_points=failed_points,
        valid_mask=arrays.get("valid_mask"),
        completed=bool(decode_scalar(group.attrs.get("completed", True))),
        metadata=metadata,
    )


def _write_kieffer_contribution(
    group: h5py.Group,
    contribution: KiefferThermodynamicContribution,
) -> None:
    """Write the Kieffer component on the sampled QHA volume grid."""
    acoustic = group.create_group("kieffer_sampled_contribution")
    for name in (
        "cutoff_frequencies_hz",
        "effective_velocities_km_s",
        "zero_point_energy",
        "thermal_energy",
        "entropy",
        "vibrational_free_energy",
        "isochoric_heat_capacity",
    ):
        write_array_dataset(acoustic, name, getattr(contribution, name))
    metadata = acoustic.create_group("metadata")
    write_mapping(metadata, contribution.metadata)


def _read_kieffer_contribution(
    group: h5py.Group,
) -> KiefferThermodynamicContribution | None:
    """Read an optional sampled-volume Kieffer component."""
    acoustic = group.get("kieffer_sampled_contribution")
    if not isinstance(acoustic, h5py.Group):
        return None
    names = (
        "cutoff_frequencies_hz",
        "effective_velocities_km_s",
        "zero_point_energy",
        "thermal_energy",
        "entropy",
        "vibrational_free_energy",
        "isochoric_heat_capacity",
    )
    values = {name: np.asarray(read_node(acoustic[name])) for name in names}
    metadata = acoustic.get("metadata")
    return KiefferThermodynamicContribution(
        cutoff_frequencies_hz=values["cutoff_frequencies_hz"],
        effective_velocities_km_s=values["effective_velocities_km_s"],
        zero_point_energy=values["zero_point_energy"],
        thermal_energy=values["thermal_energy"],
        entropy=values["entropy"],
        vibrational_free_energy=values["vibrational_free_energy"],
        isochoric_heat_capacity=values["isochoric_heat_capacity"],
        metadata=(
            read_group_mapping(metadata) if isinstance(metadata, h5py.Group) else {}
        ),
    )


def migrate_schema_1_0_thermodynamic_units(
    result: QHAResult,
    schema_version: str,
    input_data: InputData | None,
) -> str | None:
    """Convert schema 1.0 thermal derivatives to true native units.

    Quantas 2.0 schema 1.0 stored entropy and heat capacities with numerical
    values equivalent to milli-Hartree, milli-electronvolt, or milli-Rydberg
    per kelvin while the metadata reported the corresponding unprefixed energy
    unit. Schema 1.1 stores the same quantities in true native energy units per
    cell and kelvin.

    Parameters
    ----------
    result : QHAResult
        QHA result reconstructed from HDF5.
    schema_version : str
        HDF5 schema version.
    input_data : InputData or None
        Stored input metadata used to reconstruct the normalization.

    Returns
    -------
    str or None
        Migration warning when a conversion was applied, otherwise ``None``.
    """
    metadata = result.metadata if isinstance(result.metadata, dict) else {}
    if metadata.get("thermodynamic_unit_convention"):
        return None
    if str(schema_version).strip() not in {"1", "1.0"}:
        return None

    arrays = (
        "entropy",
        "isochoric_heat_capacity",
        "isobaric_heat_capacity",
        "heat_capacity_difference",
    )
    converted = False
    for name in arrays:
        value = getattr(result, name, None)
        if value is None:
            continue
        setattr(result, name, np.asarray(value, dtype=np.float64) * 1.0e-3)
        converted = True
    if not converted:
        return None

    units = metadata.setdefault("units", {})
    energy_unit = str(units.get("energy", "Ha"))
    units["entropy"] = f"{energy_unit} cell^-1 K^-1"
    units["heat_capacity"] = f"{energy_unit} cell^-1 K^-1"
    metadata["thermodynamic_unit_convention"] = "native_energy_per_cell_per_kelvin"

    raw_input = input_data.data if input_data is not None else {}
    try:
        formula_units = int(raw_input.get("formula_units", 1))
    except (TypeError, ValueError):
        formula_units = 1
    if formula_units <= 0:
        formula_units = 1
    try:
        natoms = int(raw_input.get("natoms", 0))
    except (TypeError, ValueError):
        natoms = 0
    metadata.setdefault(
        "normalization",
        {
            "native_basis": "cell",
            "formula_units_per_cell": formula_units,
            "natoms_per_cell": natoms,
            "natoms_per_formula_unit": (
                float(natoms) / float(formula_units) if natoms > 0 else 0.0
            ),
            "molar_basis": "formula_unit",
        },
    )
    metadata["unit_migration"] = {
        "from_schema": str(schema_version),
        "scale_factor": 1.0e-3,
        "assumed_formula_units_per_cell": formula_units,
    }
    return (
        "Converted schema 1.0 entropy and heat-capacity arrays from the "
        "historical implicit milli-unit scale to true native energy units "
        "per cell and kelvin."
    )


def _read_uncertainties(h5: h5py.File) -> dict[str, np.ndarray]:
    """Read uncertainty arrays.

    Parameters
    ----------
    h5 : h5py.File
        Open HDF5 file.

    Returns
    -------
    dict
        Uncertainty arrays keyed by property name.
    """
    group = None
    if "results" in h5 and "uncertainties" in h5["results"]:
        group = h5["results/uncertainties"]
    elif "uncertainties" in h5:
        group = h5["uncertainties"]
    if group is None:
        return {}
    return {
        key: np.asarray(read_node(item))
        for key, item in group.items()
        if isinstance(item, h5py.Dataset)
    }


def _read_diagnostics(h5: h5py.File) -> tuple[list[QHAFitRecord], list[QHAFailedPoint]]:
    """Read structured diagnostics.

    Parameters
    ----------
    h5 : h5py.File
        Open HDF5 file.

    Returns
    -------
    tuple
        Fit records and failed points.
    """
    if "diagnostics" not in h5:
        return [], []
    group = h5["diagnostics"]
    fit_records: list[QHAFitRecord] = []
    failed_points: list[QHAFailedPoint] = []

    if "fit_records" in group:
        for key in sorted(group["fit_records"], key=numeric_sort_key):
            record = read_group_mapping(group["fit_records"][key])
            fit_records.append(_fit_record_from_dict(record))

    if "failed_points" in group:
        for key in sorted(group["failed_points"], key=numeric_sort_key):
            point = read_group_mapping(group["failed_points"][key])
            failed_points.append(_failed_point_from_dict(point))

    return fit_records, failed_points


def read_historical_qha_warnings(h5: h5py.File) -> list[str]:
    """Read warning messages from the log group.

    Parameters
    ----------
    h5 : h5py.File
        Open HDF5 file.

    Returns
    -------
    list of str
        Stored warning messages.
    """
    warnings_group = h5.get("log/warnings")
    if warnings_group is None:
        return []
    return [
        str(decode_scalar(read_node(warnings_group[key])))
        for key in sorted(warnings_group, key=numeric_sort_key)
    ]


def _fit_record_from_dict(record: dict[str, Any]) -> QHAFitRecord:
    """Build a fit diagnostic record from a dictionary.

    Parameters
    ----------
    record : dict
        Serialized diagnostic record.

    Returns
    -------
    QHAFitRecord
        Reconstructed diagnostic record.
    """
    fit_data = record.get("fit")
    fit = None if fit_data in (None, {}) else _fit_result_from_dict(fit_data)
    return QHAFitRecord(
        quantity=str(record.get("quantity", "")),
        method=str(record.get("method", "")),
        temperature=_none_or_float(record.get("temperature")),
        pressure=_none_or_float(record.get("pressure")),
        fit=fit,
        success=bool(record.get("success", True)),
        message=str(record.get("message", "")),
        metadata=dict(record.get("metadata", {}) or {}),
    )


def _failed_point_from_dict(point: dict[str, Any]) -> QHAFailedPoint:
    """Build a failed-point record from a dictionary.

    Parameters
    ----------
    point : dict
        Serialized failed-point record.

    Returns
    -------
    QHAFailedPoint
        Reconstructed failed-point record.
    """
    return QHAFailedPoint(
        temperature=float(point.get("temperature", np.nan)),
        pressure=float(point.get("pressure", np.nan)),
        stage=str(point.get("stage", "")),
        message=str(point.get("message", "")),
        diagnostics=dict(point.get("diagnostics", {}) or {}),
    )


def _fit_result_from_dict(data: dict[str, Any]) -> FitResult:
    """Build a fit result from a dictionary.

    Parameters
    ----------
    data : dict
        Serialized fit result.

    Returns
    -------
    FitResult
        Reconstructed fit result.
    """
    return FitResult(
        success=bool(data.get("success", False)),
        status=FitStatus(str(data.get("status", FitStatus.FAILED.value))),
        quality=FitQuality(str(data.get("quality", FitQuality.FAILED.value))),
        message=str(data.get("message", "")),
        parameters=_optional_array(data.get("parameters")),
        covariance=_optional_array(data.get("covariance")),
        errors=_optional_array(data.get("errors")),
        fitted=_optional_array(data.get("fitted")),
        residuals=_optional_array(data.get("residuals")),
        rmse=_none_or_float(data.get("rmse")),
        mae=_none_or_float(data.get("mae")),
        max_abs_error=_none_or_float(data.get("max_abs_error")),
        r_squared=_none_or_float(data.get("r_squared")),
        n_points=int(data.get("n_points", 0) or 0),
        n_parameters=int(data.get("n_parameters", 0) or 0),
        dof=int(data.get("dof", 0) or 0),
        condition_number=_none_or_float(data.get("condition_number")),
        warnings=[str(item) for item in data.get("warnings", [])],
        metadata=dict(data.get("metadata", {}) or {}),
    )


def _optional_array(value: Any) -> np.ndarray | None:
    """Return an optional NumPy array.

    Parameters
    ----------
    value : Any
        Value to convert.

    Returns
    -------
    ndarray or None
        Converted array or ``None``.
    """
    if value is None:
        return None

    if isinstance(value, dict):
        if not value:
            return np.asarray([], dtype=np.float64)
        keys = sorted(value, key=numeric_sort_key)
        values = [value[key] for key in keys]
    elif isinstance(value, list):
        values = value
    else:
        return np.asarray(value, dtype=np.float64)

    arrays: list[np.ndarray] = []
    for item in values:
        array = _optional_array(item)
        if array is not None:
            arrays.append(array)
    if not arrays:
        return np.asarray([], dtype=np.float64)
    try:
        return np.stack(arrays)
    except ValueError:
        return np.asarray(arrays, dtype=np.float64)


def _none_or_float(value: Any) -> float | None:
    """Return ``None`` or a finite float-like value.

    Parameters
    ----------
    value : Any
        Value to convert.

    Returns
    -------
    float or None
        Converted float, or ``None``.
    """
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
