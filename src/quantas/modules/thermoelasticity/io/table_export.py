# -*- coding: utf-8 -*-

"""Compact human-oriented tables for thermoelastic grids and profiles."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Literal

import numpy as np

from quantas.models import ReportTable
from quantas.modules.thermoelasticity.models import (
    ThermoelasticProfileResult,
    ThermoelasticResult,
)
from quantas.modules.thermoelasticity.components import component_indices

TableFormat = Literal["csv", "text"]
TensorSelection = Literal["isothermal", "adiabatic", "both"]


def thermoelastic_point_table(
    result: ThermoelasticResult,
    *,
    tensor_condition: Literal["isothermal", "adiabatic"] = "adiabatic",
) -> ReportTable:
    """Return a compact human-readable table for one reconstructed state."""
    if result.temperature.size != 1 or result.pressure.size != 1:
        raise ValueError("point table requires exactly one pressure-temperature state")
    if tensor_condition == "isothermal":
        tensor = result.stiffness_isothermal
        sigma = result.sigma_stiffness_isothermal
        suffix = "T"
    elif tensor_condition == "adiabatic":
        tensor = result.stiffness_adiabatic
        sigma = result.sigma_stiffness_adiabatic
        suffix = "S"
    else:
        raise ValueError("tensor_condition must be isothermal or adiabatic")
    if tensor is None:
        raise ValueError(f"{tensor_condition} stiffness is unavailable")
    rows: list[list[Any]] = [
        ["Pressure", f"{float(result.pressure[0]):.4f} GPa"],
        ["Temperature", f"{float(result.temperature[0]):.2f} K"],
        ["Volume", f"{float(result.equilibrium_volume[0, 0]):.4f} Å³"],
        ["Density", f"{float(result.density[0, 0]):.2f} kg m^-3"],
    ]
    for label in result.independent_labels:
        i, j = component_indices(label)
        value = float(tensor[0, 0, i, j])
        if sigma is None or not np.isfinite(sigma[0, 0, i, j]):
            text = f"{value:.4f} GPa"
        else:
            error = float(sigma[0, 0, i, j])
            text = f"{value:.4f} +/- {error:.4f} GPa"
        rows.append([f"{label}_{suffix}", text])
    rows.extend(
        [
            [
                "Mechanically stable",
                ""
                if result.stability is None
                else bool(result.stability.stable_mask[0, 0]),
            ],
            [
                "QHA extrapolated",
                False
                if result.qha_extrapolation_mask is None
                else bool(result.qha_extrapolation_mask[0, 0]),
            ],
            ["Elastic extrapolated", bool(result.extrapolation_mask[0, 0])],
        ]
    )
    return ReportTable(
        f"Thermoelastic point ({tensor_condition})",
        ["Property", "Value"],
        rows,
    )


def write_thermoelastic_grid_table(
    result: ThermoelasticResult,
    filename: str | Path,
    *,
    tensor_condition: TensorSelection = "both",
    file_format: TableFormat = "csv",
    include_uncertainties: bool = True,
    overwrite: bool = False,
) -> Path:
    """Write one wide row per pressure-temperature grid state.

    Only symmetry-independent stiffness components are exported.  Isothermal
    and adiabatic fields occupy separate columns when ``tensor_condition`` is
    ``"both"``; rows are never duplicated.
    """
    if result.stiffness_isothermal is None:
        raise ValueError("archive does not contain a reconstructed P-T grid")
    records = [
        _grid_record(result, it, ip, tensor_condition, include_uncertainties)
        for it in range(result.temperature.size)
        for ip in range(result.pressure.size)
    ]
    return _write_records(
        records, filename, file_format=file_format, overwrite=overwrite
    )


def write_thermoelastic_profile_table(
    result: ThermoelasticResult,
    filename: str | Path,
    *,
    profile_name: str | None = None,
    tensor_condition: TensorSelection = "both",
    file_format: TableFormat = "csv",
    include_uncertainties: bool = True,
    overwrite: bool = False,
) -> Path:
    """Write one wide row per depth along one archived geological profile.

    Parameters
    ----------
    result : ThermoelasticResult
        Thermoelastic payload containing one or more evaluated profiles.
    filename : str or Path
        Destination CSV or deterministic tab-separated text file.
    profile_name : str or None, optional
        Profile to export.  Omission is accepted only when exactly one profile
        is present.
    tensor_condition : {"isothermal", "adiabatic", "both"}, optional
        Stiffness fields written as columns.
    file_format : {"csv", "text"}, optional
        Portable comma-separated or tab-separated representation.
    include_uncertainties : bool, optional
        Include one-sigma columns when available.
    overwrite : bool, optional
        Replace an existing destination.

    Returns
    -------
    Path
        Written path.
    """
    name, profile = _select_profile(result, profile_name)
    records = [
        _profile_record(
            result,
            name,
            profile,
            index,
            tensor_condition,
            include_uncertainties,
        )
        for index in range(profile.depth.size)
    ]
    return _write_records(
        records, filename, file_format=file_format, overwrite=overwrite
    )


def _grid_record(
    result: ThermoelasticResult,
    it: int,
    ip: int,
    condition: TensorSelection,
    uncertainties: bool,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "pressure_GPa": float(result.pressure[ip]),
        "temperature_K": float(result.temperature[it]),
        "volume_A3": float(result.equilibrium_volume[it, ip]),
        "density_kg_m3": float(result.density[it, ip]),
    }
    _add_components(
        record,
        result.independent_labels,
        result.stiffness_isothermal,
        result.sigma_stiffness_isothermal,
        (it, ip),
        condition,
        uncertainties,
        adiabatic=result.stiffness_adiabatic,
        sigma_adiabatic=result.sigma_stiffness_adiabatic,
    )
    record.update(
        {
            "qha_extrapolated": (
                ""
                if result.qha_extrapolation_mask is None
                else bool(result.qha_extrapolation_mask[it, ip])
            ),
            "elastic_extrapolated": bool(result.extrapolation_mask[it, ip]),
            "adiabatic_valid": (
                ""
                if result.adiabatic_valid_mask is None
                else bool(result.adiabatic_valid_mask[it, ip])
            ),
            "mechanically_stable": (
                ""
                if result.stability is None
                else bool(result.stability.stable_mask[it, ip])
            ),
            "minimum_stiffness_eigenvalue_GPa": (
                ""
                if result.stability is None
                else float(result.stability.minimum_eigenvalue[it, ip])
            ),
        }
    )
    return record


def _profile_record(
    result: ThermoelasticResult,
    profile_name: str,
    profile: ThermoelasticProfileResult,
    index: int,
    condition: TensorSelection,
    uncertainties: bool,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "profile": profile_name,
        "depth_km": float(profile.depth[index]),
        "pressure_GPa": float(profile.pressure[index]),
        "temperature_K": float(profile.temperature[index]),
        "volume_A3": float(profile.volume[index]),
        "density_kg_m3": float(profile.density[index]),
    }
    _add_components(
        record,
        result.independent_labels,
        profile.stiffness_isothermal,
        profile.sigma_stiffness_isothermal,
        (index,),
        condition,
        uncertainties,
        adiabatic=profile.stiffness_adiabatic,
        sigma_adiabatic=profile.sigma_stiffness_adiabatic,
    )
    record.update(
        {
            "qha_extrapolated": bool(profile.qha_extrapolation_mask[index]),
            "elastic_extrapolated": bool(profile.elastic_extrapolation_mask[index]),
            "adiabatic_valid": (
                ""
                if profile.adiabatic_valid_mask is None
                else bool(profile.adiabatic_valid_mask[index])
            ),
            "mechanically_stable": (
                ""
                if profile.stability is None
                else bool(profile.stability.stable_mask[index])
            ),
            "minimum_stiffness_eigenvalue_GPa": (
                ""
                if profile.stability is None
                else float(profile.stability.minimum_eigenvalue[index])
            ),
        }
    )
    return record


def _add_components(
    record: dict[str, Any],
    labels: tuple[str, ...],
    isothermal: np.ndarray | None,
    sigma_isothermal: np.ndarray | None,
    index: tuple[int, ...],
    condition: TensorSelection,
    uncertainties: bool,
    *,
    adiabatic: np.ndarray | None,
    sigma_adiabatic: np.ndarray | None,
) -> None:
    if condition not in {"isothermal", "adiabatic", "both"}:
        raise ValueError("invalid tensor condition")
    selections = []
    if condition in {"isothermal", "both"}:
        if isothermal is None:
            raise ValueError("isothermal stiffness is unavailable")
        selections.append(("T", isothermal, sigma_isothermal))
    if condition in {"adiabatic", "both"}:
        if adiabatic is None:
            if condition == "adiabatic":
                raise ValueError("adiabatic stiffness is unavailable")
        else:
            selections.append(("S", adiabatic, sigma_adiabatic))
    for suffix, stiffness, sigma in selections:
        for label in labels:
            i, j = component_indices(label)
            record[f"{label}_{suffix}_GPa"] = float(stiffness[index + (i, j)])
            if uncertainties and sigma is not None:
                record[f"sigma_{label}_{suffix}_GPa"] = float(sigma[index + (i, j)])


def _select_profile(
    result: ThermoelasticResult,
    profile_name: str | None,
) -> tuple[str, ThermoelasticProfileResult]:
    if not result.profiles:
        raise ValueError("archive does not contain evaluated depth profiles")
    if profile_name is None:
        if len(result.profiles) != 1:
            available = ", ".join(sorted(result.profiles))
            raise ValueError(
                "profile name is required because the archive contains: " + available
            )
        profile_name = next(iter(result.profiles))
    try:
        return profile_name, result.profiles[profile_name]
    except KeyError as exc:
        raise ValueError(f"profile not present in archive: {profile_name}") from exc


def _write_records(
    records: list[dict[str, Any]],
    filename: str | Path,
    *,
    file_format: TableFormat,
    overwrite: bool,
) -> Path:
    if not records:
        raise ValueError("no thermoelastic table rows were selected")
    path = Path(filename)
    if path.exists() and not overwrite:
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    delimiter = "," if file_format == "csv" else "\t"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(records[0]), delimiter=delimiter
        )
        writer.writeheader()
        writer.writerows(records)
    return path


__all__ = [
    "thermoelastic_point_table",
    "write_thermoelastic_grid_table",
    "write_thermoelastic_profile_table",
]
