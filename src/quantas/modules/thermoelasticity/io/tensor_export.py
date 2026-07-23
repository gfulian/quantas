# -*- coding: utf-8 -*-

"""Deterministic text and CSV exports of reconstructed thermoelastic tensors."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Iterable, Iterator, Literal

import numpy as np

from quantas.models import ReportTable, ResultData
from quantas.modules.thermoelasticity.models import (
    ThermoelasticProfileResult,
    ThermoelasticResult,
)
from quantas.modules.thermoelasticity.postfit import grid_step


TensorExportCondition = Literal["isothermal", "adiabatic", "both"]
_UPPER_VOIGT = tuple((i, j) for i in range(6) for j in range(i, 6))


def thermoelastic_grid_info_table(result: ThermoelasticResult) -> ReportTable:
    """Return the archived pressure-temperature coverage summary."""
    temperature_step = grid_step(result.temperature)
    pressure_step = grid_step(result.pressure)
    rows: list[list[Any]] = [
        [
            "Temperature",
            float(result.temperature[0]),
            float(result.temperature[-1]),
            result.temperature.size,
            temperature_step,
            "K",
        ],
        [
            "Pressure",
            float(result.pressure[0]),
            float(result.pressure[-1]),
            result.pressure.size,
            pressure_step,
            "GPa",
        ],
        [
            "Equilibrium volume",
            float(np.nanmin(result.equilibrium_volume)),
            float(np.nanmax(result.equilibrium_volume)),
            result.equilibrium_volume.size,
            None,
            "angstrom^3",
        ],
        [
            "Adiabatic tensor",
            None,
            None,
            (
                0
                if result.adiabatic_valid_mask is None
                else int(np.count_nonzero(result.adiabatic_valid_mask))
            ),
            None,
            "valid states",
        ],
    ]
    return ReportTable(
        "Thermoelastic archive coverage",
        ["Coordinate", "Minimum", "Maximum", "Count", "Uniform step", "Unit"],
        rows,
        metadata={
            "notes": [
                "A blank step denotes a singleton or nonuniform coordinate axis.",
                "Adiabatic tensors require QHA C_V and the Cartesian thermal-"
                "expansion tensor at nonzero temperature.",
            ]
        },
    )


def write_thermoelastic_tensor_export(
    result_data: ResultData,
    filename: str | Path,
    *,
    file_format: str = "txt",
    tensor_condition: TensorExportCondition = "isothermal",
    include_uncertainties: bool = True,
    include_grid: bool = True,
    include_profiles: bool = True,
    overwrite: bool = False,
) -> Path:
    """Write reconstructed tensors from one thermoelastic result archive.

    Parameters
    ----------
    result_data : ResultData
        Result envelope containing a completed thermoelastic payload.
    filename : str or Path
        Destination path.
    file_format : {"txt", "csv"}, optional
        Export format.
    tensor_condition : {"isothermal", "adiabatic", "both"}, optional
        Thermodynamic tensor condition to export.
    include_uncertainties : bool, optional
        Include one-sigma tensor uncertainties where available.
    include_grid, include_profiles : bool, optional
        Select rectangular-grid and depth-profile sections.
    overwrite : bool, optional
        Permit replacement of an existing file.

    Returns
    -------
    Path
        Written destination.
    """
    payload = result_data.results.get("thermoelasticity")
    if not isinstance(payload, ThermoelasticResult):
        raise ValueError("result does not contain a thermoelasticity payload")
    condition = _validated_condition(tensor_condition)
    path = Path(filename)
    if path.exists() and not overwrite:
        raise FileExistsError(path)
    normalized = file_format.strip().lower()
    if normalized not in {"txt", "csv"}:
        raise ValueError("file_format must be 'txt' or 'csv'")
    if normalized == "txt":
        path.write_text(
            format_thermoelastic_tensors_text(
                payload,
                tensor_condition=condition,
                include_uncertainties=include_uncertainties,
                include_grid=include_grid,
                include_profiles=include_profiles,
            ),
            encoding="utf-8",
        )
    else:
        _write_tensor_csv(
            payload,
            path,
            tensor_condition=condition,
            include_uncertainties=include_uncertainties,
            include_grid=include_grid,
            include_profiles=include_profiles,
        )
    return path


def format_thermoelastic_tensors_text(
    result: ThermoelasticResult,
    *,
    tensor_condition: TensorExportCondition = "isothermal",
    include_uncertainties: bool = True,
    include_grid: bool = True,
    include_profiles: bool = True,
) -> str:
    """Return a deterministic plain-text tensor export."""
    conditions = _selected_conditions(result, tensor_condition)
    _validate_sections(result, include_grid, include_profiles)
    lines = [
        "# Quantas quasi-static thermoelastic tensor export",
        f"# Tensor condition: {tensor_condition}",
        "# Stiffness unit: GPa",
        "# Voigt order: 11 22 33 23 13 12",
        "",
    ]
    qha_mask = np.asarray(result.qha_extrapolation_mask, dtype=np.bool_)
    if include_grid:
        for condition in conditions:
            stiffness, sigma, valid = _grid_tensor_arrays(result, condition)
            for it, temperature in enumerate(result.temperature):
                for ip, pressure in enumerate(result.pressure):
                    lines.extend(
                        _format_state_block(
                            kind="grid",
                            identifier=f"T{it + 1:04d}_P{ip + 1:04d}",
                            tensor_condition=condition,
                            pressure=float(pressure),
                            temperature=float(temperature),
                            volume=float(result.equilibrium_volume[it, ip]),
                            density=float(result.density[it, ip]),
                            stiffness=stiffness[it, ip],
                            sigma=(
                                None
                                if sigma is None or not include_uncertainties
                                else sigma[it, ip]
                            ),
                            tensor_valid=(
                                True if valid is None else bool(valid[it, ip])
                            ),
                            qha_extrapolated=bool(qha_mask[it, ip]),
                            elastic_extrapolated=bool(
                                result.extrapolation_mask[it, ip]
                            ),
                            mechanically_stable=(
                                None
                                if result.stability is None
                                else bool(result.stability.stable_mask[it, ip])
                            ),
                            stability_indeterminate=(
                                None
                                if result.stability is None
                                else bool(result.stability.indeterminate_mask[it, ip])
                            ),
                            minimum_stiffness_eigenvalue=(
                                None
                                if result.stability is None
                                else float(result.stability.minimum_eigenvalue[it, ip])
                            ),
                            depth=None,
                            profile=None,
                        )
                    )
    if include_profiles:
        for profile_name, profile in result.profiles.items():
            for condition in conditions:
                stiffness, sigma, valid = _profile_tensor_arrays(profile, condition)
                for index in range(profile.depth.size):
                    lines.extend(
                        _format_state_block(
                            kind="profile",
                            identifier=f"{profile_name}:{index + 1:04d}",
                            tensor_condition=condition,
                            pressure=float(profile.pressure[index]),
                            temperature=float(profile.temperature[index]),
                            volume=float(profile.volume[index]),
                            density=float(profile.density[index]),
                            stiffness=stiffness[index],
                            sigma=(
                                None
                                if sigma is None or not include_uncertainties
                                else sigma[index]
                            ),
                            tensor_valid=(
                                True if valid is None else bool(valid[index])
                            ),
                            qha_extrapolated=bool(
                                profile.qha_extrapolation_mask[index]
                            ),
                            elastic_extrapolated=bool(
                                profile.elastic_extrapolation_mask[index]
                            ),
                            mechanically_stable=(
                                None
                                if profile.stability is None
                                else bool(profile.stability.stable_mask[index])
                            ),
                            stability_indeterminate=(
                                None
                                if profile.stability is None
                                else bool(profile.stability.indeterminate_mask[index])
                            ),
                            minimum_stiffness_eigenvalue=(
                                None
                                if profile.stability is None
                                else float(profile.stability.minimum_eigenvalue[index])
                            ),
                            depth=float(profile.depth[index]),
                            profile=profile_name,
                        )
                    )
    return "\n".join(lines).rstrip() + "\n"


def _format_state_block(
    *,
    kind: str,
    identifier: str,
    tensor_condition: str,
    pressure: float,
    temperature: float,
    volume: float,
    density: float,
    stiffness: np.ndarray,
    sigma: np.ndarray | None,
    tensor_valid: bool,
    qha_extrapolated: bool,
    elastic_extrapolated: bool,
    mechanically_stable: bool | None,
    stability_indeterminate: bool | None,
    minimum_stiffness_eigenvalue: float | None,
    depth: float | None,
    profile: str | None,
) -> list[str]:
    label = (
        "C_isothermal_GPa" if tensor_condition == "isothermal" else "C_adiabatic_GPa"
    )
    sigma_label = f"sigma_{label}"
    lines = [
        f"## {identifier} [{tensor_condition}]",
        f"kind = {kind}",
        f"tensor_condition = {tensor_condition}",
        f"tensor_valid = {str(tensor_valid).lower()}",
        f"profile = {'' if profile is None else profile}",
        f"depth_km = {'' if depth is None else f'{depth:.10f}'}",
        f"pressure_GPa = {pressure:.12f}",
        f"temperature_K = {temperature:.12f}",
        f"volume_A3 = {volume:.12f}",
        f"density_kg_m3 = {density:.12f}",
        f"qha_extrapolated = {str(qha_extrapolated).lower()}",
        f"elastic_extrapolated = {str(elastic_extrapolated).lower()}",
        "mechanically_stable = "
        + ("" if mechanically_stable is None else str(mechanically_stable).lower()),
        "stability_indeterminate = "
        + (
            ""
            if stability_indeterminate is None
            else str(stability_indeterminate).lower()
        ),
        "minimum_stiffness_eigenvalue_GPa = "
        + (
            ""
            if minimum_stiffness_eigenvalue is None
            or not np.isfinite(minimum_stiffness_eigenvalue)
            else f"{minimum_stiffness_eigenvalue:.12f}"
        ),
        f"{label}:",
    ]
    lines.extend(_matrix_lines(stiffness))
    if sigma is not None:
        lines.append(f"{sigma_label}:")
        lines.extend(_matrix_lines(sigma))
    lines.append("")
    return lines


def _matrix_lines(matrix: np.ndarray) -> list[str]:
    values = np.asarray(matrix, dtype=np.float64)
    if values.shape != (6, 6):
        raise ValueError("stiffness matrix must have shape (6, 6)")
    return ["  " + " ".join(f"{value:16.8f}" for value in row) for row in values]


def _write_tensor_csv(
    result: ThermoelasticResult,
    path: Path,
    *,
    tensor_condition: TensorExportCondition,
    include_uncertainties: bool,
    include_grid: bool,
    include_profiles: bool,
) -> None:
    fieldnames = [
        "kind",
        "identifier",
        "tensor_condition",
        "tensor_valid",
        "profile",
        "depth_km",
        "pressure_GPa",
        "temperature_K",
        "volume_A3",
        "density_kg_m3",
        "qha_extrapolated",
        "elastic_extrapolated",
        "mechanically_stable",
        "stability_indeterminate",
        "minimum_stiffness_eigenvalue_GPa",
    ]
    for i, j in _UPPER_VOIGT:
        label = f"C{i + 1}{j + 1}_GPa"
        fieldnames.append(label)
        if include_uncertainties:
            fieldnames.append(f"sigma_{label}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
                stream, 
                fieldnames=fieldnames,
                lineterminator="\n",
                )
        writer.writeheader()
        writer.writerows(
            _iter_export_rows(
                result,
                tensor_condition=tensor_condition,
                include_uncertainties=include_uncertainties,
                include_grid=include_grid,
                include_profiles=include_profiles,
            )
        )


def _iter_export_rows(
    result: ThermoelasticResult,
    *,
    tensor_condition: TensorExportCondition,
    include_uncertainties: bool,
    include_grid: bool,
    include_profiles: bool,
) -> Iterator[dict[str, object]]:
    conditions = _selected_conditions(result, tensor_condition)
    _validate_sections(result, include_grid, include_profiles)
    qha_mask = np.asarray(result.qha_extrapolation_mask, dtype=np.bool_)
    if include_grid:
        for condition in conditions:
            stiffness, sigma, valid = _grid_tensor_arrays(result, condition)
            for it, temperature in enumerate(result.temperature):
                for ip, pressure in enumerate(result.pressure):
                    yield _csv_row(
                        kind="grid",
                        identifier=f"T{it + 1:04d}_P{ip + 1:04d}",
                        tensor_condition=condition,
                        tensor_valid=True if valid is None else bool(valid[it, ip]),
                        profile=None,
                        depth=None,
                        pressure=float(pressure),
                        temperature=float(temperature),
                        volume=float(result.equilibrium_volume[it, ip]),
                        density=float(result.density[it, ip]),
                        stiffness=stiffness[it, ip],
                        sigma=(
                            None
                            if sigma is None or not include_uncertainties
                            else sigma[it, ip]
                        ),
                        qha_extrapolated=bool(qha_mask[it, ip]),
                        elastic_extrapolated=bool(result.extrapolation_mask[it, ip]),
                        mechanically_stable=(
                            None
                            if result.stability is None
                            else bool(result.stability.stable_mask[it, ip])
                        ),
                        stability_indeterminate=(
                            None
                            if result.stability is None
                            else bool(result.stability.indeterminate_mask[it, ip])
                        ),
                        minimum_stiffness_eigenvalue=(
                            None
                            if result.stability is None
                            else float(result.stability.minimum_eigenvalue[it, ip])
                        ),
                        include_uncertainties=include_uncertainties,
                    )
    if include_profiles:
        for name, profile in result.profiles.items():
            for condition in conditions:
                yield from _profile_csv_rows(
                    name,
                    profile,
                    condition=condition,
                    include_uncertainties=include_uncertainties,
                )


def _profile_csv_rows(
    name: str,
    profile: ThermoelasticProfileResult,
    *,
    condition: str,
    include_uncertainties: bool,
) -> Iterable[dict[str, object]]:
    stiffness, sigma, valid = _profile_tensor_arrays(profile, condition)
    for index in range(profile.depth.size):
        yield _csv_row(
            kind="profile",
            identifier=f"{name}:{index + 1:04d}",
            tensor_condition=condition,
            tensor_valid=True if valid is None else bool(valid[index]),
            profile=name,
            depth=float(profile.depth[index]),
            pressure=float(profile.pressure[index]),
            temperature=float(profile.temperature[index]),
            volume=float(profile.volume[index]),
            density=float(profile.density[index]),
            stiffness=stiffness[index],
            sigma=(
                sigma[index] if sigma is not None and include_uncertainties else None
            ),
            qha_extrapolated=bool(profile.qha_extrapolation_mask[index]),
            elastic_extrapolated=bool(profile.elastic_extrapolation_mask[index]),
            mechanically_stable=(
                None
                if profile.stability is None
                else bool(profile.stability.stable_mask[index])
            ),
            stability_indeterminate=(
                None
                if profile.stability is None
                else bool(profile.stability.indeterminate_mask[index])
            ),
            minimum_stiffness_eigenvalue=(
                None
                if profile.stability is None
                else float(profile.stability.minimum_eigenvalue[index])
            ),
            include_uncertainties=include_uncertainties,
        )


def _csv_row(
    *,
    kind: str,
    identifier: str,
    tensor_condition: str,
    tensor_valid: bool,
    profile: str | None,
    depth: float | None,
    pressure: float,
    temperature: float,
    volume: float,
    density: float,
    stiffness: np.ndarray,
    sigma: np.ndarray | None,
    qha_extrapolated: bool,
    elastic_extrapolated: bool,
    mechanically_stable: bool | None,
    stability_indeterminate: bool | None,
    minimum_stiffness_eigenvalue: float | None,
    include_uncertainties: bool,
) -> dict[str, object]:
    row: dict[str, object] = {
        "kind": kind,
        "identifier": identifier,
        "tensor_condition": tensor_condition,
        "tensor_valid": tensor_valid,
        "profile": "" if profile is None else profile,
        "depth_km": "" if depth is None else depth,
        "pressure_GPa": pressure,
        "temperature_K": temperature,
        "volume_A3": volume,
        "density_kg_m3": density,
        "qha_extrapolated": qha_extrapolated,
        "elastic_extrapolated": elastic_extrapolated,
        "mechanically_stable": (
            "" if mechanically_stable is None else mechanically_stable
        ),
        "stability_indeterminate": (
            "" if stability_indeterminate is None else stability_indeterminate
        ),
        "minimum_stiffness_eigenvalue_GPa": (
            "" if minimum_stiffness_eigenvalue is None else minimum_stiffness_eigenvalue
        ),
    }
    values = np.asarray(stiffness, dtype=np.float64)
    sigma_values = None if sigma is None else np.asarray(sigma, dtype=np.float64)
    for i, j in _UPPER_VOIGT:
        label = f"C{i + 1}{j + 1}_GPa"
        row[label] = float(values[i, j])
        if include_uncertainties:
            row[f"sigma_{label}"] = (
                "" if sigma_values is None else float(sigma_values[i, j])
            )
    return row


def _validated_condition(value: str) -> TensorExportCondition:
    normalized = value.strip().lower()
    if normalized not in {"isothermal", "adiabatic", "both"}:
        raise ValueError(
            "tensor_condition must be 'isothermal', 'adiabatic', or 'both'"
        )
    return normalized  # type: ignore[return-value]


def _selected_conditions(
    result: ThermoelasticResult,
    value: str,
) -> tuple[str, ...]:
    condition = _validated_condition(value)
    if condition == "isothermal":
        return ("isothermal",)
    if result.stiffness_adiabatic is None:
        raise ValueError("archive does not contain adiabatic stiffness tensors")
    return (
        ("adiabatic",)
        if condition == "adiabatic"
        else (
            "isothermal",
            "adiabatic",
        )
    )


def _grid_tensor_arrays(
    result: ThermoelasticResult,
    condition: str,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None]:
    if condition == "isothermal":
        if result.stiffness_isothermal is None:
            raise ValueError("archive does not contain isothermal grid tensors")
        return (
            result.stiffness_isothermal,
            result.sigma_stiffness_isothermal,
            None,
        )
    if result.stiffness_adiabatic is None:
        raise ValueError("archive does not contain adiabatic grid tensors")
    return (
        result.stiffness_adiabatic,
        result.sigma_stiffness_adiabatic,
        result.adiabatic_valid_mask,
    )


def _profile_tensor_arrays(
    result: ThermoelasticProfileResult,
    condition: str,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None]:
    if condition == "isothermal":
        return result.stiffness_isothermal, result.sigma_stiffness_isothermal, None
    if result.stiffness_adiabatic is None:
        raise ValueError(f"profile '{result.name}' does not contain adiabatic tensors")
    return (
        result.stiffness_adiabatic,
        result.sigma_stiffness_adiabatic,
        result.adiabatic_valid_mask,
    )


def _validate_sections(
    result: ThermoelasticResult,
    include_grid: bool,
    include_profiles: bool,
) -> None:
    if not include_grid and not include_profiles:
        raise ValueError("at least one export section must be selected")
    if include_profiles and not result.profiles and not include_grid:
        raise ValueError("archive does not contain depth-profile tensors")


__all__ = [
    "TensorExportCondition",
    "format_thermoelastic_tensors_text",
    "thermoelastic_grid_info_table",
    "write_thermoelastic_tensor_export",
]
