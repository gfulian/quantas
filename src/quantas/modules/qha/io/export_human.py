# -*- coding: utf-8 -*-

"""Deterministic human-readable rendering for QHA table exports."""

from __future__ import annotations

from typing import Any, TextIO

import numpy as np

from quantas.modules.qha.formatting import QHATableFormat
from quantas.modules.qha.models import QHAResult
from quantas.modules.qha.io.export_schema import (
    ENERGY_CONTRIBUTION_PROPERTIES,
    ENERGY_PROPERTIES,
    GRUNEISEN_PROPERTIES,
    THERMAL_PROPERTIES,
    THERMODYNAMIC_POTENTIAL_PROPERTIES,
    THERMOELASTIC_PROPERTIES,
    _StructuralExportData,
    _format_for_property,
    _format_number,
    _unit_for_property,
)

def _write_human_table(
    stream: TextIO,
    result: QHAResult,
    arrays: list[tuple[str, str, str, np.ndarray, np.ndarray | None]],
    structural: _StructuralExportData | None,
    temperature: np.ndarray,
    pressure: np.ndarray,
    temperature_unit: str,
    pressure_unit: str,
    table_format: QHATableFormat,
    thermal_expansion_source: np.ndarray | None,
) -> None:
    """Write compact fixed-width tables grouped by pressure.

    Parameters
    ----------
    stream : file-like object
        Open text stream.
    result : QHAResult
        QHA result object.
    arrays : list of tuple
        Prepared property arrays and optional uncertainty arrays.
    structural : _StructuralExportData or None
        Optional structural arrays.
    temperature, pressure : ndarray
        Pressure-temperature axes.
    temperature_unit, pressure_unit : str
        Axis unit labels.
    table_format : QHATableFormat
        Numerical formatting and fixed-width configuration.
    """
    stream.write("# Quantas QHA human-readable table export\n")
    stream.write("# Values are grouped by pressure.\n")
    stream.write(
        "# Positive scientific values reserve one leading blank for sign alignment.\n"
    )
    has_structural_uncertainty = structural is not None and (
        structural.lattice_parameter_uncertainties is not None
        or structural.axial_thermal_expansion_uncertainties is not None
    )
    if (
        any(sigma is not None for _, _, _, _, sigma in arrays)
        or has_structural_uncertainty
    ):
        stream.write(
            "# sigma_* columns contain one-standard-deviation uncertainties.\n"
        )
        uncertainty_method = (
            result.metadata.get("eos_workflow", {}).get("uncertainty_method")
            if isinstance(result.metadata, dict)
            else None
        )
        if uncertainty_method not in {None, "none"}:
            stream.write(f"# Uncertainty method: {uncertainty_method}\n")
        structural_method = result.metadata.get(
            "structural_thermal_expansion",
            {},
        ).get("uncertainty_method")
        if structural_method not in {None, "none"}:
            stream.write(f"# Structural uncertainty method: {structural_method}\n")
    _write_thermal_expansion_provenance(stream, result)
    stream.write("\n")

    groups = _human_property_groups(arrays)
    for ip, pressure_value in enumerate(pressure):
        pressure_text = _format_number(float(pressure_value), table_format.pressure)
        title = f"Pressure ({pressure_unit}) = {pressure_text}"
        stream.write(title)
        stream.write("\n")
        stream.write("=" * len(title))
        stream.write("\n\n")

        structural_written = False
        if structural is not None and not any(
            group_title == "Thermoelastic properties" and group_arrays
            for group_title, group_arrays in groups
        ):
            _write_human_structural_blocks(
                stream,
                structural,
                temperature,
                ip,
                temperature_unit,
                table_format,
            )
            structural_written = True

        for group_title, group_arrays in groups:
            if not group_arrays:
                continue
            stream.write(group_title)
            stream.write("\n")
            stream.write("-" * len(group_title))
            stream.write("\n")
            _write_human_pressure_block(
                stream,
                group_arrays,
                temperature,
                ip,
                temperature_unit,
                result.metadata.get("units", {}),
                table_format,
                thermal_expansion_source,
            )
            stream.write("\n")
            if (
                structural is not None
                and not structural_written
                and group_title == "Thermoelastic properties"
            ):
                _write_human_structural_blocks(
                    stream,
                    structural,
                    temperature,
                    ip,
                    temperature_unit,
                    table_format,
                )
                structural_written = True

def _human_property_groups(
    arrays: list[tuple[str, str, str, np.ndarray, np.ndarray | None]],
) -> list[tuple[str, list[tuple[str, str, str, np.ndarray, np.ndarray | None]]]]:
    """Group properties into compact human-readable table sections."""
    if len(arrays) == 1:
        return [(arrays[0][2], arrays)]

    thermoelastic = [item for item in arrays if item[0] in THERMOELASTIC_PROPERTIES]
    energy_contributions = [
        item for item in arrays if item[0] in ENERGY_CONTRIBUTION_PROPERTIES
    ]
    potentials = [
        item for item in arrays if item[0] in THERMODYNAMIC_POTENTIAL_PROPERTIES
    ]
    thermal = [item for item in arrays if item[0] in THERMAL_PROPERTIES]
    gruneisen = [item for item in arrays if item[0] in GRUNEISEN_PROPERTIES]
    known = (
        THERMOELASTIC_PROPERTIES
        | GRUNEISEN_PROPERTIES
        | ENERGY_PROPERTIES
        | THERMAL_PROPERTIES
    )
    other = [item for item in arrays if item[0] not in known]
    return [
        ("Thermoelastic properties", thermoelastic),
        ("Gruneisen parameters", gruneisen),
        ("Energy contributions", energy_contributions),
        ("Thermodynamic potentials", potentials),
        ("Entropy and heat capacities", thermal),
        ("Additional properties", other),
    ]

def _write_human_structural_blocks(
    stream: TextIO,
    structural: _StructuralExportData,
    temperature: np.ndarray,
    pressure_index: int,
    temperature_unit: str,
    table_format: QHATableFormat,
) -> None:
    """Write structural QHA blocks at one pressure value."""
    parameters = structural.lattice_parameters[:, pressure_index]
    sigma = structural.lattice_parameter_uncertainties
    sigma_values = None if sigma is None else sigma[:, pressure_index]

    length_columns = _temperature_column(
        temperature,
        temperature_unit,
        table_format,
    )
    for index, name in enumerate(("a", "b", "c")):
        length_columns.append(
            (
                name,
                "A",
                table_format.lattice_length_width,
                table_format.lattice_length,
                parameters[:, index],
                1.0,
            )
        )
        if sigma_values is not None:
            length_columns.append(
                (
                    f"sigma_{name}",
                    "A",
                    table_format.lattice_length_width,
                    table_format.lattice_length,
                    sigma_values[:, index],
                    1.0,
                )
            )
    _write_human_named_block(stream, "Equilibrium lattice lengths", length_columns)

    angle_columns = _temperature_column(
        temperature,
        temperature_unit,
        table_format,
    )
    for offset, name in enumerate(("alpha", "beta", "gamma"), start=3):
        angle_columns.append(
            (
                name,
                "degree",
                table_format.lattice_angle_width,
                table_format.lattice_angle,
                parameters[:, offset],
                1.0,
            )
        )
        if sigma_values is not None:
            angle_columns.append(
                (
                    f"sigma_{name}",
                    "degree",
                    table_format.lattice_angle_width,
                    table_format.lattice_angle,
                    sigma_values[:, offset],
                    1.0,
                )
            )
    _write_human_named_block(stream, "Equilibrium lattice angles", angle_columns)

    axial = structural.axial_thermal_expansion
    if axial is not None:
        expansion_columns = _temperature_column(
            temperature,
            temperature_unit,
            table_format,
        )
        sigma_axial = structural.axial_thermal_expansion_uncertainties
        alpha_unit = f"{temperature_unit}^-1"
        for index, name in enumerate(("alpha_a", "alpha_b", "alpha_c")):
            expansion_columns.append(
                (
                    f"{name} x 10^5",
                    alpha_unit,
                    table_format.thermal_expansion_width,
                    table_format.thermal_expansion_scaled,
                    axial[:, pressure_index, index],
                    1.0e5,
                )
            )
            if sigma_axial is not None:
                expansion_columns.append(
                    (
                        f"sigma_{name} x 10^5",
                        alpha_unit,
                        table_format.thermal_expansion_width,
                        table_format.thermal_expansion_scaled,
                        sigma_axial[:, pressure_index, index],
                        1.0e5,
                    )
                )
        if structural.volumetric_thermal_expansion is not None:
            expansion_columns.append(
                (
                    "alphaV x 10^5",
                    alpha_unit,
                    table_format.thermal_expansion_width,
                    table_format.thermal_expansion_scaled,
                    structural.volumetric_thermal_expansion[:, pressure_index],
                    1.0e5,
                )
            )
        if structural.tensor_trace is not None:
            expansion_columns.append(
                (
                    "trace(alpha) x 10^5",
                    alpha_unit,
                    table_format.thermal_expansion_width,
                    table_format.thermal_expansion_scaled,
                    structural.tensor_trace[:, pressure_index],
                    1.0e5,
                )
            )
        if structural.trace_residual is not None:
            expansion_columns.append(
                (
                    "trace-alphaV",
                    alpha_unit,
                    table_format.thermal_expansion_width,
                    table_format.thermal_expansion,
                    structural.trace_residual[:, pressure_index],
                    1.0,
                )
            )
        if structural.extrapolation_mask is not None:
            expansion_columns.append(
                (
                    "extrapolated",
                    "-",
                    12,
                    None,
                    np.where(
                        structural.extrapolation_mask[:, pressure_index],
                        "yes",
                        "no",
                    ),
                    1.0,
                )
            )
        _write_human_named_block(
            stream,
            "Axial thermal expansion",
            expansion_columns,
        )

    tensor = structural.thermal_expansion_tensor
    if tensor is not None:
        tensor_columns = _temperature_column(
            temperature,
            temperature_unit,
            table_format,
        )
        alpha_unit = f"{temperature_unit}^-1"
        components = (
            ("alpha_11", 0, 0),
            ("alpha_22", 1, 1),
            ("alpha_33", 2, 2),
            ("alpha_23", 1, 2),
            ("alpha_13", 0, 2),
            ("alpha_12", 0, 1),
        )
        for name, i, j in components:
            tensor_columns.append(
                (
                    f"{name} x 10^5",
                    alpha_unit,
                    table_format.thermal_expansion_width,
                    table_format.thermal_expansion_scaled,
                    tensor[:, pressure_index, i, j],
                    1.0e5,
                )
            )
        _write_human_named_block(
            stream,
            "Cartesian thermal-expansion tensor",
            tensor_columns,
        )

def _temperature_column(
    temperature: np.ndarray,
    temperature_unit: str,
    table_format: QHATableFormat,
) -> list[tuple[str, str, int, str | None, np.ndarray, float]]:
    """Return the standard temperature column specification."""
    return [
        (
            "T",
            temperature_unit,
            table_format.temperature_width,
            table_format.temperature,
            temperature,
            1.0,
        )
    ]

def _write_human_named_block(
    stream: TextIO,
    title: str,
    columns: list[tuple[str, str, int, str | None, np.ndarray, float]],
) -> None:
    """Write one titled fixed-width structural table."""
    stream.write(title)
    stream.write("\n")
    stream.write("-" * len(title))
    stream.write("\n")
    _write_fixed_width_columns(stream, columns)
    stream.write("\n")

def _write_human_pressure_block(
    stream: TextIO,
    arrays: list[tuple[str, str, str, np.ndarray, np.ndarray | None]],
    temperature: np.ndarray,
    pressure_index: int,
    temperature_unit: str,
    units: dict[str, Any],
    table_format: QHATableFormat,
    thermal_expansion_source: np.ndarray | None,
) -> None:
    """Write one fixed-width property block at a selected pressure."""
    columns = _temperature_column(temperature, temperature_unit, table_format)
    for attr, symbol, _, array, sigma in arrays:
        symbol = symbol
        unit = _unit_for_property(attr, units) or "unknown"
        scale = 1.0
        if attr in {
            "thermal_expansion",
            "thermal_expansion_mixed",
            "thermal_expansion_mode",
            "thermal_expansion_numerical",
        }:
            symbol = "alphaV x 10^5"
            scale = 1.0e5
        width = _human_width_for_property(attr, table_format)
        spec = _human_format_for_property(attr, table_format)
        columns.append((symbol, unit, width, spec, array[:, pressure_index], scale))
        if attr == "thermal_expansion" and thermal_expansion_source is not None:
            columns.append(
                (
                    "alphaV method",
                    "-",
                    24,
                    None,
                    thermal_expansion_source[:, pressure_index],
                    1.0,
                )
            )
        if sigma is not None:
            columns.append(
                (
                    f"sigma_{symbol}",
                    unit,
                    width,
                    spec,
                    sigma[:, pressure_index],
                    scale,
                )
            )

    _write_fixed_width_columns(stream, columns)

def _write_fixed_width_columns(
    stream: TextIO,
    columns: list[tuple[str, str, int, str | None, np.ndarray, float]],
) -> None:
    """Write one fixed-width table from prepared column specifications."""
    if not columns:
        return
    row_count = int(np.asarray(columns[0][4]).shape[0])
    for _, _, _, _, data, _ in columns[1:]:
        if int(np.asarray(data).shape[0]) != row_count:
            raise ValueError("fixed-width table columns have inconsistent lengths")

    widths = [
        max(width, len(symbol), len(f"({unit})"))
        for symbol, unit, width, _, _, _ in columns
    ]
    stream.write(
        "  ".join(
            symbol.center(widths[index])
            for index, (symbol, _, _, _, _, _) in enumerate(columns)
        )
    )
    stream.write("\n")
    stream.write(
        "  ".join(
            f"({unit})".center(widths[index])
            for index, (_, unit, _, _, _, _) in enumerate(columns)
        )
    )
    stream.write("\n")
    stream.write("  ".join("-" * width for width in widths))
    stream.write("\n")

    for row_index in range(row_count):
        values: list[str] = []
        for column_index, (_, _, _, format_spec, data, scale) in enumerate(columns):
            text: str
            if format_spec is None:
                text = str(data[row_index])
            else:
                value = float(data[row_index]) * scale
                text = _format_number(value, str(format_spec))
            values.append(text.rjust(widths[column_index]))
        stream.write("  ".join(values))
        stream.write("\n")

def _human_format_for_property(attr: str, table_format: QHATableFormat) -> str:
    """Return a sign-aware format for fixed-width human-readable tables."""
    if attr in ENERGY_PROPERTIES or attr in THERMAL_PROPERTIES:
        return table_format.signed_scientific
    if attr in {
        "thermal_expansion",
        "thermal_expansion_mixed",
        "thermal_expansion_mode",
        "thermal_expansion_numerical",
    }:
        return table_format.thermal_expansion_scaled
    return _format_for_property(attr, table_format)

def _human_width_for_property(attr: str, table_format: QHATableFormat) -> int:
    """Return the fixed human-readable width associated with a property."""
    if attr == "equilibrium_volume":
        return table_format.volume_width
    if attr in {"isothermal_bulk_modulus", "adiabatic_bulk_modulus"}:
        return table_format.bulk_modulus_width
    if attr in {
        "bulk_modulus_derivative",
        "gruneisen",
        "mode_weighted_gruneisen",
    }:
        return table_format.bulk_derivative_width
    if attr in {
        "thermal_expansion",
        "thermal_expansion_mixed",
        "thermal_expansion_mode",
        "thermal_expansion_numerical",
    }:
        return table_format.thermal_expansion_width
    if attr in ENERGY_PROPERTIES or attr in THERMAL_PROPERTIES:
        return table_format.energy_width
    return table_format.default_width

def _write_thermal_expansion_provenance(
    stream: TextIO,
    result: QHAResult,
) -> None:
    """Write thermal-expansion method provenance as export metadata."""
    metadata = result.metadata.get("thermal_expansion", {})
    if not isinstance(metadata, dict) or not metadata:
        return
    requested = metadata.get("requested_method")
    selected = metadata.get("selected_method")
    fallback = metadata.get("fallback_method")
    if requested is not None:
        stream.write(f"# Thermal-expansion requested method: {requested}\n")
    if selected is not None:
        stream.write(f"# Thermal-expansion effective method: {selected}\n")
    if fallback is not None:
        stream.write(f"# Thermal-expansion fallback method: {fallback}\n")
    counts = metadata.get("source_counts", {})
    if isinstance(counts, dict) and counts:
        summary = ", ".join(f"{name}={int(count)}" for name, count in counts.items())
        stream.write(f"# Thermal-expansion source counts: {summary}\n")

__all__ = ["_write_human_table"]
