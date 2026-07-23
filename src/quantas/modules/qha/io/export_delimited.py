# -*- coding: utf-8 -*-

"""Delimited CSV/TSV rendering for QHA table exports."""

from __future__ import annotations

import csv
from typing import Any, TextIO

import numpy as np

from quantas.modules.qha.formatting import QHATableFormat
from quantas.modules.qha.io.export_schema import (
    _StructuralExportData,
    _format_for_property,
    _format_number,
    _unit_for_property,
)

def _write_delimited_table(
    stream: TextIO,
    arrays: list[tuple[str, str, str, np.ndarray, np.ndarray | None]],
    structural: _StructuralExportData | None,
    temperature: np.ndarray,
    pressure: np.ndarray,
    temperature_unit: str,
    pressure_unit: str,
    units: dict[str, Any],
    delimiter: str,
    table_format: QHATableFormat,
    thermal_expansion_source: np.ndarray | None,
) -> None:
    """Write a flat machine-readable pressure-temperature table.

    Parameters
    ----------
    stream : file-like object
        Open text stream.
    arrays : list of tuple
        Prepared property arrays and optional uncertainty arrays.
    structural : _StructuralExportData or None
        Optional structural arrays.
    temperature, pressure : ndarray
        Pressure-temperature axes.
    temperature_unit, pressure_unit : str
        Axis unit labels.
    units : dict
        QHA result-unit metadata.
    delimiter : str
        Field delimiter.
    table_format : QHATableFormat
        Numerical formatting rules.
    """
    writer = csv.writer(stream, delimiter=delimiter, lineterminator="\n")
    headers = [f"P ({pressure_unit})", f"T ({temperature_unit})"]
    structural_inserted = False
    include_structural_alpha_v = not any(
        attr == "thermal_expansion" for attr, _, _, _, _ in arrays
    )
    for attr, symbol, _, _, sigma in arrays:
        unit = _unit_for_property(attr, units) or "unknown"
        headers.append(f"{symbol} ({unit})")
        if attr == "thermal_expansion" and thermal_expansion_source is not None:
            headers.append("alphaV method")
        if sigma is not None:
            headers.append(f"sigma_{symbol} ({unit})")
        if attr == "equilibrium_volume" and structural is not None:
            headers.extend(
                _structural_delimited_headers(
                    structural,
                    temperature_unit=temperature_unit,
                    include_volumetric_expansion=include_structural_alpha_v,
                )
            )
            structural_inserted = True
    if structural is not None and not structural_inserted:
        headers.extend(
            _structural_delimited_headers(
                structural,
                temperature_unit=temperature_unit,
                include_volumetric_expansion=include_structural_alpha_v,
            )
        )
    writer.writerow(headers)

    for ip, pressure_value in enumerate(pressure):
        for it, temperature_value in enumerate(temperature):
            row = [
                _format_number(float(pressure_value), table_format.pressure),
                _format_number(float(temperature_value), table_format.temperature),
            ]
            structural_inserted = False
            for attr, _, _, array, sigma in arrays:
                spec = _format_for_property(attr, table_format)
                row.append(_format_number(float(array[it, ip]), spec))
                if attr == "thermal_expansion" and thermal_expansion_source is not None:
                    row.append(str(thermal_expansion_source[it, ip]))
                if sigma is not None:
                    row.append(_format_number(float(sigma[it, ip]), spec))
                if attr == "equilibrium_volume" and structural is not None:
                    row.extend(
                        _structural_delimited_row(
                            structural,
                            temperature_index=it,
                            pressure_index=ip,
                            table_format=table_format,
                            include_volumetric_expansion=(include_structural_alpha_v),
                        )
                    )
                    structural_inserted = True
            if structural is not None and not structural_inserted:
                row.extend(
                    _structural_delimited_row(
                        structural,
                        temperature_index=it,
                        pressure_index=ip,
                        table_format=table_format,
                        include_volumetric_expansion=include_structural_alpha_v,
                    )
                )
            writer.writerow(row)

def _structural_delimited_headers(
    structural: _StructuralExportData,
    *,
    temperature_unit: str,
    include_volumetric_expansion: bool,
) -> list[str]:
    """Return machine-readable headers for structural QHA properties."""
    headers: list[str] = []
    parameter_names = ("a", "b", "c", "alpha", "beta", "gamma")
    parameter_units = ("A", "A", "A", "degree", "degree", "degree")
    for index, (name, unit) in enumerate(zip(parameter_names, parameter_units)):
        headers.append(f"{name} ({unit})")
        if structural.lattice_parameter_uncertainties is not None:
            headers.append(f"sigma_{name} ({unit})")

    alpha_unit = f"{temperature_unit}^-1"
    if structural.axial_thermal_expansion is not None:
        for name in ("alpha_a", "alpha_b", "alpha_c"):
            headers.append(f"{name} ({alpha_unit})")
            if structural.axial_thermal_expansion_uncertainties is not None:
                headers.append(f"sigma_{name} ({alpha_unit})")
    if (
        include_volumetric_expansion
        and structural.volumetric_thermal_expansion is not None
    ):
        headers.append(f"alphaV ({alpha_unit})")
    if structural.tensor_trace is not None:
        headers.append(f"trace(alpha) ({alpha_unit})")
    if structural.trace_residual is not None:
        headers.append(f"trace(alpha)-alphaV ({alpha_unit})")
    if structural.thermal_expansion_tensor is not None:
        for name in (
            "alpha_11",
            "alpha_22",
            "alpha_33",
            "alpha_23",
            "alpha_13",
            "alpha_12",
        ):
            headers.append(f"{name} ({alpha_unit})")
    if structural.extrapolation_mask is not None:
        headers.append("structural_extrapolated")
    return headers

def _structural_delimited_row(
    structural: _StructuralExportData,
    *,
    temperature_index: int,
    pressure_index: int,
    table_format: QHATableFormat,
    include_volumetric_expansion: bool,
) -> list[str]:
    """Return one machine-readable structural QHA row."""
    row: list[str] = []
    values = structural.lattice_parameters[temperature_index, pressure_index]
    sigma = structural.lattice_parameter_uncertainties
    for index, value in enumerate(values):
        spec = table_format.lattice_length if index < 3 else table_format.lattice_angle
        row.append(_format_number(float(value), spec))
        if sigma is not None:
            row.append(
                _format_number(
                    float(sigma[temperature_index, pressure_index, index]),
                    spec,
                )
            )

    axial = structural.axial_thermal_expansion
    sigma_axial = structural.axial_thermal_expansion_uncertainties
    if axial is not None:
        for index in range(3):
            row.append(
                _format_number(
                    float(axial[temperature_index, pressure_index, index]),
                    table_format.thermal_expansion,
                )
            )
            if sigma_axial is not None:
                row.append(
                    _format_number(
                        float(
                            sigma_axial[
                                temperature_index,
                                pressure_index,
                                index,
                            ]
                        ),
                        table_format.thermal_expansion,
                    )
                )
    if (
        include_volumetric_expansion
        and structural.volumetric_thermal_expansion is not None
    ):
        row.append(
            _format_number(
                float(
                    structural.volumetric_thermal_expansion[
                        temperature_index,
                        pressure_index,
                    ]
                ),
                table_format.thermal_expansion,
            )
        )
    if structural.tensor_trace is not None:
        row.append(
            _format_number(
                float(structural.tensor_trace[temperature_index, pressure_index]),
                table_format.thermal_expansion,
            )
        )
    if structural.trace_residual is not None:
        row.append(
            _format_number(
                float(structural.trace_residual[temperature_index, pressure_index]),
                table_format.thermal_expansion,
            )
        )
    if structural.thermal_expansion_tensor is not None:
        tensor = structural.thermal_expansion_tensor[
            temperature_index,
            pressure_index,
        ]
        for i, j in ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)):
            row.append(
                _format_number(
                    float(tensor[i, j]),
                    table_format.thermal_expansion,
                )
            )
    if structural.extrapolation_mask is not None:
        row.append(
            "1"
            if bool(
                structural.extrapolation_mask[
                    temperature_index,
                    pressure_index,
                ]
            )
            else "0"
        )
    return row

__all__ = ["_write_delimited_table"]
