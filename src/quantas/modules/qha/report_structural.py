# -*- coding: utf-8 -*-

"""Structural tables for QHA pressure-temperature results."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import numpy as np

from quantas.models.report import ReportTable as _ReportTable
from quantas.modules.qha.formatting import (
    QHATableFormat,
    canonical_unit_label,
    format_number,
    volume_unit_label,
)
from quantas.modules.qha.models import QHAResult
from quantas.modules.qha.report_common import _normalized_indices

if TYPE_CHECKING:
    from quantas.models.report import ReportTable

def structural_property_tables(
    result: QHAResult,
    *,
    table_format: QHATableFormat | None = None,
    pressure_indices: Sequence[int] | None = None,
) -> list[ReportTable]:
    """Create structural QHA tables grouped by pressure.

    Cell dimensions are reported with the equilibrium volume because they
    describe the same structural state.  Axial expansion and propagated
    one-standard-deviation uncertainties are kept in separate tables to avoid
    obscuring the principal structural values.

    Parameters
    ----------
    result : QHAResult
        QHA result containing structural-path properties.
    table_format : QHATableFormat or None, optional
        Numerical formatting rules.
    pressure_indices : sequence of int or None, optional
        Explicit pressure-grid indices. Negative indices are accepted. ``None``
        includes every pressure.

    Returns
    -------
    list of ReportTable
        Lattice, axial-expansion, and optional uncertainty tables for each
        pressure.
    """
    if result.temperature is None or result.pressure is None:
        return []
    if result.lattice_parameters is None:
        return []
    temperature = np.asarray(result.temperature, dtype=np.float64)
    pressure = np.asarray(result.pressure, dtype=np.float64)
    shape = (temperature.size, pressure.size)
    parameters = np.asarray(result.lattice_parameters, dtype=np.float64)
    if parameters.shape != shape + (6,):
        return []

    fmt = table_format or QHATableFormat()
    metadata = result.metadata if isinstance(result.metadata, dict) else {}
    units = (
        metadata.get("units", {}) if isinstance(metadata.get("units", {}), dict) else {}
    )
    temperature_unit = canonical_unit_label(str(units.get("temperature", "K")))
    pressure_unit = canonical_unit_label(str(units.get("pressure", "unknown")))
    volume_unit = volume_unit_label(str(units.get("volume", "A")))
    equilibrium_volume = (
        None
        if result.equilibrium_volume is None
        else np.asarray(result.equilibrium_volume, dtype=np.float64)
    )
    if equilibrium_volume is not None and equilibrium_volume.shape != shape:
        equilibrium_volume = None
    alpha_abc = (
        None
        if result.axial_thermal_expansion is None
        else np.asarray(result.axial_thermal_expansion, dtype=np.float64)
    )
    alpha_v = (
        None
        if result.thermal_expansion is None
        else np.asarray(result.thermal_expansion, dtype=np.float64)
    )
    alpha_tensor = (
        None
        if result.thermal_expansion_tensor is None
        else np.asarray(result.thermal_expansion_tensor, dtype=np.float64)
    )
    sigma_parameters = result.uncertainties.get("sigma_cell_parameters")
    sigma_parameters_array = (
        None
        if sigma_parameters is None
        else np.asarray(sigma_parameters, dtype=np.float64)
    )
    sigma_volume = result.uncertainties.get("sigma_VT")
    sigma_volume_array = (
        None if sigma_volume is None else np.asarray(sigma_volume, dtype=np.float64)
    )
    if alpha_abc is not None and alpha_abc.shape != shape + (3,):
        alpha_abc = None
    if alpha_v is not None and alpha_v.shape != shape:
        alpha_v = None
    if alpha_tensor is not None and alpha_tensor.shape != shape + (3, 3):
        alpha_tensor = None
    if sigma_parameters_array is not None and sigma_parameters_array.shape != shape + (
        6,
    ):
        sigma_parameters_array = None
    if sigma_volume_array is not None and sigma_volume_array.shape != shape:
        sigma_volume_array = None

    tables: list[ReportTable] = []
    for ip in _normalized_indices(pressure.size, pressure_indices, "pressure"):
        pressure_value = pressure[ip]
        pressure_text = format_number(float(pressure_value), fmt.pressure)
        pressure_suffix = f" {pressure_unit}" if pressure_unit != "unknown" else ""
        lattice_rows: list[list[Any]] = []
        for it, temperature_value in enumerate(temperature):
            values = parameters[it, ip]
            row: list[Any] = [
                format_number(float(temperature_value), fmt.temperature),
            ]
            if equilibrium_volume is not None:
                row.append(format_number(float(equilibrium_volume[it, ip]), fmt.volume))
            row.extend(
                [format_number(float(value), fmt.volume) for value in values[:3]]
            )
            row.extend([format_number(float(value), ".6f") for value in values[3:]])
            lattice_rows.append(row)
        lattice_columns = ["T"]
        lattice_units = [temperature_unit]
        if equilibrium_volume is not None:
            lattice_columns.append("V")
            lattice_units.append(volume_unit)
        lattice_columns.extend(["a", "b", "c", "alpha", "beta", "gamma"])
        lattice_units.extend(["A", "A", "A", "degree", "degree", "degree"])
        tables.append(
            _ReportTable(
                title=(f"QHA equilibrium cell at P = {pressure_text}{pressure_suffix}"),
                columns=lattice_columns,
                rows=lattice_rows,
                metadata={
                    "formatted": True,
                    "column_units": lattice_units,
                    "structural_method": metadata.get(
                        "structural_thermal_expansion",
                        {},
                    ).get("method"),
                },
            )
        )

        if alpha_abc is not None:
            expansion_rows: list[list[Any]] = []
            for it, temperature_value in enumerate(temperature):
                values = alpha_abc[it, ip]
                row = [
                    format_number(float(temperature_value), fmt.temperature),
                    *[
                        format_number(
                            float(value) * 1.0e5,
                            fmt.thermal_expansion_scaled,
                        )
                        for value in values
                    ],
                ]
                if alpha_v is not None:
                    row.extend(
                        [
                            format_number(
                                float(alpha_v[it, ip]) * 1.0e5,
                                fmt.thermal_expansion_scaled,
                            ),
                            format_number(
                                (
                                    float(np.trace(alpha_tensor[it, ip]))
                                    - float(alpha_v[it, ip])
                                    if alpha_tensor is not None
                                    else np.nan
                                )
                                * 1.0e5,
                                fmt.thermal_expansion_scaled,
                            ),
                        ]
                    )
                expansion_rows.append(row)
            columns = ["T", "alpha_a x 10^5", "alpha_b x 10^5", "alpha_c x 10^5"]
            column_units = [temperature_unit, "K^-1", "K^-1", "K^-1"]
            if alpha_v is not None:
                columns.extend(["alphaV x 10^5", "tr(alpha)-alphaV x 10^5"])
                column_units.extend(["K^-1", "K^-1"])
            tables.append(
                _ReportTable(
                    title=(
                        "QHA axial thermal expansion at P = "
                        f"{pressure_text}{pressure_suffix}"
                    ),
                    columns=columns,
                    rows=expansion_rows,
                    metadata={
                        "formatted": True,
                        "column_units": column_units,
                        "scale": 1.0e5,
                        "approximation": metadata.get(
                            "structural_thermal_expansion",
                            {},
                        ).get("approximation"),
                    },
                )
            )

        if sigma_parameters_array is not None:
            uncertainty_rows: list[list[Any]] = []
            for it, temperature_value in enumerate(temperature):
                values = sigma_parameters_array[it, ip]
                row = [format_number(float(temperature_value), fmt.temperature)]
                if sigma_volume_array is not None:
                    row.append(
                        format_number(float(sigma_volume_array[it, ip]), fmt.default)
                    )
                row.extend(
                    [format_number(float(value), fmt.default) for value in values]
                )
                uncertainty_rows.append(row)
            uncertainty_columns = ["T"]
            uncertainty_units = [temperature_unit]
            if sigma_volume_array is not None:
                uncertainty_columns.append("sigma_V")
                uncertainty_units.append(volume_unit)
            uncertainty_columns.extend(
                [
                    "sigma_a",
                    "sigma_b",
                    "sigma_c",
                    "sigma_alpha",
                    "sigma_beta",
                    "sigma_gamma",
                ]
            )
            uncertainty_units.extend(["A", "A", "A", "degree", "degree", "degree"])
            tables.append(
                _ReportTable(
                    title=(
                        "QHA equilibrium-cell uncertainties at P = "
                        f"{pressure_text}{pressure_suffix}"
                    ),
                    columns=uncertainty_columns,
                    rows=uncertainty_rows,
                    metadata={
                        "formatted": True,
                        "column_units": uncertainty_units,
                        "uncertainty": "one_standard_deviation",
                        "method": metadata.get(
                            "structural_thermal_expansion",
                            {},
                        ).get("uncertainty_method"),
                    },
                )
            )
    return tables

__all__ = ["structural_property_tables"]
