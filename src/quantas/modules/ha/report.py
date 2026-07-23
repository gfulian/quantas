# -*- coding: utf-8 -*-

"""
Neutral reporting utilities for harmonic-approximation results.

The functions defined here convert HA input, options, and result objects into
plain table containers. These tables are independent from command-line and GUI
frontends: the CLI may render them as text, while a Dash interface may convert
exactly the same data to graphical tables.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import numpy as np

from quantas.models.report import ReportTable as _ReportTable

if TYPE_CHECKING:
    from quantas.models.report import ReportTable
from quantas.modules.ha.models import HAInput, HAOptions, HAResult


THERMODYNAMIC_LABELS = {
    "zero_point_energy": ("Uzp", "Zero-point energy"),
    "thermal_energy": ("Uth", "Thermal energy"),
    "internal_energy": ("Utot", "Internal energy"),
    "entropy": ("S", "Entropy"),
    "vibrational_free_energy": ("Fvib", "Vibrational Helmholtz free energy"),
    "free_energy": ("F", "Helmholtz free energy"),
    "isochoric_heat_capacity": ("Cv", "Isochoric heat capacity"),
}


def input_table(input_data: HAInput) -> ReportTable:
    """Create a table describing the HA input data.

    Parameters
    ----------
    input_data : HAInput
        Normalized harmonic-approximation input data.

    Returns
    -------
    ReportTable
        Table containing the main input dimensions and source information.
    """
    return _ReportTable(
        title="Input data",
        columns=["Property", "Value"],
        rows=[
            ["Job name", input_data.jobname],
            ["Number of atoms", int(input_data.natoms)],
            ["Formula units per cell (Z)", int(input_data.formula_units)],
            ["Atoms per formula unit", input_data.natoms_per_formula_unit],
            ["Volumes", int(input_data.nvol)],
            ["Q-points", int(input_data.qpoints)],
            ["Modes", int(input_data.nmodes)],
            ["Total q-point weight", float(input_data.total_q_points)],
            [
                "Source",
                str(input_data.source) if input_data.source is not None else "Unknown",
            ],
        ],
        metadata={"column_formats": [None, "general"]},
    )


def options_table(options: HAOptions) -> ReportTable:
    """Create a table describing selected HA workflow options.

    Parameters
    ----------
    options : HAOptions
        Harmonic-approximation workflow options.

    Returns
    -------
    ReportTable
        Table containing temperature range, units, and workflow switches.
    """
    return _ReportTable(
        title="Selected options",
        columns=["Option", "Value"],
        rows=[
            ["Temperature minimum", float(options.temperature_min)],
            ["Temperature maximum", float(options.temperature_max)],
            ["Temperature step", float(options.temperature_step)],
            ["Temperature unit", options.temperature_unit],
            ["Energy unit", options.energy_unit],
            ["Volume unit", options.volume_unit],
            ["Frequency unit", options.frequency_unit],
            ["Calculate thermodynamics", bool(options.calculate_thermodynamics)],
            ["Numerical precision", "float64"],
        ],
        metadata={"column_formats": [None, "general"]},
    )


def static_data_table(result: HAResult) -> ReportTable:
    """Create the volume-resolved static and zero-point energy table.

    Parameters
    ----------
    result : HAResult
        Harmonic result containing the sampled volume grid and static energy.

    Returns
    -------
    ReportTable
        One row per sampled volume.  Zero-point energy is included when
        available and is represented once per volume because it is independent
        of temperature.

    Raises
    ------
    ValueError
        If volume, static-energy, or zero-point arrays are incompatible.
    """
    if result.volume is None:
        raise ValueError("HA result volumes are not available")
    if result.static_energy is None:
        raise ValueError("HA result static energies are not available")

    volume = np.asarray(result.volume, dtype=np.float64)
    static_energy = np.asarray(result.static_energy, dtype=np.float64)
    if volume.ndim != 1 or static_energy.ndim != 1:
        raise ValueError("HA volume and static energy must be one-dimensional")
    if volume.size != static_energy.size:
        raise ValueError("HA static-energy volume dimension is incompatible")

    zero_point = None
    if result.zero_point_energy is not None:
        zero_point = _temperature_independent_volume_array(
            result.zero_point_energy,
            volume.size,
            property_name="zero_point_energy",
        )

    columns = ["V", "U0"]
    if zero_point is not None:
        columns.append("Uzp")

    rows: list[list[Any]] = []
    for index, value in enumerate(volume):
        row: list[Any] = [float(value), float(static_energy[index])]
        if zero_point is not None:
            row.append(float(zero_point[index]))
        rows.append(row)

    units = _result_units(result)
    energy_unit = str(units.get("energy", "unknown"))
    volume_unit = str(units.get("volume", "unknown"))
    metadata = {
        "column_formats": ["volume", *(["energy"] * (len(columns) - 1))],
        "column_units": [volume_unit, *([energy_unit] * (len(columns) - 1))],
        "column_alignments": ["right"] * len(columns),
        "independent_variable": "volume",
        "notes": [
            "Zero-point energy is temperature independent and is listed once "
            "for each sampled volume."
        ]
        if zero_point is not None
        else [],
    }
    return _ReportTable(
        title="Static and zero-point energies",
        columns=columns,
        rows=rows,
        metadata=metadata,
    )


def thermodynamic_summary_table(result: HAResult) -> ReportTable:
    """Create a compact summary of available HA result arrays.

    Parameters
    ----------
    result : HAResult
        Harmonic-approximation result object.

    Returns
    -------
    ReportTable
        Table listing available properties and array shapes.
    """
    rows: list[list[Any]] = []
    for attr, (symbol, description) in THERMODYNAMIC_LABELS.items():
        value = getattr(result, attr)
        rows.append(
            [
                symbol,
                description,
                "yes" if value is not None else "no",
                tuple(np.asarray(value).shape) if value is not None else "-",
            ]
        )

    return _ReportTable(
        title="Thermodynamic properties",
        columns=["Key", "Description", "Available", "Shape"],
        rows=rows,
    )


def thermodynamic_table(
    result: HAResult,
    property_name: str,
    max_rows: int | None = None,
    row_indices: Sequence[int] | None = None,
) -> ReportTable:
    """Create a neutral table for one thermodynamic property.

    Temperature-dependent properties are represented as a matrix with one row
    per selected temperature and one column per sampled volume.  Temperature-
    independent properties such as zero-point energy are represented as one row
    per volume and are never duplicated across the temperature grid.

    Parameters
    ----------
    result : HAResult
        Harmonic-approximation result object.
    property_name : str
        Name of the result attribute or historical key. Accepted values include
        ``"free_energy"`` and ``"F"``.
    max_rows : int or None, optional
        Maximum number of rows to include.  Rows refer to temperatures for
        temperature-dependent properties and volumes for temperature-independent
        properties.
    row_indices : sequence of int or None, optional
        Explicit row indices. Negative indices are accepted. Selection is
        applied before ``max_rows``.

    Returns
    -------
    ReportTable
        Frontend-neutral property table.

    Raises
    ------
    KeyError
        If the requested property is unknown.
    ValueError
        If the requested property is not available or has an unsupported shape.
    """
    attr = _resolve_property_name(property_name)
    symbol, description = THERMODYNAMIC_LABELS[attr]
    values = getattr(result, attr)

    if values is None:
        raise ValueError(f"HA result property '{property_name}' is not available")
    if result.volume is None:
        raise ValueError("HA result volumes are not available")

    volume = np.asarray(result.volume, dtype=np.float64)
    if volume.ndim != 1:
        raise ValueError("HA result volumes must be one-dimensional")
    nvol = int(volume.size)

    if attr == "zero_point_energy":
        data = _temperature_independent_volume_array(
            values,
            nvol,
            property_name=attr,
        )
        indices = _normalized_row_indices(nvol, row_indices)
        if max_rows is not None:
            indices = indices[:max_rows]
        rows = [[float(volume[index]), float(data[index])] for index in indices]
        units = _result_units(result)
        return _ReportTable(
            title=description,
            columns=["V", symbol],
            rows=rows,
            metadata={
                "property": attr,
                "symbol": symbol,
                "description": description,
                "independent_variable": "volume",
                "temperature_independent": True,
                "truncated": len(indices) < nvol,
                "total_rows": nvol,
                "column_formats": ["volume", "energy"],
                "column_units": [
                    str(units.get("volume", "unknown")),
                    str(units.get("energy", "unknown")),
                ],
                "column_alignments": ["right", "right"],
            },
        )

    if result.temperature is None:
        raise ValueError("HA result temperatures are not available")
    temperature = np.asarray(result.temperature, dtype=np.float64)
    if temperature.ndim != 1:
        raise ValueError("HA result temperatures must be one-dimensional")

    data = _temperature_volume_array(
        values,
        ntemp=int(temperature.size),
        nvol=nvol,
        property_name=attr,
    )
    indices = _normalized_row_indices(data.shape[0], row_indices)
    if max_rows is not None:
        indices = indices[:max_rows]

    rows = [
        [float(temperature[index]), *[float(value) for value in data[index]]]
        for index in indices
    ]

    value_format = "energy"
    units = _result_units(result)
    property_unit = _property_unit(attr, units)
    volume_unit = str(units.get("volume", "unknown"))
    temperature_unit = str(units.get("temperature", "K"))
    metadata = {
        "property": attr,
        "symbol": symbol,
        "description": description,
        "independent_variable": "temperature",
        "temperature_independent": False,
        "truncated": len(indices) < data.shape[0],
        "total_rows": int(data.shape[0]),
        "column_formats": ["temperature", *([value_format] * nvol)],
        "column_units": [temperature_unit, *([property_unit] * nvol)],
        "column_alignments": ["right"] * (nvol + 1),
        "notes": [f"Volume columns are labelled in {volume_unit}."],
    }

    return _ReportTable(
        title=description,
        columns=["T", *_volume_labels(volume, nvol)],
        rows=rows,
        metadata=metadata,
    )


def thermodynamic_property_tables(
    result: HAResult,
    *,
    row_indices: Sequence[int] | None = None,
    include_zero_point: bool = False,
) -> list[ReportTable]:
    """Build ordered HA property tables.

    Parameters
    ----------
    result : HAResult
        Harmonic result object.
    row_indices : sequence of int or None, optional
        Temperature rows selected for temperature-dependent properties.
    include_zero_point : bool, optional
        Include a separate zero-point table.  This is normally unnecessary when
        :func:`static_data_table` is already part of the report.

    Returns
    -------
    list of ReportTable
        One table for every available property.
    """
    tables: list[ReportTable] = []
    for property_name in THERMODYNAMIC_LABELS:
        if property_name == "zero_point_energy" and not include_zero_point:
            continue
        if getattr(result, property_name) is None:
            continue
        selected_rows = None if property_name == "zero_point_energy" else row_indices
        tables.append(
            thermodynamic_table(
                result,
                property_name,
                row_indices=selected_rows,
            )
        )
    return tables


def _normalized_row_indices(
    size: int,
    row_indices: Sequence[int] | None,
) -> list[int]:
    """Return unique, validated row indices in requested order."""
    if row_indices is None:
        return list(range(size))
    normalized: list[int] = []
    for raw_index in row_indices:
        index = raw_index if raw_index >= 0 else size + raw_index
        if index < 0 or index >= size:
            raise IndexError(f"HA report row index out of range: {raw_index}")
        if index not in normalized:
            normalized.append(index)
    return normalized


def all_tables(
    input_data: HAInput | None,
    options: HAOptions | None,
    result: HAResult,
    *,
    row_indices: Sequence[int] | None = (0, -1),
) -> list[ReportTable]:
    """Build the standard set of HA report tables.

    Parameters
    ----------
    input_data : HAInput or None
        Input data used by the workflow. If ``None``, the input table is not
        included.
    options : HAOptions or None
        Workflow options. If ``None``, the options table is not included.
    result : HAResult
        Harmonic-approximation result object.
    row_indices : sequence of int or None, optional
        Temperature rows included in the standard property excerpts.  ``None``
        includes the complete temperature grid.

    Returns
    -------
    list of ReportTable
        Neutral HA report tables including volume-resolved static data and
        thermodynamic property excerpts.
    """
    tables: list[ReportTable] = []
    if input_data is not None:
        tables.append(input_table(input_data))
    if options is not None:
        tables.append(options_table(options))
    tables.append(static_data_table(result))
    tables.append(thermodynamic_summary_table(result))
    tables.extend(
        thermodynamic_property_tables(
            result,
            row_indices=row_indices,
            include_zero_point=False,
        )
    )
    return tables


def _resolve_property_name(name: str) -> str:
    """Resolve a result attribute or historical key to a HAResult attribute."""
    if name in THERMODYNAMIC_LABELS:
        return name
    for attr, (symbol, _) in THERMODYNAMIC_LABELS.items():
        if name == symbol:
            return attr
    raise KeyError(f"unknown HA thermodynamic property: {name}")


def _temperature_independent_volume_array(
    values: Any,
    nvol: int,
    *,
    property_name: str,
) -> np.ndarray:
    """Return one temperature-independent value for every volume."""
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 1 and array.size == nvol:
        return array
    if array.ndim == 2 and array.shape == (1, nvol):
        return array[0]
    raise ValueError(
        f"HA result property '{property_name}' must have shape (nvol,) or "
        "(1, nvol)"
    )


def _temperature_volume_array(
    values: Any,
    *,
    ntemp: int,
    nvol: int,
    property_name: str,
) -> np.ndarray:
    """Return one temperature-dependent HA array with shape ``(nt, nvol)``."""
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 1 and nvol == 1 and array.size == ntemp:
        return array[:, np.newaxis]
    if array.ndim == 2 and array.shape == (ntemp, nvol):
        return array
    raise ValueError(
        f"HA result property '{property_name}' must have shape "
        f"({ntemp}, {nvol})"
    )


def _result_units(result: HAResult) -> dict[str, Any]:
    """Return normalized unit metadata from a HA result."""
    units = result.metadata.get("units", {})
    return units if isinstance(units, dict) else {}


def _property_unit(attr: str, units: dict[str, Any]) -> str:
    """Return the stored unit label for one HA thermodynamic property."""
    if attr == "entropy":
        return str(units.get("entropy", "unknown"))
    if attr == "isochoric_heat_capacity":
        return str(units.get("heat_capacity", "unknown"))
    return str(units.get("energy", "unknown"))


def _volume_labels(volume: np.ndarray | None, nvol: int) -> list[str]:
    """Return column labels for volume-dependent result arrays."""
    if volume is None or volume.shape[0] != nvol:
        return [f"V{index}" for index in range(nvol)]
    return [f"V={float(value):.8g}" for value in volume]
