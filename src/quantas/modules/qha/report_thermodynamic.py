# -*- coding: utf-8 -*-

"""Thermodynamic and fit-diagnostic tables for QHA reports."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from quantas.core.math.fitting import residual_metrics
from quantas.core.physics.units import convert_energy_per_temperature
from quantas.models.report import ReportTable as _ReportTable
from quantas.models.thermodynamics import HarmonicThermodynamicResult
from quantas.modules.qha.formatting import (
    QHATableFormat,
    canonical_unit_label,
    format_number,
    property_unit_label,
)
from quantas.modules.qha.models import QHAInput, QHAOptions, QHAResult
from quantas.modules.qha.report_common import (
    QHA_PROPERTY_LABELS,
    _normalized_indices,
    _property_rows,
    _uncertainty_for,
    resolve_property_name,
)

if TYPE_CHECKING:
    from quantas.models.report import ReportTable

def property_table(
    result: QHAResult,
    property_name: str,
    *,
    include_uncertainty: bool = True,
    max_rows: int | None = None,
) -> ReportTable:
    """Create a neutral table for a single QHA property.

    Parameters
    ----------
    result : QHAResult
        Quasi-harmonic approximation result object.
    property_name : str
        Name of the result attribute or historical key. Accepted values include
        ``"equilibrium_volume"`` and ``"VT"``.
    include_uncertainty : bool, optional
        If ``True``, include a matching uncertainty column when available.
    max_rows : int or None, optional
        Maximum number of rows to include. If ``None``, all rows are included.

    Returns
    -------
    ReportTable
        Table containing pressure, temperature, property values, and optionally
        uncertainties.

    Raises
    ------
    KeyError
        If the requested property is unknown.
    ValueError
        If the requested property is not available or has an unsupported shape.
    """
    attr = resolve_property_name(property_name)
    symbol, description = QHA_PROPERTY_LABELS[attr]
    value = getattr(result, attr)
    if value is None:
        raise ValueError(f"QHA result property '{property_name}' is not available")

    data = np.asarray(value, dtype=np.float64)
    sigma = _uncertainty_for(result, attr, symbol) if include_uncertainty else None
    if sigma is not None:
        sigma = np.asarray(sigma, dtype=np.float64)
        if sigma.shape != data.shape:
            sigma = None

    rows = _property_rows(result, data, sigma, max_rows=max_rows)
    columns = ["T", "P", symbol]
    if sigma is not None:
        columns.append(f"sigma_{symbol}")

    metadata = {
        "property": attr,
        "symbol": symbol,
        "description": description,
        "shape": tuple(data.shape),
        "has_uncertainty": sigma is not None,
        "truncated": False,
        "total_rows": int(data.size),
    }
    if max_rows is not None and len(rows) < data.size:
        metadata["truncated"] = True

    return _ReportTable(
        title=description,
        columns=columns,
        rows=rows,
        metadata=metadata,
    )

def selected_property_tables(
    result: QHAResult,
    *,
    table_format: QHATableFormat | None = None,
    pressure_indices: Sequence[int] | None = None,
) -> list[ReportTable]:
    """Create standard QHA result tables grouped by pressure.

    Available one-sigma uncertainty arrays are rendered immediately after the
    corresponding property. Compact thermodynamic reports use molar
    ``Cp-Cv`` values, while the native result arrays remain unchanged.

    Parameters
    ----------
    result : QHAResult
        Quasi-harmonic approximation result object.
    table_format : QHATableFormat or None, optional
        Numerical formatting rules. Default Quantas conventions are used when
        omitted.
    pressure_indices : sequence of int or None, optional
        Explicit pressure-grid indices. Negative indices are accepted. ``None``
        includes every pressure.

    Returns
    -------
    list of ReportTable
        One table for each pressure value.
    """
    if result.temperature is None or result.pressure is None:
        return []
    if result.equilibrium_volume is None or result.isothermal_bulk_modulus is None:
        return []

    fmt = table_format or QHATableFormat()
    temperature = np.asarray(result.temperature, dtype=np.float64)
    pressure = np.asarray(result.pressure, dtype=np.float64)
    shape = (temperature.size, pressure.size)
    metadata = result.metadata if isinstance(result.metadata, dict) else {}
    units = (
        metadata.get("units", {}) if isinstance(metadata.get("units", {}), dict) else {}
    )
    temperature_unit = canonical_unit_label(str(units.get("temperature", "K")))
    pressure_unit = canonical_unit_label(str(units.get("pressure", "unknown")))

    volume = _optional_standard_array(result.equilibrium_volume, shape)
    kt = _optional_standard_array(result.isothermal_bulk_modulus, shape)
    if volume is None or kt is None:
        return []

    kp = _optional_standard_array(result.bulk_modulus_derivative, shape)
    ks = _optional_standard_array(result.adiabatic_bulk_modulus, shape)
    alpha = _optional_standard_array(result.thermal_expansion, shape)
    gamma = _optional_standard_array(result.gruneisen, shape)
    gamma_mode = _optional_standard_array(result.mode_weighted_gruneisen, shape)
    cp_cv = None
    if result.heat_capacity_difference is not None:
        candidate = _to_jmol(
            np.asarray(result.heat_capacity_difference, dtype=np.float64),
            result,
        )
        if candidate.shape == shape:
            cp_cv = candidate

    sigma_volume = _optional_standard_array(
        _uncertainty_for(result, "equilibrium_volume", "VT"), shape
    )
    sigma_kt = _optional_standard_array(
        _uncertainty_for(result, "isothermal_bulk_modulus", "KT"), shape
    )
    sigma_kp = _optional_standard_array(
        _uncertainty_for(result, "bulk_modulus_derivative", "Kp"), shape
    )
    sigma_ks = _optional_standard_array(
        _uncertainty_for(result, "adiabatic_bulk_modulus", "KS"), shape
    )
    sigma_alpha = _optional_standard_array(
        _uncertainty_for(result, "thermal_expansion", "alphaV"), shape
    )
    sigma_cp_cv = _uncertainty_for(result, "heat_capacity_difference", "Cp-Cv")
    if sigma_cp_cv is not None:
        candidate = _to_jmol(np.asarray(sigma_cp_cv, dtype=np.float64), result)
        sigma_cp_cv = candidate if candidate.shape == shape else None

    # label, unit, values, scale, format specification
    specifications: list[tuple[str, str, np.ndarray | None, float, str]] = [
        (
            "V",
            property_unit_label("equilibrium_volume", units),
            volume,
            1.0,
            fmt.volume,
        ),
    ]
    if sigma_volume is not None:
        specifications.append(
            (
                "sigma_V",
                property_unit_label("equilibrium_volume", units),
                sigma_volume,
                1.0,
                fmt.volume,
            )
        )
    specifications.append(
        (
            "KT",
            property_unit_label("isothermal_bulk_modulus", units),
            kt,
            1.0,
            fmt.bulk_modulus,
        )
    )
    if sigma_kt is not None:
        specifications.append(
            (
                "sigma_KT",
                property_unit_label("isothermal_bulk_modulus", units),
                sigma_kt,
                1.0,
                fmt.bulk_modulus,
            )
        )
    if kp is not None:
        specifications.append(("Kp", "-", kp, 1.0, fmt.bulk_derivative))
        if sigma_kp is not None:
            specifications.append(("sigma_Kp", "-", sigma_kp, 1.0, fmt.bulk_derivative))
    if ks is not None:
        specifications.append(
            (
                "KS",
                property_unit_label("adiabatic_bulk_modulus", units),
                ks,
                1.0,
                fmt.bulk_modulus,
            )
        )
        if sigma_ks is not None:
            specifications.append(
                (
                    "sigma_KS",
                    property_unit_label("adiabatic_bulk_modulus", units),
                    sigma_ks,
                    1.0,
                    fmt.bulk_modulus,
                )
            )
    if alpha is not None:
        specifications.append(
            (
                "alphaV x 10^5",
                property_unit_label("thermal_expansion", units),
                alpha,
                1.0e5,
                fmt.thermal_expansion_scaled,
            )
        )
        if sigma_alpha is not None:
            specifications.append(
                (
                    "sigma_alphaV x 10^5",
                    property_unit_label("thermal_expansion", units),
                    sigma_alpha,
                    1.0e5,
                    fmt.thermal_expansion_scaled,
                )
            )
    if gamma is not None:
        specifications.append(("gamma", "-", gamma, 1.0, fmt.bulk_derivative))
    if gamma_mode is not None:
        specifications.append(("gamma_mode", "-", gamma_mode, 1.0, fmt.bulk_derivative))
    if cp_cv is not None:
        specifications.append(("Cp-Cv", "J mol^-1 K^-1", cp_cv, 1.0, fmt.energy))
        if sigma_cp_cv is not None:
            specifications.append(
                (
                    "sigma_Cp-Cv",
                    "J mol^-1 K^-1",
                    sigma_cp_cv,
                    1.0,
                    fmt.energy,
                )
            )

    tables: list[ReportTable] = []
    for ip in _normalized_indices(pressure.size, pressure_indices, "pressure"):
        pvalue = pressure[ip]
        rows: list[list[Any]] = []
        for it, tvalue in enumerate(temperature):
            row: list[Any] = [format_number(float(tvalue), fmt.temperature)]
            for _, _, array, scale, spec in specifications:
                row.append(
                    _format_optional_standard_value(
                        array,
                        it,
                        ip,
                        scale=scale,
                        spec=spec,
                    )
                )
            rows.append(row)
        tables.append(
            _ReportTable(
                title=(
                    "QHA results at P = "
                    f"{format_number(float(pvalue), fmt.pressure)}"
                    + (f" {pressure_unit}" if pressure_unit != "unknown" else "")
                ),
                columns=["T"] + [label for label, _, _, _, _ in specifications],
                rows=rows,
                metadata={
                    "pressure": float(pvalue),
                    "formatted": True,
                    "column_units": [temperature_unit]
                    + [unit for _, unit, _, _, _ in specifications],
                    "column_alignments": ["right"] * (len(specifications) + 1),
                    "uncertainty_columns": [
                        label
                        for label, _, _, _, _ in specifications
                        if label.startswith("sigma_")
                    ],
                    "notes": ["Cp-Cv is normalized per mole of formula units."]
                    if cp_cv is not None
                    else [],
                },
            )
        )
    return tables

def _optional_standard_array(value: Any, shape: tuple[int, int]) -> np.ndarray | None:
    """Return an optional standard-report array if it has the expected shape.

    Parameters
    ----------
    value : Any
        Candidate value to validate.
    shape : tuple of int
        Expected ``(n_temperature, n_pressure)`` shape.

    Returns
    -------
    ndarray or None
        Numeric array with the requested shape, or ``None`` when the value is
        missing or incompatible.
    """
    if value is None:
        return None
    array = np.asarray(value, dtype=np.float64)
    if array.shape != shape:
        return None
    return array

def _format_optional_standard_value(
    array: np.ndarray | None,
    temperature_index: int,
    pressure_index: int,
    *,
    scale: float = 1.0,
    spec: str = ".12E",
) -> str:
    """Format an optional standard-report value.

    Parameters
    ----------
    array : ndarray or None
        Optional property array.
    temperature_index : int
        Temperature-grid index.
    pressure_index : int
        Pressure-grid index.
    scale : float, optional
        Multiplicative factor applied before rendering.
    spec : str, optional
        Python format specification.

    Returns
    -------
    str
        Formatted value or ``"-"`` when the property is unavailable.
    """
    if array is None:
        return "-"
    value = float(array[temperature_index, pressure_index]) * scale
    return format_number(value, spec)

def debug_thermodynamic_property_tables(
    result: QHAResult,
    *,
    table_format: QHATableFormat | None = None,
) -> list[ReportTable]:
    """Create detailed thermodynamic tables with optional uncertainties.

    Parameters
    ----------
    result : QHAResult
        QHA result containing pressure-temperature properties.
    table_format : QHATableFormat or None, optional
        Numerical formatting rules.

    Returns
    -------
    list of ReportTable
        One table for each pressure value. A ``sigma_*`` column follows each
        property for which an uncertainty array is available.
    """
    if result.temperature is None or result.pressure is None:
        return []
    fmt = table_format or QHATableFormat()
    properties: tuple[tuple[str, str, str], ...] = (
        ("U_zp", "zero_point_energy", "Uzp"),
        ("U_th", "thermal_energy", "Uth"),
        ("U_tot", "internal_energy", "Utot"),
        ("S", "entropy", "S"),
        ("C_V", "isochoric_heat_capacity", "Cv"),
        ("C_P", "isobaric_heat_capacity", "Cp"),
        ("C_P-C_V", "heat_capacity_difference", "Cp-Cv"),
        ("F_vib", "vibrational_free_energy", "Fvib"),
        ("F", "free_energy", "F"),
        ("H", "enthalpy", "H"),
        ("G", "gibbs_free_energy", "G"),
        ("V", "equilibrium_volume", "VT"),
        ("K_T", "isothermal_bulk_modulus", "KT"),
        ("K'", "bulk_modulus_derivative", "Kp"),
        ("K_S", "adiabatic_bulk_modulus", "KS"),
        ("alpha_V", "thermal_expansion", "alphaV"),
    )
    temperature = np.asarray(result.temperature, dtype=np.float64)
    pressure = np.asarray(result.pressure, dtype=np.float64)
    shape = (temperature.size, pressure.size)
    metadata = result.metadata if isinstance(result.metadata, dict) else {}
    units = (
        metadata.get("units", {}) if isinstance(metadata.get("units", {}), dict) else {}
    )
    temperature_unit = canonical_unit_label(str(units.get("temperature", "K")))
    pressure_unit = canonical_unit_label(str(units.get("pressure", "unknown")))

    arrays: list[tuple[str, str, np.ndarray, str]] = []
    for label, attr, symbol in properties:
        value = getattr(result, attr)
        if value is None:
            continue
        array = np.asarray(value, dtype=np.float64)
        if array.shape != shape:
            continue
        spec = _debug_format_for_property(attr, fmt)
        arrays.append((label, attr, array, spec))
        sigma = _uncertainty_for(result, attr, symbol)
        if sigma is not None:
            sigma_array = np.asarray(sigma, dtype=np.float64)
            if sigma_array.shape == shape:
                arrays.append((f"sigma_{label}", attr, sigma_array, spec))
    if not arrays:
        return []

    tables: list[ReportTable] = []
    for ip, pvalue in enumerate(pressure):
        rows: list[list[Any]] = []
        for it, tvalue in enumerate(temperature):
            rows.append(
                [format_number(float(tvalue), fmt.temperature)]
                + [
                    format_number(float(array[it, ip]), spec)
                    for _, _, array, spec in arrays
                ]
            )
        tables.append(
            _ReportTable(
                title=(
                    "Debug thermodynamic values at P = "
                    f"{format_number(float(pvalue), fmt.pressure)}"
                    + (f" {pressure_unit}" if pressure_unit != "unknown" else "")
                ),
                columns=["T"] + [label for label, _, _, _ in arrays],
                rows=rows,
                metadata={
                    "pressure": float(pvalue),
                    "formatted": True,
                    "column_units": [temperature_unit]
                    + [property_unit_label(attr, units) for _, attr, _, _ in arrays],
                    "column_alignments": ["right"] * (len(arrays) + 1),
                    "uncertainty_columns": [
                        label for label, _, _, _ in arrays if label.startswith("sigma_")
                    ],
                },
            )
        )
    return tables

def _debug_format_for_property(attr: str, fmt: QHATableFormat) -> str:
    """Return the numerical format for a detailed QHA report property."""
    if attr == "equilibrium_volume":
        return fmt.volume
    if attr in {"isothermal_bulk_modulus", "adiabatic_bulk_modulus"}:
        return fmt.bulk_modulus
    if attr == "bulk_modulus_derivative":
        return fmt.bulk_derivative
    if attr == "thermal_expansion":
        return fmt.thermal_expansion
    return fmt.energy

def _to_jmol(values: np.ndarray, result: QHAResult) -> np.ndarray:
    """Convert an energy-per-temperature array to J mol(f.u.)^-1 K^-1."""
    source_unit = "Ha cell^-1 K^-1"
    formula_units = 1.0
    metadata = result.metadata if isinstance(result.metadata, dict) else {}
    units = metadata.get("units")
    if isinstance(units, dict):
        source_unit = str(
            units.get(
                "heat_capacity",
                f"{units.get('energy', 'Ha')} cell^-1 K^-1",
            )
        )
    normalization = metadata.get("normalization")
    if isinstance(normalization, dict):
        try:
            formula_units = float(normalization.get("formula_units_per_cell", 1.0))
        except (TypeError, ValueError):
            formula_units = 1.0
    if formula_units <= 0.0:
        formula_units = 1.0
    try:
        converted = np.asarray(
            convert_energy_per_temperature(
                values,
                source_unit,
                "J mol^-1 K^-1",
            ),
            dtype=np.float64,
        )
        if "mol" not in source_unit.lower():
            converted /= formula_units
        return converted
    except (NotImplementedError, ValueError):
        return values

def phonon_frequency_fit_tables(
    input_data: QHAInput,
    options: QHAOptions,
    *,
    zero_tolerance: float = 1.0e-12,
) -> tuple[list[ReportTable], ReportTable]:
    """Create diagnostic tables for frequency-volume polynomial fits.

    Parameters
    ----------
    input_data : QHAInput
        Input data containing mode-resolved frequencies.
    options : QHAOptions
        QHA options providing the polynomial degree.
    zero_tolerance : float, optional
        Modes whose absolute frequencies are below this value for all volumes
        are skipped.

    Returns
    -------
    tuple of list of ReportTable and ReportTable
        Per-q-point debug tables and a cumulative summary table.
    """
    if input_data.frequencies is None or input_data.volume is None:
        return [], _ReportTable(
            title="Phonon frequency fit statistics",
            columns=["Metric", "Value"],
            rows=[["Fitted modes", 0], ["Skipped zero modes", 0]],
        )

    volume = np.asarray(input_data.volume, dtype=np.float64)
    frequencies = np.asarray(input_data.frequencies, dtype=np.float64)
    qcoords: NDArray[np.float64]
    if input_data.qcoords is None:
        qcoords = np.full((frequencies.shape[0], 3), np.nan, dtype=np.float64)
    else:
        qcoords = np.asarray(input_data.qcoords, dtype=np.float64)
    degree = min(int(options.frequency_degree), max(volume.size - 1, 1))

    tables: list[ReportTable] = []
    all_r2: list[float] = []
    skipped = 0
    for iq in range(frequencies.shape[0]):
        rows: list[list[Any]] = []
        for ib in range(frequencies.shape[1]):
            values = frequencies[iq, ib, :]
            if np.all(np.abs(values) <= zero_tolerance):
                skipped += 1
                continue
            r2 = _polynomial_r_squared(volume, values, degree)
            all_r2.append(r2)
            rows.append([int(ib + 1), r2])
        if rows:
            coord = qcoords[iq] if iq < qcoords.shape[0] else np.full(3, np.nan)
            tables.append(
                _ReportTable(
                    title=(
                        f"q-point #{iq + 1} - "
                        f"[{coord[0]:.8g}, {coord[1]:.8g}, {coord[2]:.8g}]"
                    ),
                    columns=["Band", "R^2"],
                    rows=rows,
                    metadata={"qpoint_index": iq, "qcoord": coord.tolist()},
                )
            )

    summary = _r2_summary_table(
        "Phonon frequency fit statistics",
        all_r2,
        skipped=skipped,
        extra_rows=[["Polynomial degree", degree]],
    )
    return tables, summary

def _polynomial_r_squared(x: np.ndarray, y: np.ndarray, degree: int) -> float:
    """Return R-squared for a polynomial fit.

    Parameters
    ----------
    x, y : ndarray
        One-dimensional fit coordinates and observations.
    degree : int
        Polynomial degree.

    Returns
    -------
    float
        Coefficient of determination.
    """
    degree = min(int(degree), max(x.size - 1, 1))
    coeffs = np.polynomial.polynomial.polyfit(x, y, deg=degree)
    fitted = np.polynomial.polynomial.polyval(x, coeffs)
    return float(residual_metrics(y, fitted)["r_squared"])

def _r2_summary_table(
    title: str,
    values: list[float],
    *,
    skipped: int = 0,
    extra_rows: list[list[Any]] | None = None,
) -> ReportTable:
    """Create a cumulative R-squared summary table.

    Parameters
    ----------
    title : str
        Table title.
    values : list of float
        R-squared values.
    skipped : int, optional
        Number of skipped records.
    extra_rows : list of list, optional
        Additional rows prepended to the summary.

    Returns
    -------
    ReportTable
        Summary table.
    """
    rows = list(extra_rows or [])
    finite = np.asarray(
        [value for value in values if np.isfinite(value)], dtype=np.float64
    )
    rows.append(["Fitted records", int(finite.size)])
    rows.append(["Skipped zero records", int(skipped)])
    if finite.size:
        rows.extend(
            [
                ["Minimum R^2", float(np.min(finite))],
                ["Mean R^2", float(np.mean(finite))],
                ["Median R^2", float(np.median(finite))],
                ["Maximum R^2", float(np.max(finite))],
            ]
        )
    else:
        rows.extend([["Minimum R^2", None], ["Mean R^2", None], ["Maximum R^2", None]])
    return _ReportTable(title=title, columns=["Metric", "Value"], rows=rows)

def thermodynamic_fit_tables(
    sampled: HarmonicThermodynamicResult,
    options: QHAOptions,
) -> tuple[ReportTable, ReportTable]:
    """Create diagnostic tables for thermodynamic property fits.

    Parameters
    ----------
    sampled : HarmonicThermodynamicResult
        Harmonic thermodynamic properties sampled on the input volume grid.
    options : QHAOptions
        QHA options providing the polynomial degree.

    Returns
    -------
    tuple of ReportTable
        Detailed temperature-wise R-squared table and cumulative summary table.
    """
    if sampled.volume is None or sampled.temperature is None:
        empty = _ReportTable(
            title="Thermodynamic fit statistics",
            columns=["Metric", "Value"],
            rows=[["Fitted records", 0]],
        )
        return empty, empty

    volume = np.asarray(sampled.volume, dtype=np.float64)
    temperature = np.asarray(sampled.temperature, dtype=np.float64)
    degree = min(int(options.energy_degree), max(volume.size - 1, 1))
    properties: tuple[tuple[str, str], ...] = (
        ("Uzp", "zero_point_energy"),
        ("Uth", "thermal_energy"),
        ("Utot", "internal_energy"),
        ("S", "entropy"),
        ("Cv", "isochoric_heat_capacity"),
        ("Fvib", "vibrational_free_energy"),
        ("F", "free_energy"),
    )

    rows: list[list[Any]] = []
    cumulative: dict[str, list[float]] = {label: [] for label, _ in properties}
    for it, tvalue in enumerate(temperature):
        row: list[Any] = [float(tvalue)]
        for label, attr in properties:
            value = getattr(sampled, attr)
            vector = _property_vector_at_temperature(
                value, it, volume.size, temperature.size
            )
            if vector is None:
                row.append(None)
                continue
            r2 = _polynomial_r_squared(volume, vector, degree)
            cumulative[label].append(r2)
            row.append(r2)
        rows.append(row)

    details = _ReportTable(
        title="Thermodynamic property fit diagnostics",
        columns=["T"] + [label for label, _ in properties],
        rows=rows,
        metadata={"degree": degree},
    )

    summary_rows: list[list[Any]] = [["Polynomial degree", degree]]
    for label, values in cumulative.items():
        finite = np.asarray(
            [value for value in values if np.isfinite(value)], dtype=np.float64
        )
        if finite.size:
            summary_rows.append(
                [
                    label,
                    float(np.min(finite)),
                    float(np.mean(finite)),
                    float(np.max(finite)),
                    int(finite.size),
                ]
            )
        else:
            summary_rows.append([label, None, None, None, 0])
    summary = _ReportTable(
        title="Thermodynamic fit statistics",
        columns=["Property", "Minimum R^2", "Mean R^2", "Maximum R^2", "Fits"],
        rows=summary_rows,
    )
    return details, summary

def _property_vector_at_temperature(
    value: np.ndarray | None,
    temperature_index: int,
    nvolumes: int,
    ntemperatures: int,
) -> np.ndarray | None:
    """Return one volume-dependent property vector at a temperature."""
    if value is None:
        return None
    array = np.asarray(value, dtype=np.float64)
    if array.ndim == 0:
        return np.full(nvolumes, float(array), dtype=np.float64)
    if array.ndim == 1:
        if array.size == nvolumes:
            return array
        if array.size == 1:
            return np.full(nvolumes, float(array[0]), dtype=np.float64)
        return None
    if array.ndim == 2:
        if array.shape == (ntemperatures, nvolumes):
            return array[temperature_index]
        if array.shape == (nvolumes, ntemperatures):
            return array[:, temperature_index]
        if array.shape[0] == 1 and array.shape[1] == nvolumes:
            return array[0]
    return None

__all__ = [
    "debug_thermodynamic_property_tables",
    "phonon_frequency_fit_tables",
    "property_table",
    "selected_property_tables",
    "thermodynamic_fit_tables",
]
