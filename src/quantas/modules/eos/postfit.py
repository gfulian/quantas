# -*- coding: utf-8 -*-

"""Reporting and CSV export for EOS diagnostics and calculated properties."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np

from quantas.models import ReportTable

from .calculator import EOSCalculationResult
from .diagnostics import EOSDiagnosticResult


def eos_calculation_summary_table(result: EOSCalculationResult) -> ReportTable:
    """Return source, model, propagation, and extrapolation information."""
    extrapolated = result.columns.get("extrapolated")
    rows: list[list[Any]] = [
        ["Record", result.record_id],
        ["Slot", result.slot.key],
        ["Calculated states", result.nrows],
        ["Model", _model_tag(result.metadata.get("model"))],
        ["Input mode", result.metadata.get("input_mode", "unknown")],
        ["Uncertainty propagation", result.metadata.get("uncertainty_method", "none")],
        [
            "Extrapolated states",
            0 if extrapolated is None else int(np.count_nonzero(extrapolated)),
        ],
    ]
    return ReportTable(
        "EOS calculator summary",
        ["Property", "Value"],
        rows,
        metadata={"notes": list(result.warnings)},
    )


def eos_calculation_table(result: EOSCalculationResult) -> ReportTable:
    """Return calculated values with adjacent one-sigma columns."""
    names: list[str] = []
    for name in result.columns:
        names.append(name)
        if name in result.uncertainties:
            names.append(f"sigma_{name}")
    rows: list[list[Any]] = []
    for index in range(result.nrows):
        row: list[Any] = []
        for name in result.columns:
            row.append(_report_value(name, result.columns[name][index]))
            if name in result.uncertainties:
                row.append(result.uncertainties[name][index])
        rows.append(row)
    units: list[str] = []
    formats: list[str | None] = []
    for name in result.columns:
        units.append(result.units.get(name, ""))
        formats.append(_property_format(name))
        if name in result.uncertainties:
            units.append(result.units.get(name, ""))
            formats.append("eos_uncertainty")
    return ReportTable(
        "EOS calculated properties",
        [_display_name(name) for name in names],
        rows,
        metadata={
            "column_units": units,
            "column_formats": formats,
            "column_alignments": ["right"] * len(names),
            "notes": [
                "Uncertainties are propagated from the fitted parameter covariance only."
            ]
            if result.uncertainties
            else [],
        },
    )


def eos_diagnostic_summary_table(result: EOSDiagnosticResult) -> ReportTable:
    """Return residual and normalized-pressure availability information."""
    normalized = result.metadata.get("normalized_pressure", {})
    rows: list[list[Any]] = [
        ["Record", result.record_id],
        ["Slot", result.slot.key],
        ["Model", result.metadata.get("model", "unknown")],
        ["Observations", result.nrows],
        ["Selected", result.metadata.get("selected_observations", "unknown")],
        ["Excluded", result.metadata.get("excluded_observations", "unknown")],
        ["Normalized pressure", bool(normalized.get("available", False))],
        ["Strain family", normalized.get("strain_family", "not available")],
    ]
    return ReportTable(
        "EOS diagnostics summary",
        ["Property", "Value"],
        rows,
        metadata={"notes": list(result.warnings)},
    )


def eos_diagnostic_table(result: EOSDiagnosticResult) -> ReportTable:
    """Return the complete residual and finite-strain diagnostic table."""
    names = list(result.columns)
    rows = [
        [_report_value(name, result.columns[name][index]) for name in names]
        for index in range(result.nrows)
    ]
    return ReportTable(
        "EOS residual and strain diagnostics",
        [_display_name(name) for name in names],
        rows,
        metadata={
            "column_units": [result.units.get(name, "") for name in names],
            "column_formats": [_property_format(name) for name in names],
            "column_alignments": ["right"] * len(names),
        },
    )


def write_eos_calculation_csv(
    result: EOSCalculationResult,
    path: str | Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Write calculated properties and propagated uncertainties to CSV."""
    columns: list[tuple[str, np.ndarray, str]] = []
    for name, values in result.columns.items():
        columns.append((name, values, result.units.get(name, "")))
        if name in result.uncertainties:
            columns.append(
                (
                    f"sigma_{name}",
                    result.uncertainties[name],
                    result.units.get(name, ""),
                )
            )
    return _write_csv(path, columns, result.nrows, overwrite=overwrite)


def write_eos_diagnostic_csv(
    result: EOSDiagnosticResult,
    path: str | Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Write residual and finite-strain diagnostics to CSV."""
    columns = [
        (name, values, result.units.get(name, ""))
        for name, values in result.columns.items()
    ]
    return _write_csv(path, columns, result.nrows, overwrite=overwrite)


def _write_csv(
    path: str | Path,
    columns: list[tuple[str, np.ndarray, str]],
    nrows: int,
    *,
    overwrite: bool,
) -> Path:
    destination = Path(path)
    if destination.exists() and not overwrite:
        raise FileExistsError(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    headers = [name if not unit else f"{name} [{unit}]" for name, _, unit in columns]
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        for index in range(nrows):
            writer.writerow([_csv_value(values[index]) for _, values, _ in columns])
    return destination


def _report_value(name: str, value: float) -> Any:
    if name in {"included", "extrapolated", "normalized_pressure_valid"}:
        return "yes" if bool(round(float(value))) else "no"
    if name in {"row", "group"}:
        return int(round(float(value)))
    return value


def _csv_value(value: float) -> str:
    number = float(value)
    return "" if np.isnan(number) else f"{number:.16g}"


def _model_tag(value: Any) -> str:
    if isinstance(value, dict):
        return str(value.get("tag", value))
    return str(value or "unknown")


def _display_name(name: str) -> str:
    replacements = {
        "sigma_": "σ(",
        "bulk_modulus": "Bulk modulus",
        "linear_modulus": "Linear modulus",
        "expansion_coefficient": "Expansion coefficient",
        "temperature_derivative": "dX/dT",
        "dK_dT_at_pressure": "dK/dT|P",
        "finite_strain": "Finite strain",
        "normalized_pressure": "Normalized pressure",
        "standardized_residual": "Standardized residual",
        "extrapolated": "Extrapolated",
        "included": "Included",
    }
    if name.startswith("sigma_"):
        return f"σ({name[6:].replace('_', ' ')})"
    if name in replacements:
        return replacements[name]
    return name.replace("_", " ").title()


def _property_format(name: str) -> str | None:
    if "pressure" in name or "modulus" in name or name == "residual":
        return "eos_pressure"
    if name == "temperature":
        return "eos_temperature"
    if name in {
        "row",
        "group",
        "included",
        "extrapolated",
        "normalized_pressure_valid",
    }:
        return None
    if name.startswith("sigma_"):
        return "eos_uncertainty"
    return "eos_structural"


__all__ = [
    "eos_calculation_summary_table",
    "eos_calculation_table",
    "eos_diagnostic_summary_table",
    "eos_diagnostic_table",
    "write_eos_calculation_csv",
    "write_eos_diagnostic_csv",
]
