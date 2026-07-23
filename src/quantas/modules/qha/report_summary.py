# -*- coding: utf-8 -*-

"""Input, workflow, and preview tables for QHA reports."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from quantas.core.math.fitting import FitResult
from quantas.core.physics.units import energy_to_pressure
from quantas.models.report import ReportTable as _ReportTable
from quantas.modules.qha.formatting import (
    QHATableFormat,
    canonical_unit_label,
    format_number,
    native_energy_unit_label,
    volume_unit_label,
)
from quantas.modules.qha.inspect import PressureVolumePreview
from quantas.modules.qha.models import QHAInput, QHAOptions, QHAResult
from quantas.modules.qha.report_common import QHA_PROPERTY_LABELS, _uncertainty_for

if TYPE_CHECKING:
    from quantas.models.report import ReportTable

def input_table(input_data: QHAInput) -> ReportTable:
    """Create a table describing QHA input data.

    Parameters
    ----------
    input_data : QHAInput
        Normalized quasi-harmonic approximation input data.

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
            ["Mode continuity", input_data.mode_continuity_status()],
            [
                "Source",
                str(input_data.source) if input_data.source is not None else "Unknown",
            ],
        ],
    )

def options_table(options: QHAOptions) -> ReportTable:
    """Create a table describing selected QHA workflow options.

    Parameters
    ----------
    options : QHAOptions
        Quasi-harmonic approximation workflow options.

    Returns
    -------
    ReportTable
        Table containing pressure-temperature ranges, units, methods, and
        diagnostic policies.
    """
    return _ReportTable(
        title="Selected options",
        columns=["Option", "Value"],
        rows=[
            ["Temperature minimum", float(options.temperature_min)],
            ["Temperature maximum", float(options.temperature_max)],
            ["Temperature step", float(options.temperature_step)],
            ["Pressure minimum", float(options.pressure_min)],
            ["Pressure maximum", float(options.pressure_max)],
            ["Pressure step", float(options.pressure_step)],
            ["Scheme", options.scheme],
            ["Minimization", options.minimization],
            ["EOS", options.eos],
            ["Energy polynomial degree", int(options.energy_degree)],
            ["Free-energy polynomial degree", int(options.free_energy_degree)],
            ["Frequency polynomial degree", int(options.frequency_degree)],
            ["Structural-path polynomial degree", int(options.structural_degree)],
            ["Calculate Gruneisen parameter", bool(options.calculate_gruneisen)],
            ["Calculate mode Gruneisen", bool(options.calculate_mode_gruneisen)],
            ["Thermal-expansion method", options.thermal_expansion_method],
            [
                "Minimum Cv / Dulong-Petit for gamma",
                float(options.gruneisen_min_cv_fraction),
            ],
            ["Polynomial derivative method", options.polynomial_derivative_method],
            ["Polynomial local-grid points", int(options.polynomial_grid_points)],
            [
                "Polynomial local-grid separation (%)",
                float(options.polynomial_grid_separation),
            ],
            ["Temperature unit", options.temperature_unit],
            ["Pressure unit", options.pressure_unit],
            ["Energy unit", options.energy_unit],
            ["Volume unit", options.volume_unit],
            ["Frequency unit", options.frequency_unit],
            ["Debug", bool(options.debug)],
            ["Fit failure policy", options.fit_failure_policy],
            ["Maximum consecutive failures", int(options.max_consecutive_failures)],
            ["Numerical precision", "float64"],
        ],
    )

def result_summary_table(result: QHAResult) -> ReportTable:
    """Create a compact summary of available QHA result arrays.

    Parameters
    ----------
    result : QHAResult
        Quasi-harmonic approximation result object.

    Returns
    -------
    ReportTable
        Table listing available properties, array shapes, and uncertainty
        availability.
    """
    rows: list[list[Any]] = []
    for attr, (symbol, description) in QHA_PROPERTY_LABELS.items():
        value = getattr(result, attr)
        sigma = _uncertainty_for(result, attr, symbol)
        rows.append(
            [
                symbol,
                description,
                "yes" if value is not None else "no",
                tuple(np.asarray(value).shape) if value is not None else "-",
                "yes" if sigma is not None else "no",
            ]
        )

    return _ReportTable(
        title="QHA properties",
        columns=["Key", "Description", "Available", "Shape", "Uncertainty"],
        rows=rows,
        metadata={
            "completed": bool(result.completed),
            "failed_points": len(result.failed_points),
            "fit_records": len(result.fit_records),
        },
    )

def thermal_expansion_provenance_table(result: QHAResult) -> ReportTable | None:
    """Create a table describing the thermal-expansion method actually used.

    Parameters
    ----------
    result : QHAResult
        Quasi-harmonic approximation result object.

    Returns
    -------
    ReportTable or None
        Provenance table, or ``None`` when no thermal-expansion source
        information is available.
    """
    metadata = result.metadata.get("thermal_expansion", {})
    source = result.thermal_expansion_source
    if not isinstance(metadata, dict):
        metadata = {}
    if not metadata and source is None:
        return None

    requested = metadata.get("requested_method", "unknown")
    selected = metadata.get("selected_method", requested)
    fallback = metadata.get("fallback_method", "none")
    source_counts = metadata.get("source_counts", {})
    if not isinstance(source_counts, dict):
        source_counts = {}
    if not source_counts and source is not None:
        codes = metadata.get(
            "source_codes",
            {
                "invalid": 0,
                "mixed_derivative": 1,
                "mode_gruneisen": 2,
                "numerical": 3,
                "numerical_fallback": 4,
            },
        )
        if isinstance(codes, dict):
            source_array = np.asarray(source)
            source_counts = {
                name: int(np.count_nonzero(source_array == int(code)))
                for name, code in codes.items()
            }

    total = int(sum(int(value) for value in source_counts.values()))
    rows: list[list[Any]] = [
        ["Requested method", requested, None, None],
        ["Effective method", selected, None, None],
        ["Fallback method", fallback, None, None],
    ]
    labels = {
        "mixed_derivative": "Mixed derivative",
        "mode_gruneisen": "Mode Gruneisen",
        "numerical": "Equilibrium-volume derivative",
        "numerical_fallback": "Numerical fallback",
        "invalid": "Unavailable",
    }
    for key in (
        "mixed_derivative",
        "mode_gruneisen",
        "numerical",
        "numerical_fallback",
        "invalid",
    ):
        count = int(source_counts.get(key, 0))
        fraction = 100.0 * count / total if total else 0.0
        rows.append([labels[key], key, count, fraction])

    return _ReportTable(
        title="Thermal-expansion provenance",
        columns=["Entry", "Method", "Points", "Fraction (%)"],
        rows=rows,
        metadata={"total_points": total},
    )

def diagnostics_table(result: QHAResult, max_rows: int | None = None) -> ReportTable:
    """Create a table with local fit diagnostics.

    Parameters
    ----------
    result : QHAResult
        Quasi-harmonic approximation result object.
    max_rows : int or None, optional
        Maximum number of diagnostic records to include.

    Returns
    -------
    ReportTable
        Table containing local fit status, quality, residual metrics, and
        diagnostic messages.
    """
    records = result.fit_records
    limit = len(records) if max_rows is None else min(max_rows, len(records))
    rows: list[list[Any]] = []
    for record in records[:limit]:
        fit = record.fit
        rows.append(
            [
                record.quantity,
                record.method,
                record.temperature,
                record.pressure,
                record.success,
                None if fit is None else fit.status.value,
                None if fit is None else fit.quality.value,
                None if fit is None else fit.rmse,
                record.message,
            ]
        )

    return _ReportTable(
        title="Fit diagnostics",
        columns=[
            "Quantity",
            "Method",
            "T",
            "P",
            "Success",
            "Status",
            "Quality",
            "RMSE",
            "Message",
        ],
        rows=rows,
        metadata={
            "truncated": limit < len(records),
            "total_rows": len(records),
        },
    )

def failed_points_table(result: QHAResult, max_rows: int | None = None) -> ReportTable:
    """Create a table with failed pressure-temperature points.

    Parameters
    ----------
    result : QHAResult
        Quasi-harmonic approximation result object.
    max_rows : int or None, optional
        Maximum number of failed points to include.

    Returns
    -------
    ReportTable
        Table containing failed pressure-temperature points and messages.
    """
    points = result.failed_points
    limit = len(points) if max_rows is None else min(max_rows, len(points))
    rows = [
        [point.temperature, point.pressure, point.stage, point.message]
        for point in points[:limit]
    ]
    return _ReportTable(
        title="Failed points",
        columns=["T", "P", "Stage", "Message"],
        rows=rows,
        metadata={
            "truncated": limit < len(points),
            "total_rows": len(points),
        },
    )

def pressure_volume_preview_table(preview: PressureVolumePreview) -> ReportTable:
    """Create a table from a pressure-volume preview.

    Parameters
    ----------
    preview : PressureVolumePreview
        Structured preview returned by :func:`quantas.modules.qha.inspect`.

    Returns
    -------
    ReportTable
        Table containing input volumes, static energies, and pressure estimates
        from polynomial and EOS fits.
    """
    fmt = QHATableFormat()
    rows = []
    for row in preview.table_rows():
        rows.append(
            [
                format_number(
                    np.nan if row["volume"] is None else float(row["volume"]),
                    fmt.volume,
                ),
                format_number(
                    np.nan if row["energy"] is None else float(row["energy"]),
                    fmt.energy,
                ),
                format_number(
                    np.nan
                    if row["pressure_polynomial"] is None
                    else float(row["pressure_polynomial"]),
                    fmt.bulk_modulus,
                ),
                format_number(
                    np.nan
                    if row["pressure_eos"] is None
                    else float(row["pressure_eos"]),
                    fmt.bulk_modulus,
                ),
            ]
        )
    volume_unit = volume_unit_label(str(preview.metadata.get("volume_unit", "A")))
    energy_unit = native_energy_unit_label(
        str(preview.metadata.get("energy_unit", "Ha"))
    )
    pressure_unit = canonical_unit_label(preview.pressure_unit)
    return _ReportTable(
        title="Pressure-volume preview",
        columns=["V", "E", "P polynomial", "P EOS"],
        rows=rows,
        metadata={
            "success": preview.success,
            "warnings": list(preview.warnings),
            "pressure_unit": pressure_unit,
            "column_units": [
                volume_unit,
                energy_unit,
                pressure_unit,
                pressure_unit,
            ],
            "column_alignments": ["right", "right", "right", "right"],
        },
    )

def preview_diagnostics_table(preview: PressureVolumePreview) -> ReportTable:
    """Create a diagnostic table for pressure-volume preview fits.

    Parameters
    ----------
    preview : PressureVolumePreview
        Structured pressure-volume preview.

    Returns
    -------
    ReportTable
        Table containing method status, quality, RMSE, and warnings.
    """
    rows: list[list[Any]] = []
    for estimate in (preview.polynomial, preview.eos):
        if estimate is None:
            continue
        rows.append(
            [
                estimate.method,
                estimate.success,
                estimate.fit.status.value,
                estimate.fit.quality.value,
                estimate.fit.rmse,
                "; ".join(estimate.warnings),
            ]
        )
    return _ReportTable(
        title="Pressure-volume fit diagnostics",
        columns=["Method", "Success", "Status", "Quality", "RMSE", "Warnings"],
        rows=rows,
        metadata={"warnings": list(preview.warnings)},
    )

def preview_parameters_table(preview: PressureVolumePreview) -> ReportTable:
    """Create a table with fitted and implied pressure-volume parameters.

    Parameters
    ----------
    preview : PressureVolumePreview
        Structured pressure-volume preview.

    Returns
    -------
    ReportTable
        Table containing parameter values, parameter source and one-sigma
        errors for the polynomial and EOS fits used in the preview.  The EOS
        table reports the complete physical parameter set, including values
        implied by the selected EOS order.  Bulk-modulus parameters are
        rendered in the selected pressure scale.
    """
    rows: list[list[Any]] = []
    for estimate in (preview.polynomial, preview.eos):
        if estimate is None or estimate.fit.parameters is None:
            continue

        if estimate.method == "eos":
            parameter_rows = _resolved_eos_parameter_rows(estimate.fit)
        else:
            names = _fit_parameter_names(estimate.method, estimate.fit)
            errors = estimate.fit.errors
            parameter_rows = []
            for index, value in enumerate(estimate.fit.parameters):
                name = names[index] if index < len(names) else f"p{index}"
                error = None
                if errors is not None and index < len(errors):
                    error = float(errors[index])
                parameter_rows.append((name, float(value), error, "fitted"))

        for name, value, error, source in parameter_rows:
            display_value, display_error = _convert_eos_parameter_for_report(
                estimate.method,
                name,
                value,
                error,
                preview,
            )
            rows.append(
                [
                    estimate.method,
                    estimate.metadata.get("eos", estimate.metadata.get("degree", "-")),
                    name,
                    source,
                    display_value,
                    display_error,
                    _parameter_unit(estimate.method, name, preview),
                ]
            )

    return _ReportTable(
        title="Pressure-volume fit parameters",
        columns=["Method", "Model", "Parameter", "Source", "Value", "Sigma", "Unit"],
        rows=rows,
        metadata={"warnings": list(preview.warnings)},
    )

def _resolved_eos_parameter_rows(
    fit: FitResult,
) -> list[tuple[str, float, float | None, str]]:
    """Return complete EOS parameter records from fit metadata.

    Parameters
    ----------
    fit : FitResult
        EOS energy-fit result.

    Returns
    -------
    list of tuple
        ``(name, value, sigma, source)`` records.  When the fit covariance is
        available, uncertainties of implied parameters are obtained from the
        analytical complete-parameter covariance transformation.
    """
    metadata = fit.metadata or {}
    resolved_names = metadata.get("resolved_parameter_order")
    resolved_values = metadata.get("resolved_parameters")
    sources = metadata.get("parameter_sources")
    if not isinstance(resolved_names, list) or not isinstance(resolved_values, list):
        names = _fit_parameter_names("eos", fit)
        errors = fit.errors
        if fit.parameters is None:
            return []
        return [
            (
                names[index] if index < len(names) else f"p{index}",
                float(value),
                None
                if errors is None or index >= len(errors)
                else float(errors[index]),
                "fitted",
            )
            for index, value in enumerate(fit.parameters)
        ]

    resolved_errors = metadata.get("resolved_errors")
    sigma_by_name: dict[str, float] = {}
    if isinstance(resolved_errors, list) and len(resolved_errors) == len(
        resolved_names
    ):
        sigma_by_name = {
            str(name): float(error)
            for name, error in zip(resolved_names, resolved_errors, strict=True)
        }
    else:
        free_names = _fit_parameter_names("eos", fit)
        free_errors = fit.errors
        if free_errors is not None:
            sigma_by_name = {
                name: float(free_errors[index])
                for index, name in enumerate(free_names)
                if index < len(free_errors)
            }
    source_by_name = sources if isinstance(sources, dict) else {}
    rows: list[tuple[str, float, float | None, str]] = []
    for name, value in zip(resolved_names, resolved_values, strict=True):
        parameter = str(name)
        source = str(source_by_name.get(parameter, "fitted"))
        rows.append((parameter, float(value), sigma_by_name.get(parameter), source))
    return rows

def _convert_eos_parameter_for_report(
    method: str,
    name: str,
    value: float,
    error: float | None,
    preview: PressureVolumePreview,
) -> tuple[float, float | None]:
    """Convert EOS parameters from energy-density to pressure units.

    Parameters
    ----------
    method : str
        Fitting method.
    name : str
        Parameter name.
    value : float
        Parameter value in fit units.
    error : float or None
        One-sigma error in fit units.
    preview : PressureVolumePreview
        Preview containing the selected units.

    Returns
    -------
    tuple of float and float or None
        Converted value and uncertainty.
    """
    if method != "eos":
        return float(value), None if error is None else float(error)

    energy_unit = str(preview.metadata.get("energy_unit", "Ha"))
    volume_unit = str(preview.metadata.get("volume_unit", "A"))
    upper = name.upper()
    if upper == "K0":
        converted = float(
            energy_to_pressure(value, energy_unit, volume_unit, preview.pressure_unit)
        )
        converted_error = None
        if error is not None:
            converted_error = float(
                energy_to_pressure(
                    error, energy_unit, volume_unit, preview.pressure_unit
                )
            )
        return converted, converted_error
    if upper in {"KPP", "K''"}:
        pressure_factor = float(
            energy_to_pressure(1.0, energy_unit, volume_unit, preview.pressure_unit)
        )
        return (
            float(value) / pressure_factor,
            None if error is None else float(error) / pressure_factor,
        )
    return float(value), None if error is None else float(error)

def _parameter_unit(method: str, name: str, preview: PressureVolumePreview) -> str:
    """Return the display unit for a pressure-volume fit parameter.

    Parameters
    ----------
    method : str
        Fitting method.
    name : str
        Parameter name.
    preview : PressureVolumePreview
        Pressure-volume preview containing unit metadata.

    Returns
    -------
    str
        Unit label used in rendered tables.
    """
    energy = native_energy_unit_label(
        str(preview.metadata.get("energy_unit", "energy"))
    )
    volume = volume_unit_label(str(preview.metadata.get("volume_unit", "length")))
    if method == "eos":
        upper = name.upper()
        if upper == "E0":
            return energy
        if upper == "K0":
            return preview.pressure_unit
        if upper == "V0":
            return volume
        if upper in {"KP", "K'"}:
            return "-"
        if upper in {"KPP", "K''"}:
            return f"{preview.pressure_unit}^-1"
    return "model unit"

def _fit_parameter_names(method: str, fit: FitResult) -> list[str]:
    """Return parameter names for a pressure-volume fit.

    Parameters
    ----------
    method : str
        Fitting method.
    fit : FitResult
        Fit result containing optional parameter metadata.

    Returns
    -------
    list of str
        Parameter labels.
    """
    metadata = fit.metadata or {}
    order = metadata.get("parameter_order")
    if isinstance(order, list):
        return [str(item) for item in order]
    if method == "eos":
        return ["E0", "K0", "KP", "V0"]
    if method == "polynomial":
        n = 0 if fit.parameters is None else len(fit.parameters)
        return [f"a{i}" for i in range(n)]
    n = 0 if fit.parameters is None else len(fit.parameters)
    return [f"p{i}" for i in range(n)]

__all__ = [
    "diagnostics_table",
    "failed_points_table",
    "input_table",
    "options_table",
    "pressure_volume_preview_table",
    "preview_diagnostics_table",
    "preview_parameters_table",
    "result_summary_table",
    "thermal_expansion_provenance_table",
]
