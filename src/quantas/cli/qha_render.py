# -*- coding: utf-8 -*-

"""Text renderers for the quasi-harmonic approximation command-line interface.

The functions defined here convert neutral QHA report tables and structured QHA
containers into plain text. They do not print directly, so their output can be
used by Click commands, report files, tests, or other text frontends.
"""

from __future__ import annotations

from collections.abc import Sequence

from quantas.renderers.tables import render_table, render_tables


from quantas.modules.qha.inspect import PressureVolumePreview
from quantas.models import ReportTable
from quantas.models.thermodynamics import HarmonicThermodynamicResult
from quantas.modules.qha.models import QHAInput, QHAOptions, QHAResult
from quantas.modules.qha.report import (
    all_tables,
    debug_thermodynamic_property_tables,
    diagnostics_table,
    input_table,
    phonon_frequency_fit_tables,
    options_table,
    pressure_volume_preview_table,
    preview_diagnostics_table,
    preview_parameters_table,
    result_summary_table,
    selected_property_tables,
    structural_property_tables,
    thermodynamic_fit_tables,
)


def render_qha_cli_report(
    input_data: QHAInput | None,
    options: QHAOptions | None,
    result: QHAResult,
    *,
    include_diagnostics: bool = False,
    max_diagnostics: int | None = 20,
) -> str:
    """Render a compact QHA command-line report.

    Parameters
    ----------
    input_data : QHAInput or None
        Input data used by the workflow, when available.
    options : QHAOptions or None
        Options used by the workflow, when available.
    result : QHAResult
        Quasi-harmonic approximation result object.
    include_diagnostics : bool, optional
        If ``True``, include local fit diagnostics.
    max_diagnostics : int or None, optional
        Maximum number of diagnostic rows to render.

    Returns
    -------
    str
        Formatted report text.
    """
    tables = all_tables(input_data, options, result)
    tables.extend(selected_property_tables(result))
    tables.extend(structural_property_tables(result))
    if include_diagnostics and result.fit_records:
        tables.append(diagnostics_table(result, max_rows=max_diagnostics))
    return render_tables(tables)


def preview_report_tables(
    preview: PressureVolumePreview,
    *,
    include_diagnostics: bool = True,
) -> list[ReportTable]:
    """Return neutral pressure-volume preview tables."""
    tables = [pressure_volume_preview_table(preview)]
    if include_diagnostics:
        tables.append(preview_parameters_table(preview))
        tables.append(preview_diagnostics_table(preview))
    return tables


def phonon_frequency_fit_report_tables(
    input_data: QHAInput,
    options: QHAOptions,
    *,
    include_debug: bool = False,
) -> list[ReportTable]:
    """Return neutral frequency-volume fit report tables."""
    debug_tables, summary = phonon_frequency_fit_tables(input_data, options)
    tables = list(debug_tables) if include_debug else []
    tables.append(summary)
    return tables


def thermodynamic_fit_report_tables(
    sampled: HarmonicThermodynamicResult,
    options: QHAOptions,
    *,
    include_debug: bool = False,
) -> list[ReportTable]:
    """Return neutral thermodynamic fit report tables."""
    details, summary = thermodynamic_fit_tables(sampled, options)
    tables = [details] if include_debug else []
    tables.append(summary)
    return tables


def final_qha_result_tables(
    result: QHAResult,
    *,
    include_debug: bool = False,
    pressure_indices: Sequence[int] | None = None,
) -> list[ReportTable]:
    """Return neutral final QHA result tables."""
    tables = list(selected_property_tables(result, pressure_indices=pressure_indices))
    tables.extend(structural_property_tables(result, pressure_indices=pressure_indices))
    if include_debug:
        tables.extend(debug_thermodynamic_property_tables(result))
    return tables


def render_preview_report(
    preview: PressureVolumePreview,
    *,
    include_diagnostics: bool = True,
) -> str:
    """Render a pressure-volume preview report.

    Parameters
    ----------
    preview : PressureVolumePreview
        Structured pressure-volume preview.
    include_diagnostics : bool, optional
        If ``True``, include fit diagnostics for the preview fits.

    Returns
    -------
    str
        Formatted preview report text.
    """
    tables = preview_report_tables(
        preview,
        include_diagnostics=include_diagnostics,
    )
    text = "\n".join(render_table(table) for table in tables)
    if preview.warnings:
        warnings = "\n".join(f"WARNING: {warning}" for warning in preview.warnings)
        text = f"{text}\n{warnings}\n"
    return text


def render_phonon_frequency_fit_report(
    input_data: QHAInput,
    options: QHAOptions,
    *,
    include_debug: bool = False,
) -> str:
    """Render frequency-volume fit diagnostics.

    Parameters
    ----------
    input_data : QHAInput
        QHA input data containing mode-resolved frequencies.
    options : QHAOptions
        QHA options controlling the polynomial degree.
    include_debug : bool, optional
        If ``True``, render per-q-point and per-band R-squared values.

    Returns
    -------
    str
        Formatted fit report.
    """
    return render_tables(
        phonon_frequency_fit_report_tables(
            input_data,
            options,
            include_debug=include_debug,
        )
    )


def render_thermodynamic_fit_report(
    sampled: HarmonicThermodynamicResult,
    options: QHAOptions,
    *,
    include_debug: bool = False,
) -> str:
    """Render volume-dependent thermodynamic fit diagnostics.

    Parameters
    ----------
    sampled : HarmonicThermodynamicResult
        Harmonic thermodynamic properties sampled on the input volume grid.
    options : QHAOptions
        QHA options controlling the polynomial degree.
    include_debug : bool, optional
        If ``True``, include per-temperature R-squared values.

    Returns
    -------
    str
        Formatted fit report.
    """
    return render_tables(
        thermodynamic_fit_report_tables(
            sampled,
            options,
            include_debug=include_debug,
        )
    )


def render_final_qha_results(
    result: QHAResult,
    *,
    include_debug: bool = False,
) -> str:
    """Render final QHA result tables.

    Parameters
    ----------
    result : QHAResult
        QHA result object containing pressure-temperature properties.
    include_debug : bool, optional
        If ``True``, append detailed thermodynamic values at each pressure.

    Returns
    -------
    str
        Formatted final result tables.
    """
    tables = final_qha_result_tables(result, include_debug=include_debug)
    if not tables:
        return ""
    return render_tables(tables)


def render_input_data(input_data: QHAInput) -> str:
    """Render QHA input dimensions.

    Parameters
    ----------
    input_data : QHAInput
        Normalized QHA input data.

    Returns
    -------
    str
        Formatted input-data block.
    """
    return render_table(input_table(input_data))


def render_settings(options: QHAOptions) -> str:
    """Render selected QHA workflow options.

    Parameters
    ----------
    options : QHAOptions
        QHA workflow options.

    Returns
    -------
    str
        Formatted settings block.
    """
    return render_table(options_table(options))


def render_result_summary(result: QHAResult) -> str:
    """Render a compact QHA result summary.

    Parameters
    ----------
    result : QHAResult
        QHA result object.

    Returns
    -------
    str
        Formatted result-summary block.
    """
    return render_table(result_summary_table(result))


__all__ = [
    "final_qha_result_tables",
    "phonon_frequency_fit_report_tables",
    "preview_report_tables",
    "render_final_qha_results",
    "render_input_data",
    "render_settings",
    "render_phonon_frequency_fit_report",
    "render_preview_report",
    "render_qha_cli_report",
    "render_result_summary",
    "render_table",
    "render_thermodynamic_fit_report",
    "thermodynamic_fit_report_tables",
]
