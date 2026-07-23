# -*- coding: utf-8 -*-

"""Stable facade for neutral quasi-harmonic report tables."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quantas.modules.qha.models import QHAInput, QHAOptions, QHAResult
from quantas.modules.qha.report_common import QHA_PROPERTY_LABELS, resolve_property_name
from quantas.modules.qha.report_structural import structural_property_tables
from quantas.modules.qha.report_summary import (
    diagnostics_table,
    failed_points_table,
    input_table,
    options_table,
    pressure_volume_preview_table,
    preview_diagnostics_table,
    preview_parameters_table,
    result_summary_table,
    thermal_expansion_provenance_table,
)
from quantas.modules.qha.report_thermodynamic import (
    debug_thermodynamic_property_tables,
    phonon_frequency_fit_tables,
    property_table,
    selected_property_tables,
    thermodynamic_fit_tables,
)

if TYPE_CHECKING:
    from quantas.models.report import ReportTable

def all_tables(
    input_data: QHAInput | None,
    options: QHAOptions | None,
    result: QHAResult,
) -> list[ReportTable]:
    """Build the standard set of QHA report tables.

    Parameters
    ----------
    input_data : QHAInput or None
        Input data used by the workflow. If ``None``, the input table is not
        included.
    options : QHAOptions or None
        Workflow options. If ``None``, the options table is not included.
    result : QHAResult
        Quasi-harmonic approximation result object.

    Returns
    -------
    list of ReportTable
        Neutral QHA report tables.
    """
    tables: list[ReportTable] = []
    if input_data is not None:
        tables.append(input_table(input_data))
    if options is not None:
        tables.append(options_table(options))
    tables.append(result_summary_table(result))
    provenance = thermal_expansion_provenance_table(result)
    if provenance is not None:
        tables.append(provenance)
    if result.failed_points:
        tables.append(failed_points_table(result))
    return tables

__all__ = [
    "QHA_PROPERTY_LABELS",
    "all_tables",
    "debug_thermodynamic_property_tables",
    "diagnostics_table",
    "failed_points_table",
    "input_table",
    "options_table",
    "phonon_frequency_fit_tables",
    "pressure_volume_preview_table",
    "preview_diagnostics_table",
    "preview_parameters_table",
    "property_table",
    "resolve_property_name",
    "result_summary_table",
    "selected_property_tables",
    "structural_property_tables",
    "thermal_expansion_provenance_table",
    "thermodynamic_fit_tables",
]
