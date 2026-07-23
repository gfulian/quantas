# -*- coding: utf-8 -*-

"""Text report adapters for the harmonic-approximation CLI."""

from __future__ import annotations

from quantas.modules.ha.models import HAInput, HAOptions, HAResult
from quantas.modules.ha.report import (
    input_table,
    options_table,
    thermodynamic_summary_table,
)
from quantas.renderers.tables import render_table, render_tables


def render_ha_cli_report(
    input_data: HAInput | None,
    options: HAOptions | None,
    result: HAResult,
) -> str:
    """Render the standard HA report with the shared text backend.

    Parameters
    ----------
    input_data : HAInput or None
        Input data used by the workflow.
    options : HAOptions or None
        Scientific workflow options.
    result : HAResult
        Harmonic result payload.

    Returns
    -------
    str
        Plain-text report.
    """
    tables = []
    if options is not None:
        tables.append(options_table(options))
    if input_data is not None:
        tables.append(input_table(input_data))
    tables.append(thermodynamic_summary_table(result))
    return render_tables(tables)


def render_settings(options: HAOptions) -> str:
    """Render selected HA settings."""
    return render_table(options_table(options))


def render_input_data(input_data: HAInput) -> str:
    """Render the normalized HA input summary."""
    return render_table(input_table(input_data))


def render_thermodynamic_summary(result: HAResult) -> str:
    """Render the available harmonic properties."""
    return render_table(thermodynamic_summary_table(result))


__all__ = [
    "render_ha_cli_report",
    "render_input_data",
    "render_settings",
    "render_table",
    "render_thermodynamic_summary",
]
