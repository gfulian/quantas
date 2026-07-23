# -*- coding: utf-8 -*-

"""Text report adapters for the seismic CLI."""

from __future__ import annotations

from quantas.models import ResultData
from quantas.api.seismic import build_report
from quantas.renderers.tables import render_tables


def render_seismic_cli_report(result: ResultData) -> str:
    """Render the complete seismic report with the shared text backend.

    Parameters
    ----------
    result : ResultData
        Complete Quantas seismic result.

    Returns
    -------
    str
        Plain-text report suitable for terminal and file output.
    """
    return render_tables(build_report(result))


__all__ = ["render_seismic_cli_report"]
