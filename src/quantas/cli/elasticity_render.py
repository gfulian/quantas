# -*- coding: utf-8 -*-

"""Text report adapters for the elasticity CLI."""

from __future__ import annotations

from quantas.modules.elasticity.models import (
    ElasticityInput,
    ElasticityOptions,
    ElasticityResult,
)
from quantas.modules.elasticity.report import build_elasticity_report
from quantas.renderers.tables import render_table, render_tables


def render_elasticity_cli_report(
    input_data: ElasticityInput,
    options: ElasticityOptions,
    result: ElasticityResult,
) -> str:
    """Render the complete elasticity report with the shared text backend."""
    return render_tables(build_elasticity_report(input_data, options, result))


__all__ = ["render_elasticity_cli_report", "render_table"]
