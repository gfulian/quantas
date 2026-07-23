# -*- coding: utf-8 -*-

"""Thermoelastic plotting command group."""

from __future__ import annotations

import click

from quantas.cli.output import CLIOutput
from quantas.cli.thermoelastic_plot_compare import compare_plot
from quantas.cli.thermoelastic_plot_domain import domain_plot
from quantas.cli.thermoelastic_plot_fit import fit_plot
from quantas.cli.thermoelastic_plot_profile import profile_plot
from quantas.cli.thermoelastic_plot_pt import pt_plot
from quantas.models import ReportTable

_PLOT_DESCRIPTIONS = (
    (
        "fit",
        "Observed Cij(V), finite-strain fit, confidence band, and residuals.",
    ),
    (
        "pt",
        "Cij(P,T), one-sigma uncertainty, or relative-uncertainty contour maps.",
    ),
    (
        "profile",
        "Absolute or relative elastic components along an archived depth profile.",
    ),
    (
        "compare",
        "Direct comparison of isothermal C^T and adiabatic C^S at fixed P or T.",
    ),
    (
        "domain",
        "QHA/elastic P-T coverage, extrapolation masks, and geological paths.",
    ),
)


@click.group(name="plot", invoke_without_command=True)
@click.option(
    "--list-plots",
    is_flag=True,
    default=False,
    help="List thermoelastic plot families and exit.",
)
@click.pass_context
def plot(ctx: click.Context, list_plots: bool) -> None:
    """Generate thermoelastic diagnostic and scientific figures."""
    if list_plots:
        output = CLIOutput(show_progress=False)
        try:
            output.table(
                ReportTable(
                    "Thermoelastic plot families",
                    ["Subcommand", "Purpose"],
                    [list(item) for item in _PLOT_DESCRIPTIONS],
                )
            )
            output.save()
        finally:
            output.close()
        return
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


plot.add_command(compare_plot)
plot.add_command(fit_plot)
plot.add_command(pt_plot)
plot.add_command(profile_plot)
plot.add_command(domain_plot)


__all__ = ["plot"]
