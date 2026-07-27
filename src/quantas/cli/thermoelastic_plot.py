# -*- coding: utf-8 -*-

"""Thermoelastic plotting command group."""

from __future__ import annotations

from pathlib import Path

import click

from quantas.cli.output import CLIOutput
from quantas.cli.thermoelastic_plot_compare import compare_plot
from quantas.cli.thermoelastic_plot_domain import domain_plot
from quantas.cli.thermoelastic_plot_fit import fit_plot
from quantas.cli.thermoelastic_plot_profile import profile_plot
from quantas.cli.thermoelastic_plot_pt import pt_plot
from quantas.api.common import ReportTable
from quantas.api.thermoelasticity import describe_plots, read_result

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
@click.option(
    "--archive",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Native thermoelastic result used to list only plot families actually "
        "available in that archive. Valid only with --list-plots."
    ),
)
@click.pass_context
def plot(ctx: click.Context, list_plots: bool, archive: Path | None) -> None:
    """Generate thermoelastic diagnostic and scientific figures."""
    if list_plots:
        output = CLIOutput(show_progress=False)
        try:
            if archive is None:
                output.table(
                    ReportTable(
                        "Thermoelastic plot families",
                        ["Subcommand", "Purpose"],
                        [list(item) for item in _PLOT_DESCRIPTIONS],
                    )
                )
            else:
                inventory = describe_plots(read_result(archive))
                output.table(
                    ReportTable(
                        "Available thermoelastic plot families",
                        ["Subcommand", "Name", "Kind", "Properties"],
                        [
                            [
                                item.key,
                                item.name,
                                item.plot_kind,
                                ", ".join(item.property_keys) or "context-defined",
                            ]
                            for item in inventory.representations
                        ],
                    )
                )
                if inventory.warnings:
                    output.table(
                        ReportTable(
                            "Plot discovery warnings",
                            ["Warning"],
                            [[warning] for warning in inventory.warnings],
                        )
                    )
            output.save()
        finally:
            output.close()
        return
    if archive is not None:
        raise click.UsageError("--archive is valid only with --list-plots")
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


plot.add_command(compare_plot)
plot.add_command(fit_plot)
plot.add_command(pt_plot)
plot.add_command(profile_plot)
plot.add_command(domain_plot)


__all__ = ["plot"]
