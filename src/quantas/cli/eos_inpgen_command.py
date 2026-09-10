# -*- coding: utf-8 -*-

"""CLI adapter for Energy EOS input generation."""

from __future__ import annotations

from pathlib import Path
from typing import cast

import click

from quantas.api import eos as eos_api
from quantas.cli.contracts import OUTPUT_GROUP
from quantas.cli.grouped_options import GroupedCommand, grouped_option
from quantas.cli.messages import confirm, quantas_finish, quantas_title
from quantas.cli.output import CLIOutput
from quantas.io.path import ensure_suffix
from quantas.models import ReportTable


@click.command(name="inpgen", cls=GroupedCommand)
@click.argument(
    "source",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@grouped_option(
    "-o",
    "--output",
    "outfile",
    group=OUTPUT_GROUP,
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Destination Energy EOS text input.",
)
@grouped_option(
    "-f",
    "--force",
    group=OUTPUT_GROUP,
    is_flag=True,
    default=False,
    help="Replace an existing Energy EOS input without prompting.",
)
@grouped_option(
    "--interface",
    group="Input selection",
    type=click.Choice(["crystal"], case_sensitive=False),
    default="crystal",
    show_default=True,
    help="Electronic-structure output interface.",
)
@grouped_option(
    "--list",
    "is_list",
    group="Input selection",
    is_flag=True,
    default=False,
    help=(
        "Interpret SOURCE as a text file listing backend outputs. Each output "
        "may contribute one or several Energy EOS states."
    ),
)
@grouped_option(
    "--jobname",
    group="Metadata",
    default="Quantas Energy EOS input",
    show_default=True,
    help="Human-readable dataset title.",
)
def inpgen(
    source: Path,
    outfile: Path,
    force: bool,
    interface: str,
    is_list: bool,
    jobname: str,
) -> None:
    """Generate an Energy EOS dataset from electronic-structure output.

    CRYSTAL sources may contain one static/optimized state or a complete native
    multi-volume EOS series; ``--list`` flattens compatible sources into one
    volume-sorted table containing cell metrics and total energies.
    """
    destination = ensure_suffix(outfile, ".dat")
    if destination.exists() and not force and not confirm(
        f"Output '{destination}' already exists. Replace it?",
        default=False,
    ):
        raise click.Abort()

    try:
        written = eos_api.create_input(
            source,
            destination,
            interface=cast(eos_api.InputInterface, interface.lower()),
            is_list=is_list,
            jobname=jobname,
        )
        dataset = eos_api.read_input(written)
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc

    volume = dataset.column("volume")
    energy = dataset.column("energy")
    output = CLIOutput(show_progress=False)
    output.message(quantas_title(), bold=True)
    output.table(
        ReportTable(
            "Energy EOS input summary",
            ["Property", "Value"],
            [
                ["Output", str(written)],
                ["Interface", interface.lower()],
                ["States", dataset.npoints],
                ["Volume minimum (Å³)", float(volume.min())],
                ["Volume maximum (Å³)", float(volume.max())],
                ["Energy unit", dataset.units.get("energy", "Ha")],
                ["Minimum-energy volume (Å³)", float(volume[int(energy.argmin())])],
            ],
        )
    )
    output.message(quantas_finish())
    output.save()


__all__ = ["inpgen"]
