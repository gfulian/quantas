# -*- coding: utf-8 -*-

"""Input-generation commands for quasi-static thermoelasticity."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import click

from quantas.cli.contracts import OUTPUT_GROUP
from quantas.cli.grouped_options import GroupedCommand, grouped_option
from quantas.cli.messages import quantas_finish, quantas_title
from quantas.cli.output import CLIOutput
from quantas.cli.thermoelastic_common import require_output_replacement
from quantas.models import ReportTable
from quantas.api.thermoelasticity import (
    create_input as create_thermoelastic_input,
    read_input as read_thermoelastic_input,
    write_profile_template as write_thermoelastic_profile_template,
)


@click.command(name="inpgen", cls=GroupedCommand)
@click.argument(
    "sources",
    nargs=-1,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@grouped_option(
    "-o",
    "--output",
    "outfile",
    group=OUTPUT_GROUP,
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Destination thermoelastic YAML input.",
)
@grouped_option(
    "-f",
    "--force",
    group=OUTPUT_GROUP,
    is_flag=True,
    default=False,
    help="Replace an existing YAML file without prompting.",
)
@grouped_option(
    "--list",
    "is_list",
    group="Input selection",
    is_flag=True,
    default=False,
    help="Interpret the single SOURCE argument as a text file listing CRYSTAL outputs.",
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
    "--jobname",
    group="Metadata",
    default="Quantas quasi-static thermoelastic input",
    show_default=True,
    help="Human-readable workflow description.",
)
@grouped_option(
    "--reference",
    group="Reference state",
    type=click.IntRange(min=0),
    default=None,
    help=(
        "Reference point index after sorting by increasing volume. "
        "Default: closest to zero pressure."
    ),
)
@grouped_option(
    "--symprec",
    group="Symmetry controls",
    type=click.FloatRange(min=0.0, min_open=True),
    default=1.0e-5,
    show_default=True,
    help="spglib Cartesian symmetry tolerance in angstrom.",
)
@grouped_option(
    "--angle-tolerance",
    group="Symmetry controls",
    type=float,
    default=-1.0,
    show_default=True,
    help="spglib angular tolerance in degrees; -1 uses its internal choice.",
)
@grouped_option(
    "--elastic-tolerance",
    group="Symmetry controls",
    type=click.FloatRange(min=0.0),
    default=1.0e-3,
    show_default=True,
    help="Absolute tolerance in GPa for elastic-symmetry detection.",
)
@grouped_option(
    "--pressure-tolerance",
    group="CRYSTAL validation",
    type=click.FloatRange(min=0.0),
    default=5.0e-2,
    show_default=True,
    help=(
        "Maximum difference in GPa between PRESSURE and the corrected elastic pressure."
    ),
)
@grouped_option(
    "--structure-correspondence-tolerance",
    group="CRYSTAL validation",
    type=click.FloatRange(min=0.0),
    default=5.0e-1,
    show_default=True,
    help=(
        "Maximum ordered-atom displacement in angstrom along the sampled "
        "structural path before input generation fails."
    ),
)
def inpgen(
    sources: tuple[Path, ...],
    outfile: Path,
    force: bool,
    is_list: bool,
    interface: str,
    jobname: str,
    reference: int | None,
    symprec: float,
    angle_tolerance: float,
    elastic_tolerance: float,
    pressure_tolerance: float,
    structure_correspondence_tolerance: float,
) -> None:
    """Generate readable QSA input from CRYSTAL SOEC output files."""
    if not sources:
        raise click.UsageError("provide at least one CRYSTAL output or one list file")
    if is_list and len(sources) != 1:
        raise click.UsageError("--list requires exactly one SOURCE list file")
    require_output_replacement(outfile, force)
    source_value: Path | Sequence[Path] = sources[0] if is_list else sources
    try:
        destination = create_thermoelastic_input(
            source_value,
            outfile,
            interface=interface,
            is_list=is_list,
            jobname=jobname,
            reference=reference,
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            elastic_tolerance=elastic_tolerance,
            pressure_tolerance=pressure_tolerance,
            structure_correspondence_tolerance=(structure_correspondence_tolerance),
        )
        parsed = read_thermoelastic_input(destination)
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc
    series = parsed.elastic_series
    output = CLIOutput(show_progress=False)
    output.message(quantas_title(), bold=True)
    output.table(
        ReportTable(
            "Thermoelastic input summary",
            ["Property", "Value"],
            [
                ["Output", str(destination)],
                ["Elastic points", series.npoints],
                ["Volume minimum (Å³)", series.volume_bounds[0]],
                ["Volume maximum (Å³)", series.volume_bounds[1]],
                ["Reference index", series.reference_index],
                ["Space group", series.symmetry.international_symbol],
                ["Elastic symmetry", series.elastic_symmetry],
                ["Pre-stress convention", "CRYSTAL PRESSURE / Wallace"],
                [
                    "Frame normalization",
                    series.metadata.get("frame_normalization", {}).get(
                        "method", "unavailable"
                    ),
                ],
                [
                    "Maximum removed rotation (degrees)",
                    series.metadata.get("frame_normalization", {}).get(
                        "maximum_removed_rotation_degrees"
                    ),
                ],
            ],
        )
    )
    output.message(quantas_finish())
    output.save()


@click.command(name="profile-template", cls=GroupedCommand)
@click.argument(
    "outfile",
    type=click.Path(dir_okay=False, path_type=Path),
    required=False,
    default=Path("thermoelastic_profile.yaml"),
)
@grouped_option(
    "-f",
    "--force",
    group=OUTPUT_GROUP,
    is_flag=True,
    default=False,
    help="Replace an existing template without prompting.",
)
def profile_template(outfile: Path, force: bool) -> None:
    """Write an editable geothermobarometric profile YAML template."""
    destination = outfile.with_suffix(".yaml")
    require_output_replacement(destination, force)
    try:
        written = write_thermoelastic_profile_template(destination, overwrite=True)
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"Thermoelastic profile template written to: {written}")
    click.echo(
        "Review the lithospheric parameters, then use "
        "'quantas thermoelasticity analysis profile FIT.hdf5 "
        f"--profile-spec {written}'."
    )


__all__ = ["inpgen", "profile_template"]
