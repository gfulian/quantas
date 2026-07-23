# -*- coding: utf-8 -*-

"""Compact thermoelastic table commands."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import click
import numpy as np

from quantas.cli.contracts import OUTPUT_GROUP
from quantas.cli.grouped_options import GroupedCommand, grouped_option
from quantas.cli.output import CLIOutput
from quantas.cli.thermoelastic_common import (
    approve_output_replacement,
    thermoelastic_payload,
)
from quantas.api.thermoelasticity import (
    analyze_grid as analyze_thermoelastic_result,
    point_table as thermoelastic_point_table,
    read_result as read_thermoelastic_hdf5,
    write_grid_table as write_thermoelastic_grid_table,
    write_profile_table as write_thermoelastic_profile_table,
)


def _table_output_options(function):
    decorators = (
        grouped_option(
            "-o",
            "--output",
            group=OUTPUT_GROUP,
            type=click.Path(dir_okay=False, path_type=Path),
            required=True,
            help="Destination CSV or deterministic text table.",
        ),
        grouped_option(
            "--format",
            "file_format",
            group=OUTPUT_GROUP,
            type=click.Choice(["csv", "text"], case_sensitive=False),
            default="csv",
            show_default=True,
            help="Output table format.",
        ),
        grouped_option(
            "--tensor-condition",
            group="Scientific selection",
            type=click.Choice(
                ["isothermal", "adiabatic", "both"], case_sensitive=False
            ),
            default="both",
            show_default=True,
            help="Select isothermal, adiabatic, or both stiffness tensors.",
        ),
        grouped_option(
            "--uncertainties/--no-uncertainties",
            group="Scientific selection",
            default=True,
            show_default=True,
            help="Include or omit one-sigma tensor uncertainties.",
        ),
        grouped_option(
            "-f",
            "--force",
            group=OUTPUT_GROUP,
            is_flag=True,
            default=False,
            help="Replace OUTPUT without prompting.",
        ),
    )
    for decorator in reversed(decorators):
        function = decorator(function)
    return function


@click.group(name="table")
def table() -> None:
    """Print or export compact scientific tables with independent Cij only."""


@table.command(name="point", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("pressure", type=float)
@click.argument("temperature", type=float)
@grouped_option(
    "--tensor-condition",
    group="Scientific selection",
    type=click.Choice(["isothermal", "adiabatic"], case_sensitive=False),
    default="adiabatic",
    show_default=True,
    help="Select the isothermal or adiabatic stiffness tensor.",
)
@grouped_option(
    "--extrapolation",
    group="Scientific selection",
    type=click.Choice(["fail", "warn", "allow"], case_sensitive=False),
    default="fail",
    show_default=True,
    help="Policy for a state outside QHA or elastic-volume support.",
)
def point_table(
    archive: Path,
    pressure: float,
    temperature: float,
    tensor_condition: str,
    extrapolation: str,
) -> None:
    """Print one pressure-temperature state as a compact table."""
    analyzed = analyze_thermoelastic_result(
        read_thermoelastic_hdf5(archive),
        pressure=np.asarray([pressure], dtype=np.float64),
        temperature=np.asarray([temperature], dtype=np.float64),
        extrapolation_policy=extrapolation.lower(),
    )
    output = CLIOutput(show_progress=False)
    try:
        output.table(
            thermoelastic_point_table(
                thermoelastic_payload(analyzed),
                tensor_condition=cast(
                    Literal["isothermal", "adiabatic"],
                    tensor_condition.lower(),
                ),
            )
        )
        output.save()
    finally:
        output.close()


@table.command(name="grid", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@_table_output_options
def grid_table(
    archive: Path,
    output: Path,
    file_format: str,
    tensor_condition: str,
    uncertainties: bool,
    force: bool,
) -> None:
    """Export an archived P-T grid in wide row-oriented form."""
    path = _normalize_path(output, file_format)
    if not approve_output_replacement(path, force):
        click.echo("Table not saved")
        return
    written = write_thermoelastic_grid_table(
        thermoelastic_payload(read_thermoelastic_hdf5(archive)),
        path,
        tensor_condition=cast(
            Literal["isothermal", "adiabatic", "both"], tensor_condition.lower()
        ),
        file_format=cast(Literal["csv", "text"], file_format.lower()),
        include_uncertainties=uncertainties,
        overwrite=True,
    )
    click.echo(f"Wide thermoelastic grid table written to: {written}")


@table.command(name="profile", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@grouped_option(
    "--profile-name",
    group="Profile selection",
    default=None,
    help="Optional only when the archive contains exactly one profile.",
)
@_table_output_options
def profile_table(
    archive: Path,
    profile_name: str | None,
    output: Path,
    file_format: str,
    tensor_condition: str,
    uncertainties: bool,
    force: bool,
) -> None:
    """Export one archived geological profile in wide row-oriented form."""
    path = _normalize_path(output, file_format)
    if not approve_output_replacement(path, force):
        click.echo("Table not saved")
        return
    written = write_thermoelastic_profile_table(
        thermoelastic_payload(read_thermoelastic_hdf5(archive)),
        path,
        profile_name=profile_name,
        tensor_condition=cast(
            Literal["isothermal", "adiabatic", "both"], tensor_condition.lower()
        ),
        file_format=cast(Literal["csv", "text"], file_format.lower()),
        include_uncertainties=uncertainties,
        overwrite=True,
    )
    click.echo(f"Wide thermoelastic profile table written to: {written}")


def _normalize_path(path: Path, file_format: str) -> Path:
    suffix = ".csv" if file_format.lower() == "csv" else ".dat"
    return path if path.suffix.lower() in {".csv", ".dat"} else path.with_suffix(suffix)


__all__ = ["table"]
