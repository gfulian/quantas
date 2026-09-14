# -*- coding: utf-8 -*-

"""Shared CLI adapter for Kieffer enrichment of HA/QHA YAML inputs."""

from __future__ import annotations

from pathlib import Path

import click

from quantas.api import ha as ha_api
from quantas.api import qha as qha_api
from quantas.cli.contracts import NUMERICAL_GROUP, OUTPUT_GROUP, SCIENTIFIC_GROUP
from quantas.cli.grouped_options import GroupedCommand, grouped_option
from quantas.cli.messages import confirm, echo_error
from quantas.cli.eos_model_type import ENERGY_EOS_MODEL


@click.command(
    name="add-kieffer",
    cls=GroupedCommand,
    help=(
        "Create a new HA/QHA YAML input enriched with Kieffer acoustic "
        "cutoffs calculated from elastic outputs."
    ),
)
@click.argument(
    "filename",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.argument(
    "elastic_outputs",
    nargs=-1,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@grouped_option(
    "--interface",
    type=click.Choice(["crystal"]),
    default="crystal",
    show_default=True,
    group=SCIENTIFIC_GROUP,
    help="Interface used to read the elastic output files.",
)
@grouped_option(
    "--elastic-list",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    group=SCIENTIFIC_GROUP,
    help="Text file listing elastic outputs, one path per line.",
)
@grouped_option(
    "--pressure-source",
    type=click.Choice(
        [
            "auto",
            "output-stress",
            "manual",
            "energy-eos",
            "energy-polynomial",
        ]
    ),
    default="auto",
    show_default=True,
    group=SCIENTIFIC_GROUP,
    help="Source of hydrostatic pressure for raw elastic tensors.",
)
@grouped_option(
    "--pressure",
    "manual_pressures",
    type=float,
    multiple=True,
    group=SCIENTIFIC_GROUP,
    help="Manual pressure in GPa; repeat once per elastic output.",
)
@grouped_option(
    "--eos",
    type=ENERGY_EOS_MODEL,
    default="BM3",
    show_default=True,
    group=SCIENTIFIC_GROUP,
    help=("Energy EOS used with --pressure-source energy-eos. Run "
        "'quantas eos show-models --domain ev' to list compatible models."),
)
@grouped_option(
    "--degree",
    "polynomial_degree",
    type=click.IntRange(min=2),
    default=3,
    show_default=True,
    group=NUMERICAL_GROUP,
    help="E(V) polynomial degree used with energy-polynomial pressure.",
)
@grouped_option(
    "--mu-order",
    type=click.IntRange(min=2),
    default=12,
    show_default=True,
    group=NUMERICAL_GROUP,
    help="Gauss-Legendre order in cos(theta) before refinement.",
)
@grouped_option(
    "--phi-order",
    type=click.IntRange(min=4),
    default=24,
    show_default=True,
    group=NUMERICAL_GROUP,
    help="Periodic azimuthal quadrature order before refinement.",
)
@grouped_option(
    "--refinement-factor",
    type=click.IntRange(min=2),
    default=2,
    show_default=True,
    group=NUMERICAL_GROUP,
    help="Integer refinement applied to both directional orders.",
)
@grouped_option(
    "--batch-size",
    type=click.IntRange(min=1),
    default=512,
    show_default=True,
    group=NUMERICAL_GROUP,
    help="Maximum Christoffel directions evaluated per batch.",
)
@grouped_option(
    "-o",
    "--output",
    "outfile",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    group=OUTPUT_GROUP,
    help="Output YAML file. Default: INPUT stem + '-kieffer.yaml'.",
)
@grouped_option(
    "-f",
    "--force",
    is_flag=True,
    default=False,
    group=OUTPUT_GROUP,
    help="Overwrite an existing output file without prompting.",
)
@click.pass_context
def add_kieffer(
    ctx: click.Context,
    filename: Path,
    elastic_outputs: tuple[Path, ...],
    interface: str,
    elastic_list: Path | None,
    pressure_source: str,
    manual_pressures: tuple[float, ...],
    eos: str,
    polynomial_degree: int,
    mu_order: int,
    phi_order: int,
    refinement_factor: int,
    batch_size: int,
    outfile: Path | None,
    force: bool,
) -> None:
    """Add Kieffer cutoffs to a copy of a Quantas phonon YAML input."""
    workflow = "qha" if ctx.parent and ctx.parent.info_name == "qha" else "ha"
    destination = outfile or filename.with_name(f"{filename.stem}-kieffer.yaml")
    if destination.resolve() == filename.resolve():
        raise click.UsageError("output must differ from the source phonon input")
    if elastic_list is not None and elastic_outputs:
        raise click.UsageError(
            "ELASTIC_OUTPUTS and --elastic-list are mutually exclusive"
        )
    outputs = (
        _read_elastic_list(elastic_list)
        if elastic_list is not None
        else elastic_outputs
    )
    if not outputs:
        raise click.UsageError(
            "provide one or more ELASTIC_OUTPUTS or use --elastic-list"
        )
    policy = pressure_source.replace("-", "_")
    if policy == "manual" and not manual_pressures:
        raise click.UsageError("manual pressure source requires --pressure values")
    if policy != "manual" and manual_pressures:
        raise click.UsageError("--pressure requires --pressure-source manual")
    if workflow == "ha" and policy in {"energy_eos", "energy_polynomial"}:
        raise click.UsageError(
            "energy-derived pressure requires a multi-volume QHA input"
        )
    if destination.exists() and not force:
        if not confirm(
            f"Output file {destination} exists. Overwrite it?", default=False
        ):
            return

    try:
        if workflow == "qha":
            output = qha_api.add_kieffer_input(
                filename,
                destination,
                outputs,
                interface=interface,
                pressure_policy=policy,
                manual_pressures_gpa=(manual_pressures if policy == "manual" else None),
                eos=eos,
                polynomial_degree=polynomial_degree,
                mu_order=mu_order,
                phi_order=phi_order,
                refinement_factor=refinement_factor,
                batch_size=batch_size,
            )
        else:
            output = ha_api.add_kieffer_input(
                filename,
                destination,
                outputs,
                interface=interface,
                pressure_policy=policy,
                manual_pressures_gpa=(manual_pressures if policy == "manual" else None),
                mu_order=mu_order,
                phi_order=phi_order,
                refinement_factor=refinement_factor,
                batch_size=batch_size,
            )
    except Exception as exc:
        echo_error(str(exc))
        raise click.Abort() from exc
    click.echo(f"Kieffer-enriched input written to: {output}")


def _read_elastic_list(filename: Path) -> tuple[Path, ...]:
    """Read non-empty, non-comment paths relative to a list file."""
    outputs: list[Path] = []
    for raw_line in filename.read_text(encoding="utf-8").splitlines():
        value = raw_line.strip()
        if not value or value.startswith("#"):
            continue
        path = Path(value)
        if not path.is_absolute():
            path = filename.parent / path
        if not path.is_file():
            raise click.UsageError(f"elastic output does not exist: {path}")
        outputs.append(path)
    return tuple(outputs)


__all__ = ["add_kieffer"]
