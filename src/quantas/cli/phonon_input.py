# -*- coding: utf-8 -*-

"""Shared command-line adapter for HA/QHA phonon input generation.

The HA and QHA command groups register the same Click command object from this
module.  The numerical and parsing work remains delegated to the common
HA/QHA input generator in :mod:`quantas.modules.ha.api`; no reader, structure,
or YAML logic is duplicated in the command-line layer.
"""

from __future__ import annotations

from pathlib import Path

import click

from quantas.cli.grouped_options import GroupedCommand
from quantas.cli.messages import confirm, echo, echo_error, echo_warning
from quantas.io.phonons import PhononInputFileReader
from quantas.api.ha import create_input as create_ha_input


@click.command(
    name="inpgen",
    cls=GroupedCommand,
    help=(
        "Generate a Quantas HA/QHA YAML input from one phonon output, a "
        "CRYSTAL QHA output, or a list of volume-dependent phonon outputs."
    ),
)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "-o",
    "--output",
    "outfile",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Output Quantas YAML file. Default: input file base name + '.yaml'.",
)
@click.option(
    "--interface",
    type=click.Choice(["crystal", "crystal-qha", "phonopy"], case_sensitive=False),
    default="crystal",
    show_default=True,
    help="Interface used to read the QM phonon output.",
)
@click.option(
    "--list",
    "is_list",
    is_flag=True,
    default=False,
    help="Treat FILENAME as a text file listing QM output files.",
)
@click.option(
    "--reference",
    type=int,
    default=0,
    show_default=True,
    help="Reference file index for multi-file input generation.",
)
@click.option(
    "--jobname",
    default=None,
    help=(
        "Job description written to the YAML input. If omitted, it is "
        "requested interactively."
    ),
)
@click.option(
    "-Z",
    "--formula-units",
    type=click.IntRange(min=1),
    default=1,
    show_default=True,
    help="Number of chemical formula units in the normalization cell.",
)
@click.pass_context
def phonon_inpgen(
    ctx: click.Context,
    filename: Path,
    outfile: Path | None,
    interface: str,
    is_list: bool,
    reference: int,
    jobname: str | None,
    formula_units: int,
) -> None:
    """Generate a shared Quantas HA/QHA YAML input from phonon outputs.

    Parameters
    ----------
    ctx : click.Context
        Active Click context.  Its parent command determines whether the
        interactive default job description refers to HA or QHA.
    filename : pathlib.Path
        Quantum-mechanical output file or text file listing output files.
    outfile : pathlib.Path or None
        Destination YAML path.  If ``None``, the source suffix is replaced by
        ``.yaml``.
    interface : str
        Input interface: ``crystal``, ``crystal-qha``, or ``phonopy``.
    is_list : bool
        Interpret ``filename`` as a file list when ``True``.
    reference : int
        Reference item for multi-file structural ordering and metadata.
    jobname : str or None
        Description stored in the generated input.  If omitted, prompt the
        user with a workflow-specific default.
    formula_units : int
        Number of formula units in the thermodynamic normalization cell.

    Raises
    ------
    click.Abort
        If parsing, structural reconstruction, or YAML generation fails.
    """
    workflow = _parent_workflow_name(ctx)
    destination = outfile if outfile is not None else filename.with_suffix(".yaml")

    if jobname is None:
        jobname = click.prompt(
            "Provide a small description of the current system",
            default=f"Quantas {workflow.upper()} input",
            show_default=True,
        )

    if destination.exists():
        overwrite = confirm(
            f"Output file {destination} exists. Would you like to overwrite it?",
            default=False,
        )
        if not overwrite:
            echo("Input file not written")
            return

    try:
        output = create_ha_input(
            filename,
            destination,
            interface=interface.lower(),
            is_list=is_list,
            reference=reference,
            jobname=jobname,
            formula_units=formula_units,
        )
        summary = _read_generated_structure_summary(output)
    except Exception as exc:
        echo_error(str(exc))
        raise click.Abort() from exc

    echo(f"Input file written to: {output}")
    if summary is None:
        echo_warning(
            "No structural volume path was written. HA/QHA thermodynamics remain "
            "available, but lattice parameters and linear/tensorial thermal "
            "expansion cannot be calculated from this input."
        )
    else:
        natoms, nvol, symmetry = summary
        suffix = f", symmetry {symmetry}" if symmetry else ""
        echo(
            "Structural volume path written: "
            f"{nvol} structure(s), {natoms} atom(s) in the primitive "
            f"representation{suffix}."
        )


def _parent_workflow_name(ctx: click.Context) -> str:
    """Return the top-level workflow name owning the shared command.

    Parameters
    ----------
    ctx : click.Context
        Context of the shared ``inpgen`` command.

    Returns
    -------
    str
        ``"qha"`` when registered under the QHA group, otherwise ``"ha"``.
    """
    if ctx.parent is not None and ctx.parent.info_name == "qha":
        return "qha"
    return "ha"


def _read_generated_structure_summary(
    filename: str | Path,
) -> tuple[int, int, str | None] | None:
    """Read compact structural metadata from a generated phonon YAML file.

    Parameters
    ----------
    filename : str or pathlib.Path
        Generated Quantas HA/QHA YAML file.

    Returns
    -------
    tuple of int, int, and str or None, or None
        Primitive atom count, number of structural states, and optional
        international space-group symbol.  ``None`` is returned for valid
        historical inputs without a ``structure`` block.

    Raises
    ------
    ValueError
        If the newly generated file cannot be read back consistently.
    """
    reader = PhononInputFileReader(filename)
    if not reader.completed:
        raise ValueError(reader.error or "Generated phonon input is not readable")

    # Read the summary from the raw optional YAML block rather than from the
    # ``PhononInputFileReader.structure`` convenience property.  This keeps
    # the shared CLI command compatible with installations in which the input
    # generator has already been updated but the neutral reader still belongs
    # to an earlier Quantas snapshot.  The numerical QHA reader will still
    # require the full structural models before it can consume the block.
    raw_data = getattr(reader, "data", None)
    if raw_data is None:
        raw_data = getattr(reader, "_data", None)
    if not isinstance(raw_data, dict):
        raise ValueError("Generated phonon input does not contain a YAML mapping")

    structure = raw_data.get("structure")
    if not isinstance(structure, dict):
        return None

    atomic_numbers = structure.get("atomic_numbers")
    if not isinstance(atomic_numbers, list):
        reference = structure.get("reference", {})
        atomic_numbers = (
            reference.get("atomic_numbers") if isinstance(reference, dict) else None
        )
    if not isinstance(atomic_numbers, list) or not atomic_numbers:
        raise ValueError("Generated structure block has no primitive atomic numbers")

    volume_series = structure.get("volume_series")
    if not isinstance(volume_series, dict):
        raise ValueError("Generated structure block has no structural volume series")
    lattices = volume_series.get("lattice")
    if not isinstance(lattices, list) or not lattices:
        raise ValueError("Generated structure block has no lattice matrices")

    symbol = None
    symmetry = structure.get("symmetry")
    if isinstance(symmetry, dict):
        raw_symbol = symmetry.get("international_symbol")
        if raw_symbol not in (None, ""):
            symbol = str(raw_symbol)

    return len(atomic_numbers), len(lattices), symbol
