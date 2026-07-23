# -*- coding: utf-8 -*-

"""EOS specification-template command."""

from __future__ import annotations

from pathlib import Path

import click

from quantas.cli.contracts import OUTPUT_GROUP

from quantas.cli.eos_common import GroupedCommand, grouped_option
from quantas.cli.output import print_terminal_message
from quantas.api.eos import (
    SPEC_TEMPLATE_FILENAME as EOS_SPEC_TEMPLATE_FILENAME,
    write_spec_template as write_eos_spec_template,
)


@click.command(name="spec-template", cls=GroupedCommand)
@click.argument(
    "output",
    required=False,
    default=EOS_SPEC_TEMPLATE_FILENAME,
    type=click.Path(dir_okay=False, path_type=Path),
)
@grouped_option(
    "-f",
    "--force",
    group=OUTPUT_GROUP,
    is_flag=True,
    default=False,
    help="Replace an existing template file. The command never prompts.",
)
def spec_template(output: Path, force: bool) -> None:
    """Generate a complete commented ``QUANTAS EOS SPEC 1`` template.

    OUTPUT defaults to ``eos.spec``.  The suffix has no semantic meaning.
    The generated file contains one active minimal job and commented examples
    for every supported section, model family, solver, and constraint form.
    """
    try:
        destination = write_eos_spec_template(output, overwrite=force)
    except FileExistsError as exc:
        raise click.UsageError(
            f"template file exists: {output}; use --force to replace it"
        ) from exc
    except OSError as exc:
        raise click.ClickException(str(exc)) from exc
    print_terminal_message(f"EOS specification template written to: {destination}")


__all__ = ["spec_template"]
