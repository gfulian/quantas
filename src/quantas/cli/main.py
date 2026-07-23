# -*- coding: utf-8 -*-

"""Top-level Click entry point for the Quantas command-line interface."""

from __future__ import annotations

import click

from quantas.cli.general import MostSimilarCommandGroup, print_version
from quantas.cli.registry import register_commands
from quantas.cli.reference_help import apply_reference_help


@click.group(
    cls=MostSimilarCommandGroup,
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-v",
    "--version",
    is_flag=True,
    callback=print_version,
    expose_value=False,
    is_eager=True,
    help="Show the software version and exit.",
)
def main():
    """Quantas command-line interface."""
    pass


register_commands(main)
apply_reference_help(main, ("quantas",), recursive=False)
