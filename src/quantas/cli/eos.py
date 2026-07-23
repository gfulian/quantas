# -*- coding: utf-8 -*-

"""Command registry for equation-of-state workflows.

Individual command implementations live in dedicated modules.  This historical
import path remains the stable Click group entry point.
"""

from __future__ import annotations

import click

from quantas.cli.reference_help import apply_reference_help

from quantas.cli.eos_postfit_commands import calculate, diagnose, plot
from quantas.cli.eos_run_command import run
from quantas.cli.eos_spec_command import spec_template


@click.group(name="eos")
def eos() -> None:
    """Equation-of-state analysis and persistent EOS archives."""
    return


eos.add_command(spec_template)
eos.add_command(run)
eos.add_command(diagnose)
eos.add_command(plot)
eos.add_command(calculate)


__all__ = ["eos"]

apply_reference_help(eos, ('eos',))
