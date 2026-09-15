# -*- coding: utf-8 -*-

"""Registry of top-level Quantas command groups."""

from __future__ import annotations

from quantas.cli.completion_command import completion
from quantas.cli.elasticity import elasticity
from quantas.cli.eos import eos
from quantas.cli.ha import ha
from quantas.cli.qha import qha
from quantas.cli.seismic import seismic
from quantas.cli.thermoelastic import thermoelasticity


COMMANDS = [completion, elasticity, eos, ha, qha, seismic, thermoelasticity]


def register_commands(cli_group) -> None:
    """Register all top-level Quantas commands on one Click group.

    Parameters
    ----------
    cli_group : click.Group
        Root command group that receives the supported CLI command families.

    Notes
    -----
    Registration binds frontend adapters only; numerical workflow discovery lives
    in :mod:`quantas.api.registry`."""
    for command in COMMANDS:
        cli_group.add_command(command)
