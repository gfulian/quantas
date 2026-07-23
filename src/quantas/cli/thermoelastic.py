# -*- coding: utf-8 -*-

"""Command-line entry point for quasi-static thermoelasticity."""

from __future__ import annotations

import click

from quantas.cli.reference_help import apply_reference_help

from quantas.cli.thermoelastic_analysis import analysis
from quantas.cli.thermoelastic_export import export
from quantas.cli.thermoelastic_inpgen import inpgen, profile_template
from quantas.cli.thermoelastic_inspect import inspect_archive
from quantas.cli.thermoelastic_plot import plot
from quantas.cli.thermoelastic_run import run
from quantas.cli.thermoelastic_table import table


@click.group(name="thermoelasticity")
def thermoelasticity() -> None:
    """Quasi-static QHA--elasticity coupling and Cij(P,T) analysis."""


thermoelasticity.add_command(inpgen)
thermoelasticity.add_command(profile_template)
thermoelasticity.add_command(run)
thermoelasticity.add_command(analysis)
thermoelasticity.add_command(table)
thermoelasticity.add_command(export)
thermoelasticity.add_command(inspect_archive)
thermoelasticity.add_command(plot)


__all__ = ["thermoelasticity"]

apply_reference_help(thermoelasticity, ('thermoelasticity',))
