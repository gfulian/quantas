# -*- coding: utf-8 -*-

"""Shared CLI helpers for the thermoelastic workflow."""

from __future__ import annotations

from pathlib import Path

import click

from quantas.cli.messages import confirm
from quantas.models import ResultData
from quantas.api.thermoelasticity import Result as ThermoelasticResult


def thermoelastic_payload(result_data: ResultData) -> ThermoelasticResult:
    """Return the thermoelastic payload from a generic result envelope."""
    payload = result_data.results.get("thermoelasticity")
    if not isinstance(payload, ThermoelasticResult):
        raise click.ClickException("archive lacks a thermoelasticity payload")
    return payload


def approve_output_replacement(path: Path, force: bool) -> bool:
    """Return whether an output may be created or replaced."""
    return (
        not path.exists()
        or force
        or confirm(
            f"Output '{path}' already exists. Replace it?",
            default=False,
        )
    )


def require_output_replacement(path: Path, force: bool) -> None:
    """Abort cleanly when an existing output is not approved for replacement."""
    if not approve_output_replacement(path, force):
        raise click.Abort()


__all__ = [
    "approve_output_replacement",
    "require_output_replacement",
    "thermoelastic_payload",
]
