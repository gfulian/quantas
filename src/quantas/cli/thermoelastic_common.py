# -*- coding: utf-8 -*-

"""Shared CLI helpers for the thermoelastic workflow."""

from __future__ import annotations

from pathlib import Path

import click

from quantas.cli.messages import confirm
from quantas.models import ResultData
from quantas.api.thermoelasticity import Result as ThermoelasticResult


def thermoelastic_payload(result_data: ResultData) -> ThermoelasticResult:
    """Extract the thermoelastic payload from a generic result envelope.

    Parameters
    ----------
    result_data : ResultData
        Generic Quantas result envelope read through the public API.

    Returns
    -------
    ThermoelasticResult
        Typed thermoelastic scientific payload.

    Raises
    ------
    click.ClickException
        If the envelope does not contain a thermoelastic result."""
    payload = result_data.results.get("thermoelasticity")
    if not isinstance(payload, ThermoelasticResult):
        raise click.ClickException("archive lacks a thermoelasticity payload")
    return payload


def approve_output_replacement(path: Path, force: bool) -> bool:
    """Return whether one output path may be created or replaced.

    Parameters
    ----------
    path : Path
        Proposed output destination.
    force : bool
        Replace existing output without prompting when ``True``.

    Returns
    -------
    bool
        ``True`` when the path is absent, replacement is forced, or the user
        explicitly approves replacement."""
    return (
        not path.exists()
        or force
        or confirm(
            f"Output '{path}' already exists. Replace it?",
            default=False,
        )
    )


def require_output_replacement(path: Path, force: bool) -> None:
    """Require approval before replacing an existing output path.

    Parameters
    ----------
    path : Path
        Proposed output destination.
    force : bool
        Replace existing output without prompting when ``True``.

    Raises
    ------
    click.Abort
        If the path exists and replacement is not approved."""
    if not approve_output_replacement(path, force):
        raise click.Abort()


__all__ = [
    "approve_output_replacement",
    "require_output_replacement",
    "thermoelastic_payload",
]
