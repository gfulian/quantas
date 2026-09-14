# -*- coding: utf-8 -*-

"""Generate shell-completion registration scripts for Quantas."""

from __future__ import annotations

import click
from click.shell_completion import get_completion_class

# Importing the module registers the PowerShell backend with Click.
from quantas.cli.shell_completion import PowerShellComplete  # noqa: F401


@click.command(name="completion")
@click.argument(
    "shell",
    type=click.Choice(["powershell", "bash", "zsh", "fish"], case_sensitive=False),
)
@click.pass_context
def completion(ctx: click.Context, shell: str) -> None:
    """Print the shell-completion registration script for SHELL.

    The output is intended to be evaluated by the selected shell.  For
    PowerShell, for example, run ``quantas completion powershell | Out-String |
    Invoke-Expression`` in the current session or place that line in
    ``$PROFILE`` for persistent completion.
    """
    root = ctx.find_root()
    cls = get_completion_class(shell.lower())
    if cls is None:  # pragma: no cover - choices and registration guard this
        raise click.ClickException(f"shell completion is unavailable for {shell!r}")
    completer = cls(
        root.command,
        {},
        root.info_name or "quantas",
        "_QUANTAS_COMPLETE",
    )
    click.echo(completer.source())


__all__ = ["completion"]
