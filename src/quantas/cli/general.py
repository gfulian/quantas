# -*- coding: utf-8 -*-
"""Friendly command discovery and version handling for the Quantas CLI."""

from __future__ import annotations

import difflib
from typing import Any

import click

from .messages import echo_error, echo_highlight, quantas_error, quantas_title


class MostSimilarCommandGroup(click.Group):
    """Suggest nearby command names when a user mistypes a subcommand."""

    def get_command(
        self,
        ctx: click.Context,
        cmd_name: str,
    ) -> click.Command | None:
        """Resolve a command or report a short list of likely alternatives.

        Parameters
        ----------
        ctx : click.Context
            Active Click context.
        cmd_name : str
            Command name supplied by the user.

        Returns
        -------
        click.Command or None
            Matching command. Click terminates the invocation after reporting
            alternatives when no exact match exists.
        """
        command = click.Group.get_command(self, ctx, cmd_name)
        if command is not None:
            return command

        matches = difflib.get_close_matches(
            cmd_name,
            self.list_commands(ctx),
            cutoff=0.5,
        )
        if not matches:
            matches = [
                candidate
                for candidate in sorted(self.list_commands(ctx))
                if candidate.startswith(cmd_name)
            ][:3]

        echo_error(quantas_error())
        if matches:
            choices = "\n".join(f"\t{match}" for match in sorted(matches))
            ctx.fail(
                f"'{cmd_name}' is not a Quantas command.\n\n"
                f"The most similar commands are:\n{choices}"
            )
        ctx.fail(
            f"'{cmd_name}' is not a Quantas command.\n\n"
            "No similar commands were found."
        )
        return None


def print_version(
    ctx: click.Context,
    param: click.Parameter,
    value: Any,
) -> None:
    """Print the Quantas banner for an eager ``--version`` option.

    Parameters
    ----------
    ctx : click.Context
        Active Click context.
    param : click.Parameter
        Version option invoking the callback.
    value : Any
        Parsed option value.
    """
    del param
    if not value or ctx.resilient_parsing:
        return
    echo_highlight(quantas_title(), bold=True)
    ctx.exit()
