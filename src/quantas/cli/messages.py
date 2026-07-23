# -*- coding: utf-8 -*-

"""Shared terminal messages and prompts for the Quantas Click interface."""

from __future__ import annotations

from pathlib import Path

import click

from quantas._metadata import AUTHOR_TEXT, COPYRIGHT_TEXT
from quantas._version import __version__
from quantas.cli.output import print_terminal_message
from quantas.core.events import EventLevel


def quantas_title() -> str:
    """Return the Quantas software title banner."""
    banner = r"""
________                       __
\_____  \  __ _______    _____/  |______    ______
 /  / \  \|  |  \__  \  /    \   __\__  \  /  ___/
/   \_/.  \  |  // __ \|   |  \  |  / __ \_\___ \
\_____\ \_/____/(____  /___|  /__| (____  /____  >
       \__>          \/     \/          \/     \/ """
    return (
        f"{banner}\n"
        f"{'v':>45}{__version__}\n"
        f"Authors: {AUTHOR_TEXT}\n"
        f"{COPYRIGHT_TEXT}\n"
    )


def quantas_error() -> str:
    """Return the large Quantas error banner."""
    return (
        r"""
_______________________________ ________ __________._.
\_   _____/\______   \______   \\_____  \\______   \ |
 |    __)_  |       _/|       _/ /   |   \|       _/ |
 |        \ |    |   \|    |   \/    |    \    |   \\|
/_______  / |____|_  /|____|_  /\_______  /____|_  /__
        \/         \/        \/         \/       \/ \/"""
        + "\n"
    )


def quantas_warning() -> str:
    """Return the large Quantas warning banner."""
    return (
        r"""
 __      __                     .__
/  \    /  \_____ _______  ____ |__| ____    ____
\   \/\/   /\__  \\_  __ \/    \|  |/    \  / ___\
 \        /  / __ \|  | \/   |  \  |   |  \/ /_/  >
  \__/\  /  (____  /__|  |___|  /__|___|  /\___  /
       \/        \/           \/        \//_____/"""
        + "\n"
    )


def quantas_finish() -> str:
    """Return the Quantas closing banner."""
    return (
        r"""
___________.__       .__       .__
\_   _____/|__| ____ |__| _____|  |__
 |    __)  |  |/    \|  |/  ___/  |  \
 |     \   |  |   |  \  |\___ \|   Y  \
 \___  /   |__|___|  /__/____  >___|  /
     \/            \/        \/     \/"""
        + "\nThank you for using this software!\n"
    )


def init_logfile(logfile: str | Path | None = None) -> None:
    """Create or truncate a CLI log file.

    Parameters
    ----------
    logfile : str, Path, or None, optional
        Destination log file. No file is created when omitted.
    """
    if logfile is not None:
        Path(logfile).write_text("", encoding="utf-8")


def echo(
    msg: str,
    logfile: str | Path | None = None,
    silent: bool = False,
    bold: bool = False,
) -> None:
    """Write a standard message to the terminal and optional log file."""
    print_terminal_message(msg, bold=bold, silent=silent)
    _append_log(logfile, msg)


def echo_error(
    msg: str,
    logfile: str | Path | None = None,
    silent: bool = False,
    bold: bool = False,
) -> None:
    """Write an error message to the terminal and optional log file."""
    print_terminal_message(msg, level=EventLevel.ERROR, bold=bold, silent=silent)
    _append_log(logfile, msg)


def echo_warning(
    msg: str,
    logfile: str | Path | None = None,
    silent: bool = False,
    bold: bool = False,
) -> None:
    """Write a warning message to the terminal and optional log file."""
    print_terminal_message(msg, level=EventLevel.WARNING, bold=bold, silent=silent)
    _append_log(logfile, msg)


def echo_highlight(
    msg: str,
    logfile: str | Path | None = None,
    silent: bool = False,
    bold: bool = True,
) -> None:
    """Write a highlighted message to the terminal and optional log file."""
    print_terminal_message(msg, bold=bold, silent=silent)
    _append_log(logfile, msg)


def confirm(
    question: str,
    default: bool = False,
    suffix: str = ": ",
    show_default: bool = True,
) -> bool:
    """Prompt the user for a yes-or-no confirmation.

    Parameters
    ----------
    question : str
        Prompt text.
    default : bool, optional
        Default answer selected by pressing Return.
    suffix : str, optional
        Prompt suffix.
    show_default : bool, optional
        Whether the default answer is shown.

    Returns
    -------
    bool
        User response.
    """
    return click.confirm(
        click.style(question, bold=True),
        default=default,
        prompt_suffix=suffix,
        show_default=show_default,
    )


def prompt(
    question: str,
    default: str | None = None,
    show_default: bool = True,
) -> str:
    """Prompt the user for a text value.

    Parameters
    ----------
    question : str
        Prompt text.
    default : str or None, optional
        Default response.
    show_default : bool, optional
        Whether the default response is shown.

    Returns
    -------
    str
        User response.
    """
    return str(
        click.prompt(
            click.style(question, bold=True),
            default=default,
            show_default=show_default,
        )
    )


def _append_log(logfile: str | Path | None, message: str) -> None:
    """Append one message to a log file when configured."""
    if logfile is None:
        return
    with Path(logfile).open("a", encoding="utf-8") as stream:
        stream.write(f"{message}\n")
