# -*- coding: utf-8 -*-

"""Shared terminal and report output services for the Quantas CLI.

The terminal path uses Rich, while report files always use the deterministic
plain-text table renderer.  This module deliberately contains no scientific
logic and no workflow-specific event interpretation.
"""

from __future__ import annotations

import os
from pathlib import Path
import sys
from typing import Sequence, TextIO

from rich.console import Console
from rich.text import Text

from quantas.core.events import Event, EventLevel
from quantas.models import ReportTable
from quantas.renderers.tables import build_rich_renderable, render_table

from .progress import RichProgressDisplay


def create_console(*, file: TextIO | None = None) -> Console:
    """Create a Rich console bound directly to a standard text stream.

    Parameters
    ----------
    file : text stream or None, optional
        Explicit output stream. When omitted, :data:`sys.stdout` is used
        directly. This keeps Click responsible only for command parsing and
        lets Rich perform its own terminal and Windows capability detection.

    Returns
    -------
    rich.console.Console
        Console configured for automatic terminal capability detection, with
        a conservative plain-text fallback for redirected or explicitly
        non-ANSI output.
    """
    stream = sys.stdout if file is None else file
    plain_output = _requires_plain_output(stream)
    return Console(
        file=stream,
        color_system=None if plain_output else "auto",
        force_terminal=False if plain_output else None,
        no_color=plain_output,
        highlight=False,
        markup=False,
        soft_wrap=plain_output,
    )


def _requires_plain_output(stream: TextIO) -> bool:
    """Return whether Rich styling must be disabled for one stream.

    Parameters
    ----------
    stream : text stream
        Candidate terminal or redirected output stream.

    Returns
    -------
    bool
        ``True`` for redirected output, ``TERM=dumb``, or an explicit
        ``QUANTAS_NO_ANSI``/``NO_COLOR`` request.
    """
    if os.environ.get("QUANTAS_NO_ANSI") or os.environ.get("NO_COLOR"):
        return True
    if os.environ.get("TERM", "").lower() == "dumb":
        return True
    isatty = getattr(stream, "isatty", None)
    if not callable(isatty):
        return True
    try:
        return not bool(isatty())
    except (OSError, ValueError):
        return True


class CLIOutput:
    """Render terminal output with Rich and collect a plain-text report.

    Parameters
    ----------
    report_file : str, Path or None, optional
        Destination written by :meth:`save`.
    silent : bool, optional
        If ``True``, terminal rendering is disabled while report collection
        remains active.
    console : rich.console.Console or None, optional
        Explicit console, primarily useful for tests.
    show_progress : bool, optional
        Enable terminal-only live progress rendering when supported.
    """

    def __init__(
        self,
        report_file: str | Path | None = None,
        *,
        silent: bool = False,
        console: Console | None = None,
        show_progress: bool = True,
    ) -> None:
        self.report_file = Path(report_file) if report_file is not None else None
        self.silent = bool(silent)
        self.console = create_console() if console is None else console
        self.progress_display = RichProgressDisplay(
            self.console,
            enabled=show_progress and not self.silent,
        )
        self._chunks: list[tuple[str, str]] = []
        self._terminal_last_kind: str | None = None

    def __enter__(self) -> "CLIOutput":
        """Return this output service for use as a context manager."""
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        """Always release live terminal resources on context exit."""
        self.close()

    def message(
        self,
        message: str,
        *,
        level: EventLevel = EventLevel.INFO,
        persist: bool = True,
        bold: bool = False,
    ) -> None:
        """Render one message and optionally append it to the plain report."""
        if not message:
            return
        plain = _plain_message(message, level)
        if persist:
            self._append_chunk("message", plain)
        if self.silent:
            return
        if level == EventLevel.ERROR:
            self.progress_display.stop()
        self._prepare_terminal_block("message")
        self.console.print(_terminal_message(message, level, bold=bold))
        self._terminal_last_kind = "message"

    def text_block(self, text: str, *, persist: bool = True) -> None:
        """Render preformatted text without interpreting Rich markup."""
        if not text:
            return
        if persist:
            self._append_chunk("block", text)
        if not self.silent:
            self._prepare_terminal_block("block")
            self.console.print(Text(text))
            self._terminal_last_kind = "block"

    def table(self, table: ReportTable, *, persist: bool = True) -> None:
        """Render one neutral table to terminal and plain-text report."""
        if persist:
            self._append_chunk("table", render_table(table).rstrip())
        if not self.silent:
            self._prepare_terminal_block("table")
            self.console.print(build_rich_renderable(table))
            self._terminal_last_kind = "table"

    def tables(
        self,
        tables: Sequence[ReportTable],
        *,
        persist: bool = True,
    ) -> None:
        """Render an ordered sequence of neutral tables."""
        for table in tables:
            self.table(table, persist=persist)

    def progress(self, event: Event) -> None:
        """Update the terminal-only live progress display."""
        if self.silent:
            return
        self.progress_display.update(event)

    def finish_progress(self) -> None:
        """Stop and clear the current live progress display."""
        self.progress_display.finish()

    def close(self) -> None:
        """Release live terminal resources without forcing completion."""
        self.progress_display.stop()

    def text(self) -> str:
        """Return the deterministic plain-text report collected so far."""
        if not self._chunks:
            return ""
        rendered = self._chunks[0][1]
        previous_kind = self._chunks[0][0]
        for kind, chunk in self._chunks[1:]:
            is_block_boundary = "table" in {previous_kind, kind} or "block" in {
                previous_kind,
                kind,
            }
            separator = "\n\n" if is_block_boundary else "\n"
            rendered += separator + chunk
            previous_kind = kind
        return rendered

    def save(self) -> None:
        """Stop progress and write the collected plain-text report."""
        self.close()
        if self.report_file is not None:
            self.report_file.write_text(self.text(), encoding="utf-8")

    def _append_chunk(self, kind: str, text: str) -> None:
        """Append one typed report block without terminal decoration."""
        if text:
            self._chunks.append((kind, text))

    def _prepare_terminal_block(self, kind: str) -> None:
        """Insert visual whitespace around tables and preformatted blocks."""
        if self._terminal_last_kind is None:
            return
        if kind in {"table", "block"} or self._terminal_last_kind in {"table", "block"}:
            self.console.print()


def print_terminal_message(
    message: str,
    *,
    level: EventLevel = EventLevel.INFO,
    bold: bool = False,
    silent: bool = False,
    console: Console | None = None,
) -> None:
    """Print one standalone terminal message using the shared style policy."""
    if silent:
        return
    target = create_console() if console is None else console
    target.print(_terminal_message(message, level, bold=bold))


def _plain_message(message: str, level: EventLevel) -> str:
    """Return the portable report representation of one message."""
    prefixes = {
        EventLevel.WARNING: "WARNING: ",
        EventLevel.ERROR: "ERROR: ",
        EventLevel.DEBUG: "DEBUG: ",
        EventLevel.PROGRESS: "Progress: ",
    }
    return prefixes.get(level, "") + message


def _terminal_message(message: str, level: EventLevel, *, bold: bool) -> Text:
    """Return a Rich text object following the Quantas terminal style policy."""
    style = ""
    prefix = ""
    if level == EventLevel.WARNING:
        style = "yellow"
        prefix = "WARNING: "
    elif level == EventLevel.ERROR:
        style = "red"
        prefix = "ERROR: "
    elif level == EventLevel.DEBUG:
        prefix = "DEBUG: "
    elif level == EventLevel.PROGRESS:
        prefix = "Progress: "
    if bold:
        style = f"bold {style}".strip()
    return Text(prefix + message, style=style)


__all__ = [
    "CLIOutput",
    "create_console",
    "print_terminal_message",
]
