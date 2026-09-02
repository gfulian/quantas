# -*- coding: utf-8 -*-

"""Terminal observer for HA/QHA phonon input generation."""

from __future__ import annotations

from rich.console import Console

from quantas.cli.output import CLIOutput
from quantas.core.events import Event, EventLevel
from quantas.models import ReportTable


class PhononInputTextObserver:
    """Render frontend-neutral phonon input-generation events.

    Input generation intentionally creates no persistent text report. The
    generated YAML file is the scientific artifact; this observer only renders
    operational diagnostics to the terminal.

    Parameters
    ----------
    silent : bool, optional
        Suppress terminal output when ``True``.
    debug : bool, optional
        Render detailed debug messages and mode-tracking tables when ``True``.
    console : rich.console.Console or None, optional
        Explicit Rich console, primarily useful for tests.
    """

    def __init__(
        self,
        *,
        silent: bool = False,
        debug: bool = False,
        console: Console | None = None,
    ) -> None:
        self.debug = bool(debug)
        self.output = CLIOutput(
            report_file=None,
            silent=silent,
            console=console,
            show_progress=False,
        )

    def __call__(self, event: Event) -> None:
        """Receive and render one input-generation event."""
        table = event.data.get("table")
        if event.level == EventLevel.DEBUG:
            if not self.debug:
                return
            if isinstance(table, ReportTable):
                self.output.table(table, persist=False)
            elif event.message:
                self.output.message(
                    event.message,
                    level=EventLevel.DEBUG,
                    persist=False,
                )
            return

        if event.level == EventLevel.RESULT:
            if isinstance(table, ReportTable):
                self.output.table(table, persist=False)
            elif event.message:
                self.output.message(event.message, persist=False)
            return

        if event.level in {EventLevel.INFO, EventLevel.WARNING, EventLevel.ERROR}:
            self.output.message(event.message, level=event.level, persist=False)

    def close(self) -> None:
        """Release terminal resources owned by the output service."""
        self.output.close()


__all__ = ["PhononInputTextObserver"]
