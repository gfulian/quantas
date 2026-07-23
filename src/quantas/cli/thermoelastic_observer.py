# -*- coding: utf-8 -*-

"""CLI observer for quasi-static thermoelastic workflows."""

from __future__ import annotations

from pathlib import Path

from rich.console import Console

from quantas.cli.contracts import ReportVerbosity, parse_verbosity
from quantas.cli.output import CLIOutput
from quantas.core.events import Event, EventLevel


class ThermoelasticTextObserver:
    """Render workflow events and collect a deterministic plain-text report."""

    def __init__(
        self,
        report_file: str | Path | None = None,
        *,
        silent: bool = False,
        show_progress: bool = True,
        verbosity: str | ReportVerbosity = ReportVerbosity.STANDARD,
        console: Console | None = None,
    ) -> None:
        self.show_progress = bool(show_progress)
        self.verbosity = parse_verbosity(verbosity)
        self.output = CLIOutput(
            report_file=report_file,
            silent=silent,
            console=console,
            show_progress=show_progress,
        )

    def __call__(self, event: Event) -> None:
        """Render one Quantas event."""
        if event.level == EventLevel.PROGRESS:
            if self.show_progress:
                self.output.progress(event)
            return
        if event.level == EventLevel.DEBUG and not self.verbosity.includes_debug:
            return
        if event.level == EventLevel.RESULT:
            if self.verbosity.includes_debug:
                self.output.message(event.message, level=EventLevel.DEBUG)
            return
        self.output.message(event.message, level=event.level)

    def save(self) -> None:
        """Write the configured report file."""
        self.output.save()

    def close(self) -> None:
        """Release terminal progress resources."""
        self.output.close()

    def text(self) -> str:
        """Return the collected plain-text report."""
        return self.output.text()


__all__ = ["ThermoelasticTextObserver"]
