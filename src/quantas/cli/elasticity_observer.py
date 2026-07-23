# -*- coding: utf-8 -*-

"""Command-line observer for elasticity workflow events."""

from __future__ import annotations

from pathlib import Path

from rich.console import Console

from quantas.cli.contracts import ReportVerbosity, parse_verbosity
from quantas.cli.output import CLIOutput
from quantas.core.events import Event, EventLevel
from quantas.modules.elasticity.report import (
    averages_table,
    compliance_table,
    input_table,
    options_table,
    stability_table,
    source_stiffness_table,
    stiffness_table,
    tensor_rotation_table,
    variations_table,
)


class ElasticityTextObserver:
    """Render elasticity events with Rich and collect a plain-text report.

    Live progress is terminal-only and never enters the persistent report.
    """

    def __init__(
        self,
        report_file: str | Path | None = None,
        silent: bool = False,
        show_progress: bool = True,
        verbosity: str | ReportVerbosity = ReportVerbosity.STANDARD,
        console: Console | None = None,
    ) -> None:
        self.verbosity = parse_verbosity(verbosity)
        self.output = CLIOutput(
            report_file=report_file,
            silent=silent,
            console=console,
            show_progress=show_progress,
        )

    @property
    def report_file(self) -> Path | None:
        """Return the configured report file."""
        return self.output.report_file

    @property
    def silent(self) -> bool:
        """Return whether terminal output is disabled."""
        return self.output.silent

    def __call__(self, event: Event) -> None:
        """Receive and render one Quantas event."""
        if event.level in {EventLevel.INFO, EventLevel.WARNING, EventLevel.ERROR}:
            self.output.message(event.message, level=event.level)
        elif event.level == EventLevel.PROGRESS:
            self.output.progress(event)
        elif event.level == EventLevel.DEBUG and self.verbosity.includes_debug:
            self.output.message(event.message, level=EventLevel.DEBUG)
        elif event.level == EventLevel.RESULT:
            self._handle_result_event(event)

    def save(self) -> None:
        """Save the deterministic plain-text report."""
        self.output.save()

    def close(self) -> None:
        """Release terminal progress resources."""
        self.output.close()

    def text(self) -> str:
        """Return the deterministic plain-text report."""
        return self.output.text()

    def _handle_result_event(self, event: Event) -> None:
        """Render a structured elasticity result event."""
        kind = event.data.get("kind")
        if kind == "settings":
            self.output.table(options_table(event.data["options"]))
        elif kind == "rotation":
            self.output.table(source_stiffness_table(event.data["source_stiffness"]))
            self.output.table(tensor_rotation_table(event.data["rotation"]))
        elif kind == "input":
            input_data = event.data["input"]
            result = event.data["result"]
            self.output.table(input_table(input_data, result))
            self.output.table(stiffness_table(result))
            if self.verbosity.includes_extended:
                self.output.table(compliance_table(result))
        elif kind == "averages":
            self.output.table(averages_table(event.data["result"]))
        elif kind == "stability":
            self.output.table(stability_table(event.data["result"]))
        elif kind == "variations":
            self.output.table(variations_table(event.data["result"]))
