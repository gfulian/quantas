# -*- coding: utf-8 -*-

"""Observers used by the harmonic-approximation command-line interface."""

from __future__ import annotations

from pathlib import Path

from rich.console import Console

from quantas.cli.contracts import ReportVerbosity, parse_verbosity
from quantas.cli.output import CLIOutput
from quantas.core.events import Event, EventLevel
from quantas.modules.ha.models import HAInput, HAResult
from quantas.modules.ha.report import (
    input_table,
    options_table,
    static_data_table,
    thermodynamic_property_tables,
    thermodynamic_summary_table,
)


class HATextObserver:
    """Render HA events with Rich and collect a plain-text report.

    Live progress is terminal-only and is never written to report files or
    embedded in HDF5 workflow logs.
    """

    def __init__(
        self,
        report_file: str | Path | None = None,
        silent: bool = False,
        show_progress: bool = True,
        verbosity: str | ReportVerbosity = ReportVerbosity.STANDARD,
        include_timing: bool = False,
        console: Console | None = None,
    ) -> None:
        self.show_progress = show_progress
        self.verbosity = parse_verbosity(verbosity)
        self.include_timing = bool(include_timing)
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
        """Receive and render one HA event."""
        if event.level in {EventLevel.INFO, EventLevel.WARNING, EventLevel.ERROR}:
            self.output.message(event.message, level=event.level)
        elif event.level == EventLevel.PROGRESS and self.show_progress:
            self.output.progress(event)
        elif event.level == EventLevel.DEBUG and (
            self.verbosity.includes_debug or self.include_timing
        ):
            self.output.message(self._render_debug(event), level=EventLevel.DEBUG)
        elif event.level == EventLevel.RESULT:
            self._handle_result_event(event)

    def save(self) -> None:
        """Save the collected plain-text report."""
        self.output.save()

    def close(self) -> None:
        """Release terminal progress resources."""
        self.output.close()

    def text(self) -> str:
        """Return the collected plain-text report."""
        return self.output.text()

    def _handle_result_event(self, event: Event) -> None:
        """Render a structured HA result event."""
        kind = event.data.get("kind")
        if kind == "settings":
            self.output.table(options_table(event.data["options"]))
        elif kind in {"input", "input_summary"}:
            input_data = event.data.get("input")
            if isinstance(input_data, HAInput):
                self.output.table(input_table(input_data))
            else:
                self.output.message(
                    "Input summary: "
                    f"{event.data.get('natoms', '?')} atoms, "
                    f"{event.data.get('qpoints', '?')} q-points, "
                    f"{event.data.get('nvol', '?')} volumes"
                )
        elif kind == "thermodynamics":
            result = event.data["result"]
            if isinstance(result, HAResult):
                self.output.table(static_data_table(result))
                self.output.table(thermodynamic_summary_table(result))
                row_indices = None if self.verbosity.includes_extended else (0, -1)
                for table in thermodynamic_property_tables(
                    result,
                    row_indices=row_indices,
                    include_zero_point=False,
                ):
                    self.output.table(table)
        elif kind == "thermodynamic_backend" and (
            self.verbosity.includes_debug or self.include_timing
        ):
            backend = event.data.get("backend", "unknown")
            self.output.message(
                f"Thermodynamic backend: {backend}",
                level=EventLevel.DEBUG,
            )

    @staticmethod
    def _render_debug(event: Event) -> str:
        """Render selected debug events."""
        if event.data.get("kind") == "timing":
            label = event.data.get("label", event.message)
            elapsed = float(event.data.get("elapsed_seconds", 0.0))
            backend = event.data.get("backend", "unknown")
            return f"Timing: {label} [{backend}] = {elapsed:.6f} s"
        return event.message
