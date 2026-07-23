# -*- coding: utf-8 -*-

"""Observers used by the quasi-harmonic approximation command-line interface."""

from __future__ import annotations

from pathlib import Path

from rich.console import Console

from quantas.cli.contracts import ReportVerbosity, parse_verbosity
from quantas.cli.output import CLIOutput
from quantas.cli.qha_render import (
    final_qha_result_tables,
    phonon_frequency_fit_report_tables,
    preview_report_tables,
    thermodynamic_fit_report_tables,
)
from quantas.core.events import Event, EventLevel
from quantas.modules.qha.models import QHAInput, QHAResult
from quantas.modules.qha.report import input_table, options_table, result_summary_table


class QHATextObserver:
    """Render QHA events with Rich and collect a plain-text report.

    Live progress is terminal-only and is never persisted in the report or
    embedded in HDF5 workflow logs.
    """

    def __init__(
        self,
        report_file: str | Path | None = None,
        silent: bool = False,
        show_progress: bool = True,
        verbosity: str | ReportVerbosity = ReportVerbosity.STANDARD,
        max_debug_records: int | None = None,
        console: Console | None = None,
    ) -> None:
        self.show_progress = show_progress
        self.verbosity = parse_verbosity(verbosity)
        self.max_debug_records = max_debug_records
        self._debug_records = 0
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
        """Receive and render one QHA event."""
        if event.level == EventLevel.INFO:
            self.output.message(event.message)
        elif event.level == EventLevel.WARNING:
            self.output.message(
                self._state_aware_message(event),
                level=EventLevel.WARNING,
            )
            if self.show_progress and event.progress is not None:
                self.output.progress(event)
        elif event.level == EventLevel.ERROR:
            self.output.message(
                self._state_aware_message(event),
                level=EventLevel.ERROR,
            )
        elif event.level == EventLevel.PROGRESS:
            if (
                self.verbosity.includes_debug
                and event.data.get("kind") == "qha_point_completed"
            ):
                text = self._render_point_completed(event)
                if text:
                    self.output.text_block(text)
            elif self.show_progress:
                self.output.progress(event)
        elif event.level == EventLevel.DEBUG and self.verbosity.includes_debug:
            text = self._render_debug(event)
            if text:
                self.output.text_block(text)
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

    @staticmethod
    def _state_aware_message(event: Event) -> str:
        """Return an event message with pressure-temperature context."""
        pressure = event.data.get("pressure")
        temperature = event.data.get("temperature")
        pressure_unit = event.data.get("pressure_unit", "GPa")
        temperature_unit = event.data.get("temperature_unit", "K")

        if pressure is not None and temperature is not None:
            if "pressure-temperature state" in event.message:
                return (
                    f"{event.message}: {float(pressure):g} {pressure_unit} - "
                    f"{float(temperature):g} {temperature_unit}"
                )
            return (
                f"{event.message} at pressure-temperature conditions: "
                f"{float(pressure):g} {pressure_unit} - "
                f"{float(temperature):g} {temperature_unit}"
            )
        if pressure is not None:
            return f"{event.message} at {float(pressure):g} {pressure_unit}"
        if temperature is not None:
            return f"{event.message} at {float(temperature):g} {temperature_unit}"
        return event.message

    def _handle_result_event(self, event: Event) -> None:
        """Render a structured QHA result event."""
        kind = event.data.get("kind")
        if kind == "settings":
            self.output.table(options_table(event.data["options"]))
        elif kind == "input":
            input_data = event.data.get("input")
            if isinstance(input_data, QHAInput):
                self.output.table(input_table(input_data))
        elif kind == "input_summary":
            self.output.message(
                "Input summary: "
                f"{event.data.get('nvol', '?')} volumes, "
                f"{event.data.get('qpoints', '?')} q-points, "
                f"{event.data.get('nmodes', '?')} modes, "
                f"mode continuity={event.data.get('mode_continuity', 'unknown')}"
            )
        elif kind == "static_data" and self.verbosity.includes_debug:
            self.output.message(
                "Static data: "
                f"volume unit={event.data.get('volume_unit', '-')}, "
                f"energy unit={event.data.get('energy_unit', '-')}"
            )
        elif kind == "pressure_volume_preview":
            preview = event.data["preview"]
            self.output.tables(preview_report_tables(preview, include_diagnostics=True))
            for warning in preview.warnings:
                self.output.message(warning, level=EventLevel.WARNING)
        elif kind == "frequency_fit_report":
            self.output.tables(
                phonon_frequency_fit_report_tables(
                    event.data["input"],
                    event.data["options"],
                    include_debug=self.verbosity.includes_debug,
                )
            )
        elif kind == "thermodynamic_fit_report":
            self.output.tables(
                thermodynamic_fit_report_tables(
                    event.data["sampled_thermodynamics"],
                    event.data["options"],
                    include_debug=self.verbosity.includes_debug,
                )
            )
        elif kind == "gruneisen_summary":
            metadata = event.data.get("metadata", {})
            fit_stats = metadata.get("mode_fit_r_squared", {})
            self.output.message(
                "Mode Gruneisen summary: "
                f"continuity={metadata.get('mode_continuity', 'unknown')}, "
                f"usable modes={metadata.get('n_usable_modes', '?')}/"
                f"{metadata.get('n_modes', '?')}, "
                f"minimum R^2={fit_stats.get('minimum', '-')}"
            )
        elif kind == "workflow_summary":
            workflow = event.data.get("workflow", {})
            self.output.message(
                "Workflow summary: "
                f"completed={workflow.get('completed', True)}, "
                f"failed={workflow.get('failed_points', 0)}, "
                f"warnings={workflow.get('warnings', 0)}"
            )
        elif kind == "final_results":
            self.output.tables(
                final_qha_result_tables(
                    event.data["result"],
                    include_debug=self.verbosity.includes_debug,
                    pressure_indices=(
                        None if self.verbosity.includes_extended else (0, -1)
                    ),
                )
            )
        elif kind == "qha_point_completed" and self.verbosity.includes_debug:
            text = self._render_point_completed(event)
            if text:
                self.output.text_block(text)
        elif kind == "qha_completed":
            self.output.message("QHA workflow completed")
        elif kind == "thermodynamics":
            result = event.data["result"]
            if isinstance(result, QHAResult):
                self.output.table(result_summary_table(result))

    def _render_debug(self, event: Event) -> str:
        """Render selected debug events."""
        if event.data.get("kind") in {"qha_point_started", "qha_started"}:
            return ""
        if event.data.get("kind") == "qha_fit_record":
            if (
                self.max_debug_records is not None
                and self._debug_records >= self.max_debug_records
            ):
                return ""
            self._debug_records += 1
            return self._render_fit_record(event)
        return event.message

    def _render_fit_record(self, event: Event) -> str:
        """Render a local fit diagnostic event."""
        fit = event.data.get("fit", {})
        metadata = event.data.get("metadata", {})
        if not isinstance(metadata, dict):
            metadata = {}
        if not isinstance(fit, dict):
            fit = {}
        volume = metadata.get("volume", event.data.get("volume"))
        volume_sigma = metadata.get("volume_sigma")
        kt = metadata.get("bulk_modulus")
        kp = metadata.get("bulk_modulus_derivative")
        range_status = metadata.get("range_status", "-")
        lines = [
            (
                "Minimize volume at "
                f"P = {event.data.get('pressure', '-')} and "
                f"T = {event.data.get('temperature', '-')}"
            ),
            f"  calculated volume: {_format_number(volume)} +/- "
            f"{_format_number(volume_sigma)}",
            "  fitting parameters:",
        ]
        params = fit.get("parameters")
        errors = fit.get("errors")
        names = (
            fit.get("metadata", {}).get("parameter_order")
            if isinstance(fit.get("metadata"), dict)
            else None
        )
        if not isinstance(names, list):
            names = _default_parameter_names(event.data.get("method", "-"), params)
        if isinstance(params, list):
            for index, value in enumerate(params):
                error = None
                if isinstance(errors, list) and index < len(errors):
                    error = errors[index]
                name = names[index] if index < len(names) else f"p{index}"
                lines.append(
                    f"    {name} = {_format_number(value)} +/- {_format_number(error)}"
                )
        else:
            lines.append("    -")
        lines.extend(
            [
                "  diagnostics:",
                f"    success = {event.data.get('success', '-')}",
                f"    status = {fit.get('status', '-')}",
                f"    quality = {fit.get('quality', '-')}",
                f"    RMSE = {_format_number(fit.get('rmse'))}",
                f"    R^2 = {_format_number(fit.get('r_squared'))}",
                f"    range = {range_status}",
                f"    KT = {_format_number(kt)}",
                f"    Kp = {_format_number(kp)}",
                "",
            ]
        )
        return "\n".join(lines)

    @staticmethod
    def _render_point_completed(event: Event) -> str:
        """Render a completed pressure-temperature point when requested."""
        return ""


def _format_number(value: object) -> str:
    """Return a compact string for a diagnostic value."""
    if value is None:
        return "-"
    try:
        return f"{float(str(value)):.12e}"
    except (TypeError, ValueError):
        return str(value)


def _default_parameter_names(method: object, params: object) -> list[str]:
    """Return default parameter names for a fit record."""
    nparams = len(params) if isinstance(params, list) else 0
    if method == "eos" and nparams == 4:
        return ["E0", "K0", "KP", "V0"]
    if method in {"poly", "polynomial"}:
        return [f"c{i}" for i in range(nparams)]
    return [f"p{i}" for i in range(nparams)]
