"""Characterization tests for Rich terminal and plain-report separation."""

from __future__ import annotations

from io import StringIO
import sys

import pytest
from rich.console import Console

from quantas.cli.output import CLIOutput, create_console
from quantas.core.events import Event, EventLevel
from quantas.models import ReportTable
from quantas.renderers.tables import build_rich_table


pytestmark = pytest.mark.cli


def _console(stream: StringIO) -> Console:
    return Console(
        file=stream,
        force_terminal=True,
        color_system="standard",
        width=80,
        highlight=False,
        markup=False,
    )


def test_rich_table_uses_neutral_values_and_bold_headers() -> None:
    table = ReportTable(
        title="Elastic properties",
        columns=["Quantity", "Value"],
        rows=[["Young modulus", 123.456789]],
        metadata={
            "column_formats": [None, ".4f"],
            "column_units": ["", "GPa"],
            "column_alignments": ["left", "right"],
        },
    )
    rendered = build_rich_table(table)
    assert rendered.title.plain == "Elastic properties"
    assert rendered.columns[1].header == "Value\n(GPa)"
    assert rendered.columns[1].justify == "right"


def test_terminal_uses_rich_but_report_remains_plain_text(tmp_path) -> None:
    stream = StringIO()
    report = tmp_path / "report.txt"
    output = CLIOutput(
        report_file=report,
        console=_console(stream),
    )
    output.message("check this result", level=EventLevel.WARNING)
    output.message("calculation failed", level=EventLevel.ERROR)
    output.table(
        ReportTable(
            title="Energy",
            columns=["Quantity", "Value"],
            rows=[["E", -123.456789]],
            metadata={"column_formats": [None, "energy_ha"]},
        )
    )
    output.save()

    terminal = stream.getvalue()
    plain = report.read_text(encoding="utf-8")
    assert "\x1b[" in terminal
    assert "WARNING: check this result" in terminal
    assert "ERROR: calculation failed" in terminal
    assert "Energy" in terminal
    assert "WARNING: check this result" in plain
    assert "ERROR: calculation failed" in plain
    assert "-1.234567890000E+02" in plain
    assert "\x1b[" not in plain
    assert not {"┏", "┓", "┗", "┛", "━", "┃"}.intersection(plain)


def test_progress_messages_are_terminal_only(tmp_path) -> None:
    stream = StringIO()
    report = tmp_path / "report.txt"
    output = CLIOutput(report_file=report, console=_console(stream))
    output.message("grid (2/4, 50.0%)", level=EventLevel.PROGRESS, persist=False)
    output.message("completed")
    output.save()

    assert "Progress: grid (2/4, 50.0%)" in stream.getvalue()
    assert "Progress:" not in output.text()
    assert report.read_text(encoding="utf-8") == "completed"


class _TTYStringIO(StringIO):
    """String stream advertising terminal capabilities for renderer tests."""

    def isatty(self) -> bool:
        """Return ``True`` to emulate an interactive terminal."""
        return True


def test_default_console_uses_real_sys_stdout(monkeypatch) -> None:
    """The Rich adapter bypasses Click stream wrappers."""
    stream = StringIO()
    monkeypatch.setattr(sys, "stdout", stream)
    console = create_console()
    assert console.file is stream


def test_redirected_stream_disables_all_ansi_sequences() -> None:
    """Non-TTY output is plain even when Rich Text contains styles."""
    stream = StringIO()
    console = create_console(file=stream)
    console.print("plain redirected output", style="bold red")
    rendered = stream.getvalue()
    assert "plain redirected output" in rendered
    assert "\x1b[" not in rendered


def test_explicit_no_ansi_disables_styles_on_tty(monkeypatch) -> None:
    """QUANTAS_NO_ANSI provides a conservative terminal fallback."""
    stream = _TTYStringIO()
    monkeypatch.setenv("QUANTAS_NO_ANSI", "1")
    console = create_console(file=stream)
    console.print("plain terminal output", style="bold yellow")
    rendered = stream.getvalue()
    assert "plain terminal output" in rendered
    assert "\x1b[" not in rendered


def test_default_tty_console_emits_rich_styles(monkeypatch) -> None:
    """A capable TTY retains Rich styling and color support."""
    stream = _TTYStringIO()
    monkeypatch.delenv("QUANTAS_NO_ANSI", raising=False)
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")
    console = create_console(file=stream)
    console.print("styled terminal output", style="bold yellow")
    rendered = stream.getvalue()
    assert "styled terminal output" in rendered
    assert "\x1b[" in rendered


def test_progress_warning_do_not_pollute_report() -> None:
    """Warnings print above a live task while progress remains terminal-only."""
    stream = _TTYStringIO()
    console = Console(
        file=stream,
        force_terminal=True,
        color_system="standard",
        width=100,
        highlight=False,
        markup=False,
    )
    output = CLIOutput(console=console, show_progress=True)
    output.progress(
        Event(
            "grid",
            level=EventLevel.PROGRESS,
            progress=0.5,
            data={"label": "pressure-temperature grid", "current": 2, "total": 4},
        )
    )
    output.message("local minimization failed", level=EventLevel.WARNING)
    output.progress(
        Event(
            "grid",
            level=EventLevel.PROGRESS,
            progress=1.0,
            data={"label": "pressure-temperature grid", "current": 4, "total": 4},
        )
    )
    output.close()

    terminal = stream.getvalue()
    assert "pressure-temperature grid" in terminal
    assert "50%" in terminal
    assert "WARNING: local minimization failed" in terminal
    assert "Progress:" not in output.text()
    assert output.text() == "WARNING: local minimization failed"


def test_no_progress_disables_live_renderer_on_tty() -> None:
    """The debugging-oriented no-progress mode emits no live control output."""
    stream = _TTYStringIO()
    console = Console(
        file=stream,
        force_terminal=True,
        color_system="standard",
        width=80,
        highlight=False,
        markup=False,
    )
    output = CLIOutput(console=console, show_progress=False)
    output.progress(
        Event(
            "grid",
            level=EventLevel.PROGRESS,
            progress=0.5,
            data={"current": 1, "total": 2},
        )
    )
    output.close()
    assert stream.getvalue() == ""


def test_tables_are_separated_from_surrounding_messages() -> None:
    """Plain reports use one empty line around each table block."""
    output = CLIOutput(silent=True)
    output.message("Before table")
    output.table(
        ReportTable(
            title="Options",
            columns=["Option", "Value"],
            rows=[["Mode", "freq"]],
        )
    )
    output.message("After table")
    text = output.text()
    assert "Before table\n\nOptions" in text
    assert "freq\n\nAfter table" in text


def test_all_run_commands_enable_progress_by_default() -> None:
    """Every scientific run command exposes an opt-out live progress flag."""
    from click.testing import CliRunner

    from quantas.cli.elasticity import elasticity
    from quantas.cli.ha import ha
    from quantas.cli.qha import qha
    from quantas.cli.seismic import seismic

    runner = CliRunner()
    for command in (elasticity, ha, qha, seismic):
        result = runner.invoke(command, ["run", "--help"])
        assert result.exit_code == 0
        normalized_help = " ".join(result.output.split())
        assert "--progress / --no-progress" in normalized_help
        assert "[default: progress]" in normalized_help


def test_cli_output_context_manager_stops_live_progress() -> None:
    """Context exit releases a live Rich task even when an error occurs."""
    stream = _TTYStringIO()
    console = Console(
        file=stream,
        force_terminal=True,
        color_system="standard",
        width=80,
        highlight=False,
        markup=False,
    )
    output = CLIOutput(console=console, show_progress=True)

    with pytest.raises(RuntimeError):
        with output:
            output.progress(
                Event(
                    "grid",
                    level=EventLevel.PROGRESS,
                    progress=0.5,
                    data={"current": 1, "total": 2},
                )
            )
            raise RuntimeError("stop")

    assert output.progress_display.active is False
