"""CLI contract tests for the fixed double-precision policy."""

from __future__ import annotations

from click.testing import CliRunner

from quantas.cli.main import main


_RUN_COMMANDS = (
    ("elasticity",),
    ("seismic",),
    ("ha",),
    ("qha",),
)


def test_scientific_runs_do_not_expose_precision_switches() -> None:
    runner = CliRunner()
    for command in _RUN_COMMANDS:
        response = runner.invoke(main, [*command, "run", "--help"])
        assert response.exit_code == 0, response.output
        assert "--storage-precision" not in response.output
        assert "--working-precision" not in response.output
