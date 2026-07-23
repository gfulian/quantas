"""Tests for the deterministic staged test-suite runner."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pytest


pytestmark = pytest.mark.infrastructure


def _load_runner() -> ModuleType:
    """Load ``tools/run_tests.py`` without requiring ``tools`` as a package."""
    path = Path(__file__).resolve().parents[2] / "tools" / "run_tests.py"
    spec = importlib.util.spec_from_file_location("quantas_test_runner", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_all_target_expands_to_isolated_stages() -> None:
    runner = _load_runner()
    assert runner._selected_suites("all") == (
        "core",
        "elasticity",
        "seismic",
        "ha",
        "qha",
        "eos",
        "thermoelasticity",
        "cli",
        "elasticity-plotting",
        "seismic-plotting",
        "ha-plotting",
        "qha-plotting",
        "eos-plotting",
        "thermoelasticity-plotting",
        "examples",
    )
    assert runner._selected_suites("library") == runner.LIBRARY_STAGES
    assert runner._selected_suites("qha") == ("qha",)


def test_staged_runner_executes_every_stage_and_forwards_pytest_args(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = _load_runner()
    commands: list[list[str]] = []

    def fake_run(command: list[str], **kwargs: object) -> int:
        commands.append(command)
        environment = kwargs["environment"]
        assert isinstance(environment, dict)
        assert environment["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] == "1"
        assert environment["MPLBACKEND"] == "Agg"
        assert environment["OPENBLAS_NUM_THREADS"] == "1"
        assert environment["OMP_NUM_THREADS"] == "1"
        assert "MPLCONFIGDIR" in environment
        assert Path(environment["MPLCONFIGDIR"]).is_dir()
        assert kwargs["timeout"] == 600.0
        return 0

    monkeypatch.setattr(runner, "_spawn_stage", fake_run)
    status = runner.main(["all", "-v", "--", "-q", "--tb=short"])

    assert status == 0
    assert [command[4] for command in commands] == [
        "core",
        "elasticity",
        "seismic",
        "ha",
        "qha",
        "eos",
        "thermoelasticity",
        "cli",
        "elasticity-plotting",
        "seismic-plotting",
        "ha-plotting",
        "qha-plotting",
        "eos-plotting",
        "thermoelasticity-plotting",
        "examples",
    ]
    assert all("-v" in command for command in commands)
    assert all("-q" in command and "--tb=short" in command for command in commands)


def test_staged_runner_returns_first_failure_but_finishes_by_default(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = _load_runner()
    statuses = iter([0, 3, *([0] * 13)])
    calls: list[str] = []

    def fake_run(command: list[str], **_: object) -> int:
        calls.append(command[4])
        return next(statuses)

    monkeypatch.setattr(runner, "_spawn_stage", fake_run)
    status = runner.main(["all"])

    assert status == 3
    assert calls == [
        "core",
        "elasticity",
        "seismic",
        "ha",
        "qha",
        "eos",
        "thermoelasticity",
        "cli",
        "elasticity-plotting",
        "seismic-plotting",
        "ha-plotting",
        "qha-plotting",
        "eos-plotting",
        "thermoelasticity-plotting",
        "examples",
    ]


def test_fail_fast_stops_after_first_failed_stage(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = _load_runner()
    statuses = iter([0, 2])
    calls: list[str] = []

    def fake_run(command: list[str], **_: object) -> int:
        calls.append(command[4])
        return next(statuses)

    monkeypatch.setattr(runner, "_spawn_stage", fake_run)
    status = runner.main(["all", "--fail-fast"])

    assert status == 2
    assert calls == ["core", "elasticity"]


def test_stage_timeout_returns_shell_timeout_status(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = _load_runner()

    def fake_run(*args: object, **kwargs: object) -> int:
        return 124

    monkeypatch.setattr(runner, "_spawn_stage", fake_run)
    status = runner.main(["core", "--stage-timeout", "1"])

    assert status == 124
