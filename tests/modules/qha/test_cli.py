from __future__ import annotations

import numpy as np

from click.testing import CliRunner

from quantas.core.events import Event, EventLevel
from quantas.cli.qha import qha
from quantas.cli.qha_observer import QHATextObserver
from quantas.cli.qha_render import (
    render_preview_report,
    render_qha_cli_report,
    render_table,
)
from quantas.core.math.fitting import FitQuality, FitResult, FitStatus
from quantas.modules.qha.inspect import PressureEstimate, PressureVolumePreview
from quantas.modules.qha.models import QHAInput, QHAOptions, QHAResult
from quantas.modules.qha.report import input_table


def _input() -> QHAInput:
    return QHAInput(
        jobname="MgO",
        natoms=2,
        qpoints=1,
        volume=np.array([18.0, 19.0, 20.0]),
        energy=np.array([-10.0, -10.2, -10.0]),
        frequencies=np.ones((1, 6, 3)),
        weights=np.array([1.0]),
    )


def _result() -> QHAResult:
    return QHAResult(
        jobname="MgO",
        temperature=np.array([300.0]),
        pressure=np.array([0.0]),
        volume=np.array([18.0, 19.0, 20.0]),
        static_energy=np.array([-10.0, -10.2, -10.0]),
        equilibrium_volume=np.array([[19.0]]),
        isothermal_bulk_modulus=np.array([[160.0]]),
        metadata={"workflow": {"completed": True, "failed_points": 0, "warnings": 0}},
    )


def test_render_table_formats_neutral_table() -> None:
    text = render_table(input_table(_input()))

    assert "Input data" in text
    assert "Job name" in text
    assert "MgO" in text


def test_render_qha_cli_report_contains_standard_sections() -> None:
    text = render_qha_cli_report(_input(), QHAOptions(), _result())

    assert "Input data" in text
    assert "Selected options" in text
    assert "QHA properties" in text


def test_render_preview_report_includes_warnings_and_diagnostics() -> None:
    fit = FitResult(
        success=True,
        status=FitStatus.SUCCESS,
        quality=FitQuality.GOOD,
        parameters=np.array([1.0, 2.0]),
        rmse=0.01,
    )
    estimate = PressureEstimate(
        method="poly",
        pressure=np.array([1.0, 0.0]),
        fit=fit,
        unit="GPa",
        warnings=["near boundary"],
    )
    preview = PressureVolumePreview(
        volume=np.array([18.0, 19.0]),
        energy=np.array([-10.0, -10.2]),
        pressure_unit="GPa",
        polynomial=estimate,
        eos=estimate,
        warnings=["global warning"],
    )

    text = render_preview_report(preview)

    assert "Pressure-volume preview" in text
    assert "Pressure-volume fit diagnostics" in text
    assert "WARNING: global warning" in text


def test_qha_text_observer_collects_result_and_debug_events() -> None:
    observer = QHATextObserver(silent=True, verbosity="debug")
    options = QHAOptions(debug=True)

    observer(
        Event(
            "settings",
            level=EventLevel.RESULT,
            data={"kind": "settings", "options": options},
        )
    )
    observer(
        Event(
            "fit",
            level=EventLevel.DEBUG,
            data={
                "kind": "qha_fit_record",
                "quantity": "VT",
                "method": "poly",
                "temperature": 300.0,
                "pressure": 0.0,
                "success": True,
                "fit": {"status": "success", "quality": "good", "rmse": 0.0},
            },
        )
    )

    text = observer.text()
    assert "Selected options" in text
    assert "Minimize volume at P = 0.0 and T = 300.0" in text


def test_qha_group_exposes_expected_commands() -> None:
    runner = CliRunner()
    result = runner.invoke(qha, ["--help"])

    assert result.exit_code == 0
    assert "inspect" in result.output
    assert "run" in result.output
    assert "export" in result.output


def test_qha_export_help_documents_structural_selection() -> None:
    result = CliRunner().invoke(qha, ["export", "--help"])

    assert result.exit_code == 0
    assert "'structure'" in result.output
    assert "anisotropic thermal expansion" in result.output


def test_qha_progress_is_not_persisted_in_plain_report() -> None:
    observer = QHATextObserver(silent=True, show_progress=True)
    observer(
        Event(
            "pressure-temperature point completed",
            level=EventLevel.PROGRESS,
            progress=0.5,
            data={
                "kind": "qha_point_completed",
                "label": "pressure-temperature grid",
                "current": 2,
                "total": 4,
            },
        )
    )

    assert "2/4, 50.0%" not in observer.text()


def test_qha_text_observer_renders_warning_state_and_units() -> None:
    observer = QHATextObserver(silent=True)
    observer(
        Event(
            "the free-energy EOS fit failed",
            level=EventLevel.WARNING,
            data={
                "kind": "qha_point_failed",
                "pressure": 2.0,
                "temperature": 1940.0,
                "pressure_unit": "GPa",
                "temperature_unit": "K",
            },
        )
    )

    assert observer.text() == (
        "WARNING: the free-energy EOS fit failed at pressure-temperature "
        "conditions: 2 GPa - 1940 K"
    )


def test_qha_text_observer_avoids_repeating_state_wording() -> None:
    observer = QHATextObserver(silent=True)
    observer(
        Event(
            "the pressure-temperature state produced diagnostic warnings",
            level=EventLevel.WARNING,
            data={
                "pressure": 2.0,
                "temperature": 1940.0,
                "pressure_unit": "GPa",
                "temperature_unit": "K",
            },
        )
    )

    assert observer.text() == (
        "WARNING: the pressure-temperature state produced diagnostic warnings: "
        "2 GPa - 1940 K"
    )


def test_qha_run_help_exposes_polynomial_derivative_options() -> None:
    runner = CliRunner()
    result = runner.invoke(qha, ["run", "--help"])

    assert result.exit_code == 0
    assert "--poly-derivatives" in result.output
    assert "--poly-grid-points" in result.output
    assert "--poly-grid-separation" in result.output
    assert "--mode-continuity" in result.output
    assert "--mode-gruneisen" in result.output
    assert "--gruneisen-min-cv-fraction" in result.output


def test_qha_inspect_failure_reports_original_error_without_name_error(
    tmp_path, monkeypatch
) -> None:
    """Input inspection failures must not reference an uninitialized observer."""
    filename = tmp_path / "broken.yaml"
    filename.write_text("job: broken\n", encoding="utf-8")

    def fail_inspection(*args, **kwargs):
        raise ValueError("invalid QHA input")

    monkeypatch.setattr("quantas.cli.qha.inspect_qha_input", fail_inspection)
    result = CliRunner().invoke(qha, ["inspect", str(filename)])

    assert result.exit_code != 0
    assert "invalid QHA input" in result.output
    assert "observer" not in result.output.lower()


def test_qha_run_closes_observer_when_input_override_fails(
    tmp_path, monkeypatch
) -> None:
    """Any pre-run QHA failure must release terminal progress resources."""
    filename = tmp_path / "broken.yaml"
    filename.write_text("job: broken\n", encoding="utf-8")
    state = {"closed": False}

    class FakeObserver:
        def __init__(self, *args, **kwargs) -> None:
            pass

        def __call__(self, event) -> None:
            pass

        def close(self) -> None:
            state["closed"] = True

        def save(self) -> None:
            raise AssertionError("save must not run after a failed calculation")

        def text(self) -> str:
            return ""

    def fail_read(*args, **kwargs):
        raise ValueError("cannot override mode continuity")

    monkeypatch.setattr("quantas.cli.qha.QHATextObserver", FakeObserver)
    monkeypatch.setattr("quantas.cli.qha.read_qha_input", fail_read)
    result = CliRunner().invoke(
        qha,
        ["run", str(filename), "--mode-continuity", "verified"],
    )

    assert result.exit_code != 0
    assert "cannot override mode continuity" in result.output
    assert state["closed"] is True


def test_qha_inpgen_help_exposes_shared_generator_options() -> None:
    """The QHA group exposes the shared HA/QHA input generator."""
    result = CliRunner().invoke(qha, ["inpgen", "--help"])

    assert result.exit_code == 0
    assert "--interface" in result.output
    assert "crystal-qha" in result.output
    assert "--list" in result.output
    assert "--reference" in result.output
    assert "--formula-units" in result.output


def test_qha_plot_help_exposes_exact_grid_sections() -> None:
    result = CliRunner().invoke(qha, ["plot", "--help"])

    assert result.exit_code == 0
    assert "--axis" in result.output
    assert "--pressure" in result.output
    assert "--temperature" in result.output
    assert "Exact native" in result.output


def test_qha_plot_forwards_exact_grid_sections_to_public_api(
    tmp_path, monkeypatch
) -> None:
    from types import SimpleNamespace

    filename = tmp_path / "qha.h5"
    filename.write_bytes(b"test")
    captured = {}

    monkeypatch.setattr("quantas.cli.qha.read_qha_hdf5", lambda path: object())

    def fake_build(result, *, properties, options):
        captured["properties"] = properties
        captured["options"] = options
        return object()

    monkeypatch.setattr("quantas.cli.qha.build_qha_plots", fake_build)
    monkeypatch.setattr(
        "quantas.cli.qha.render_plot_collection",
        lambda collection, options: SimpleNamespace(warnings=[], artifacts=[]),
    )

    result = CliRunner().invoke(
        qha,
        [
            "plot",
            str(filename),
            "--property",
            "KT",
            "--axis",
            "pressure",
            "--temperature",
            "300",
            "--temperature",
            "600",
            "--2d",
        ],
    )

    assert result.exit_code == 0
    assert captured["properties"] == ("KT",)
    assert captured["options"].curve_axis == "pressure"
    assert captured["options"].selected_temperatures == (300.0, 600.0)
    assert captured["options"].include_contours is True
