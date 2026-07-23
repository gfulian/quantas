from __future__ import annotations

import numpy as np
import pytest

from quantas.modules.qha.models import QHAInput, QHAOptions
from quantas.modules.qha.workflow import (
    QHAWorkflowEvent,
    QHAWorkflowState,
    run_volume_minimization_workflow,
)


def make_input() -> QHAInput:
    volume = np.linspace(9.0, 11.0, 7)
    energy = 0.5 * (volume - 10.0) ** 2 - 5.0
    return QHAInput(jobname="workflow-test", volume=volume, energy=energy)


def test_workflow_state_reports_progress() -> None:
    state = QHAWorkflowState(total_points=4, processed_points=1)

    assert state.progress == pytest.approx(0.25)
    assert state.as_dict()["total_points"] == 4


def test_workflow_runs_volume_minimization_grid() -> None:
    options = QHAOptions(
        temperature_min=300.0,
        temperature_max=400.0,
        temperature_step=100.0,
        pressure_min=0.0,
        pressure_max=1.0,
        pressure_step=1.0,
        minimization="poly",
        free_energy_degree=2,
    )

    result = run_volume_minimization_workflow(make_input(), options)

    assert result.completed is True
    assert result.valid_mask.shape == (2, 2)
    assert result.valid_mask.all()
    assert result.metadata["workflow"]["total_points"] == 4
    assert result.metadata["workflow"]["successful_points"] == 4
    assert result.equilibrium_volume[0, 0] == pytest.approx(10.0, abs=1.0e-8)


def test_workflow_emits_structured_events() -> None:
    events: list[QHAWorkflowEvent] = []
    options = QHAOptions(
        temperature_min=300.0,
        temperature_max=300.0,
        pressure_min=0.0,
        pressure_max=0.0,
        minimization="poly",
        free_energy_degree=2,
        debug=True,
    )

    result = run_volume_minimization_workflow(
        make_input(),
        options,
        callback=events.append,
    )

    assert result.completed is True
    assert [event.kind for event in events] == [
        "started",
        "fit_record",
        "point_started",
        "fit_record",
        "point_completed",
        "completed",
    ]
    assert events[-1].progress == pytest.approx(1.0)
    assert events[1].data["quantity"] == "F"
    assert events[1].pressure is None
    assert events[3].data["quantity"] == "KT/Kp"


def test_workflow_stops_after_consecutive_failures() -> None:
    input_data = QHAInput(
        jobname="bad-workflow",
        volume=np.array([1.0, 2.0, 3.0]),
        energy=np.array([1.0, np.nan, 3.0]),
    )
    options = QHAOptions(
        temperature_min=300.0,
        temperature_max=301.0,
        temperature_step=1.0,
        pressure_min=0.0,
        pressure_max=2.0,
        pressure_step=1.0,
        minimization="poly",
        free_energy_degree=2,
        fit_failure_policy="stop",
        max_consecutive_failures=2,
        debug=True,
    )

    result = run_volume_minimization_workflow(input_data, options)

    assert result.completed is False
    assert len(result.failed_points) == 6
    assert result.metadata["workflow"]["processed_points"] == 6
    assert result.metadata["workflow"]["stopped"] is True
    assert result.metadata["stop"]["max_consecutive_failures"] == 2
    assert result.metadata["stop"]["failure_scope"] == "temperature"


def test_workflow_failure_policy_raise_propagates_runtime_error() -> None:
    input_data = QHAInput(
        volume=np.array([1.0, 2.0, 3.0]),
        energy=np.array([1.0, np.nan, 3.0]),
    )
    options = QHAOptions(
        minimization="poly",
        fit_failure_policy="raise",
        free_energy_degree=2,
    )

    with pytest.raises(RuntimeError):
        run_volume_minimization_workflow(input_data, options)
