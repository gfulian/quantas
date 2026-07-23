"""Tests for the HA calculator."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.events import EventLevel, ListObserver
from quantas.models import ResultData
from quantas.modules.ha.calculator import HACalculator
from quantas.modules.ha.io.reader import HAInputFileReader
from quantas.modules.ha.models import HAInput, HAOptions

DATA = Path(__file__).parent / "data/mgo_b3lyp_qha.yaml"


def _minimal_input() -> HAInput:
    frequencies = np.ones((1, 6, 1), dtype=np.float64) * 100.0
    return HAInput(
        jobname="Minimal HA",
        natoms=2,
        supercell=np.eye(3),
        qpoints=1,
        volume=np.array([18.0], dtype=np.float64),
        energy=np.array([-100.0], dtype=np.float64),
        frequencies=frequencies,
        weights=np.array([1.0], dtype=np.float64),
    )


def test_ha_calculator_returns_result_data() -> None:
    """HACalculator returns a generic Quantas result object."""
    calculator = HACalculator(
        _minimal_input(),
        options=HAOptions(
            temperature_min=0.0,
            temperature_max=100.0,
            temperature_step=100.0,
        ),
    )

    result = calculator.execute()

    assert isinstance(result, ResultData)
    assert result.metadata.module == "ha"
    assert result.metadata.method == "harmonic"
    assert result.metadata.schema_version == "2.0"
    assert result.results["ha"].jobname == "Minimal HA"
    assert result.results["ha"].temperature.shape == (2,)
    assert result.results["ha"].free_energy.shape == (2, 1)
    assert result.results["ha"].has_thermodynamic_data()
    assert (
        result.results["ha"].metadata["thermodynamic_unit_convention"]
        == "native_energy_per_cell_per_kelvin"
    )
    assert calculator.completed is True


def test_ha_calculator_emits_structured_events() -> None:
    """HACalculator emits INFO, PROGRESS and RESULT events."""
    observer = ListObserver()
    calculator = HACalculator(
        _minimal_input(),
        options=HAOptions(
            temperature_min=0.0,
            temperature_max=100.0,
            temperature_step=100.0,
        ),
        observer=observer,
    )

    calculator.execute()

    levels = [event.level for event in observer.events]
    assert EventLevel.INFO in levels
    assert EventLevel.PROGRESS in levels
    assert EventLevel.RESULT in levels

    progress_events = [
        event for event in observer.events if event.level == EventLevel.PROGRESS
    ]
    assert len(progress_events) == 6
    assert progress_events[-1].progress == pytest.approx(1.0)
    assert progress_events[-1].data["current"] == progress_events[-1].data["total"]

    result_events = [
        event
        for event in observer.events
        if event.level == EventLevel.RESULT and event.data.get("kind")
    ]
    assert {event.data["kind"] for event in result_events} >= {
        "settings",
        "input_summary",
        "thermodynamics",
    }


def test_ha_calculator_rejects_invalid_input() -> None:
    """Invalid input data are rejected during prepare()."""
    invalid = _minimal_input()
    invalid.weights = np.array([0.0], dtype=np.float64)
    calculator = HACalculator(invalid)

    with pytest.raises(ValueError, match="sum of q-point weights"):
        calculator.execute()


def test_ha_calculator_runs_on_mgo_yaml_example() -> None:
    """HACalculator runs on the MgO YAML file used by HA/QHA."""
    reader = HAInputFileReader(DATA)
    ha_input = reader.to_input()
    options = HAOptions(
        temperature_min=0.0,
        temperature_max=300.0,
        temperature_step=300.0,
    )
    result = HACalculator(ha_input, options=options).execute()

    assert result.metadata.module == "ha"
    assert result.results["ha"].volume.size == ha_input.nvol
    assert result.results["ha"].free_energy.shape == (2, ha_input.nvol)
    assert np.all(np.isfinite(result.results["ha"].free_energy))
