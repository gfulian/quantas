from __future__ import annotations

import numpy as np
import pytest

from quantas.core.events import EventLevel, ListObserver
from quantas.models import ResultData
from quantas.modules.qha.calculator import QHACalculator
from quantas.modules.qha.models import QHAInput, QHAOptions, QHAResult


def make_input() -> QHAInput:
    volume = np.linspace(9.0, 11.0, 7)
    energy = 0.5 * (volume - 10.0) ** 2 - 5.0
    frequencies = np.ones((1, 1, volume.size), dtype=float) * 100.0
    return QHAInput(
        jobname="calculator-test",
        volume=volume,
        energy=energy,
        qpoints=1,
        weights=np.array([1.0]),
        frequencies=frequencies,
    )


def make_options(**kwargs) -> QHAOptions:
    values = dict(
        temperature_min=300.0,
        temperature_max=300.0,
        pressure_min=0.0,
        pressure_max=0.0,
        minimization="poly",
        free_energy_degree=2,
    )
    values.update(kwargs)
    return QHAOptions(**values)


def test_calculator_execute_returns_result_data() -> None:
    observer = ListObserver()
    calculator = QHACalculator(make_input(), make_options(), observer=observer)

    result = calculator.execute()

    assert isinstance(result, ResultData)
    assert result.metadata.module == "qha"
    assert result.metadata.method == "quasi-harmonic"
    assert result.metadata.schema_version == "2.0"
    assert isinstance(result.results["qha"], QHAResult)
    assert (
        result.results["qha"].metadata["thermodynamic_unit_convention"]
        == "native_energy_per_cell_per_kelvin"
    )
    assert result.results["qha"].completed is True
    assert result.results["qha"].equilibrium_volume[0, 0] == pytest.approx(
        10.0, abs=1.0e-8
    )
    assert calculator.completed is True
    assert any(event.message == "QHA workflow completed" for event in observer.events)


def test_calculator_emits_debug_fit_records_when_requested() -> None:
    observer = ListObserver()
    calculator = QHACalculator(
        make_input(),
        make_options(debug=True),
        observer=observer,
    )

    result = calculator.execute()

    debug_events = [
        event for event in observer.events if event.level is EventLevel.DEBUG
    ]
    fit_events = [
        event for event in debug_events if event.data.get("kind") == "qha_fit_record"
    ]

    assert result.results["qha"].fit_records
    assert fit_events
    assert fit_events[0].data["quantity"] == "F"


def test_calculator_stores_warnings_for_stopped_workflow() -> None:
    input_data = QHAInput(
        jobname="bad-calculator",
        volume=np.array([1.0, 2.0, 3.0]),
        energy=np.array([1.0, np.nan, 3.0]),
    )
    options = make_options(
        temperature_max=301.0,
        temperature_step=1.0,
        pressure_min=0.0,
        pressure_max=1.0,
        pressure_step=1.0,
        fit_failure_policy="stop",
        max_consecutive_failures=2,
        debug=True,
    )
    observer = ListObserver()
    calculator = QHACalculator(input_data, options, observer=observer)

    result = calculator.execute()

    assert result.results["qha"].completed is False
    assert result.results["qha"].failed_points
    assert result.warnings
    assert any(event.level is EventLevel.WARNING for event in observer.events)


def test_calculator_prepare_validates_options() -> None:
    calculator = QHACalculator(
        make_input(),
        make_options(temperature_step=0.0),
    )

    with pytest.raises(ValueError):
        calculator.prepare()
