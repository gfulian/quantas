from __future__ import annotations

import numpy as np
import pytest

from quantas.modules.qha.analysis import (
    analyze_volume_minimization,
    initialize_result,
    prepare_free_energy_grid,
)
from quantas.modules.qha.models import QHAInput, QHAOptions


def make_input() -> QHAInput:
    volume = np.linspace(9.0, 11.0, 7)
    energy = 0.5 * (volume - 10.0) ** 2 - 5.0
    return QHAInput(
        jobname="analysis-test",
        volume=volume,
        energy=energy,
    )


def test_initialize_result_sets_grids_and_metadata() -> None:
    options = QHAOptions(
        temperature_min=300.0,
        temperature_max=400.0,
        temperature_step=100.0,
        pressure_min=0.0,
        pressure_max=1.0,
        pressure_step=1.0,
    )
    result = initialize_result(make_input(), options)

    assert result.jobname == "analysis-test"
    assert result.temperature.tolist() == [300.0, 400.0]
    assert result.pressure.tolist() == [0.0, 1.0]
    assert result.equilibrium_volume.shape == (2, 2)
    assert result.valid_mask.shape == (2, 2)
    assert result.metadata["minimization"] == "poly"


def test_prepare_free_energy_grid_broadcasts_static_curve() -> None:
    grid = prepare_free_energy_grid([1.0, 2.0, 3.0], ntemperatures=2, nvolumes=3)

    assert grid.shape == (2, 3)
    np.testing.assert_allclose(grid[0], [1.0, 2.0, 3.0])
    np.testing.assert_allclose(grid[1], [1.0, 2.0, 3.0])


def test_prepare_free_energy_grid_rejects_wrong_shape() -> None:
    with pytest.raises(ValueError, match="shape"):
        prepare_free_energy_grid(np.ones((2, 4)), ntemperatures=2, nvolumes=3)


def test_polynomial_volume_minimization_grid() -> None:
    options = QHAOptions(
        temperature_min=300.0,
        temperature_max=300.0,
        temperature_step=1.0,
        pressure_min=0.0,
        pressure_max=0.0,
        pressure_step=1.0,
        minimization="poly",
        free_energy_degree=2,
        debug=True,
    )
    result = analyze_volume_minimization(make_input(), options)

    assert result.completed is True
    assert result.failed_points == []
    assert result.valid_mask.shape == (1, 1)
    assert bool(result.valid_mask[0, 0]) is True
    assert result.equilibrium_volume[0, 0] == pytest.approx(10.0, abs=1.0e-8)
    assert len(result.fit_records) == 2
    assert result.fit_records[0].quantity == "F"
    assert result.fit_records[0].pressure is None
    assert result.fit_records[1].quantity == "KT/Kp"


def test_eos_volume_minimization_grid() -> None:
    volume = np.linspace(9.5, 10.5, 9)
    energy = -7.0 + 0.08 * (volume - 10.0) ** 2
    input_data = QHAInput(jobname="eos-analysis", volume=volume, energy=energy)
    options = QHAOptions(
        temperature_min=300.0,
        temperature_max=300.0,
        pressure_min=0.0,
        pressure_max=0.0,
        minimization="eos",
        eos="birch-murnaghan",
        debug=True,
    )

    result = analyze_volume_minimization(input_data, options)

    assert result.completed is True
    assert bool(result.valid_mask[0, 0]) is True
    assert result.equilibrium_volume[0, 0] == pytest.approx(10.0, abs=1.0e-2)
    assert len(result.fit_records) == 1
    assert result.fit_records[0].method == "eos"


def test_failed_fits_are_recorded_and_stop_after_limit() -> None:
    input_data = QHAInput(
        jobname="bad-analysis",
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

    result = analyze_volume_minimization(input_data, options)

    assert result.completed is False
    assert len(result.failed_points) == 6
    assert result.valid_mask.sum() == 0
    assert result.metadata["stop"]["max_consecutive_failures"] == 2
    assert result.metadata["stop"]["failure_scope"] == "temperature"


def test_failure_policy_raise_propagates_runtime_error() -> None:
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
        analyze_volume_minimization(input_data, options)
