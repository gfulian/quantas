"""Public phase-solver API and acoustic-mode semantics."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ElasticMedium
from quantas.core.physics.seismic import (
    MODE_DESCRIPTIONS,
    MODE_INDEX,
    MODE_ORDER,
    MODE_SYMBOLS,
    ChristoffelSolver,
    WaveMode,
)


@pytest.mark.physics
@pytest.mark.seismic
def test_wave_mode_order_symbols_and_descriptions_are_unambiguous() -> None:
    assert MODE_ORDER == (WaveMode.V_S2, WaveMode.V_S1, WaveMode.V_P)
    assert MODE_INDEX == {
        WaveMode.V_S2: 0,
        WaveMode.V_S1: 1,
        WaveMode.V_P: 2,
    }
    assert MODE_SYMBOLS == {
        WaveMode.V_S2: "V_S2",
        WaveMode.V_S1: "V_S1",
        WaveMode.V_P: "V_P",
    }
    assert MODE_DESCRIPTIONS[WaveMode.V_S2] == "slow quasi-shear wave"
    assert MODE_DESCRIPTIONS[WaveMode.V_S1] == "fast quasi-shear wave"
    assert MODE_DESCRIPTIONS[WaveMode.V_P] == "quasi-longitudinal wave"


@pytest.mark.physics
@pytest.mark.seismic
def test_directional_result_exposes_named_mode_results(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    result = solver.solve_direction([1.0, 2.0, 3.0])

    slow = result.for_mode(WaveMode.V_S2)
    fast = result.for_mode("v_s1")
    primary = result.for_mode(WaveMode.V_P)

    assert slow.mode is WaveMode.V_S2
    assert fast.mode is WaveMode.V_S1
    assert primary.mode is WaveMode.V_P
    assert slow.phase_speed == pytest.approx(result.phase_speeds[0])
    assert fast.phase_speed == pytest.approx(result.phase_speeds[1])
    assert primary.phase_speed == pytest.approx(result.phase_speeds[2])
    assert slow.phase_speed <= fast.phase_speed <= primary.phase_speed
    np.testing.assert_array_equal(slow.polarization, result.polarizations[0])
    assert slow.valid and fast.valid and primary.valid
    assert not slow.degenerate

    with pytest.raises(ValueError):
        result.for_mode("primary")


@pytest.mark.physics
@pytest.mark.seismic
def test_solver_and_results_own_read_only_numerical_arrays(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    tensor = ElasticTensor(stiffness)
    medium = ElasticMedium(tensor, density)
    solver = ChristoffelSolver(medium)
    result = solver.solve_direction([1.0, 2.0, 3.0])

    arrays = (
        solver.reduced_stiffness,
        solver.christoffel_hessian,
        result.direction,
        result.christoffel,
        result.eigenvalues,
        result.polarizations,
        result.phase_speeds,
        result.eigenvalue_gaps,
        result.relative_eigenvalue_gaps,
        result.valid_mask,
        result.clamped_mask,
        result.degeneracy_mask,
        result.invalid_mask,
    )
    assert all(not array.flags.writeable for array in arrays)
    with pytest.raises(ValueError, match="read-only"):
        result.phase_speeds[0] = 0.0
