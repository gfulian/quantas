"""Degeneracy and invalid-eigenvalue diagnostics of the phase solver."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ElasticMedium
from quantas.core.physics.seismic import ChristoffelSolver, WaveMode


def _isotropic_stiffness(lame: float, shear: float) -> np.ndarray:
    stiffness = np.zeros((6, 6), dtype=float)
    stiffness[:3, :3] = lame
    np.fill_diagonal(stiffness[:3, :3], lame + 2.0 * shear)
    np.fill_diagonal(stiffness[3:, 3:], shear)
    return stiffness


@pytest.mark.physics
@pytest.mark.seismic
def test_isotropic_shear_degeneracy_is_valid() -> None:
    stiffness = _isotropic_stiffness(lame=70.0, shear=50.0)
    density = 3000.0
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))

    for direction in ([1.0, 0.0, 0.0], [1.0, 2.0, 3.0]):
        result = solver.solve_direction(direction)
        expected_s = np.sqrt(1000.0 * 50.0 / density)
        expected_p = np.sqrt(1000.0 * 170.0 / density)
        np.testing.assert_allclose(
            result.phase_speeds,
            [expected_s, expected_s, expected_p],
            rtol=1.0e-12,
        )
        np.testing.assert_array_equal(result.valid_mask, [True, True, True])
        np.testing.assert_array_equal(result.clamped_mask, [False, False, False])
        np.testing.assert_array_equal(result.degeneracy_mask, [True, True, False])
        assert result.has_degeneracy
        assert result.is_fully_valid
        assert result.for_mode(WaveMode.V_S2).degenerate
        assert result.for_mode(WaveMode.V_S1).degenerate
        assert not result.for_mode(WaveMode.V_P).degenerate

        q = result.direction
        p_polarization = result.for_mode(WaveMode.V_P).polarization
        assert abs(np.dot(p_polarization, q)) == pytest.approx(1.0)
        np.testing.assert_allclose(
            result.polarizations[:2] @ q, np.zeros(2), atol=1.0e-12
        )


@pytest.mark.physics
@pytest.mark.seismic
def test_small_negative_eigenvalue_is_clamped_and_reported() -> None:
    stiffness = np.diag([-1.0e-12, 30.0, 40.0, 30.0, 20.0, 10.0])
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), 1000.0))
    result = solver.solve_direction([1.0, 0.0, 0.0])

    np.testing.assert_array_equal(result.valid_mask, [True, True, True])
    np.testing.assert_array_equal(result.clamped_mask, [True, False, False])
    assert result.phase_speeds[0] == pytest.approx(0.0)
    assert result.eigenvalues[0] < 0.0
    assert result.eigenvalues[0] >= -result.eigenvalue_threshold


@pytest.mark.physics
@pytest.mark.seismic
def test_negative_eigenvalue_is_rejected() -> None:
    stiffness = np.diag([-1.0, 30.0, 40.0, 30.0, 20.0, 10.0])
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), 1000.0))
    result = solver.solve_direction([1.0, 0.0, 0.0])

    np.testing.assert_array_equal(result.valid_mask, [False, True, True])
    np.testing.assert_array_equal(result.invalid_mask, [True, False, False])
    np.testing.assert_array_equal(result.clamped_mask, [False, False, False])
    assert np.isnan(result.phase_speeds[0])
    assert not result.is_fully_valid
    assert not result.for_mode(WaveMode.V_S2).valid
