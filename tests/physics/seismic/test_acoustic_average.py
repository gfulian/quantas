"""Acoustic velocity averages and Kieffer cutoff validation."""

from __future__ import annotations

import math

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import (
    ChristoffelSolver,
    ElasticMedium,
    average_acoustic_phase_velocities,
    kieffer_cutoff_frequencies,
    solve_phase_directions,
)


def _isotropic_stiffness(lame: float, shear: float) -> np.ndarray:
    """Return an isotropic stiffness matrix in GPa."""
    stiffness = np.zeros((6, 6), dtype=np.float64)
    stiffness[:3, :3] = lame
    np.fill_diagonal(stiffness[:3, :3], lame + 2.0 * shear)
    np.fill_diagonal(stiffness[3:, 3:], shear)
    return stiffness


def _isotropic_solver() -> ChristoffelSolver:
    """Return a stable isotropic reference solver."""
    return ChristoffelSolver(
        ElasticMedium(ElasticTensor(_isotropic_stiffness(70.0, 50.0)), 3000.0)
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_arbitrary_direction_solver_normalizes_rows() -> None:
    solver = _isotropic_solver()
    directions = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 2.0]])
    result = solve_phase_directions(solver, directions, batch_size=1)

    assert result.n_points == 2
    np.testing.assert_allclose(np.linalg.norm(result.directions, axis=1), 1.0)
    pointwise = solver.solve_direction(directions[1])
    np.testing.assert_allclose(result.phase_speeds[1], pointwise.phase_speeds)


@pytest.mark.physics
@pytest.mark.seismic
def test_isotropic_average_matches_analytical_velocities() -> None:
    solver = _isotropic_solver()
    result = average_acoustic_phase_velocities(
        solver, mu_order=4, phi_order=8, refinement_factor=2
    )
    shear = math.sqrt(1000.0 * 50.0 / 3000.0)
    longitudinal = math.sqrt(1000.0 * 170.0 / 3000.0)

    np.testing.assert_allclose(
        result.effective_velocities,
        [shear, shear, longitudinal],
        rtol=2.0e-15,
    )
    np.testing.assert_allclose(result.relative_errors, 0.0, atol=2.0e-15)
    assert result.direction_count == 128
    assert result.degenerate_direction_count == 128


@pytest.mark.physics
@pytest.mark.seismic
def test_anisotropic_average_converges_under_refinement(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    refined = average_acoustic_phase_velocities(
        solver, mu_order=24, phi_order=48, refinement_factor=2
    )

    assert np.all(refined.relative_errors < 2.0e-4)


@pytest.mark.physics
@pytest.mark.seismic
def test_cutoff_conversion_matches_defining_equations() -> None:
    velocities = np.array([4.0, 5.0, 8.0], dtype=np.float64)
    volume = 250.0
    result = kieffer_cutoff_frequencies(velocities, volume)
    radius = (6.0 * np.pi**2 / (volume * 1.0e-30)) ** (1.0 / 3.0)

    assert result.brillouin_radius_m_inv == pytest.approx(radius)
    np.testing.assert_allclose(
        result.angular_frequencies_rad_s,
        (2.0 / np.pi) * velocities * 1.0e3 * radius,
    )
    np.testing.assert_allclose(
        result.frequencies_hz,
        velocities * 1.0e3 * radius / np.pi**2,
    )
    np.testing.assert_allclose(
        result.wavenumbers_cm1,
        result.frequencies_hz / (299_792_458.0 * 100.0),
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_acoustic_inputs_are_rejected_explicitly() -> None:
    solver = _isotropic_solver()
    for directions in (
        np.ones(3),
        np.empty((0, 3)),
        np.array([[0.0, 0.0, 0.0]]),
        np.array([[np.nan, 0.0, 1.0]]),
    ):
        with pytest.raises(ValueError):
            solve_phase_directions(solver, directions)

    with pytest.raises(ValueError):
        average_acoustic_phase_velocities(solver, mu_order=1)
    with pytest.raises(ValueError):
        average_acoustic_phase_velocities(solver, refinement_factor=0)
    with pytest.raises(ValueError):
        kieffer_cutoff_frequencies(np.ones(2), 100.0)
    with pytest.raises(ValueError):
        kieffer_cutoff_frequencies(np.ones(3), 0.0)
