"""Analytical group velocities, ray directions and power-flow angles."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ElasticMedium
from quantas.core.physics.seismic import ChristoffelSolver, WaveMode
from tests.reference.seismic_reference import SeismicFormulaReference


def _isotropic_stiffness(lame: float, shear: float) -> np.ndarray:
    stiffness = np.zeros((6, 6), dtype=float)
    stiffness[:3, :3] = lame
    np.fill_diagonal(stiffness[:3, :3], lame + 2.0 * shear)
    np.fill_diagonal(stiffness[3:, 3:], shear)
    return stiffness


def _raw_eigenvalues(
    solver: ChristoffelSolver,
    direction: np.ndarray,
) -> np.ndarray:
    matrix = np.einsum(
        "j,ijkl,l->ik",
        direction,
        solver.reduced_stiffness,
        direction,
        optimize=True,
    )
    matrix = 0.5 * (matrix + matrix.T)
    return np.linalg.eigvalsh(matrix)


@pytest.mark.physics
@pytest.mark.seismic
@pytest.mark.baseline
def test_group_solver_matches_characterized_directional_reference(
    hydroxylapatite_data: tuple[np.ndarray, float],
    hydroxylapatite_reference: SeismicFormulaReference,
    generic_directions: np.ndarray,
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))

    for direction in generic_directions:
        actual = solver.solve_group_direction(direction)
        expected = hydroxylapatite_reference.solve(direction)
        np.testing.assert_allclose(
            actual.eigenvalue_gradients,
            expected.eigenvalue_gradients,
            rtol=1.0e-13,
            atol=3.0e-14,
        )
        np.testing.assert_allclose(
            actual.group_velocities,
            expected.group_velocities,
            rtol=1.0e-13,
            atol=5.0e-15,
        )
        np.testing.assert_allclose(
            actual.group_speeds,
            expected.group_speeds,
            rtol=1.0e-13,
            atol=5.0e-15,
        )
        np.testing.assert_allclose(
            actual.ray_directions,
            expected.group_directions,
            rtol=1.0e-13,
            atol=5.0e-15,
        )
        np.testing.assert_allclose(
            actual.power_flow_angles,
            expected.power_flow_angles,
            rtol=0.0,
            atol=1.0e-8,
        )
        np.testing.assert_array_equal(actual.valid_mask, True)
        np.testing.assert_array_equal(actual.resolved_mask, True)


@pytest.mark.physics
@pytest.mark.seismic
def test_eigenvalue_gradients_match_centered_finite_differences(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    direction = np.asarray([1.0, 2.0, 3.0], dtype=float)
    direction /= np.linalg.norm(direction)
    analytical = solver.solve_group_direction(direction).eigenvalue_gradients
    step = 2.0e-5

    numerical = np.empty((3, 3), dtype=float)
    for coordinate in range(3):
        offset = np.zeros(3, dtype=float)
        offset[coordinate] = step
        numerical[:, coordinate] = (
            _raw_eigenvalues(solver, direction + offset)
            - _raw_eigenvalues(solver, direction - offset)
        ) / (2.0 * step)

    np.testing.assert_allclose(
        analytical,
        numerical,
        rtol=2.0e-8,
        atol=2.0e-7,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_group_velocity_satisfies_radial_identity_and_angle_relation(
    hydroxylapatite_data: tuple[np.ndarray, float],
    generic_directions: np.ndarray,
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))

    for direction in generic_directions:
        result = solver.solve_group_direction(direction)
        radial_components = np.einsum(
            "mi,i->m",
            result.group_velocities,
            result.direction,
            optimize=True,
        )
        np.testing.assert_allclose(
            radial_components,
            result.phase_speeds,
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        assert np.all(result.group_speeds >= result.phase_speeds)
        np.testing.assert_allclose(
            np.cos(result.power_flow_angles),
            result.phase_speeds / result.group_speeds,
            rtol=1.0e-12,
            atol=1.0e-12,
        )


@pytest.mark.physics
@pytest.mark.seismic
def test_isotropic_group_velocity_is_radial_for_all_modes() -> None:
    density = 3000.0
    solver = ChristoffelSolver(
        ElasticMedium(
            ElasticTensor(_isotropic_stiffness(lame=70.0, shear=50.0)), density
        )
    )
    result = solver.solve_group_direction([1.0, 2.0, 3.0])

    expected = result.phase_speeds[:, np.newaxis] * result.direction
    np.testing.assert_allclose(result.group_velocities, expected, atol=2.0e-15)
    np.testing.assert_allclose(result.group_speeds, result.phase_speeds)
    np.testing.assert_allclose(
        result.ray_directions,
        np.broadcast_to(result.direction, (3, 3)),
        atol=2.0e-15,
    )
    np.testing.assert_allclose(result.power_flow_angles, 0.0)
    np.testing.assert_array_equal(result.valid_mask, [True, True, True])
    np.testing.assert_array_equal(result.resolved_mask, [False, False, True])
    assert not result.for_mode(WaveMode.V_S2).resolved
    assert not result.for_mode(WaveMode.V_S1).resolved
    assert result.for_mode(WaveMode.V_P).resolved


@pytest.mark.physics
@pytest.mark.seismic
def test_zero_or_invalid_phase_speed_masks_group_quantities() -> None:
    stiffness = np.diag([-1.0e-12, 30.0, 40.0, 30.0, 20.0, 10.0])
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), 1000.0))
    clamped = solver.solve_group_direction([1.0, 0.0, 0.0])

    assert clamped.phase.clamped_mask[0]
    assert not clamped.valid_mask[0]
    assert np.isnan(clamped.group_velocities[0]).all()
    assert np.isnan(clamped.group_speeds[0])
    assert np.isnan(clamped.ray_directions[0]).all()
    assert np.isnan(clamped.power_flow_angles[0])

    invalid_stiffness = np.diag([-1.0, 30.0, 40.0, 30.0, 20.0, 10.0])
    invalid_solver = ChristoffelSolver(
        ElasticMedium(ElasticTensor(invalid_stiffness), 1000.0)
    )
    invalid = invalid_solver.solve_group_direction([1.0, 0.0, 0.0])
    assert not invalid.valid_mask[0]
    assert np.isnan(invalid.group_velocities[0]).all()


@pytest.mark.physics
@pytest.mark.seismic
def test_group_result_arrays_are_read_only(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    result = ChristoffelSolver(
        ElasticMedium(ElasticTensor(stiffness), density)
    ).solve_group_direction([1.0, 2.0, 3.0])

    arrays = (
        result.eigenvalue_gradients,
        result.group_velocities,
        result.group_speeds,
        result.ray_directions,
        result.power_flow_angles,
        result.valid_mask,
        result.resolved_mask,
    )
    assert all(not array.flags.writeable for array in arrays)
    with pytest.raises(ValueError, match="read-only"):
        result.group_speeds[0] = 0.0


@pytest.mark.physics
@pytest.mark.seismic
def test_group_quantities_follow_density_scaling(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    direction = np.asarray([1.0, 2.0, 3.0], dtype=float)
    reference = ChristoffelSolver(
        ElasticMedium(ElasticTensor(stiffness), density)
    ).solve_group_direction(direction)
    scaled = ChristoffelSolver(
        ElasticMedium(ElasticTensor(stiffness), 4.0 * density)
    ).solve_group_direction(direction)

    np.testing.assert_allclose(
        scaled.eigenvalue_gradients,
        reference.eigenvalue_gradients / 4.0,
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        scaled.group_velocities,
        reference.group_velocities / 2.0,
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        scaled.group_speeds,
        reference.group_speeds / 2.0,
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(scaled.ray_directions, reference.ray_directions)
    np.testing.assert_allclose(
        scaled.power_flow_angles,
        reference.power_flow_angles,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_group_field_has_expected_antipodal_symmetry(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    direction = np.asarray([1.0, -2.0, 3.0], dtype=float)
    positive = solver.solve_group_direction(direction)
    negative = solver.solve_group_direction(-direction)

    np.testing.assert_allclose(negative.phase_speeds, positive.phase_speeds)
    np.testing.assert_allclose(
        negative.eigenvalue_gradients,
        -positive.eigenvalue_gradients,
    )
    np.testing.assert_allclose(
        negative.group_velocities,
        -positive.group_velocities,
    )
    np.testing.assert_allclose(negative.group_speeds, positive.group_speeds)
    np.testing.assert_allclose(
        negative.ray_directions,
        -positive.ray_directions,
    )
    np.testing.assert_allclose(
        negative.power_flow_angles,
        positive.power_flow_angles,
    )
