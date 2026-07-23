"""Analytical eigenvalue Hessians and acoustic enhancement factors."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ElasticMedium
from quantas.core.physics.seismic import (
    ChristoffelSolver,
    WaveMode,
    cofactor_matrix,
)
from tests.reference.seismic_reference import SeismicFormulaReference


def _isotropic_stiffness(lame: float, shear: float) -> np.ndarray:
    """Return an isotropic stiffness matrix in Voigt notation."""
    stiffness = np.zeros((6, 6), dtype=float)
    stiffness[:3, :3] = lame
    np.fill_diagonal(stiffness[:3, :3], lame + 2.0 * shear)
    np.fill_diagonal(stiffness[3:, 3:], shear)
    return stiffness


def _raw_eigenvalues(
    solver: ChristoffelSolver,
    direction: np.ndarray,
) -> np.ndarray:
    """Evaluate Christoffel eigenvalues without normalizing the direction."""
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
def test_enhancement_solver_matches_characterized_directional_reference(
    hydroxylapatite_data: tuple[np.ndarray, float],
    hydroxylapatite_reference: SeismicFormulaReference,
    generic_directions: np.ndarray,
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))

    for direction in generic_directions:
        actual = solver.solve_enhancement_direction(direction)
        expected = hydroxylapatite_reference.solve(direction)
        np.testing.assert_allclose(
            actual.eigenvalue_hessians,
            expected.eigenvalue_hessians,
            rtol=2.0e-13,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            actual.enhancement,
            expected.enhancement,
            rtol=5.0e-13,
            atol=5.0e-14,
        )
        np.testing.assert_allclose(
            actual.log10_enhancement,
            np.log10(expected.enhancement),
            rtol=5.0e-13,
            atol=5.0e-14,
        )
        np.testing.assert_array_equal(actual.valid_mask, True)
        np.testing.assert_array_equal(actual.resolved_mask, True)
        np.testing.assert_array_equal(actual.finite_mask, True)
        np.testing.assert_array_equal(actual.caustic_candidate_mask, False)


@pytest.mark.physics
@pytest.mark.seismic
def test_eigenvalue_hessians_match_centered_finite_differences(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    direction = np.asarray([1.0, 2.0, 3.0], dtype=float)
    direction /= np.linalg.norm(direction)
    analytical = solver.solve_enhancement_direction(direction).eigenvalue_hessians
    step = 2.0e-5
    basis = np.eye(3)

    numerical = np.empty((3, 3, 3), dtype=float)
    for first in range(3):
        for second in range(3):
            numerical[:, first, second] = (
                _raw_eigenvalues(
                    solver,
                    direction + step * basis[first] + step * basis[second],
                )
                - _raw_eigenvalues(
                    solver,
                    direction + step * basis[first] - step * basis[second],
                )
                - _raw_eigenvalues(
                    solver,
                    direction - step * basis[first] + step * basis[second],
                )
                + _raw_eigenvalues(
                    solver,
                    direction - step * basis[first] - step * basis[second],
                )
            ) / (4.0 * step * step)

    np.testing.assert_allclose(
        analytical,
        numerical,
        rtol=3.0e-5,
        atol=3.0e-4,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_ray_direction_gradients_match_centered_finite_differences(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    direction = np.asarray([1.0, 2.0, 3.0], dtype=float)
    direction /= np.linalg.norm(direction)
    analytical = solver.solve_enhancement_direction(direction).ray_direction_gradients
    step = 1.0e-6

    numerical = np.empty((3, 3, 3), dtype=float)
    for coordinate in range(3):
        offset = np.zeros(3, dtype=float)
        offset[coordinate] = step
        plus = solver.solve_group_direction(direction + offset).ray_directions
        minus = solver.solve_group_direction(direction - offset).ray_directions
        numerical[:, :, coordinate] = (plus - minus) / (2.0 * step)

    np.testing.assert_allclose(
        analytical,
        numerical,
        rtol=5.0e-7,
        atol=5.0e-8,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_eigenvalue_hessian_obeys_quadratic_homogeneity(
    hydroxylapatite_data: tuple[np.ndarray, float],
    generic_directions: np.ndarray,
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))

    for direction in generic_directions:
        result = solver.solve_enhancement_direction(direction)
        hessian_times_direction = np.einsum(
            "mij,j->mi",
            result.eigenvalue_hessians,
            result.direction,
            optimize=True,
        )
        np.testing.assert_allclose(
            hessian_times_direction,
            result.group.eigenvalue_gradients,
            rtol=2.0e-12,
            atol=2.0e-12,
        )
        radial_curvature = np.einsum(
            "i,mij,j->m",
            result.direction,
            result.eigenvalue_hessians,
            result.direction,
            optimize=True,
        )
        np.testing.assert_allclose(
            radial_curvature,
            2.0 * result.group.phase.eigenvalues,
            rtol=2.0e-12,
            atol=2.0e-12,
        )
        radial_ray_gradient = np.einsum(
            "mij,j->mi",
            result.ray_direction_gradients,
            result.direction,
            optimize=True,
        )
        np.testing.assert_allclose(
            radial_ray_gradient,
            0.0,
            atol=3.0e-14,
        )
        tangent_constraint = np.einsum(
            "mi,mij->mj",
            result.group.ray_directions,
            result.ray_direction_gradients,
            optimize=True,
        )
        np.testing.assert_allclose(
            tangent_constraint,
            0.0,
            atol=3.0e-14,
        )
        np.testing.assert_allclose(
            result.enhancement,
            1.0 / result.area_factors,
            rtol=2.0e-14,
            atol=2.0e-14,
        )


@pytest.mark.physics
@pytest.mark.seismic
def test_isotropic_p_mode_has_unit_enhancement() -> None:
    direction = np.asarray([1.0, 2.0, 3.0], dtype=float)
    direction /= np.linalg.norm(direction)
    solver = ChristoffelSolver(
        ElasticMedium(
            ElasticTensor(_isotropic_stiffness(lame=70.0, shear=50.0)), 3000.0
        )
    )
    result = solver.solve_enhancement_direction(direction)
    p_mode = result.for_mode(WaveMode.V_P)

    expected_gradient = np.eye(3) - np.outer(direction, direction)
    np.testing.assert_allclose(
        p_mode.ray_direction_gradient,
        expected_gradient,
        atol=2.0e-15,
    )
    assert p_mode.area_factor == pytest.approx(1.0, abs=2.0e-15)
    assert p_mode.enhancement == pytest.approx(1.0, abs=2.0e-15)
    assert p_mode.log10_enhancement == pytest.approx(0.0, abs=1.0e-15)
    assert p_mode.valid
    assert p_mode.resolved
    assert p_mode.finite
    assert not p_mode.caustic_candidate

    np.testing.assert_array_equal(result.valid_mask, [False, False, True])
    assert np.isnan(result.enhancement[:2]).all()
    assert np.isnan(result.eigenvalue_hessians[:2]).all()


@pytest.mark.physics
@pytest.mark.seismic
def test_enhancement_is_invariant_under_density_scaling(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    direction = np.asarray([1.0, -3.0, 2.0], dtype=float)
    reference = ChristoffelSolver(
        ElasticMedium(ElasticTensor(stiffness), density)
    ).solve_enhancement_direction(direction)
    scaled = ChristoffelSolver(
        ElasticMedium(ElasticTensor(stiffness), 4.0 * density)
    ).solve_enhancement_direction(direction)

    np.testing.assert_allclose(
        scaled.eigenvalue_hessians,
        reference.eigenvalue_hessians / 4.0,
        rtol=2.0e-12,
        atol=2.0e-12,
    )
    np.testing.assert_allclose(
        scaled.ray_direction_gradients,
        reference.ray_direction_gradients,
        rtol=2.0e-12,
        atol=2.0e-12,
    )
    np.testing.assert_allclose(scaled.area_factors, reference.area_factors)
    np.testing.assert_allclose(scaled.enhancement, reference.enhancement)
    np.testing.assert_allclose(
        scaled.log10_enhancement,
        reference.log10_enhancement,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_enhancement_has_expected_antipodal_symmetry(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    direction = np.asarray([1.0, -2.0, 3.0], dtype=float)
    positive = solver.solve_enhancement_direction(direction)
    negative = solver.solve_enhancement_direction(-direction)

    np.testing.assert_allclose(
        negative.eigenvalue_hessians,
        positive.eigenvalue_hessians,
    )
    np.testing.assert_allclose(
        negative.ray_direction_gradients,
        positive.ray_direction_gradients,
    )
    np.testing.assert_allclose(negative.area_factors, positive.area_factors)
    np.testing.assert_allclose(negative.enhancement, positive.enhancement)
    np.testing.assert_allclose(
        negative.log10_enhancement,
        positive.log10_enhancement,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_caustic_threshold_flags_without_modifying_finite_values(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    direction = np.asarray([1.0, 2.0, 3.0], dtype=float)
    reference = ChristoffelSolver(
        ElasticMedium(ElasticTensor(stiffness), density)
    ).solve_enhancement_direction(direction)
    flagged = ChristoffelSolver(
        ElasticMedium(ElasticTensor(stiffness), density),
        caustic_atol=100.0,
        caustic_rtol=0.0,
    ).solve_enhancement_direction(direction)

    np.testing.assert_allclose(flagged.enhancement, reference.enhancement)
    np.testing.assert_allclose(
        flagged.log10_enhancement,
        reference.log10_enhancement,
    )
    np.testing.assert_array_equal(flagged.caustic_candidate_mask, True)
    np.testing.assert_array_equal(flagged.finite_mask, True)
    assert flagged.has_caustic_candidate


@pytest.mark.physics
@pytest.mark.seismic
def test_invalid_or_degenerate_modes_do_not_report_curvature() -> None:
    isotropic = ChristoffelSolver(
        ElasticMedium(
            ElasticTensor(_isotropic_stiffness(lame=70.0, shear=50.0)), 3000.0
        )
    ).solve_enhancement_direction([1.0, 2.0, 3.0])
    np.testing.assert_array_equal(isotropic.resolved_mask, [False, False, True])
    np.testing.assert_array_equal(isotropic.valid_mask, [False, False, True])

    stiffness = np.diag([-1.0, 30.0, 40.0, 30.0, 20.0, 10.0])
    invalid = ChristoffelSolver(
        ElasticMedium(ElasticTensor(stiffness), 1000.0)
    ).solve_enhancement_direction([1.0, 0.0, 0.0])
    assert not invalid.valid_mask[0]
    assert np.isnan(invalid.enhancement[0])
    assert np.isnan(invalid.log10_enhancement[0])


@pytest.mark.physics
@pytest.mark.seismic
def test_cofactor_matrix_supports_regular_and_singular_inputs() -> None:
    regular = np.asarray(
        [[2.0, 1.0, 3.0], [0.0, 4.0, -1.0], [1.0, 2.0, 5.0]],
        dtype=float,
    )
    expected = np.linalg.det(regular) * np.linalg.inv(regular).T
    np.testing.assert_allclose(cofactor_matrix(regular), expected)

    singular = np.diag([1.0, 1.0, 0.0])
    np.testing.assert_array_equal(
        cofactor_matrix(singular),
        np.diag([0.0, 0.0, 1.0]),
    )
    with pytest.raises(ValueError, match="shape"):
        cofactor_matrix(np.eye(2))


@pytest.mark.physics
@pytest.mark.seismic
def test_enhancement_result_arrays_are_read_only(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    result = ChristoffelSolver(
        ElasticMedium(ElasticTensor(stiffness), density)
    ).solve_enhancement_direction([1.0, 2.0, 3.0])

    arrays = (
        result.eigenvalue_hessians,
        result.ray_direction_gradients,
        result.area_factors,
        result.caustic_thresholds,
        result.enhancement,
        result.log10_enhancement,
        result.valid_mask,
        result.resolved_mask,
        result.finite_mask,
        result.caustic_candidate_mask,
    )
    assert all(not array.flags.writeable for array in arrays)
    with pytest.raises(ValueError, match="read-only"):
        result.enhancement[0] = 0.0
