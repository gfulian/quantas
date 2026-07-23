"""Directional regression tests for the production phase solver."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ElasticMedium
from quantas.core.physics.seismic import ChristoffelSolver
from tests.reference.seismic_reference import SeismicFormulaReference


@pytest.mark.physics
@pytest.mark.seismic
@pytest.mark.baseline
def test_phase_solver_matches_characterized_directional_reference(
    hydroxylapatite_data: tuple[np.ndarray, float],
    hydroxylapatite_reference: SeismicFormulaReference,
    generic_directions: np.ndarray,
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))

    np.testing.assert_allclose(
        solver.christoffel_hessian,
        hydroxylapatite_reference.christoffel_hessian,
        rtol=0.0,
        atol=0.0,
    )
    for direction in generic_directions:
        actual = solver.solve_direction(direction)
        expected = hydroxylapatite_reference.solve(direction)
        np.testing.assert_allclose(
            actual.christoffel,
            expected.christoffel,
            rtol=1.0e-13,
            atol=1.0e-13,
        )
        np.testing.assert_allclose(actual.eigenvalues, expected.eigenvalues)
        np.testing.assert_allclose(actual.phase_speeds, expected.phase_speeds)
        np.testing.assert_allclose(
            solver.christoffel_gradient(direction),
            expected.christoffel_gradient,
            rtol=0.0,
            atol=0.0,
        )
        overlaps = np.abs(
            np.einsum(
                "mi,mi->m",
                actual.polarizations,
                expected.polarizations,
            )
        )
        np.testing.assert_allclose(overlaps, np.ones(3))
        np.testing.assert_array_equal(actual.valid_mask, np.ones(3, dtype=bool))
        np.testing.assert_array_equal(actual.clamped_mask, np.zeros(3, dtype=bool))
        np.testing.assert_array_equal(
            actual.degeneracy_mask,
            np.zeros(3, dtype=bool),
        )


@pytest.mark.physics
@pytest.mark.seismic
def test_christoffel_eigenpairs_are_symmetric_and_orthonormal(
    hydroxylapatite_data: tuple[np.ndarray, float],
    generic_directions: np.ndarray,
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))

    for direction in generic_directions:
        result = solver.solve_direction(direction)
        np.testing.assert_allclose(result.christoffel, result.christoffel.T)
        np.testing.assert_allclose(
            result.polarizations @ result.polarizations.T,
            np.eye(3),
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        for index in range(3):
            np.testing.assert_allclose(
                result.christoffel @ result.polarizations[index],
                result.eigenvalues[index] * result.polarizations[index],
                rtol=1.0e-12,
                atol=1.0e-12,
            )


@pytest.mark.physics
@pytest.mark.seismic
def test_christoffel_gradient_and_hessian_match_finite_differences(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    direction = np.asarray([1.0, 2.0, 3.0], dtype=float)
    direction /= np.linalg.norm(direction)
    step = 1.0e-5

    numerical_gradient = np.empty((3, 3, 3), dtype=float)
    numerical_hessian = np.empty((3, 3, 3, 3), dtype=float)
    for coordinate in range(3):
        offset = np.zeros(3, dtype=float)
        offset[coordinate] = step
        plus = _unconstrained_christoffel(solver, direction + offset)
        minus = _unconstrained_christoffel(solver, direction - offset)
        numerical_gradient[coordinate] = (plus - minus) / (2.0 * step)

        plus_gradient = _unconstrained_gradient(solver, direction + offset)
        minus_gradient = _unconstrained_gradient(solver, direction - offset)
        numerical_hessian[coordinate] = (plus_gradient - minus_gradient) / (2.0 * step)

    np.testing.assert_allclose(
        solver.christoffel_gradient(direction),
        numerical_gradient,
        rtol=2.0e-10,
        atol=2.0e-9,
    )
    np.testing.assert_allclose(
        solver.christoffel_hessian,
        numerical_hessian,
        rtol=2.0e-10,
        atol=2.0e-9,
    )


def _unconstrained_christoffel(
    solver: ChristoffelSolver,
    direction: np.ndarray,
) -> np.ndarray:
    return np.einsum(
        "j,ijkl,l->ik",
        direction,
        solver.reduced_stiffness,
        direction,
        optimize=True,
    )


def _unconstrained_gradient(
    solver: ChristoffelSolver,
    direction: np.ndarray,
) -> np.ndarray:
    symmetrized = solver.reduced_stiffness + np.transpose(
        solver.reduced_stiffness,
        (0, 2, 1, 3),
    )
    return np.einsum("k,iakl->ail", direction, symmetrized, optimize=True)
