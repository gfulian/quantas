"""Finite-difference checks of the analytical SEISMIC derivatives."""

from __future__ import annotations

import numpy as np
import pytest

from tests.reference.seismic_reference import SeismicFormulaReference


def _raw_christoffel(reference: SeismicFormulaReference, q: np.ndarray) -> np.ndarray:
    """Evaluate the quadratic Christoffel form without normalizing q."""
    return np.dot(q, np.dot(q, reference.reduced_stiffness))


def _raw_eigenvalues(reference: SeismicFormulaReference, q: np.ndarray) -> np.ndarray:
    """Evaluate ascending Christoffel eigenvalues without normalizing q."""
    return np.linalg.eigvalsh(_raw_christoffel(reference, q))


@pytest.mark.physics
@pytest.mark.seismic
def test_christoffel_gradient_and_hessian_match_finite_differences(
    hydroxylapatite_reference: SeismicFormulaReference,
) -> None:
    """Analytical matrix derivatives reproduce centered finite differences."""
    q = np.asarray([1.0, 2.0, 3.0], dtype=float)
    q /= np.linalg.norm(q)
    analytical_gradient = hydroxylapatite_reference.christoffel_gradient(q)
    analytical_hessian = hydroxylapatite_reference.christoffel_hessian
    step = 1.0e-5
    basis = np.eye(3)

    numerical_gradient = np.empty((3, 3, 3), dtype=float)
    numerical_hessian = np.empty((3, 3, 3, 3), dtype=float)
    for a in range(3):
        numerical_gradient[a] = (
            _raw_christoffel(hydroxylapatite_reference, q + step * basis[a])
            - _raw_christoffel(hydroxylapatite_reference, q - step * basis[a])
        ) / (2.0 * step)
        for b in range(3):
            numerical_hessian[a, b] = (
                _raw_christoffel(
                    hydroxylapatite_reference,
                    q + step * basis[a] + step * basis[b],
                )
                - _raw_christoffel(
                    hydroxylapatite_reference,
                    q + step * basis[a] - step * basis[b],
                )
                - _raw_christoffel(
                    hydroxylapatite_reference,
                    q - step * basis[a] + step * basis[b],
                )
                + _raw_christoffel(
                    hydroxylapatite_reference,
                    q - step * basis[a] - step * basis[b],
                )
            ) / (4.0 * step * step)

    np.testing.assert_allclose(
        analytical_gradient,
        numerical_gradient,
        rtol=2.0e-10,
        atol=2.0e-9,
    )
    np.testing.assert_allclose(
        analytical_hessian,
        numerical_hessian,
        rtol=2.0e-6,
        atol=2.0e-5,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_eigenvalue_gradient_and_hessian_match_finite_differences(
    hydroxylapatite_reference: SeismicFormulaReference,
) -> None:
    """Analytical eigenvalue derivatives agree away from degeneracies."""
    q = np.asarray([1.0, 2.0, 3.0], dtype=float)
    q /= np.linalg.norm(q)
    result = hydroxylapatite_reference.solve(q)
    step = 2.0e-5
    basis = np.eye(3)

    numerical_gradient = np.empty((3, 3), dtype=float)
    numerical_hessian = np.empty((3, 3, 3), dtype=float)
    for a in range(3):
        numerical_gradient[:, a] = (
            _raw_eigenvalues(hydroxylapatite_reference, q + step * basis[a])
            - _raw_eigenvalues(hydroxylapatite_reference, q - step * basis[a])
        ) / (2.0 * step)
        for b in range(3):
            numerical_hessian[:, a, b] = (
                _raw_eigenvalues(
                    hydroxylapatite_reference,
                    q + step * basis[a] + step * basis[b],
                )
                - _raw_eigenvalues(
                    hydroxylapatite_reference,
                    q + step * basis[a] - step * basis[b],
                )
                - _raw_eigenvalues(
                    hydroxylapatite_reference,
                    q - step * basis[a] + step * basis[b],
                )
                + _raw_eigenvalues(
                    hydroxylapatite_reference,
                    q - step * basis[a] - step * basis[b],
                )
            ) / (4.0 * step * step)

    np.testing.assert_allclose(
        result.eigenvalue_gradients,
        numerical_gradient,
        rtol=2.0e-8,
        atol=2.0e-7,
    )
    np.testing.assert_allclose(
        result.eigenvalue_hessians,
        numerical_hessian,
        rtol=3.0e-5,
        atol=3.0e-4,
    )
