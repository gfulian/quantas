"""Rotation covariance tests for future elastic-frame rotations."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ElasticMedium
from quantas.core.physics.seismic import ChristoffelSolver
from tests.reference.tensor_rotation import rotate_voigt_stiffness, rotation_matrix


@pytest.mark.physics
@pytest.mark.seismic
def test_christoffel_rotation_covariance(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    """Rotating material and wave normal preserves speeds and rotates modes."""
    stiffness, density = hydroxylapatite_data
    rotation = rotation_matrix((1.0, 2.0, -1.0), 0.413)
    rotated_stiffness = rotate_voigt_stiffness(stiffness, rotation)

    original_solver = ChristoffelSolver(
        ElasticMedium(ElasticTensor(stiffness), density)
    )
    rotated_solver = ChristoffelSolver(
        ElasticMedium(ElasticTensor(rotated_stiffness), density)
    )
    direction = np.asarray([1.0, -2.0, 3.0], dtype=float)
    direction /= np.linalg.norm(direction)
    rotated_direction = rotation @ direction

    original = original_solver.solve_direction(direction)
    rotated = rotated_solver.solve_direction(rotated_direction)

    np.testing.assert_allclose(
        rotated.phase_speeds,
        original.phase_speeds,
        rtol=3.0e-13,
        atol=3.0e-13,
    )
    np.testing.assert_allclose(
        rotated.christoffel,
        rotation @ original.christoffel @ rotation.T,
        rtol=3.0e-13,
        atol=3.0e-13,
    )
    expected_polarizations = original.polarizations @ rotation.T
    overlaps = np.abs(
        np.einsum("mi,mi->m", rotated.polarizations, expected_polarizations)
    )
    np.testing.assert_allclose(overlaps, np.ones(3), atol=3.0e-13)
