"""Validation of the experimental batched NumPy SEISMIC kernel."""

from __future__ import annotations

import numpy as np
import pytest

from tests.reference.seismic_reference import SeismicFormulaReference
from tests.reference.seismic_vectorized import solve_batched


@pytest.mark.physics
@pytest.mark.seismic
@pytest.mark.parametrize("calculate_enhancement", [False, True])
def test_batched_kernel_matches_scalar_reference(
    hydroxylapatite_data: tuple[np.ndarray, float],
    generic_directions: np.ndarray,
    calculate_enhancement: bool,
) -> None:
    """Batched phase, group, and enhancement values match the pointwise reference."""
    stiffness, density = hydroxylapatite_data
    reference = SeismicFormulaReference(stiffness, density)
    batched = solve_batched(
        stiffness,
        density,
        generic_directions,
        calculate_enhancement=calculate_enhancement,
    )

    for index, direction in enumerate(generic_directions):
        scalar = reference.solve(direction)
        np.testing.assert_allclose(batched.christoffel[index], scalar.christoffel)
        np.testing.assert_allclose(batched.eigenvalues[index], scalar.eigenvalues)
        np.testing.assert_allclose(batched.phase_speeds[index], scalar.phase_speeds)
        np.testing.assert_allclose(
            batched.eigenvalue_gradients[index], scalar.eigenvalue_gradients
        )
        np.testing.assert_allclose(
            batched.group_velocities[index], scalar.group_velocities
        )
        np.testing.assert_allclose(batched.group_speeds[index], scalar.group_speeds)
        np.testing.assert_allclose(
            batched.group_directions[index], scalar.group_directions
        )
        np.testing.assert_allclose(
            batched.power_flow_angles[index], scalar.power_flow_angles
        )

        overlaps = np.abs(
            np.einsum(
                "mi,mi->m",
                batched.polarizations[index],
                scalar.polarizations,
            )
        )
        np.testing.assert_allclose(overlaps, np.ones(3))

        if calculate_enhancement:
            assert batched.eigenvalue_hessians is not None
            assert batched.enhancement is not None
            np.testing.assert_allclose(
                batched.eigenvalue_hessians[index],
                scalar.eigenvalue_hessians,
                rtol=1.0e-11,
                atol=1.0e-11,
            )
            np.testing.assert_allclose(
                batched.enhancement[index],
                scalar.enhancement,
                rtol=1.0e-10,
                atol=1.0e-10,
            )
        else:
            assert batched.eigenvalue_hessians is None
            assert batched.enhancement is None
