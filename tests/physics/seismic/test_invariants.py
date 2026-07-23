"""Scientific invariants of the historical Christoffel implementation."""

from __future__ import annotations

import numpy as np
import pytest

from tests.reference.seismic_reference import SeismicFormulaReference


@pytest.mark.physics
@pytest.mark.seismic
def test_christoffel_eigenpairs_and_group_identity(
    hydroxylapatite_reference: SeismicFormulaReference,
    generic_directions: np.ndarray,
) -> None:
    """Eigenpairs are orthonormal and vg dot n equals the phase speed."""
    for direction in generic_directions:
        result = hydroxylapatite_reference.solve(direction)
        np.testing.assert_allclose(result.christoffel, result.christoffel.T)
        np.testing.assert_allclose(
            result.polarizations @ result.polarizations.T,
            np.eye(3),
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        for mode in range(3):
            np.testing.assert_allclose(
                result.christoffel @ result.polarizations[mode],
                result.eigenvalues[mode] * result.polarizations[mode],
                rtol=1.0e-12,
                atol=1.0e-12,
            )
        np.testing.assert_allclose(
            result.group_velocities @ result.direction,
            result.phase_speeds,
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        assert np.all(result.phase_speeds[:-1] <= result.phase_speeds[1:])


@pytest.mark.physics
@pytest.mark.seismic
def test_antipodal_symmetry_of_phase_group_and_enhancement(
    hydroxylapatite_reference: SeismicFormulaReference,
    generic_directions: np.ndarray,
) -> None:
    """Wave-normal reversal preserves speeds and reverses group vectors."""
    for direction in generic_directions:
        positive = hydroxylapatite_reference.solve(direction)
        negative = hydroxylapatite_reference.solve(-direction)
        np.testing.assert_allclose(negative.christoffel, positive.christoffel)
        np.testing.assert_allclose(negative.eigenvalues, positive.eigenvalues)
        np.testing.assert_allclose(negative.phase_speeds, positive.phase_speeds)
        np.testing.assert_allclose(
            negative.group_velocities, -positive.group_velocities
        )
        np.testing.assert_allclose(negative.group_speeds, positive.group_speeds)
        np.testing.assert_allclose(
            negative.group_directions, -positive.group_directions
        )
        np.testing.assert_allclose(
            negative.power_flow_angles, positive.power_flow_angles
        )
        np.testing.assert_allclose(
            negative.enhancement, positive.enhancement, rtol=1.0e-11
        )
        overlap = np.abs(
            np.einsum("mi,mi->m", negative.polarizations, positive.polarizations)
        )
        np.testing.assert_allclose(overlap, np.ones(3))


@pytest.mark.physics
@pytest.mark.seismic
def test_isotropic_speeds_and_power_flow() -> (
    None
):
    """The isotropic analytical limit is reproduced, including S degeneracy."""
    lame = 70.0
    shear = 50.0
    density = 3000.0
    stiffness = np.zeros((6, 6), dtype=float)
    stiffness[:3, :3] = lame
    np.fill_diagonal(stiffness[:3, :3], lame + 2.0 * shear)
    np.fill_diagonal(stiffness[3:, 3:], shear)
    reference = SeismicFormulaReference(stiffness, density)
    expected = np.asarray(
        [
            np.sqrt(1000.0 * shear / density),
            np.sqrt(1000.0 * shear / density),
            np.sqrt(1000.0 * (lame + 2.0 * shear) / density),
        ]
    )
    for direction in (
        np.asarray([1.0, 0.0, 0.0]),
        np.asarray([1.0, 2.0, 3.0]),
        np.asarray([-2.0, 1.0, 4.0]),
    ):
        result = reference.solve(direction)
        np.testing.assert_allclose(result.phase_speeds, expected, rtol=1.0e-12)
        np.testing.assert_allclose(result.group_speeds, expected, rtol=1.0e-12)
        np.testing.assert_allclose(result.power_flow_angles, np.zeros(3), atol=1.0e-12)
        p_polarization = result.polarizations[2]
        assert abs(np.dot(p_polarization, result.direction)) == pytest.approx(1.0)
        s_projection = result.polarizations[:2] @ result.direction
        np.testing.assert_allclose(s_projection, np.zeros(2), atol=1.0e-12)
