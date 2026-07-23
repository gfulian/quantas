"""Units, mode names, and known Quantas 0.9 SEISMIC conventions."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from tests.reference.seismic_reference import (
    REFERENCE_MODE_NAMES,
    MODE_ORDER,
    SeismicFormulaReference,
    ReferenceWaveMode,
    original_cartesian_angles,
    reference_stiffness_tensor,
    physical_cartesian_angles,
)


def _isotropic_stiffness(lame: float, shear: float) -> np.ndarray:
    """Build an isotropic stiffness matrix in Voigt notation."""
    stiffness = np.zeros((6, 6), dtype=float)
    stiffness[:3, :3] = lame
    np.fill_diagonal(stiffness[:3, :3], lame + 2.0 * shear)
    np.fill_diagonal(stiffness[3:, 3:], shear)
    return stiffness


@pytest.mark.physics
@pytest.mark.seismic
def test_official_mode_order_matches_ascending_christoffel_eigenvalues() -> None:
    """Ascending eigenvalues map to V_S2, V_S1, and V_P."""
    assert MODE_ORDER == (
        ReferenceWaveMode.V_S2,
        ReferenceWaveMode.V_S1,
        ReferenceWaveMode.V_P,
    )
    assert REFERENCE_MODE_NAMES == {
        ReferenceWaveMode.V_S2: "slow_secondary",
        ReferenceWaveMode.V_S1: "fast_secondary",
        ReferenceWaveMode.V_P: "primary",
    }


@pytest.mark.physics
@pytest.mark.seismic
def test_reference_voigt_mapping_matches_elastic_tensor(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    """The independent reference mapping agrees with the shared ElasticTensor."""
    stiffness, density = hydroxylapatite_data
    historical = reference_stiffness_tensor(stiffness)
    current = ElasticTensor(stiffness).stiffness_tensor
    np.testing.assert_allclose(historical, current, rtol=0.0, atol=0.0)


@pytest.mark.physics
@pytest.mark.seismic
def test_density_factor_matches_a_full_si_calculation() -> None:
    """The factor 1000/rho converts GPa directly to (km/s)^2."""
    stiffness_gpa = _isotropic_stiffness(lame=70.0, shear=50.0)
    density = 3000.0
    direction = np.asarray([1.0, 2.0, 3.0])
    direction /= np.linalg.norm(direction)

    reference = SeismicFormulaReference(stiffness_gpa, density)
    result = reference.solve(direction)

    stiffness_pa = reference.reduced_stiffness * density / 1000.0 * 1.0e9
    gamma_si = np.dot(direction, np.dot(direction, stiffness_pa)) / density
    speeds_si_km_s = np.sqrt(np.linalg.eigvalsh(gamma_si)) / 1000.0

    np.testing.assert_allclose(
        result.phase_speeds,
        speeds_si_km_s,
        rtol=1.0e-13,
        atol=1.0e-13,
    )
    expected_s = np.sqrt(1000.0 * 50.0 / density)
    expected_p = np.sqrt(1000.0 * (70.0 + 2.0 * 50.0) / density)
    np.testing.assert_allclose(
        result.phase_speeds,
        [expected_s, expected_s, expected_p],
        rtol=1.0e-13,
        atol=1.0e-13,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_density_scaling_and_enhancement_invariance(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    """Changing density scales velocities but not directions or focusing."""
    stiffness, density = hydroxylapatite_data
    direction = np.asarray([1.0, 2.0, 3.0])
    first = SeismicFormulaReference(stiffness, density).solve(direction)
    second = SeismicFormulaReference(stiffness, 4.0 * density).solve(direction)

    np.testing.assert_allclose(second.eigenvalues, first.eigenvalues / 4.0)
    np.testing.assert_allclose(second.phase_speeds, first.phase_speeds / 2.0)
    np.testing.assert_allclose(second.group_velocities, first.group_velocities / 2.0)
    np.testing.assert_allclose(second.group_speeds, first.group_speeds / 2.0)
    np.testing.assert_allclose(second.group_directions, first.group_directions)
    np.testing.assert_allclose(second.power_flow_angles, first.power_flow_angles)
    np.testing.assert_allclose(second.enhancement, first.enhancement, rtol=1.0e-11)


@pytest.mark.physics
@pytest.mark.seismic
def test_original_azimuth_differs_from_physical_convention() -> None:
    """Freeze the known original atan2(x, y) azimuth defect."""
    direction = np.asarray([1.0, 0.0, 0.0])
    _, original_phi = original_cartesian_angles(direction)
    _, physical_phi = physical_cartesian_angles(direction)
    assert original_phi == pytest.approx(np.pi / 2.0)
    assert physical_phi == pytest.approx(0.0)
