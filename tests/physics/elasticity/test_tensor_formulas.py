"""Frozen characterization of tensor conventions and original Quantas formulas."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import (
    ElasticTensor,
    linear_compressibility,
    poisson_ratio,
    shear_modulus,
    young_modulus,
)
from tests.reference.elasticity_reference import ElasticityFormulaReference
from tests.reference.tensor_rotation import (
    cartesian_stiffness_to_voigt,
    rotate_cartesian_stiffness,
    rotate_voigt_stiffness,
    rotation_matrix,
    voigt_stiffness_to_cartesian,
)

_VOIGT_PAIRS = (
    (0, 0),
    (1, 1),
    (2, 2),
    (1, 2),
    (0, 2),
    (0, 1),
)


def _isotropic_stiffness(lame: float, shear: float) -> np.ndarray:
    """Return an isotropic stiffness matrix in the Quantas Voigt ordering."""
    matrix = np.zeros((6, 6), dtype=float)
    matrix[:3, :3] = lame
    np.fill_diagonal(matrix[:3, :3], lame + 2.0 * shear)
    np.fill_diagonal(matrix[3:, 3:], shear)
    return matrix


def _cartesian_compliance_to_voigt(tensor: np.ndarray) -> np.ndarray:
    """Reconstruct engineering compliance coefficients from ``S_ijkl``."""
    matrix = np.empty((6, 6), dtype=float)
    factors = (1.0, 1.0, 1.0, 2.0, 2.0, 2.0)
    for p, (i, j) in enumerate(_VOIGT_PAIRS):
        for q, (k, ell) in enumerate(_VOIGT_PAIRS):
            matrix[p, q] = factors[p] * factors[q] * tensor[i, j, k, ell]
    return matrix


@pytest.mark.physics
@pytest.mark.elasticity
def test_stiffness_voigt_cartesian_round_trip_is_exact(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """The documented ``11,22,33,23,13,12`` ordering is fixed explicitly."""
    independent = voigt_stiffness_to_cartesian(hydroxylapatite_tensor.stiffness)
    np.testing.assert_array_equal(
        hydroxylapatite_tensor.stiffness_tensor,
        independent,
    )
    np.testing.assert_array_equal(
        cartesian_stiffness_to_voigt(hydroxylapatite_tensor.stiffness_tensor),
        hydroxylapatite_tensor.stiffness,
    )


@pytest.mark.physics
@pytest.mark.elasticity
def test_compliance_round_trip_preserves_shear_factors(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """Compliance conversion applies factors two and four to shear entries."""
    tensor = hydroxylapatite_tensor
    reconstructed = _cartesian_compliance_to_voigt(tensor.compliance_tensor)
    np.testing.assert_allclose(reconstructed, tensor.compliance, atol=0.0, rtol=0.0)

    assert tensor.compliance_tensor[1, 2, 1, 2] == pytest.approx(
        tensor.compliance[3, 3] / 4.0
    )
    assert tensor.compliance_tensor[0, 0, 1, 2] == pytest.approx(
        tensor.compliance[0, 3] / 2.0
    )
    assert tensor.stiffness_tensor[1, 2, 1, 2] == pytest.approx(tensor.stiffness[3, 3])
    assert tensor.stiffness_tensor[0, 0, 1, 2] == pytest.approx(tensor.stiffness[0, 3])


@pytest.mark.physics
@pytest.mark.elasticity
@pytest.mark.parametrize("name", ["stiffness_tensor", "compliance_tensor"])
def test_cartesian_tensors_obey_minor_and_major_symmetry(
    hydroxylapatite_tensor: ElasticTensor,
    name: str,
) -> None:
    """Both fourth-rank tensors obey the elastic minor and major symmetries."""
    tensor = getattr(hydroxylapatite_tensor, name)
    np.testing.assert_allclose(tensor, np.swapaxes(tensor, 0, 1), atol=0.0, rtol=0.0)
    np.testing.assert_allclose(tensor, np.swapaxes(tensor, 2, 3), atol=0.0, rtol=0.0)
    np.testing.assert_allclose(
        tensor,
        np.transpose(tensor, (2, 3, 0, 1)),
        atol=1.0e-18,
        rtol=3.0e-16,
    )


@pytest.mark.physics
@pytest.mark.elasticity
def test_directional_operations_match_soec_formulas(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """Pointwise public operations remain equivalent to Quantas 0.9 SOEC."""
    reference = ElasticityFormulaReference(hydroxylapatite_tensor.stiffness)
    angles = (
        (0.321, 0.654, 0.987),
        (1.123, 0.789, 0.456),
        (2.234, 4.321, 1.234),
        (0.5 * np.pi, 1.5 * np.pi, 0.25 * np.pi),
    )
    for theta, phi, chi in angles:
        np.testing.assert_allclose(
            young_modulus(hydroxylapatite_tensor, (theta, phi)),
            reference.young_modulus(theta, phi),
            rtol=2.0e-14,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            linear_compressibility(hydroxylapatite_tensor, (theta, phi)),
            reference.linear_compressibility(theta, phi),
            rtol=2.0e-14,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            shear_modulus(hydroxylapatite_tensor, (theta, phi, chi)),
            reference.shear_modulus(theta, phi, chi),
            rtol=2.0e-14,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            poisson_ratio(hydroxylapatite_tensor, (theta, phi, chi)),
            reference.poisson_ratio(theta, phi, chi),
            rtol=2.0e-14,
            atol=2.0e-14,
        )


@pytest.mark.physics
@pytest.mark.elasticity
def test_isotropic_tensor_is_rotation_invariant() -> None:
    """A fourth-rank isotropic tensor is unchanged by a proper rotation."""
    tensor = voigt_stiffness_to_cartesian(_isotropic_stiffness(70.0, 50.0))
    rotation = rotation_matrix((1.0, -2.0, 3.0), 0.731)
    rotated = rotate_cartesian_stiffness(tensor, rotation)
    np.testing.assert_allclose(rotated, tensor, rtol=2.0e-15, atol=3.0e-14)


@pytest.mark.physics
@pytest.mark.elasticity
def test_rotation_can_exchange_trigonal_c14_and_c15_conventions() -> None:
    """A basal-plane rotation can move a trigonal coefficient between slots.

    This characterization uses an active right-handed rotation of ``+30°``
    around Cartesian ``z``.  Under this convention a tensor with ``C14 = 10``
    and ``C15 = 0`` is represented after rotation with ``C14 = 0`` and
    ``C15 = -10``.  The test captures why frame metadata is scientifically
    necessary before reporting symmetry-dependent matrix entries.
    """
    c11, c12, c13, c33, c44, c14 = 200.0, 80.0, 70.0, 180.0, 60.0, 10.0
    c66 = 0.5 * (c11 - c12)
    stiffness = np.asarray(
        [
            [c11, c12, c13, c14, 0.0, 0.0],
            [c12, c11, c13, -c14, 0.0, 0.0],
            [c13, c13, c33, 0.0, 0.0, 0.0],
            [c14, -c14, 0.0, c44, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, c44, c14],
            [0.0, 0.0, 0.0, 0.0, c14, c66],
        ],
        dtype=float,
    )
    rotation = rotation_matrix((0.0, 0.0, 1.0), np.deg2rad(30.0))
    rotated = rotate_voigt_stiffness(stiffness, rotation)

    assert rotated[0, 3] == pytest.approx(0.0, abs=5.0e-14)
    assert rotated[0, 4] == pytest.approx(-c14, abs=5.0e-14)
    assert rotated[3, 5] == pytest.approx(c14, abs=5.0e-14)
    assert rotated[4, 5] == pytest.approx(0.0, abs=5.0e-14)
