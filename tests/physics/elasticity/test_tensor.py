"""Tests for the shared elastic-tensor representation."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor


def test_tensor_stores_voigt_and_cartesian_representations(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    tensor = hydroxylapatite_tensor
    assert tensor.stiffness.shape == (6, 6)
    assert tensor.compliance.shape == (6, 6)
    assert tensor.stiffness_tensor.shape == (3, 3, 3, 3)
    assert tensor.compliance_tensor.shape == (3, 3, 3, 3)
    np.testing.assert_allclose(
        tensor.stiffness @ tensor.compliance,
        np.eye(6),
        atol=1.0e-10,
    )


def test_tensor_rejects_invalid_shape_values_symmetry_and_singularity() -> None:
    with pytest.raises(ValueError, match="shape"):
        ElasticTensor(np.eye(3))

    non_finite = np.eye(6)
    non_finite[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        ElasticTensor(non_finite)

    asymmetric = np.eye(6)
    asymmetric[0, 1] = 1.0
    with pytest.raises(ValueError, match="symmetric"):
        ElasticTensor(asymmetric)

    with pytest.raises(ValueError, match="singular"):
        ElasticTensor(np.zeros((6, 6)))


def test_tensor_owns_immutable_numerical_arrays() -> None:
    """External mutation cannot invalidate cached tensor representations."""
    stiffness = np.eye(6)
    tensor = ElasticTensor(stiffness)
    stiffness[0, 0] = 2.0

    assert tensor.stiffness[0, 0] == 1.0
    with pytest.raises(ValueError, match="read-only"):
        tensor.stiffness[0, 0] = 3.0


def test_tensor_eigenvalues_match_characterized_reference(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    np.testing.assert_allclose(
        hydroxylapatite_tensor.eigenvalues,
        [39.687, 39.687, 61.007, 116.82176234, 122.015, 358.23723766],
        rtol=1.0e-10,
        atol=1.0e-10,
    )
