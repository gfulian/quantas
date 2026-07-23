"""Tests for elastic symmetry detection and specialization."""

from __future__ import annotations

import numpy as np

from quantas.core.physics.elasticity import (
    ElasticTensor,
    OrthorhombicElasticTensor,
    detect_elastic_symmetry,
    is_elastic_symmetry,
    specialize_elastic_tensor,
)


def test_hydroxylapatite_is_detected_as_hexagonal(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    matrix = hydroxylapatite_tensor.stiffness
    assert detect_elastic_symmetry(matrix) == "hexagonal"
    assert is_elastic_symmetry(matrix, "hexagonal") is True
    assert is_elastic_symmetry(matrix, "cubic") is False
    np.testing.assert_allclose(
        matrix[5, 5], (matrix[0, 0] - matrix[0, 1]) / 2.0, atol=1e-3
    )


def test_orthorhombic_tensor_is_specialized_without_aliases() -> None:
    matrix = np.diag([180.0, 170.0, 160.0, 60.0, 55.0, 50.0])
    matrix[0, 1] = matrix[1, 0] = 70.0
    matrix[0, 2] = matrix[2, 0] = 65.0
    matrix[1, 2] = matrix[2, 1] = 62.0
    tensor = ElasticTensor(matrix)
    assert detect_elastic_symmetry(matrix) == "orthorhombic"
    specialized = specialize_elastic_tensor(tensor, "orthorhombic")
    assert isinstance(specialized, OrthorhombicElasticTensor)
    np.testing.assert_allclose(specialized.stiffness, tensor.stiffness)
