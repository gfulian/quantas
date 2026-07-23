"""Tests for high-level elasticity analysis helpers."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import OrthorhombicElasticTensor
from quantas.modules.elasticity.analysis import (
    calculate_basic_properties,
    create_elastic_tensor,
    create_principal_plane_grids,
    specialize_tensor,
    validate_input,
)
from quantas.modules.elasticity.models import ElasticityInput, ElasticityOptions


def test_validate_input_rejects_invalid_matrix_and_options() -> None:
    options = ElasticityOptions()
    with pytest.raises(ValueError, match="missing"):
        validate_input(ElasticityInput(stiffness=None), options)
    with pytest.raises(ValueError, match=r"shape \(6, 6\)"):
        validate_input(ElasticityInput(stiffness=np.eye(3)), options)
    asymmetric = np.eye(6)
    asymmetric[0, 1] = 1.0
    with pytest.raises(ValueError, match="symmetric"):
        validate_input(ElasticityInput(stiffness=asymmetric), options)
    non_finite = np.eye(6)
    non_finite[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        validate_input(ElasticityInput(stiffness=non_finite), options)
    with pytest.raises(ValueError, match="at least 2"):
        validate_input(
            ElasticityInput(stiffness=np.eye(6)),
            ElasticityOptions(ntheta=1),
        )


def test_create_elastic_tensor_reports_singular_matrix() -> None:
    with pytest.raises(ValueError, match="matrix is singular"):
        create_elastic_tensor(ElasticityInput(stiffness=np.zeros((6, 6))))


def test_polar_planes_are_closed_and_follow_cartesian_planes() -> None:
    planes = create_principal_plane_grids(5)
    expected = np.linspace(0.0, 2.0 * np.pi, 5)
    np.testing.assert_allclose(planes["xy"]["theta"], np.pi / 2.0)
    np.testing.assert_allclose(planes["xy"]["phi"], expected)
    np.testing.assert_allclose(planes["xz"]["theta"], expected)
    np.testing.assert_allclose(planes["xz"]["phi"], 0.0)
    np.testing.assert_allclose(planes["yz"]["theta"], expected)
    np.testing.assert_allclose(planes["yz"]["phi"], np.pi / 2.0)


def test_workflow_uses_orthorhombic_specialization() -> None:
    matrix = np.diag([150.0, 160.0, 170.0, 50.0, 60.0, 70.0])
    input_data = ElasticityInput(stiffness=matrix)
    tensor = create_elastic_tensor(input_data)
    result = calculate_basic_properties(tensor, input_data)
    specialized = specialize_tensor(tensor, result)
    assert result.crystal_system == "orthorhombic"
    assert isinstance(specialized, OrthorhombicElasticTensor)
