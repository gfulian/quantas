"""Tests for production tensor conventions and user-defined rotations."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.geometry import (
    TensorRotation,
    axis_angle_rotation,
    tensor_frame_mapping,
    validate_rotation_matrix,
    xyz_rotation_matrix,
)
from quantas.core.physics.elasticity import (
    ElasticTensor,
    StiffnessSymmetryCriterion,
    cartesian_compliance_to_voigt,
    cartesian_stiffness_to_voigt,
    rotate_voigt_stiffness,
    stiffness_symmetry_residual,
    validate_stiffness_matrix,
    voigt_compliance_to_cartesian,
    voigt_stiffness_to_cartesian,
)
from tests.reference.tensor_rotation import (
    rotate_voigt_stiffness as reference_rotate_voigt_stiffness,
)
from tests.reference.tensor_rotation import (
    voigt_stiffness_to_cartesian as reference_voigt_stiffness_to_cartesian,
)


@pytest.mark.physics
@pytest.mark.elasticity
def test_voigt_conversion_matches_reference_oracle(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """The extracted production conversion preserves the frozen reference convention."""
    expected = reference_voigt_stiffness_to_cartesian(hydroxylapatite_tensor.stiffness)
    actual = voigt_stiffness_to_cartesian(hydroxylapatite_tensor.stiffness)
    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(
        cartesian_stiffness_to_voigt(actual),
        hydroxylapatite_tensor.stiffness,
    )


@pytest.mark.physics
@pytest.mark.elasticity
def test_compliance_conversion_round_trip_is_exact(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """Engineering shear factors survive a full compliance round trip."""
    cartesian = voigt_compliance_to_cartesian(hydroxylapatite_tensor.compliance)
    reconstructed = cartesian_compliance_to_voigt(cartesian)
    np.testing.assert_allclose(
        reconstructed,
        hydroxylapatite_tensor.compliance,
        atol=0.0,
        rtol=0.0,
    )


@pytest.mark.physics
@pytest.mark.elasticity
def test_validation_preserves_frobenius_criterion() -> None:
    """The extracted validator preserves the frozen symmetry threshold."""
    matrix = np.eye(6)
    matrix[0, 1] = 4.0e-4
    matrix[1, 0] = 0.0

    assert stiffness_symmetry_residual(matrix) == pytest.approx(np.sqrt(2.0) * 4.0e-4)
    np.testing.assert_array_equal(validate_stiffness_matrix(matrix), matrix)

    matrix[0, 1] = 8.0e-4
    with pytest.raises(ValueError, match="symmetric"):
        validate_stiffness_matrix(matrix)


@pytest.mark.physics
@pytest.mark.elasticity
def test_central_validation_can_preserve_elementwise_seismic_criterion() -> None:
    """The shared validator explicitly represents the stricter SEISMIC rule."""
    matrix = np.eye(6)
    matrix[0, 1] = 8.0e-9
    validate_stiffness_matrix(
        matrix,
        symmetry_tolerance=1.0e-8,
        symmetry_criterion=StiffnessSymmetryCriterion.ELEMENTWISE,
        copy=False,
    )

    matrix[0, 1] = 1.2e-8
    with pytest.raises(ValueError, match="symmetric"):
        validate_stiffness_matrix(
            matrix,
            symmetry_tolerance=1.0e-8,
            symmetry_criterion="elementwise",
            copy=False,
        )


@pytest.mark.physics
@pytest.mark.elasticity
def test_rotation_matches_reference_oracle(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """The public rotation operation matches the independent rank-four oracle."""
    rotation = axis_angle_rotation((1.0, -2.0, 3.0), 0.731)
    expected = reference_rotate_voigt_stiffness(
        hydroxylapatite_tensor.stiffness,
        rotation,
    )
    actual = rotate_voigt_stiffness(hydroxylapatite_tensor.stiffness, rotation)
    np.testing.assert_allclose(actual, expected, atol=3.0e-14, rtol=2.0e-15)
    np.testing.assert_allclose(
        hydroxylapatite_tensor.rotate(rotation).stiffness,
        expected,
        atol=3.0e-14,
        rtol=2.0e-15,
    )


@pytest.mark.physics
@pytest.mark.elasticity
def test_rotation_round_trip_recovers_original_components(
    hydroxylapatite_tensor: ElasticTensor,
) -> None:
    """A frame transformation followed by its inverse is lossless numerically."""
    rotation = axis_angle_rotation((0.4, -0.7, 0.2), 1.137)
    rotated = hydroxylapatite_tensor.rotate(rotation)
    recovered = rotated.rotate(rotation.T)
    np.testing.assert_allclose(
        recovered.stiffness,
        hydroxylapatite_tensor.stiffness,
        atol=2.0e-13,
        rtol=2.0e-15,
    )


@pytest.mark.physics
@pytest.mark.elasticity
def test_rotation_validation_rejects_left_handed_frame() -> None:
    """A reflection cannot be silently accepted as a tensor rotation."""
    reflection = np.diag([1.0, 1.0, -1.0])
    with pytest.raises(ValueError, match="right-handed"):
        validate_rotation_matrix(reflection)


@pytest.mark.physics
@pytest.mark.elasticity
def test_xyz_rotation_uses_fixed_source_axes_in_documented_order() -> None:
    """The XYZ helper is exactly Rz @ Ry @ Rx."""
    expected = (
        axis_angle_rotation((0.0, 0.0, 1.0), np.deg2rad(30.0))
        @ axis_angle_rotation((0.0, 1.0, 0.0), np.deg2rad(20.0))
        @ axis_angle_rotation((1.0, 0.0, 0.0), np.deg2rad(10.0))
    )
    actual = xyz_rotation_matrix((10.0, 20.0, 30.0))
    np.testing.assert_allclose(actual, expected, atol=2.0e-15, rtol=0.0)


@pytest.mark.physics
@pytest.mark.elasticity
def test_tensor_rotation_preserves_user_provenance() -> None:
    """A rotation specification owns an immutable validated matrix."""
    rotation = TensorRotation.from_xyz(0.0, 0.0, 30.0)
    assert rotation.kind.value == "xyz"
    assert rotation.angles == (0.0, 0.0, 30.0)
    assert rotation.angle_unit == "degree"
    assert rotation.matrix.flags.writeable is False

    metadata = tensor_frame_mapping(rotation)
    assert metadata["source_frame"] == "source"
    assert metadata["analysis_frame"] == "rotated"
    assert metadata["transformed"] is True
    np.testing.assert_allclose(metadata["component_transform"], rotation.matrix)
    assert tensor_frame_mapping(None) == {
        "source_frame": "source",
        "analysis_frame": "source",
        "transformed": False,
    }
