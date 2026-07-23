"""Sampled acoustic enhancement fields and mode-specific views."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ElasticMedium
from quantas.core.physics.seismic import (
    ChristoffelSolver,
    DirectionalEnhancementResult,
    EnhancementFieldResult,
    EnhancementModeField,
    WaveMode,
    build_enhancement_field,
)


@pytest.mark.physics
@pytest.mark.seismic
def test_enhancement_results_form_read_only_field(
    hydroxylapatite_data: tuple[np.ndarray, float],
    generic_directions: np.ndarray,
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    directional = [
        solver.solve_enhancement_direction(direction)
        for direction in generic_directions
    ]

    field = build_enhancement_field(directional)

    assert isinstance(field, EnhancementFieldResult)
    assert field.n_points == 4
    assert field.directions.shape == (4, 3)
    assert field.eigenvalue_hessians.shape == (4, 3, 3, 3)
    assert field.ray_direction_gradients.shape == (4, 3, 3, 3)
    assert field.area_factors.shape == (4, 3)
    assert field.enhancement.shape == (4, 3)
    assert field.log10_enhancement.shape == (4, 3)
    np.testing.assert_allclose(field.directions, generic_directions)
    assert not field.has_caustic_candidate

    fast = field.for_mode(WaveMode.V_S1)
    assert isinstance(fast, EnhancementModeField)
    assert fast.mode is WaveMode.V_S1
    assert fast.group.mode is WaveMode.V_S1
    np.testing.assert_array_equal(
        fast.eigenvalue_hessians,
        field.eigenvalue_hessians[:, 1, :, :],
    )
    np.testing.assert_array_equal(
        fast.log10_enhancement,
        field.log10_enhancement[:, 1],
    )

    arrays = (
        field.eigenvalue_hessians,
        field.ray_direction_gradients,
        field.area_factors,
        field.caustic_thresholds,
        field.enhancement,
        field.log10_enhancement,
        field.valid_mask,
        field.resolved_mask,
        field.finite_mask,
        field.caustic_candidate_mask,
        fast.enhancement,
        fast.group.group_speeds,
    )
    assert all(not array.flags.writeable for array in arrays)
    with pytest.raises(ValueError, match="read-only"):
        field.enhancement[0, 0] = 0.0


@pytest.mark.physics
@pytest.mark.seismic
def test_enhancement_field_preserves_directional_results(
    hydroxylapatite_data: tuple[np.ndarray, float],
    generic_directions: np.ndarray,
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    directional = tuple(
        solver.solve_enhancement_direction(direction)
        for direction in generic_directions
    )
    field = build_enhancement_field(directional)

    np.testing.assert_allclose(
        field.eigenvalue_hessians,
        np.stack([item.eigenvalue_hessians for item in directional]),
    )
    np.testing.assert_allclose(
        field.ray_direction_gradients,
        np.stack([item.ray_direction_gradients for item in directional]),
    )
    np.testing.assert_allclose(
        field.enhancement,
        np.stack([item.enhancement for item in directional]),
    )
    np.testing.assert_allclose(
        field.log10_enhancement,
        np.stack([item.log10_enhancement for item in directional]),
    )
    np.testing.assert_allclose(
        field.group.group_velocities,
        np.stack([item.group.group_velocities for item in directional]),
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_enhancement_field_rejects_invalid_entries() -> None:
    with pytest.raises(ValueError, match="At least one"):
        build_enhancement_field([])
    with pytest.raises(TypeError, match="DirectionalEnhancementResult"):
        build_enhancement_field([object()])  # type: ignore[list-item]


@pytest.mark.physics
@pytest.mark.seismic
def test_enhancement_result_type_is_public() -> None:
    assert DirectionalEnhancementResult.__name__ == "DirectionalEnhancementResult"
