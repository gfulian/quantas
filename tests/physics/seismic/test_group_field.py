"""Sampled group fields and mode-specific views."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ElasticMedium
from quantas.core.physics.seismic import (
    ChristoffelSolver,
    DirectionalGroupResult,
    GroupFieldResult,
    GroupModeField,
    WaveMode,
    build_group_field,
)


@pytest.mark.physics
@pytest.mark.seismic
def test_directional_group_results_stack_into_read_only_field(
    hydroxylapatite_data: tuple[np.ndarray, float],
    generic_directions: np.ndarray,
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    directional = [
        solver.solve_group_direction(direction) for direction in generic_directions
    ]

    field = build_group_field(directional)

    assert isinstance(field, GroupFieldResult)
    assert field.n_points == 4
    assert field.directions.shape == (4, 3)
    assert field.eigenvalue_gradients.shape == (4, 3, 3)
    assert field.group_velocities.shape == (4, 3, 3)
    assert field.group_speeds.shape == (4, 3)
    assert field.ray_directions.shape == (4, 3, 3)
    assert field.power_flow_angles.shape == (4, 3)
    np.testing.assert_allclose(field.directions, generic_directions)

    fast = field.for_mode(WaveMode.V_S1)
    assert isinstance(fast, GroupModeField)
    assert fast.mode is WaveMode.V_S1
    assert fast.phase.mode is WaveMode.V_S1
    np.testing.assert_array_equal(
        fast.group_velocities,
        field.group_velocities[:, 1, :],
    )
    np.testing.assert_array_equal(
        fast.phase.phase_speeds,
        field.phase.phase_speeds[:, 1],
    )

    arrays = (
        field.eigenvalue_gradients,
        field.group_velocities,
        field.group_speeds,
        field.ray_directions,
        field.power_flow_angles,
        field.valid_mask,
        field.resolved_mask,
        fast.group_velocities,
        fast.phase.phase_speeds,
    )
    assert all(not array.flags.writeable for array in arrays)
    with pytest.raises(ValueError, match="read-only"):
        field.group_speeds[0, 0] = 0.0


@pytest.mark.physics
@pytest.mark.seismic
def test_group_field_preserves_directional_results(
    hydroxylapatite_data: tuple[np.ndarray, float],
    generic_directions: np.ndarray,
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    directional = tuple(
        solver.solve_group_direction(direction) for direction in generic_directions
    )
    field = build_group_field(directional)

    np.testing.assert_allclose(
        field.eigenvalue_gradients,
        np.stack([item.eigenvalue_gradients for item in directional]),
    )
    np.testing.assert_allclose(
        field.group_velocities,
        np.stack([item.group_velocities for item in directional]),
    )
    np.testing.assert_allclose(
        field.group_speeds,
        np.stack([item.group_speeds for item in directional]),
    )
    np.testing.assert_allclose(
        field.ray_directions,
        np.stack([item.ray_directions for item in directional]),
    )
    np.testing.assert_allclose(
        field.power_flow_angles,
        np.stack([item.power_flow_angles for item in directional]),
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_group_field_builder_rejects_empty_or_invalid_entries() -> None:
    with pytest.raises(ValueError, match="At least one"):
        build_group_field([])
    with pytest.raises(TypeError, match="DirectionalGroupResult"):
        build_group_field([object()])  # type: ignore[list-item]


@pytest.mark.physics
@pytest.mark.seismic
def test_group_result_type_is_public() -> None:
    assert DirectionalGroupResult.__name__ == "DirectionalGroupResult"
