"""Tracking diagnostics across a crystallographic acoustic axis."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ElasticMedium
from quantas.core.physics.seismic import (
    ChristoffelSolver,
    ModePair,
    build_phase_field,
    track_polarizations,
)


@pytest.mark.physics
@pytest.mark.seismic
def test_hexagonal_axis_preserves_polarization_ambiguity(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    offsets = np.linspace(-0.1, 0.1, 21)
    directions = np.column_stack(
        [offsets, np.zeros_like(offsets), np.ones_like(offsets)]
    )

    field = build_phase_field(solver.solve_direction(q) for q in directions)
    tracked = track_polarizations(field)
    center = len(offsets) // 2

    np.testing.assert_array_equal(
        field.candidate_indices(ModePair.V_S2_V_S1),
        [center],
    )
    assert field.pair_eigenvalue_gaps[center, 0] == pytest.approx(0.0)
    assert field.phase_speeds[center, 0] == pytest.approx(field.phase_speeds[center, 1])
    np.testing.assert_array_equal(tracked.branch_mode_indices[center, :2], -1)
    np.testing.assert_array_equal(tracked.resolved_mask[center, :2], False)
    assert tracked.shear_subspace_rotation_mask[center]
    assert tracked.shear_permutation_ambiguous_mask[center]

    outside_axis = np.ones(field.n_points, dtype=bool)
    outside_axis[center] = False
    assert np.all(tracked.resolved_mask[outside_axis, :2])
    assert np.nanmin(tracked.continuity_scores[:, :2]) > 0.99
    assert np.all(field.phase_speeds[:, 0] <= field.phase_speeds[:, 1])
    assert np.all(field.phase_speeds[:, 1] <= field.phase_speeds[:, 2])
