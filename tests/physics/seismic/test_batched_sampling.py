"""Batched sampling of phase, group and enhancement fields."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.geometry import Hemisphere, create_spherical_grid
from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ElasticMedium
from quantas.core.physics.seismic import (
    ChristoffelSolver,
    SamplingLevel,
    SeismicFieldResult,
    sample_seismic_field,
)


@pytest.mark.physics
@pytest.mark.seismic
def test_batched_enhancement_sampling_matches_pointwise_solver(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    grid = create_spherical_grid(5, 7, hemisphere=Hemisphere.UPPER)

    sampled = sample_seismic_field(
        solver,
        grid,
        level=SamplingLevel.ENHANCEMENT,
        batch_size=4,
    )
    pointwise = [
        solver.solve_enhancement_direction(direction)
        for direction in grid.flat_directions
    ]

    np.testing.assert_allclose(
        sampled.phase.eigenvalues,
        np.stack([item.group.phase.eigenvalues for item in pointwise]),
        rtol=2.0e-13,
        atol=2.0e-13,
    )
    np.testing.assert_allclose(
        sampled.phase.phase_speeds,
        np.stack([item.group.phase.phase_speeds for item in pointwise]),
        rtol=2.0e-13,
        atol=2.0e-13,
    )
    assert sampled.group is not None
    np.testing.assert_allclose(
        sampled.group.group_velocities,
        np.stack([item.group.group_velocities for item in pointwise]),
        rtol=5.0e-13,
        atol=5.0e-13,
    )
    assert sampled.enhancement is not None
    np.testing.assert_allclose(
        sampled.enhancement.log10_enhancement,
        np.stack([item.log10_enhancement for item in pointwise]),
        rtol=2.0e-11,
        atol=2.0e-11,
        equal_nan=True,
    )
    np.testing.assert_array_equal(
        sampled.enhancement.caustic_candidate_mask,
        np.stack([item.caustic_candidate_mask for item in pointwise]),
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_batched_sampler_preserves_the_isotropic_limit() -> None:
    bulk = 100.0
    shear = 50.0
    c11 = bulk + 4.0 * shear / 3.0
    c12 = bulk - 2.0 * shear / 3.0
    stiffness = np.asarray(
        [
            [c11, c12, c12, 0.0, 0.0, 0.0],
            [c12, c11, c12, 0.0, 0.0, 0.0],
            [c12, c12, c11, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, shear, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, shear, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, shear],
        ],
        dtype=float,
    )
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), 3000.0))
    grid = create_spherical_grid(4, 7, hemisphere="upper")

    sampled = sample_seismic_field(solver, grid, batch_size=5)

    np.testing.assert_allclose(
        sampled.phase.phase_speeds[:, 0],
        sampled.phase.phase_speeds[:, 1],
        atol=1.0e-13,
    )
    assert np.all(sampled.phase.degeneracy_mask[:, :2])
    assert sampled.group is not None
    assert np.all(~sampled.group.resolved_mask[:, :2])
    assert sampled.enhancement is not None
    assert np.all(~sampled.enhancement.valid_mask[:, :2])
    np.testing.assert_allclose(
        sampled.enhancement.enhancement[:, 2],
        1.0,
        atol=2.0e-13,
    )
    np.testing.assert_allclose(
        sampled.enhancement.log10_enhancement[:, 2],
        0.0,
        atol=2.0e-13,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_sampling_levels_build_only_requested_quantities(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    grid = create_spherical_grid(3, 5, hemisphere="upper")

    phase = sample_seismic_field(solver, grid, level="phase", batch_size=6)
    group = sample_seismic_field(solver, grid, level="group", batch_size=6)

    assert phase.level is SamplingLevel.PHASE
    assert phase.group is None
    assert phase.enhancement is None
    assert group.level is SamplingLevel.GROUP
    assert group.group is not None
    assert group.enhancement is None


@pytest.mark.physics
@pytest.mark.seismic
def test_batch_size_does_not_change_numerical_results(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    grid = create_spherical_grid(6, 9, hemisphere="upper")

    small = sample_seismic_field(solver, grid, batch_size=1)
    large = sample_seismic_field(solver, grid, batch_size=1000)

    np.testing.assert_allclose(small.phase.phase_speeds, large.phase.phase_speeds)
    assert small.group is not None and large.group is not None
    np.testing.assert_allclose(small.group.group_speeds, large.group.group_speeds)
    assert small.enhancement is not None and large.enhancement is not None
    np.testing.assert_allclose(
        small.enhancement.log10_enhancement,
        large.enhancement.log10_enhancement,
        equal_nan=True,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_sampling_reports_completed_directions_per_batch(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    grid = create_spherical_grid(4, 5, hemisphere="upper")
    calls: list[tuple[int, int]] = []

    sample_seismic_field(
        solver,
        grid,
        level="phase",
        batch_size=6,
        progress_callback=lambda current, total: calls.append((current, total)),
    )

    assert calls == [(6, 20), (12, 20), (18, 20), (20, 20)]


@pytest.mark.physics
@pytest.mark.seismic
def test_sampler_uses_batched_kernel_and_can_track_axes(
    hydroxylapatite_data: tuple[np.ndarray, float],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    grid = create_spherical_grid(4, 7, hemisphere="upper")

    def fail(*args: object, **kwargs: object) -> None:
        raise AssertionError("The pointwise solver must not be called by the sampler.")

    monkeypatch.setattr(solver, "solve_direction", fail)
    monkeypatch.setattr(solver, "solve_group_direction", fail)
    monkeypatch.setattr(solver, "solve_enhancement_direction", fail)

    sampled = sample_seismic_field(
        solver,
        grid,
        level="phase",
        batch_size=5,
        track_polarization_axes=True,
    )

    assert isinstance(sampled, SeismicFieldResult)
    assert sampled.tracking is not None
    assert sampled.tracking.polarizations.shape == (grid.size, 3, 3)
    resolved = sampled.tracking.resolved_mask
    points, branches = np.nonzero(resolved)
    local_modes = sampled.tracking.branch_mode_indices[points, branches]
    overlaps = np.abs(
        np.einsum(
            "ni,ni->n",
            sampled.tracking.polarizations[points, branches],
            sampled.phase.polarizations[points, local_modes],
        )
    )
    np.testing.assert_allclose(overlaps, 1.0, atol=1.0e-12)


@pytest.mark.physics
@pytest.mark.seismic
def test_sampled_fields_are_read_only(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    grid = create_spherical_grid(3, 5, hemisphere="upper")
    sampled = sample_seismic_field(solver, grid, batch_size=4)

    arrays = [sampled.phase.phase_speeds, sampled.phase.polarizations]
    assert sampled.group is not None
    arrays.extend([sampled.group.group_speeds, sampled.group.ray_directions])
    assert sampled.enhancement is not None
    arrays.extend(
        [
            sampled.enhancement.log10_enhancement,
            sampled.enhancement.caustic_candidate_mask,
        ]
    )
    assert all(not array.flags.writeable for array in arrays)


@pytest.mark.physics
@pytest.mark.seismic
@pytest.mark.parametrize("batch_size", [0, -1, True, 1.5])
def test_sampler_rejects_invalid_batch_sizes(
    hydroxylapatite_data: tuple[np.ndarray, float],
    batch_size: object,
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    grid = create_spherical_grid(3, 5, hemisphere="upper")
    with pytest.raises(ValueError, match="batch_size"):
        sample_seismic_field(solver, grid, batch_size=batch_size)  # type: ignore[arg-type]
