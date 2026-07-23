"""Axial sign alignment and continuity tracking of polarization fields."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ElasticMedium
from quantas.core.physics.seismic import (
    BRANCH_INDEX,
    BRANCH_ORDER,
    ChristoffelSolver,
    PolarizationBranch,
    WaveMode,
    align_axial_vector,
    build_phase_field,
    track_polarizations,
)


def _read_only(array: np.ndarray) -> np.ndarray:
    copied = np.array(array, copy=True)
    copied.setflags(write=False)
    return copied


def _isotropic_stiffness(lame: float, shear: float) -> np.ndarray:
    stiffness = np.zeros((6, 6), dtype=float)
    stiffness[:3, :3] = lame
    np.fill_diagonal(stiffness[:3, :3], lame + 2.0 * shear)
    np.fill_diagonal(stiffness[3:, 3:], shear)
    return stiffness


@pytest.mark.physics
@pytest.mark.seismic
def test_polarization_branch_order_is_separate_from_wave_mode_order() -> None:
    assert BRANCH_ORDER == (
        PolarizationBranch.SHEAR_A,
        PolarizationBranch.SHEAR_B,
        PolarizationBranch.P,
    )
    assert BRANCH_INDEX == {
        PolarizationBranch.SHEAR_A: 0,
        PolarizationBranch.SHEAR_B: 1,
        PolarizationBranch.P: 2,
    }


@pytest.mark.physics
@pytest.mark.seismic
def test_align_axial_vector_selects_sign_without_changing_magnitude() -> None:
    reference = np.asarray([1.0, 2.0, 3.0])
    vector = np.asarray([-2.0, -4.0, -6.0])
    aligned = align_axial_vector(reference, vector)

    assert np.dot(reference, aligned) >= 0.0
    assert np.linalg.norm(aligned) == pytest.approx(np.linalg.norm(vector))
    np.testing.assert_array_equal(aligned, -vector)
    assert not aligned.flags.writeable

    for invalid in ([1.0, 2.0], [0.0, 0.0, 0.0], [np.nan, 0.0, 0.0]):
        with pytest.raises(ValueError):
            align_axial_vector(reference, invalid)


@pytest.mark.physics
@pytest.mark.seismic
def test_tracking_removes_arbitrary_axial_sign_reversals(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    result = solver.solve_direction([1.0, 2.0, 3.0])
    results = [result, result, result]

    modified = []
    for index, result in enumerate(results):
        polarizations = np.array(result.polarizations, copy=True)
        if index == 1:
            polarizations *= -1.0
        modified.append(replace(result, polarizations=_read_only(polarizations)))

    field = build_phase_field(modified)
    tracked = track_polarizations(field)

    overlaps = np.einsum(
        "nij,nij->ni",
        tracked.polarizations[:-1],
        tracked.polarizations[1:],
    )
    assert np.all(overlaps >= 0.0)
    np.testing.assert_array_equal(tracked.sign_flip_mask[1], [True, True, True])
    assert not np.any(tracked.segment_start_mask[1:])
    assert np.all(tracked.resolved_mask)


@pytest.mark.physics
@pytest.mark.seismic
def test_shear_tracking_can_exchange_local_modes_without_renaming_them(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    first = solver.solve_direction([1.0, 2.0, 3.0])
    second = solver.solve_direction([1.01, 2.0, 3.0])

    swapped = np.array(second.polarizations, copy=True)
    swapped[:2] = swapped[[1, 0]]
    second_swapped = replace(second, polarizations=_read_only(swapped))
    field = build_phase_field([first, second_swapped])
    original_speeds = np.array(field.phase_speeds, copy=True)

    tracked = track_polarizations(field)

    np.testing.assert_array_equal(tracked.branch_mode_indices[1], [1, 0, 2])
    assert tracked.shear_swap_mask[1]
    assert tracked.mode_at(1, PolarizationBranch.SHEAR_A) is WaveMode.V_S1
    assert tracked.mode_at(1, PolarizationBranch.SHEAR_B) is WaveMode.V_S2
    np.testing.assert_array_equal(field.phase_speeds, original_speeds)
    assert (
        field.for_mode(WaveMode.V_S2).phase_speeds[1]
        <= field.for_mode(WaveMode.V_S1).phase_speeds[1]
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_degenerate_shear_subspace_alignment() -> None:
    stiffness = _isotropic_stiffness(lame=70.0, shear=50.0)
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), 3000.0))
    directions = (
        [1.0, 0.0, 0.0],
        [1.0, 0.1, 0.0],
        [1.0, 0.2, 0.1],
    )
    field = build_phase_field(solver.solve_direction(q) for q in directions)
    tracked = track_polarizations(field)

    np.testing.assert_array_equal(tracked.branch_mode_indices[:, :2], -1)
    np.testing.assert_array_equal(tracked.resolved_mask[:, :2], False)
    np.testing.assert_array_equal(
        tracked.shear_subspace_rotation_mask,
        [False, True, True],
    )
    np.testing.assert_array_equal(
        tracked.shear_permutation_ambiguous_mask,
        [True, True, True],
    )
    assert tracked.mode_at(1, PolarizationBranch.SHEAR_A) is None
    assert np.all(np.isfinite(tracked.polarizations))
    np.testing.assert_allclose(
        np.einsum("nij,nkj->nik", tracked.polarizations, tracked.polarizations),
        np.broadcast_to(np.eye(3), (3, 3, 3)),
        atol=1.0e-12,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_invalid_modes_start_new_tracking_segments() -> None:
    stiffness = np.diag([-1.0, 30.0, 40.0, 30.0, 20.0, 10.0])
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), 1000.0))
    field = build_phase_field(
        [
            solver.solve_direction([1.0, 0.0, 0.0]),
            solver.solve_direction([0.0, 1.0, 0.0]),
        ]
    )
    tracked = track_polarizations(field)

    assert not field.valid_mask[0, 0]
    np.testing.assert_array_equal(tracked.segment_start_mask[0], True)
    np.testing.assert_array_equal(tracked.segment_start_mask[1, :2], True)
    assert np.isnan(tracked.continuity_scores[1, 0])


@pytest.mark.physics
@pytest.mark.seismic
def test_triple_degeneracy_starts_new_segments_for_all_branches() -> None:
    shear = 50.0
    stiffness = _isotropic_stiffness(lame=-shear, shear=shear)
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), 2500.0))
    field = build_phase_field(
        solver.solve_direction(q) for q in ([1.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    )
    tracked = track_polarizations(field)

    np.testing.assert_array_equal(field.triple_degeneracy_mask, True)
    np.testing.assert_array_equal(tracked.branch_mode_indices, -1)
    np.testing.assert_array_equal(tracked.resolved_mask, False)
    np.testing.assert_array_equal(tracked.segment_start_mask, True)
    assert not np.any(tracked.shear_subspace_rotation_mask)


@pytest.mark.physics
@pytest.mark.seismic
def test_tracking_results_and_branch_views_are_read_only(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    field = build_phase_field(
        solver.solve_direction(q) for q in ([1.0, 2.0, 3.0], [1.1, 2.0, 3.0])
    )
    tracked = track_polarizations(field)
    branch = tracked.for_branch(PolarizationBranch.SHEAR_A)

    arrays = (
        tracked.polarizations,
        tracked.branch_mode_indices,
        tracked.sign_flip_mask,
        tracked.continuity_scores,
        tracked.resolved_mask,
        branch.polarizations,
        branch.local_mode_indices,
    )
    assert all(not array.flags.writeable for array in arrays)
    with pytest.raises(ValueError, match="read-only"):
        tracked.polarizations[0, 0, 0] = 0.0
    with pytest.raises(IndexError):
        tracked.mode_at(10, PolarizationBranch.P)
    with pytest.raises(ValueError):
        tracked.for_branch("v_s1")


@pytest.mark.physics
@pytest.mark.seismic
def test_tracking_validates_inputs(
    hydroxylapatite_data: tuple[np.ndarray, float],
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    field = build_phase_field([solver.solve_direction([1.0, 2.0, 3.0])])

    with pytest.raises(TypeError, match="PhaseFieldResult"):
        track_polarizations(object())  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="permutation_tolerance"):
        track_polarizations(field, permutation_tolerance=-1.0)
