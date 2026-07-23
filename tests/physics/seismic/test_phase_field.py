"""Sampled phase-field models and mode-pair diagnostics."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ElasticMedium
from quantas.core.physics.seismic import (
    MODE_PAIR_INDEX,
    MODE_PAIR_ORDER,
    MODE_PAIR_SYMBOLS,
    ChristoffelSolver,
    ModePair,
    WaveMode,
    build_phase_field,
)


def _isotropic_stiffness(lame: float, shear: float) -> np.ndarray:
    stiffness = np.zeros((6, 6), dtype=float)
    stiffness[:3, :3] = lame
    np.fill_diagonal(stiffness[:3, :3], lame + 2.0 * shear)
    np.fill_diagonal(stiffness[3:, 3:], shear)
    return stiffness


@pytest.mark.physics
@pytest.mark.seismic
def test_mode_pair_order_and_symbols_are_explicit() -> None:
    assert MODE_PAIR_ORDER == (
        ModePair.V_S2_V_S1,
        ModePair.V_S1_V_P,
    )
    assert MODE_PAIR_INDEX == {
        ModePair.V_S2_V_S1: 0,
        ModePair.V_S1_V_P: 1,
    }
    assert MODE_PAIR_SYMBOLS == {
        ModePair.V_S2_V_S1: "V_S2-V_S1",
        ModePair.V_S1_V_P: "V_S1-V_P",
    }


@pytest.mark.physics
@pytest.mark.seismic
def test_directional_results_stack_into_read_only_phase_field(
    hydroxylapatite_data: tuple[np.ndarray, float],
    generic_directions: np.ndarray,
) -> None:
    stiffness, density = hydroxylapatite_data
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), density))
    directional = [
        solver.solve_direction(direction) for direction in generic_directions
    ]

    field = build_phase_field(directional)

    assert field.n_points == len(generic_directions)
    assert field.directions.shape == (4, 3)
    assert field.eigenvalues.shape == (4, 3)
    assert field.phase_speeds.shape == (4, 3)
    assert field.polarizations.shape == (4, 3, 3)
    assert field.pair_eigenvalue_gaps.shape == (4, 2)
    assert field.pair_degeneracy_mask.shape == (4, 2)
    np.testing.assert_allclose(field.directions, generic_directions)

    fast = field.for_mode(WaveMode.V_S1)
    np.testing.assert_array_equal(fast.phase_speeds, field.phase_speeds[:, 1])
    np.testing.assert_array_equal(fast.polarizations, field.polarizations[:, 1])

    shear_pair = field.for_pair(ModePair.V_S2_V_S1)
    np.testing.assert_array_equal(
        shear_pair.eigenvalue_gaps,
        field.pair_eigenvalue_gaps[:, 0],
    )

    arrays = (
        field.directions,
        field.eigenvalues,
        field.phase_speeds,
        field.polarizations,
        field.pair_degeneracy_mask,
        fast.phase_speeds,
        shear_pair.degeneracy_mask,
    )
    assert all(not array.flags.writeable for array in arrays)
    with pytest.raises(ValueError, match="read-only"):
        field.phase_speeds[0, 0] = 0.0


@pytest.mark.physics
@pytest.mark.seismic
def test_isotropic_field_reports_shear_axes() -> None:
    stiffness = _isotropic_stiffness(lame=70.0, shear=50.0)
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), 3000.0))
    directions = np.asarray(
        [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 2.0, 3.0]],
        dtype=float,
    )
    field = build_phase_field(solver.solve_direction(q) for q in directions)

    np.testing.assert_array_equal(field.shear_axis_candidate_mask, True)
    np.testing.assert_array_equal(field.upper_axis_candidate_mask, False)
    np.testing.assert_array_equal(field.acoustic_axis_candidate_mask, True)
    np.testing.assert_array_equal(field.triple_degeneracy_mask, False)
    np.testing.assert_array_equal(field.candidate_indices(), [0, 1, 2])
    np.testing.assert_array_equal(
        field.candidate_indices(ModePair.V_S2_V_S1),
        [0, 1, 2],
    )
    assert field.candidate_indices(ModePair.V_S1_V_P).size == 0
    np.testing.assert_allclose(
        field.candidate_directions(ModePair.V_S2_V_S1),
        field.directions,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_triple_degeneracy_is_distinguished_from_shear_degeneracy() -> None:
    shear = 50.0
    stiffness = _isotropic_stiffness(lame=-shear, shear=shear)
    solver = ChristoffelSolver(ElasticMedium(ElasticTensor(stiffness), 2500.0))
    field = build_phase_field([solver.solve_direction([1.0, 2.0, 3.0])])

    np.testing.assert_array_equal(field.pair_degeneracy_mask, [[True, True]])
    np.testing.assert_array_equal(field.triple_degeneracy_mask, [True])
    np.testing.assert_array_equal(field.acoustic_axis_candidate_mask, [True])


@pytest.mark.physics
@pytest.mark.seismic
def test_phase_field_builder_rejects_empty_or_invalid_entries() -> None:
    with pytest.raises(ValueError, match="At least one"):
        build_phase_field([])
    with pytest.raises(TypeError, match="DirectionalPhaseResult"):
        build_phase_field([object()])  # type: ignore[list-item]
