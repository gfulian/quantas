"""Tests for backend-neutral phonon-mode continuity tracking."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from quantas.core.numerics.phonon_tracking import track_phonon_modes
from quantas.interfaces.crystal.phonon_modes import CrystalPhononModeParser


def test_mode_tracking_recovers_permuted_branches() -> None:
    """Global overlap assignment should restore deliberately permuted modes."""
    reference = np.eye(3, dtype=np.complex128).reshape(3, 1, 3)
    permutation = np.array([1, 2, 0])
    source = reference[permutation]
    frequencies = np.array(
        [
            [[100.0, 200.0, 300.0]],
            [[200.0, 300.0, 100.0]],
        ]
    )
    eigenvectors = np.array([[reference], [source]])

    result = track_phonon_modes(
        frequencies,
        eigenvectors,
        [10.0, 11.0],
        reference_index=0,
    )

    assert result.verified
    np.testing.assert_allclose(result.frequencies[1, 0], [100.0, 200.0, 300.0])
    np.testing.assert_array_equal(result.permutations[1, 0], [2, 0, 1])


def test_mode_tracking_is_invariant_to_complex_mode_phase() -> None:
    """Arbitrary eigenvector phases must not alter mode identity."""
    reference = np.eye(3, dtype=np.complex128).reshape(3, 1, 3)
    phases = np.exp(1j * np.array([0.4, -1.1, 2.2]))
    source = reference * phases[:, None, None]
    frequencies = np.array(
        [
            [[100.0, 200.0, 300.0]],
            [[101.0, 201.0, 301.0]],
        ]
    )
    eigenvectors = np.array([[reference], [source]])

    result = track_phonon_modes(frequencies, eigenvectors, [10.0, 11.0])

    assert result.verified
    np.testing.assert_allclose(result.overlaps[1, 0], 1.0, atol=1.0e-12)


def test_mode_tracking_resolves_rotated_degenerate_subspace() -> None:
    """A rotated basis inside a degenerate eigenspace should remain verified."""
    reference = np.eye(3, dtype=np.complex128).reshape(3, 1, 3)
    rotation = np.array(
        [
            [1.0, 1.0],
            [-1.0, 1.0],
        ],
        dtype=np.complex128,
    ) / np.sqrt(2.0)
    source = reference.copy()
    source[:2] = (rotation @ reference[:2].reshape(2, 3)).reshape(2, 1, 3)
    frequencies = np.array(
        [
            [[100.0, 100.0, 300.0]],
            [[101.0, 101.0, 301.0]],
        ]
    )
    eigenvectors = np.array([[reference], [source]])

    result = track_phonon_modes(frequencies, eigenvectors, [10.0, 11.0])

    assert result.verified
    assert result.degenerate_subspaces == 1
    assert result.ambiguous_assignments == 0
    np.testing.assert_allclose(
        result.minimum_subspace_singular_value,
        1.0,
        atol=1.0e-12,
    )


def test_ambiguous_assignment_is_a_caution_not_an_automatic_failure() -> None:
    """A strong but non-unique overlap should remain usable with a caution."""
    reference = np.eye(3, dtype=np.complex128).reshape(3, 1, 3)
    source = reference.copy()
    source[0, 0] = np.array([0.8, 0.6, 0.0])
    source[1, 0] = np.array([-0.6, 0.8, 0.0])
    frequencies = np.array(
        [
            [[100.0, 200.0, 300.0]],
            [[101.0, 201.0, 301.0]],
        ]
    )
    eigenvectors = np.array([[reference], [source]])

    result = track_phonon_modes(frequencies, eigenvectors, [10.0, 11.0])

    assert result.verified
    assert result.ambiguous_assignments > 0
    assert result.caution_assignments > 0
    assert result.unresolved_assignments == 0


def test_low_overlap_without_independent_fit_support_is_unreliable() -> None:
    """Two sampled states cannot use a trivial fit to rescue a weak overlap."""
    reference = np.eye(3, dtype=np.complex128).reshape(3, 1, 3)
    source = reference.copy()
    source[0, 0] = np.array([0.8, 0.6, 0.0])
    source[1, 0] = np.array([-0.6, 0.8, 0.0])
    frequencies = np.array(
        [
            [[100.0, 200.0, 300.0]],
            [[101.0, 201.0, 301.0]],
        ]
    )
    eigenvectors = np.array([[reference], [source]])

    result = track_phonon_modes(
        frequencies,
        eigenvectors,
        [10.0, 11.0],
        minimum_overlap=0.9,
    )

    assert result.status == "unreliable"
    assert result.fit.degree is None
    assert result.unresolved_assignments == 2


def test_frequency_fit_can_support_weak_overlap_after_tracking() -> None:
    """A smooth multi-volume branch may demote a weak overlap to a caution."""
    reference = np.eye(3, dtype=np.complex128).reshape(3, 1, 3)
    mixed = reference.copy()
    mixed[0, 0] = np.array([0.8, 0.6, 0.0])
    mixed[1, 0] = np.array([-0.6, 0.8, 0.0])
    eigenvectors = np.array(
        [[reference], [reference], [reference], [reference], [mixed]]
    )
    frequencies = np.array(
        [
            [[100.0, 200.0, 300.0]],
            [[101.0, 202.0, 303.0]],
            [[102.0, 204.0, 306.0]],
            [[103.0, 206.0, 309.0]],
            [[104.0, 208.0, 312.0]],
        ]
    )

    result = track_phonon_modes(
        frequencies,
        eigenvectors,
        [10.0, 11.0, 12.0, 13.0, 14.0],
        minimum_overlap=0.9,
    )

    assert result.verified
    assert result.fit.degree == 3
    assert result.fit.residual_degrees_of_freedom == 1
    assert result.fit.predictive_degree == 2
    assert result.fit.predictive_residual_degrees_of_freedom == 1
    assert result.low_overlap_assignments == 2
    assert result.caution_assignments >= 2
    assert result.unresolved_assignments == 0



def test_leave_one_out_rejects_bad_weak_overlap() -> None:
    """A low overlap must not be rescued by a fit containing its own bad point."""
    reference = np.eye(3, dtype=np.complex128).reshape(3, 1, 3)
    mixed = reference.copy()
    mixed[0, 0] = np.array([0.8, 0.6, 0.0])
    mixed[1, 0] = np.array([-0.6, 0.8, 0.0])
    eigenvectors = np.array(
        [[reference], [reference], [reference], [reference], [mixed]]
    )
    frequencies = np.array(
        [
            [[100.0, 200.0, 300.0]],
            [[101.0, 202.0, 303.0]],
            [[102.0, 204.0, 306.0]],
            [[103.0, 206.0, 309.0]],
            [[107.0, 214.0, 312.0]],
        ]
    )

    result = track_phonon_modes(
        frequencies,
        eigenvectors,
        [10.0, 11.0, 12.0, 13.0, 14.0],
        minimum_overlap=0.9,
    )

    assert np.all(result.fit.supported)
    assert result.fit.predictive_degree == 2
    assert result.fit.predictive_residual_degrees_of_freedom == 1
    assert result.fit.predictive_residuals[4, 0, 0] > 2.0
    assert not result.fit.predictive_supported[4, 0, 0]
    assert result.status == "unreliable"
    assert result.unresolved_assignments == 2


def test_local_pair_assignments_do_not_depend_on_reference_state() -> None:
    """Changing branch labels must not alter adjacent-volume assignments."""
    basis = np.eye(3, dtype=np.complex128).reshape(3, 1, 3)
    state0 = basis
    state1 = basis[[1, 2, 0]]
    state2 = basis[[2, 0, 1]]
    eigenvectors = np.array([[state0], [state1], [state2]])
    frequencies = np.array(
        [
            [[100.0, 200.0, 300.0]],
            [[200.0, 300.0, 100.0]],
            [[300.0, 100.0, 200.0]],
        ]
    )

    low_reference = track_phonon_modes(
        frequencies,
        eigenvectors,
        [10.0, 11.0, 12.0],
        reference_index=0,
    )
    high_reference = track_phonon_modes(
        frequencies,
        eigenvectors,
        [10.0, 11.0, 12.0],
        reference_index=2,
    )

    assert low_reference.volume_order == high_reference.volume_order == (0, 1, 2)
    for first, second in zip(low_reference.steps, high_reference.steps, strict=True):
        assert (first.predecessor_index, first.source_index, first.qpoint_index) == (
            second.predecessor_index,
            second.source_index,
            second.qpoint_index,
        )
        np.testing.assert_array_equal(first.permutation, second.permutation)
        np.testing.assert_allclose(first.overlaps, second.overlaps)


def test_mode_tracking_traverses_unsorted_volumes_from_reference() -> None:
    """Tracking should use adjacent volumes rather than input file order."""
    reference = np.eye(3, dtype=np.complex128).reshape(3, 1, 3)
    eigenvectors = np.array([[reference], [reference], [reference]])
    frequencies = np.array(
        [
            [[110.0, 210.0, 310.0]],
            [[90.0, 190.0, 290.0]],
            [[100.0, 200.0, 300.0]],
        ]
    )

    result = track_phonon_modes(
        frequencies,
        eigenvectors,
        [11.0, 9.0, 10.0],
        reference_index=2,
    )

    assert result.volume_order == (1, 2, 0)
    assert result.traversal == (2, 1, 0)
    assert result.verified


def test_dolomite_real_dataset_is_verified_with_cautions():
    root = Path(__file__).resolve().parents[3]
    files = sorted(
        (root / "examples/qha/crystal-phonons").glob(
            "dol_pbe0_crystal_pdisp_p0[1-7].out"
        )
    )
    assert len(files) == 7
    volume_pattern = re.compile(
        r"PRIMITIVE CELL .*?VOLUME=\s*([0-9.+\-DEde]+)"
    )
    volumes = []
    frequencies = []
    eigenvectors = []
    for filename in files:
        lines = filename.read_text(encoding="utf-8", errors="replace").splitlines()
        volume = next(
            float(match.group(1).replace("D", "E").replace("d", "e"))
            for line in lines
            if (match := volume_pattern.search(line)) is not None
        )
        modes = CrystalPhononModeParser(filename).parse(30)
        assert modes is not None
        volumes.append(volume)
        frequencies.append(modes.frequencies)
        eigenvectors.append(modes.eigenvectors)

    result = track_phonon_modes(
        np.asarray(frequencies),
        np.asarray(eigenvectors),
        np.asarray(volumes),
        reference_index=0,
    )

    assert result.status == "verified"
    assert result.caution_assignments == 274
    assert result.unresolved_assignments == 0
    assert result.low_overlap_assignments == 6
    assert result.ambiguous_assignments == 274
    assert result.reordered_assignments == 420
    assert result.local_reordered_assignments == 261
    assert result.fit.predictive_degree == 3
    assert result.fit.predictive_residual_degrees_of_freedom == 2
    np.testing.assert_allclose(
        result.minimum_overlap,
        0.4209753995491722,
        rtol=0.0,
        atol=1.0e-12,
    )
    assert result.fit.degree == 3
    assert result.fit.residual_degrees_of_freedom == 3
    assert np.count_nonzero(result.fit.supported) == 27 * 30
