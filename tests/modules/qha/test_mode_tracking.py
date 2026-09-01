
"""Tests for backend-neutral phonon-mode continuity tracking."""

from __future__ import annotations

import numpy as np

from quantas.core.numerics.phonon_tracking import track_phonon_modes



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


def test_mode_tracking_marks_unresolved_mixture_unreliable() -> None:
    """Non-degenerate modes with competing overlaps should not be verified."""
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

    assert result.status == "unreliable"
    assert result.ambiguous_assignments > 0


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

    assert result.traversal == (2, 1, 0)
    assert result.verified
