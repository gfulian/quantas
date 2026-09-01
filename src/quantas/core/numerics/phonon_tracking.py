# -*- coding: utf-8 -*-

"""Backend-neutral phonon-mode continuity tracking across sampled states."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import linear_sum_assignment


ModeTrackingStatus = Literal["verified", "unreliable"]


@dataclass(frozen=True, slots=True)
class PhononModeTrackingStep:
    """Diagnostics for one adjacent-state phonon-mode assignment.

    Parameters
    ----------
    predecessor_index, source_index : int
        Original state indices connected by this tracking step.
    qpoint_index : int
        Zero-based q-point index.
    permutation : ndarray
        Mapping from tracked branch index to the raw mode index in the source
        state.
    overlaps : ndarray
        Absolute scalar products for the selected mode assignments.
    ambiguous_mask : ndarray
        Mask of non-degenerate assignments whose scalar-product separation is
        smaller than the configured ambiguity margin.
    low_overlap_mask : ndarray
        Mask of non-degenerate assignments below the minimum overlap.
    degenerate_subspaces : tuple of tuple of int
        Tracked branch indices treated as degenerate eigenspaces.
    subspace_min_singular_values : tuple of float
        Minimum singular value for each degenerate-subspace overlap matrix.
    """

    predecessor_index: int
    source_index: int
    qpoint_index: int
    permutation: NDArray[np.int64]
    overlaps: NDArray[np.float64]
    ambiguous_mask: NDArray[np.bool_]
    low_overlap_mask: NDArray[np.bool_]
    degenerate_subspaces: tuple[tuple[int, ...], ...]
    subspace_min_singular_values: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class PhononModeTrackingResult:
    """Mode-continuous frequencies and diagnostics across sampled states.

    Parameters
    ----------
    frequencies : ndarray
        Reordered frequencies with shape ``(states, qpoints, modes)``.
    permutations : ndarray
        Mapping from tracked branches to raw source modes, with shape
        ``(states, qpoints, modes)``.
    overlaps : ndarray
        Selected absolute scalar products. The reference state contains NaN.
    status : {"verified", "unreliable"}
        Overall mode-continuity assessment.
    reference_index : int
        Original state index used to initialize branch labels.
    traversal : tuple of int
        Original state indices in the order in which tracking was performed.
    steps : tuple of PhononModeTrackingStep
        Per-q-point assignment diagnostics.
    """

    frequencies: NDArray[np.float64]
    permutations: NDArray[np.int64]
    overlaps: NDArray[np.float64]
    status: ModeTrackingStatus
    reference_index: int
    traversal: tuple[int, ...]
    steps: tuple[PhononModeTrackingStep, ...]

    @property
    def verified(self) -> bool:
        """Return whether all tracked assignments passed continuity checks."""
        return self.status == "verified"

    @property
    def ambiguous_assignments(self) -> int:
        """Return the number of unresolved scalar-product ambiguities."""
        return int(sum(np.count_nonzero(step.ambiguous_mask) for step in self.steps))

    @property
    def low_overlap_assignments(self) -> int:
        """Return the number of unresolved low-overlap assignments."""
        return int(sum(np.count_nonzero(step.low_overlap_mask) for step in self.steps))

    @property
    def reordered_assignments(self) -> int:
        """Return the number of selected raw mode indices differing from identity."""
        return int(
            sum(
                np.count_nonzero(
                    step.permutation != np.arange(step.permutation.size)
                )
                for step in self.steps
            )
        )

    @property
    def degenerate_subspaces(self) -> int:
        """Return the number of degenerate subspaces treated during tracking."""
        return int(sum(len(step.degenerate_subspaces) for step in self.steps))

    @property
    def minimum_overlap(self) -> float:
        """Return the minimum selected scalar product across tracked steps."""
        values = np.concatenate([step.overlaps for step in self.steps])
        return float(np.nanmin(values)) if values.size else float("nan")

    @property
    def minimum_subspace_singular_value(self) -> float:
        """Return the weakest resolved degenerate-subspace singular value."""
        values = [
            value
            for step in self.steps
            for value in step.subspace_min_singular_values
        ]
        return float(min(values)) if values else float("nan")


def track_phonon_modes(
    frequencies: ArrayLike,
    eigenvectors: ArrayLike,
    volumes: ArrayLike,
    *,
    reference_index: int = 0,
    ambiguity_margin: float = 0.4,
    minimum_overlap: float = 0.5,
    degeneracy_atol: float = 5.0e-2,
    degeneracy_rtol: float = 1.0e-6,
    minimum_subspace_overlap: float = 0.8,
) -> PhononModeTrackingResult:
    """Track phonon branches across volumes by eigenvector scalar products.

    Tracking starts from ``reference_index`` and proceeds through adjacent
    volumes independently toward compression and expansion. At each q point a
    global linear assignment maximizes absolute eigenvector overlaps. Rotated
    bases inside numerically degenerate eigenspaces are aligned by a unitary
    Procrustes transformation before the next volume is processed.

    Parameters
    ----------
    frequencies : array_like
        Frequencies with shape ``(states, qpoints, modes)``.
    eigenvectors : array_like
        Unit-norm complex eigenvectors with shape
        ``(states, qpoints, modes, atoms, 3)``.
    volumes : array_like
        One volume per state. Values need not be supplied in sorted order.
    reference_index : int, optional
        State fixing the branch labels.
    ambiguity_margin : float, optional
        Minimum separation between a selected scalar product and competing
        row/column products. ``0.4`` follows the diagnostic criterion described
        for CRYSTAL QHA mode-continuity analysis.
    minimum_overlap : float, optional
        Minimum acceptable scalar product for a non-degenerate assignment.
    degeneracy_atol, degeneracy_rtol : float, optional
        Absolute and relative frequency tolerances used to identify numerical
        degeneracies.
    minimum_subspace_overlap : float, optional
        Minimum singular value accepted when matching degenerate eigenspaces.

    Returns
    -------
    PhononModeTrackingResult
        Reordered frequencies, branch permutations, and continuity diagnostics.

    Raises
    ------
    ValueError
        If dimensions, state counts, norms, reference index, volumes, or
        tolerances are invalid.
    """
    freq = np.asarray(frequencies, dtype=np.float64)
    eig = np.asarray(eigenvectors, dtype=np.complex128)
    volume = np.asarray(volumes, dtype=np.float64)
    _validate_inputs(freq, eig, volume, reference_index)
    for name, value in (
        ("ambiguity_margin", ambiguity_margin),
        ("minimum_overlap", minimum_overlap),
        ("degeneracy_atol", degeneracy_atol),
        ("degeneracy_rtol", degeneracy_rtol),
        ("minimum_subspace_overlap", minimum_subspace_overlap),
    ):
        if not np.isfinite(value) or value < 0.0:
            raise ValueError(f"{name} must be finite and non-negative")
    if minimum_overlap > 1.0 or minimum_subspace_overlap > 1.0:
        raise ValueError("overlap thresholds cannot exceed one")

    nstates, nqpoints, nmodes = freq.shape
    sorted_indices = np.argsort(volume, kind="stable")
    sorted_position = int(np.flatnonzero(sorted_indices == reference_index)[0])
    lower = tuple(int(index) for index in sorted_indices[:sorted_position][::-1])
    upper = tuple(int(index) for index in sorted_indices[sorted_position + 1 :])
    traversal = (int(reference_index),) + lower + upper

    tracked_freq = np.empty_like(freq)
    permutations = np.empty((nstates, nqpoints, nmodes), dtype=np.int64)
    overlaps = np.full((nstates, nqpoints, nmodes), np.nan, dtype=np.float64)
    tracked_vectors = np.empty_like(eig)
    identity = np.arange(nmodes, dtype=np.int64)

    tracked_freq[reference_index] = freq[reference_index]
    tracked_vectors[reference_index] = eig[reference_index]
    permutations[reference_index] = identity[np.newaxis, :]

    steps: list[PhononModeTrackingStep] = []
    status: ModeTrackingStatus = "verified"

    for branch in (lower, upper):
        predecessor = int(reference_index)
        for source in branch:
            for qpoint in range(nqpoints):
                step, aligned = _track_pair(
                    tracked_freq[predecessor, qpoint],
                    tracked_vectors[predecessor, qpoint],
                    freq[source, qpoint],
                    eig[source, qpoint],
                    predecessor_index=predecessor,
                    source_index=source,
                    qpoint_index=qpoint,
                    ambiguity_margin=ambiguity_margin,
                    minimum_overlap=minimum_overlap,
                    degeneracy_atol=degeneracy_atol,
                    degeneracy_rtol=degeneracy_rtol,
                    minimum_subspace_overlap=minimum_subspace_overlap,
                )
                permutation = step.permutation
                tracked_freq[source, qpoint] = freq[source, qpoint, permutation]
                tracked_vectors[source, qpoint] = aligned
                permutations[source, qpoint] = permutation
                overlaps[source, qpoint] = step.overlaps
                if np.any(step.ambiguous_mask) or np.any(step.low_overlap_mask):
                    status = "unreliable"
                steps.append(step)
            predecessor = int(source)

    return PhononModeTrackingResult(
        frequencies=tracked_freq,
        permutations=permutations,
        overlaps=overlaps,
        status=status,
        reference_index=int(reference_index),
        traversal=traversal,
        steps=tuple(steps),
    )


def _track_pair(
    reference_frequencies: NDArray[np.float64],
    reference_vectors: NDArray[np.complex128],
    source_frequencies: NDArray[np.float64],
    source_vectors: NDArray[np.complex128],
    *,
    predecessor_index: int,
    source_index: int,
    qpoint_index: int,
    ambiguity_margin: float,
    minimum_overlap: float,
    degeneracy_atol: float,
    degeneracy_rtol: float,
    minimum_subspace_overlap: float,
) -> tuple[PhononModeTrackingStep, NDArray[np.complex128]]:
    """Track one q point between adjacent states."""
    nmodes = reference_frequencies.size
    reference_flat = reference_vectors.reshape(nmodes, -1)
    source_flat = source_vectors.reshape(nmodes, -1)
    matrix = np.abs(reference_flat.conj() @ source_flat.T)
    matrix = np.clip(matrix, 0.0, 1.0)
    rows, columns = linear_sum_assignment(-matrix)
    permutation = np.empty(nmodes, dtype=np.int64)
    permutation[rows] = columns
    selected = matrix[np.arange(nmodes), permutation]

    aligned = np.array(source_vectors[permutation], dtype=np.complex128, copy=True)
    components = _degenerate_components(
        reference_frequencies,
        source_frequencies[permutation],
        atol=degeneracy_atol,
        rtol=degeneracy_rtol,
    )
    resolved_mask = np.zeros(nmodes, dtype=bool)
    accepted_groups: list[tuple[int, ...]] = []
    singular_values: list[float] = []

    for component in components:
        indices = np.asarray(component, dtype=np.int64)
        current_indices = permutation[indices]
        overlap = reference_flat[indices].conj() @ source_flat[current_indices].T
        left, singular, right_transpose = np.linalg.svd(overlap)
        minimum = float(np.min(singular))
        singular_values.append(minimum)
        if minimum < minimum_subspace_overlap:
            continue
        rotation = left @ right_transpose
        aligned_flat = rotation @ source_flat[current_indices]
        aligned[indices] = aligned_flat.reshape(aligned[indices].shape)
        resolved_mask[indices] = True
        accepted_groups.append(tuple(int(index) for index in indices))

    ambiguous = np.zeros(nmodes, dtype=bool)
    low_overlap = np.zeros(nmodes, dtype=bool)
    for mode, source_mode in enumerate(permutation):
        if resolved_mask[mode]:
            continue
        selected_value = float(selected[mode])
        row_alternatives = np.delete(matrix[mode], source_mode)
        column_alternatives = np.delete(matrix[:, source_mode], mode)
        competitor = max(
            float(np.max(row_alternatives)) if row_alternatives.size else 0.0,
            float(np.max(column_alternatives)) if column_alternatives.size else 0.0,
        )
        ambiguous[mode] = selected_value - competitor < ambiguity_margin
        low_overlap[mode] = selected_value < minimum_overlap

    step = PhononModeTrackingStep(
        predecessor_index=int(predecessor_index),
        source_index=int(source_index),
        qpoint_index=int(qpoint_index),
        permutation=permutation,
        overlaps=np.asarray(selected, dtype=np.float64),
        ambiguous_mask=ambiguous,
        low_overlap_mask=low_overlap,
        degenerate_subspaces=tuple(accepted_groups),
        subspace_min_singular_values=tuple(singular_values),
    )
    return step, aligned


def _degenerate_components(
    reference: NDArray[np.float64],
    source: NDArray[np.float64],
    *,
    atol: float,
    rtol: float,
) -> tuple[tuple[int, ...], ...]:
    """Return connected mode groups degenerate in either adjacent state."""
    nmodes = reference.size
    adjacency: list[set[int]] = [set() for _ in range(nmodes)]
    for values in (reference, source):
        order = np.argsort(values, kind="stable")
        for left, right in zip(order[:-1], order[1:], strict=True):
            if _frequency_degenerate(
                values[left],
                values[right],
                atol=atol,
                rtol=rtol,
            ):
                first = int(left)
                second = int(right)
                adjacency[first].add(second)
                adjacency[second].add(first)

    components: list[tuple[int, ...]] = []
    visited: set[int] = set()
    for start in range(nmodes):
        if start in visited or not adjacency[start]:
            continue
        stack = [start]
        visited.add(start)
        component: list[int] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbour in adjacency[current]:
                if neighbour not in visited:
                    visited.add(neighbour)
                    stack.append(neighbour)
        components.append(tuple(sorted(component)))
    return tuple(components)


def _frequency_degenerate(
    first: float,
    second: float,
    *,
    atol: float,
    rtol: float,
) -> bool:
    """Return whether two frequencies belong to the same numerical manifold."""
    scale = max(abs(float(first)), abs(float(second)), 1.0)
    return abs(float(first) - float(second)) <= atol + rtol * scale


def _validate_inputs(
    frequencies: NDArray[np.float64],
    eigenvectors: NDArray[np.complex128],
    volumes: NDArray[np.float64],
    reference_index: int,
) -> None:
    """Validate array layouts and normalized eigenvectors."""
    if frequencies.ndim != 3:
        raise ValueError("frequencies must have shape (states, qpoints, modes)")
    if eigenvectors.ndim != 5 or eigenvectors.shape[-1] != 3:
        raise ValueError(
            "eigenvectors must have shape (states, qpoints, modes, atoms, 3)"
        )
    if eigenvectors.shape[:3] != frequencies.shape:
        raise ValueError("frequency and eigenvector dimensions are inconsistent")
    if volumes.ndim != 1 or volumes.shape[0] != frequencies.shape[0]:
        raise ValueError("volumes must contain one value per phonon state")
    if not 0 <= int(reference_index) < frequencies.shape[0]:
        raise ValueError("reference_index is outside the sampled states")
    if not np.all(np.isfinite(frequencies)) or not np.all(np.isfinite(volumes)):
        raise ValueError("frequencies and volumes must be finite")
    if np.unique(volumes).size != volumes.size:
        raise ValueError("sampled phonon volumes must be unique")
    norms = np.linalg.norm(eigenvectors.reshape(*frequencies.shape, -1), axis=-1)
    if not np.all(np.isfinite(norms)) or not np.allclose(
        norms, 1.0, rtol=1.0e-10, atol=1.0e-12
    ):
        raise ValueError("phonon eigenvectors must be normalized to unit norm")
