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
class PhononModeFitDiagnostics:
    """Polynomial diagnostics for tracked phonon-frequency branches.

    Parameters
    ----------
    degree : int or None
        Polynomial degree used for every branch. ``None`` means that too few
        volume points were available for an independent continuity check.
    residual_degrees_of_freedom : int
        Number of residual degrees of freedom left after fitting one branch.
    r_squared : ndarray
        Coefficient of determination for each ``(qpoint, branch)``.
    rmse : ndarray
        Root-mean-square residual in the same units as the input frequencies.
    max_residual : ndarray
        Maximum absolute residual in the same units as the input frequencies.
    supported : ndarray
        Boolean mask identifying branches whose global fitted frequency path is
        smooth according to the configured descriptive thresholds.  This mask
        is diagnostic only and does not rescue low-overlap assignments.
    predictive_degree : int or None
        Polynomial degree used by leave-one-out predictions. ``None`` means
        that too few volume points were available for an independent
        prediction.
    predictive_residual_degrees_of_freedom : int
        Residual degrees of freedom in each leave-one-out training fit.
    predictive_residuals : ndarray
        Absolute leave-one-out prediction residuals with shape
        ``(states, qpoints, branches)``.
    predictive_tolerances : ndarray
        State-wise acceptance limits for leave-one-out residuals, in the same
        units as the input frequencies.
    predictive_supported : ndarray
        Boolean mask identifying state/branch values independently supported
        by leave-one-out prediction.
    """

    degree: int | None
    residual_degrees_of_freedom: int
    r_squared: NDArray[np.float64]
    rmse: NDArray[np.float64]
    max_residual: NDArray[np.float64]
    supported: NDArray[np.bool_]
    predictive_degree: int | None
    predictive_residual_degrees_of_freedom: int
    predictive_residuals: NDArray[np.float64]
    predictive_tolerances: NDArray[np.float64]
    predictive_supported: NDArray[np.bool_]


@dataclass(frozen=True, slots=True)
class PhononModeTrackingStep:
    """Diagnostics for one adjacent-volume phonon-mode assignment.

    The step always connects two states that are adjacent after sorting by
    volume. The local assignment is therefore independent of the state chosen
    to label the final phonon branches.

    Parameters
    ----------
    predecessor_index, source_index : int
        Original state indices for the lower- and higher-volume states.
    qpoint_index : int
        Zero-based q-point index.
    permutation : ndarray
        Mapping from the raw mode index in the lower-volume state to the raw
        mode index in the higher-volume state.
    branch_indices : ndarray
        Mapping from each raw lower-volume mode to the final branch label fixed
        by ``reference_index``.
    overlaps : ndarray
        Absolute scalar products for the selected local assignments.
    competitor_overlaps : ndarray
        Best competing scalar product for each selected assignment.
    overlap_gaps : ndarray
        Difference between selected and best competing scalar products.
    ambiguous_mask : ndarray
        Assignments whose scalar-product separation is smaller than the
        configured ambiguity margin.
    low_overlap_mask : ndarray
        Non-degenerate assignments below the minimum overlap threshold.
    caution_mask : ndarray
        Assignments retained as scientifically usable but deserving inspection.
    unresolved_mask : ndarray
        Assignments for which continuity could not be established reliably.
    degenerate_mask : ndarray
        Modes treated as members of a degenerate eigenspace.
    degenerate_subspaces : tuple of tuple of int
        Lower-volume raw mode indices for successfully matched degenerate
        eigenspaces.
    unresolved_degenerate_subspaces : tuple of tuple of int
        Lower-volume raw mode indices for degenerate eigenspaces that failed
        the configured subspace-overlap criterion.
    subspace_min_singular_values : tuple of float
        Minimum singular value for every detected degenerate-subspace overlap.
    """

    predecessor_index: int
    source_index: int
    qpoint_index: int
    permutation: NDArray[np.int64]
    branch_indices: NDArray[np.int64]
    overlaps: NDArray[np.float64]
    competitor_overlaps: NDArray[np.float64]
    overlap_gaps: NDArray[np.float64]
    ambiguous_mask: NDArray[np.bool_]
    low_overlap_mask: NDArray[np.bool_]
    caution_mask: NDArray[np.bool_]
    unresolved_mask: NDArray[np.bool_]
    degenerate_mask: NDArray[np.bool_]
    degenerate_subspaces: tuple[tuple[int, ...], ...]
    unresolved_degenerate_subspaces: tuple[tuple[int, ...], ...]
    subspace_min_singular_values: tuple[float, ...]
    degenerate_subspace_min_singular_values: tuple[float, ...]
    unresolved_degenerate_subspace_min_singular_values: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class PhononModeTrackingResult:
    """Mode-continuous frequencies and diagnostics across sampled states.

    Parameters
    ----------
    frequencies : ndarray
        Reordered frequencies with shape ``(states, qpoints, modes)``.
    permutations : ndarray
        Mapping from final branch labels to raw source modes, with shape
        ``(states, qpoints, modes)``.
    overlaps : ndarray
        Selected adjacent-volume scalar products expressed in final branch
        order. The reference state contains NaN.
    status : {"verified", "unreliable"}
        Overall mode-continuity assessment. Cautions do not make a dataset
        unreliable; unresolved assignments do.
    reference_index : int
        Original state index used only to label the final phonon branches.
    volume_order : tuple of int
        Original state indices sorted from smallest to largest volume.
    traversal : tuple of int
        Reference-centered ordering retained for compact provenance and
        backward compatibility with the first tracking implementation.
    steps : tuple of PhononModeTrackingStep
        Local adjacent-volume assignment diagnostics.
    fit : PhononModeFitDiagnostics
        Frequency-path diagnostics evaluated after eigenvector tracking.
    """

    frequencies: NDArray[np.float64]
    permutations: NDArray[np.int64]
    overlaps: NDArray[np.float64]
    status: ModeTrackingStatus
    reference_index: int
    volume_order: tuple[int, ...]
    traversal: tuple[int, ...]
    steps: tuple[PhononModeTrackingStep, ...]
    fit: PhononModeFitDiagnostics

    @property
    def verified(self) -> bool:
        """Return whether no mode assignment remains unresolved."""
        return self.status == "verified"

    @property
    def ambiguous_assignments(self) -> int:
        """Return the number of scalar-product ambiguities detected."""
        return int(sum(np.count_nonzero(step.ambiguous_mask) for step in self.steps))

    @property
    def low_overlap_assignments(self) -> int:
        """Return the number of non-degenerate assignments below threshold."""
        return int(sum(np.count_nonzero(step.low_overlap_mask) for step in self.steps))

    @property
    def caution_assignments(self) -> int:
        """Return the number of usable assignments that deserve inspection."""
        return int(sum(np.count_nonzero(step.caution_mask) for step in self.steps))

    @property
    def unresolved_assignments(self) -> int:
        """Return the number of assignments whose continuity is unresolved."""
        return int(sum(np.count_nonzero(step.unresolved_mask) for step in self.steps))

    @property
    def reordered_assignments(self) -> int:
        """Return final branch-to-raw assignments differing from identity."""
        identity = np.arange(self.permutations.shape[-1], dtype=np.int64)
        return int(np.count_nonzero(self.permutations != identity))

    @property
    def local_reordered_assignments(self) -> int:
        """Return adjacent-volume raw-mode assignments differing from identity."""
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
        """Return the number of resolved degenerate subspaces."""
        return int(sum(len(step.degenerate_subspaces) for step in self.steps))

    @property
    def unresolved_degenerate_subspaces(self) -> int:
        """Return the number of degenerate subspaces that failed matching."""
        return int(
            sum(len(step.unresolved_degenerate_subspaces) for step in self.steps)
        )

    @property
    def minimum_overlap(self) -> float:
        """Return the minimum selected non-degenerate scalar product."""
        values = [
            step.overlaps[~step.degenerate_mask]
            for step in self.steps
            if np.any(~step.degenerate_mask)
        ]
        if not values:
            return float("nan")
        merged = np.concatenate(values)
        return float(np.nanmin(merged)) if merged.size else float("nan")

    @property
    def minimum_subspace_singular_value(self) -> float:
        """Return the weakest detected degenerate-subspace singular value."""
        values = [
            value
            for step in self.steps
            for value in step.subspace_min_singular_values
        ]
        return float(min(values)) if values else float("nan")


@dataclass(frozen=True, slots=True)
class _LocalPairResult:
    """Internal local assignment before final branch labels are composed."""

    predecessor_index: int
    source_index: int
    qpoint_index: int
    permutation: NDArray[np.int64]
    overlaps: NDArray[np.float64]
    competitor_overlaps: NDArray[np.float64]
    overlap_gaps: NDArray[np.float64]
    ambiguous_mask: NDArray[np.bool_]
    low_overlap_mask: NDArray[np.bool_]
    degenerate_mask: NDArray[np.bool_]
    subspace_failed_mask: NDArray[np.bool_]
    degenerate_subspaces: tuple[tuple[int, ...], ...]
    unresolved_degenerate_subspaces: tuple[tuple[int, ...], ...]
    subspace_min_singular_values: tuple[float, ...]
    degenerate_subspace_min_singular_values: tuple[float, ...]
    unresolved_degenerate_subspace_min_singular_values: tuple[float, ...]


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
    fit_min_r_squared: float = 0.98,
    fit_max_rmse: float = 2.0,
    predictive_max_residual: float = 2.0,
    predictive_max_range_fraction: float = 0.10,
) -> PhononModeTrackingResult:
    """Track phonon branches across volumes by eigenvector scalar products.

    Every adjacent-volume pair is matched independently after sorting the
    sampled states by volume. The selected ``reference_index`` is used only to
    label the connected branches; it does not alter any local overlap matrix or
    local Hungarian assignment. Numerical degeneracies are evaluated as
    eigenspaces through singular values rather than as individual eigenvectors.

    Weak individual overlaps are retained as cautions only when both endpoints
    of the local assignment are independently supported by leave-one-out
    polynomial prediction of the tracked ``nu(V)`` branch. Frequency fitting
    never participates in the initial mode assignment.

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
        State fixing the final branch labels only.
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
    fit_min_r_squared : float, optional
        Minimum coefficient of determination used by the descriptive global
        branch-fit diagnostic.  It does not decide low-overlap acceptance.
    fit_max_rmse : float, optional
        Maximum global-fit RMSE, in the input frequency unit, used only for
        descriptive branch diagnostics.
    predictive_max_residual : float, optional
        Absolute leave-one-out residual floor, in the input frequency unit,
        used when validating a low-overlap assignment.
    predictive_max_range_fraction : float, optional
        Relative leave-one-out residual limit expressed as a fraction of the
        training-branch frequency range. The effective predictive tolerance is
        the larger of this relative limit and ``predictive_max_residual``.

    Returns
    -------
    PhononModeTrackingResult
        Reordered frequencies, branch permutations, fit diagnostics, and
        continuity diagnostics.

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
        ("fit_min_r_squared", fit_min_r_squared),
        ("fit_max_rmse", fit_max_rmse),
        ("predictive_max_residual", predictive_max_residual),
        ("predictive_max_range_fraction", predictive_max_range_fraction),
    ):
        if not np.isfinite(value) or value < 0.0:
            raise ValueError(f"{name} must be finite and non-negative")
    if minimum_overlap > 1.0 or minimum_subspace_overlap > 1.0:
        raise ValueError("overlap thresholds cannot exceed one")
    if fit_min_r_squared > 1.0:
        raise ValueError("fit_min_r_squared cannot exceed one")
    if predictive_max_range_fraction > 1.0:
        raise ValueError("predictive_max_range_fraction cannot exceed one")

    nstates, nqpoints, nmodes = freq.shape
    volume_order = tuple(int(index) for index in np.argsort(volume, kind="stable"))
    sorted_position = volume_order.index(int(reference_index))
    lower = tuple(reversed(volume_order[:sorted_position]))
    upper = volume_order[sorted_position + 1 :]
    traversal = (int(reference_index),) + lower + upper

    local_steps: dict[tuple[int, int, int], _LocalPairResult] = {}
    for lower_index, upper_index in zip(volume_order[:-1], volume_order[1:], strict=True):
        for qpoint in range(nqpoints):
            local_steps[(lower_index, upper_index, qpoint)] = _track_pair(
                freq[lower_index, qpoint],
                eig[lower_index, qpoint],
                freq[upper_index, qpoint],
                eig[upper_index, qpoint],
                predecessor_index=lower_index,
                source_index=upper_index,
                qpoint_index=qpoint,
                ambiguity_margin=ambiguity_margin,
                minimum_overlap=minimum_overlap,
                degeneracy_atol=degeneracy_atol,
                degeneracy_rtol=degeneracy_rtol,
                minimum_subspace_overlap=minimum_subspace_overlap,
            )

    permutations = np.empty((nstates, nqpoints, nmodes), dtype=np.int64)
    identity = np.arange(nmodes, dtype=np.int64)
    permutations[reference_index] = identity[np.newaxis, :]

    for qpoint in range(nqpoints):
        current = permutations[reference_index, qpoint]
        for source_index in upper:
            predecessor_index = volume_order[volume_order.index(source_index) - 1]
            local = local_steps[(predecessor_index, source_index, qpoint)]
            current = local.permutation[current]
            permutations[source_index, qpoint] = current

        current = permutations[reference_index, qpoint]
        for source_index in lower:
            upper_position = volume_order.index(source_index) + 1
            predecessor_index = volume_order[upper_position]
            local = local_steps[(source_index, predecessor_index, qpoint)]
            inverse = _inverse_permutation(local.permutation)
            current = inverse[current]
            permutations[source_index, qpoint] = current

    tracked_freq = np.take_along_axis(freq, permutations, axis=2)
    predictive_targets: set[tuple[int, int, int]] = set()
    for lower_index, upper_index in zip(
        volume_order[:-1], volume_order[1:], strict=True
    ):
        for qpoint in range(nqpoints):
            local = local_steps[(lower_index, upper_index, qpoint)]
            raw_to_branch = _inverse_permutation(permutations[lower_index, qpoint])
            for raw_mode in np.flatnonzero(local.low_overlap_mask):
                branch = int(raw_to_branch[raw_mode])
                predictive_targets.add((lower_index, qpoint, branch))
                predictive_targets.add((upper_index, qpoint, branch))

    fit = _fit_frequency_branches(
        tracked_freq,
        volume,
        volume_order,
        predictive_targets=predictive_targets,
        fit_min_r_squared=fit_min_r_squared,
        fit_max_rmse=fit_max_rmse,
        predictive_max_residual=predictive_max_residual,
        predictive_max_range_fraction=predictive_max_range_fraction,
    )

    steps: list[PhononModeTrackingStep] = []
    for lower_index, upper_index in zip(volume_order[:-1], volume_order[1:], strict=True):
        for qpoint in range(nqpoints):
            local = local_steps[(lower_index, upper_index, qpoint)]
            raw_to_branch = _inverse_permutation(permutations[lower_index, qpoint])
            branch_indices = raw_to_branch.copy()
            unresolved = np.asarray(local.subspace_failed_mask, dtype=bool).copy()
            caution = np.zeros(nmodes, dtype=bool)

            for raw_mode in range(nmodes):
                if local.degenerate_mask[raw_mode]:
                    continue
                branch = int(branch_indices[raw_mode])
                if local.low_overlap_mask[raw_mode]:
                    lower_supported = fit.predictive_supported[
                        lower_index, qpoint, branch
                    ]
                    upper_supported = fit.predictive_supported[
                        upper_index, qpoint, branch
                    ]
                    if lower_supported and upper_supported:
                        caution[raw_mode] = True
                    else:
                        unresolved[raw_mode] = True
                if local.ambiguous_mask[raw_mode] and not unresolved[raw_mode]:
                    caution[raw_mode] = True
            caution &= ~unresolved

            steps.append(
                PhononModeTrackingStep(
                    predecessor_index=lower_index,
                    source_index=upper_index,
                    qpoint_index=qpoint,
                    permutation=local.permutation,
                    branch_indices=branch_indices,
                    overlaps=local.overlaps,
                    competitor_overlaps=local.competitor_overlaps,
                    overlap_gaps=local.overlap_gaps,
                    ambiguous_mask=local.ambiguous_mask,
                    low_overlap_mask=local.low_overlap_mask,
                    caution_mask=caution,
                    unresolved_mask=unresolved,
                    degenerate_mask=local.degenerate_mask,
                    degenerate_subspaces=local.degenerate_subspaces,
                    unresolved_degenerate_subspaces=(
                        local.unresolved_degenerate_subspaces
                    ),
                    subspace_min_singular_values=(
                        local.subspace_min_singular_values
                    ),
                    degenerate_subspace_min_singular_values=(
                        local.degenerate_subspace_min_singular_values
                    ),
                    unresolved_degenerate_subspace_min_singular_values=(
                        local.unresolved_degenerate_subspace_min_singular_values
                    ),
                )
            )

    status: ModeTrackingStatus = (
        "unreliable"
        if any(np.any(step.unresolved_mask) for step in steps)
        else "verified"
    )
    overlaps = _reference_centered_overlaps(
        tuple(steps),
        permutations,
        volume_order,
        reference_index,
    )
    return PhononModeTrackingResult(
        frequencies=tracked_freq,
        permutations=permutations,
        overlaps=overlaps,
        status=status,
        reference_index=int(reference_index),
        volume_order=volume_order,
        traversal=traversal,
        steps=tuple(steps),
        fit=fit,
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
) -> _LocalPairResult:
    """Track one q point between two adjacent raw volume states."""
    nmodes = reference_frequencies.size
    reference_flat = reference_vectors.reshape(nmodes, -1)
    source_flat = source_vectors.reshape(nmodes, -1)
    matrix = np.abs(reference_flat.conj() @ source_flat.T)
    matrix = np.clip(matrix, 0.0, 1.0)
    rows, columns = linear_sum_assignment(-matrix)
    permutation = np.empty(nmodes, dtype=np.int64)
    permutation[rows] = columns
    selected = matrix[np.arange(nmodes), permutation]

    components = _degenerate_components(
        reference_frequencies,
        source_frequencies[permutation],
        atol=degeneracy_atol,
        rtol=degeneracy_rtol,
    )
    degenerate_mask = np.zeros(nmodes, dtype=bool)
    subspace_failed_mask = np.zeros(nmodes, dtype=bool)
    accepted_groups: list[tuple[int, ...]] = []
    failed_groups: list[tuple[int, ...]] = []
    singular_values: list[float] = []
    accepted_singular_values: list[float] = []
    failed_singular_values: list[float] = []

    for component in components:
        indices = np.asarray(component, dtype=np.int64)
        current_indices = permutation[indices]
        overlap = reference_flat[indices].conj() @ source_flat[current_indices].T
        singular = np.linalg.svd(overlap, compute_uv=False)
        minimum = float(np.min(singular))
        singular_values.append(minimum)
        degenerate_mask[indices] = True
        group = tuple(int(index) for index in indices)
        if minimum < minimum_subspace_overlap:
            subspace_failed_mask[indices] = True
            failed_groups.append(group)
            failed_singular_values.append(minimum)
        else:
            accepted_groups.append(group)
            accepted_singular_values.append(minimum)

    competitor = np.zeros(nmodes, dtype=np.float64)
    ambiguous = np.zeros(nmodes, dtype=bool)
    low_overlap = np.zeros(nmodes, dtype=bool)
    for mode, source_mode in enumerate(permutation):
        if degenerate_mask[mode]:
            continue
        selected_value = float(selected[mode])
        row_alternatives = np.delete(matrix[mode], source_mode)
        column_alternatives = np.delete(matrix[:, source_mode], mode)
        competitor_value = max(
            float(np.max(row_alternatives)) if row_alternatives.size else 0.0,
            float(np.max(column_alternatives)) if column_alternatives.size else 0.0,
        )
        competitor[mode] = competitor_value
        ambiguous[mode] = selected_value - competitor_value < ambiguity_margin
        low_overlap[mode] = selected_value < minimum_overlap

    return _LocalPairResult(
        predecessor_index=int(predecessor_index),
        source_index=int(source_index),
        qpoint_index=int(qpoint_index),
        permutation=permutation,
        overlaps=np.asarray(selected, dtype=np.float64),
        competitor_overlaps=competitor,
        overlap_gaps=np.asarray(selected - competitor, dtype=np.float64),
        ambiguous_mask=ambiguous,
        low_overlap_mask=low_overlap,
        degenerate_mask=degenerate_mask,
        subspace_failed_mask=subspace_failed_mask,
        degenerate_subspaces=tuple(accepted_groups),
        unresolved_degenerate_subspaces=tuple(failed_groups),
        subspace_min_singular_values=tuple(singular_values),
        degenerate_subspace_min_singular_values=tuple(accepted_singular_values),
        unresolved_degenerate_subspace_min_singular_values=tuple(
            failed_singular_values
        ),
    )


def _fit_frequency_branches(
    frequencies: NDArray[np.float64],
    volumes: NDArray[np.float64],
    volume_order: tuple[int, ...],
    *,
    predictive_targets: set[tuple[int, int, int]],
    fit_min_r_squared: float,
    fit_max_rmse: float,
    predictive_max_residual: float,
    predictive_max_range_fraction: float,
) -> PhononModeFitDiagnostics:
    """Fit global and leave-one-out ``nu(V)`` continuity diagnostics."""
    nstates, nqpoints, nmodes = frequencies.shape
    shape = (nqpoints, nmodes)
    predictive_shape = (nstates, nqpoints, nmodes)
    r_squared = np.full(shape, np.nan, dtype=np.float64)
    rmse = np.full(shape, np.nan, dtype=np.float64)
    max_residual = np.full(shape, np.nan, dtype=np.float64)
    supported = np.zeros(shape, dtype=bool)
    predictive_residuals = np.full(predictive_shape, np.nan, dtype=np.float64)
    predictive_tolerances = np.full(predictive_shape, np.nan, dtype=np.float64)
    predictive_supported = np.zeros(predictive_shape, dtype=bool)

    predictive_degree: int | None = None
    predictive_residual_dof = 0
    if nstates >= 4:
        predictive_degree = min(3, nstates - 3)
        predictive_residual_dof = (nstates - 1) - (predictive_degree + 1)

    if nstates < 3:
        return PhononModeFitDiagnostics(
            degree=None,
            residual_degrees_of_freedom=0,
            r_squared=r_squared,
            rmse=rmse,
            max_residual=max_residual,
            supported=supported,
            predictive_degree=predictive_degree,
            predictive_residual_degrees_of_freedom=predictive_residual_dof,
            predictive_residuals=predictive_residuals,
            predictive_tolerances=predictive_tolerances,
            predictive_supported=predictive_supported,
        )

    degree = min(3, nstates - 2)
    residual_dof = nstates - (degree + 1)
    ordered = np.asarray(volume_order, dtype=np.int64)
    x = volumes[ordered]
    span = float(np.ptp(x))
    if span <= 0.0:
        raise ValueError("sampled phonon volumes must span a finite interval")
    x_scaled = (x - float(np.mean(x))) / span

    for qpoint in range(nqpoints):
        for branch in range(nmodes):
            values = frequencies[ordered, qpoint, branch]
            coefficients = np.polyfit(x_scaled, values, degree)
            fitted = np.polyval(coefficients, x_scaled)
            residuals = values - fitted
            ss_res = float(np.sum(residuals**2))
            centered = values - float(np.mean(values))
            ss_tot = float(np.sum(centered**2))
            current_rmse = float(np.sqrt(np.mean(residuals**2)))
            current_max = float(np.max(np.abs(residuals)))
            current_r2 = (
                float(1.0 - ss_res / ss_tot)
                if ss_tot > np.finfo(np.float64).eps
                else float("nan")
            )
            r_squared[qpoint, branch] = current_r2
            rmse[qpoint, branch] = current_rmse
            max_residual[qpoint, branch] = current_max
            supported[qpoint, branch] = (
                (np.isfinite(current_r2) and current_r2 >= fit_min_r_squared)
                or current_rmse <= fit_max_rmse
            )

            if predictive_degree is None:
                continue
            for target_position, target_state in enumerate(ordered):
                target = (int(target_state), qpoint, branch)
                if target not in predictive_targets:
                    continue
                mask = np.ones(nstates, dtype=bool)
                mask[target_position] = False
                training_x = x[mask]
                training_values = values[mask]
                training_span = float(np.ptp(training_x))
                if training_span <= 0.0:
                    continue
                training_center = float(np.mean(training_x))
                training_scaled = (training_x - training_center) / training_span
                target_scaled = (
                    float(x[target_position]) - training_center
                ) / training_span
                coefficients = np.polyfit(
                    training_scaled,
                    training_values,
                    predictive_degree,
                )
                predicted = float(np.polyval(coefficients, target_scaled))
                residual = abs(float(values[target_position]) - predicted)
                branch_range = float(np.ptp(training_values))
                tolerance = max(
                    predictive_max_residual,
                    predictive_max_range_fraction * branch_range,
                )
                predictive_residuals[target_state, qpoint, branch] = residual
                predictive_tolerances[target_state, qpoint, branch] = tolerance
                predictive_supported[target_state, qpoint, branch] = (
                    residual <= tolerance
                )

    return PhononModeFitDiagnostics(
        degree=degree,
        residual_degrees_of_freedom=residual_dof,
        r_squared=r_squared,
        rmse=rmse,
        max_residual=max_residual,
        supported=supported,
        predictive_degree=predictive_degree,
        predictive_residual_degrees_of_freedom=predictive_residual_dof,
        predictive_residuals=predictive_residuals,
        predictive_tolerances=predictive_tolerances,
        predictive_supported=predictive_supported,
    )


def _reference_centered_overlaps(
    steps: tuple[PhononModeTrackingStep, ...],
    permutations: NDArray[np.int64],
    volume_order: tuple[int, ...],
    reference_index: int,
) -> NDArray[np.float64]:
    """Express adjacent-volume overlaps in final branch order by state."""
    nstates, nqpoints, nmodes = permutations.shape
    overlaps = np.full((nstates, nqpoints, nmodes), np.nan, dtype=np.float64)
    lookup = {
        (step.predecessor_index, step.source_index, step.qpoint_index): step
        for step in steps
    }
    ref_position = volume_order.index(int(reference_index))

    for position, state in enumerate(volume_order):
        if state == reference_index:
            continue
        if position < ref_position:
            lower_index = state
            upper_index = volume_order[position + 1]
        else:
            lower_index = volume_order[position - 1]
            upper_index = state
        for qpoint in range(nqpoints):
            step = lookup[(lower_index, upper_index, qpoint)]
            raw_lower = permutations[lower_index, qpoint]
            overlaps[state, qpoint] = step.overlaps[raw_lower]
    return overlaps


def _inverse_permutation(permutation: NDArray[np.int64]) -> NDArray[np.int64]:
    """Return the inverse of one one-dimensional integer permutation."""
    inverse = np.empty_like(permutation)
    inverse[permutation] = np.arange(permutation.size, dtype=np.int64)
    return inverse


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
