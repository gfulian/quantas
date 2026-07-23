# -*- coding: utf-8 -*-

"""Continuity tracking for axial acoustic-polarization fields."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .field import PhaseFieldResult
from .modes import (
    BRANCH_INDEX,
    MODE_INDEX,
    PolarizationBranch,
    WaveMode,
)


@dataclass(frozen=True, slots=True)
class PolarizationBranchField:
    """Tracked polarization data for one continuity branch.

    Parameters
    ----------
    branch : PolarizationBranch
        Continuity branch identifier.
    polarizations : ndarray
        Axial polarization field with shape ``(n_points, 3)``.
    local_mode_indices : ndarray
        Indices of locally ordered wave modes. A value of ``-1`` indicates
        that no unique local eigenvector can be assigned.
    sign_flip_mask : ndarray
        Boolean mask identifying sign reversals applied for continuity.
    continuity_scores : ndarray
        Absolute dot products with the preceding tracked axis. Segment starts
        contain ``nan``.
    resolved_mask : ndarray
        Boolean mask identifying uniquely resolved local polarization axes.
    segment_start_mask : ndarray
        Boolean mask identifying points initialized without a predecessor.
    """

    branch: PolarizationBranch
    polarizations: NDArray[np.float64]
    local_mode_indices: NDArray[np.int64]
    sign_flip_mask: NDArray[np.bool_]
    continuity_scores: NDArray[np.float64]
    resolved_mask: NDArray[np.bool_]
    segment_start_mask: NDArray[np.bool_]


@dataclass(frozen=True, slots=True)
class PolarizationTrackingResult:
    """Continuity representation of sampled acoustic polarizations.

    Parameters
    ----------
    directions : ndarray
        Unit wave-normal directions with shape ``(n_points, 3)``.
    polarizations : ndarray
        Tracked axial polarizations with shape ``(n_points, 3, 3)`` in branch
        order ``shear_a``, ``shear_b`` and ``p``.
    branch_mode_indices : ndarray
        Mapping from each continuity branch to the locally ordered acoustic
        mode. ``-1`` denotes a non-unique assignment in a degenerate subspace.
    sign_flip_mask : ndarray
        Boolean mask for axial sign reversals.
    continuity_scores : ndarray
        Absolute dot products with preceding tracked axes.
    resolved_mask : ndarray
        Boolean mask for uniquely resolved local axes.
    segment_start_mask : ndarray
        Boolean mask for branch initialization points.
    shear_swap_mask : ndarray
        Boolean mask indicating that the two local shear modes were swapped
        to preserve branch continuity.
    shear_permutation_ambiguous_mask : ndarray
        Boolean mask for points where the two shear permutations are not
        distinguishable under the selected overlap tolerance.
    shear_subspace_rotation_mask : ndarray
        Boolean mask for points where a degenerate shear basis was rotated
        inside its two-dimensional eigenspace for continuity.
    """

    directions: NDArray[np.float64]
    polarizations: NDArray[np.float64]
    branch_mode_indices: NDArray[np.int64]
    sign_flip_mask: NDArray[np.bool_]
    continuity_scores: NDArray[np.float64]
    resolved_mask: NDArray[np.bool_]
    segment_start_mask: NDArray[np.bool_]
    shear_swap_mask: NDArray[np.bool_]
    shear_permutation_ambiguous_mask: NDArray[np.bool_]
    shear_subspace_rotation_mask: NDArray[np.bool_]

    @property
    def n_points(self) -> int:
        """Return the number of tracked directions."""
        return int(self.directions.shape[0])

    def for_branch(
        self,
        branch: PolarizationBranch | str,
    ) -> PolarizationBranchField:
        """Return tracked data for one continuity branch.

        Parameters
        ----------
        branch : PolarizationBranch or str
            Continuity branch enum or its string value.

        Returns
        -------
        PolarizationBranchField
            Read-only views of the selected branch data.

        Raises
        ------
        ValueError
            If ``branch`` is not supported.
        """
        resolved = PolarizationBranch(branch)
        index = BRANCH_INDEX[resolved]
        return PolarizationBranchField(
            branch=resolved,
            polarizations=_readonly_float_view(self.polarizations[:, index, :]),
            local_mode_indices=_readonly_int_view(self.branch_mode_indices[:, index]),
            sign_flip_mask=_readonly_bool_view(self.sign_flip_mask[:, index]),
            continuity_scores=_readonly_float_view(self.continuity_scores[:, index]),
            resolved_mask=_readonly_bool_view(self.resolved_mask[:, index]),
            segment_start_mask=_readonly_bool_view(self.segment_start_mask[:, index]),
        )

    def mode_at(
        self,
        point_index: int,
        branch: PolarizationBranch | str,
    ) -> WaveMode | None:
        """Return the local wave mode assigned to one tracked branch.

        Parameters
        ----------
        point_index : int
            Sampled point index.
        branch : PolarizationBranch or str
            Continuity branch.

        Returns
        -------
        WaveMode or None
            Locally ordered mode, or ``None`` when the assignment is not
            unique.

        Raises
        ------
        IndexError
            If ``point_index`` is outside the sampled field.
        ValueError
            If ``branch`` is not supported.
        """
        if point_index < -self.n_points or point_index >= self.n_points:
            raise IndexError("Polarization point index is out of range.")
        resolved = PolarizationBranch(branch)
        local_index = int(self.branch_mode_indices[point_index, BRANCH_INDEX[resolved]])
        if local_index < 0:
            return None
        return (WaveMode.V_S2, WaveMode.V_S1, WaveMode.V_P)[local_index]


def align_axial_vector(
    reference: ArrayLike,
    vector: ArrayLike,
) -> NDArray[np.float64]:
    """Align an axial vector with a reference by selecting its sign.

    Parameters
    ----------
    reference, vector : array_like
        Non-zero finite Cartesian vectors with shape ``(3,)``.

    Returns
    -------
    ndarray
        Read-only copy of ``vector`` whose dot product with ``reference`` is
        non-negative. The vector magnitude is preserved.

    Raises
    ------
    ValueError
        If either vector is invalid.
    """
    aligned, _ = _align_axis(
        _validate_axis(reference, "reference"),
        _validate_axis(vector, "vector"),
    )
    aligned.setflags(write=False)
    return aligned


def track_polarizations(
    field: PhaseFieldResult,
    *,
    permutation_tolerance: float = 1.0e-10,
) -> PolarizationTrackingResult:
    """Track axial polarizations along an ordered directional field.

    The locally ordered modes ``V_S2``, ``V_S1`` and ``V_P`` remain unchanged
    in ``field``. The returned ``shear_a`` and ``shear_b`` branches are
    continuity labels that may exchange their local shear-mode assignment.

    Parameters
    ----------
    field : PhaseFieldResult
        Sampled phase field ordered along a path or another deterministic
        traversal.
    permutation_tolerance : float, optional
        Absolute tolerance used when comparing the direct and exchanged shear
        overlap scores.

    Returns
    -------
    PolarizationTrackingResult
        Tracked axial field, local-mode mapping and continuity diagnostics.

    Raises
    ------
    TypeError
        If ``field`` is not a :class:`PhaseFieldResult`.
    ValueError
        If ``permutation_tolerance`` is non-finite or negative.
    """
    if not isinstance(field, PhaseFieldResult):
        raise TypeError("track_polarizations requires a PhaseFieldResult.")
    tolerance = float(permutation_tolerance)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("permutation_tolerance must be finite and non-negative.")

    n_points = field.n_points
    tracked = np.empty_like(field.polarizations)
    branch_mode_indices = np.full((n_points, 3), -1, dtype=np.int64)
    sign_flip_mask = np.zeros((n_points, 3), dtype=bool)
    continuity_scores = np.full((n_points, 3), np.nan, dtype=float)
    resolved_mask = np.zeros((n_points, 3), dtype=bool)
    segment_start_mask = np.zeros((n_points, 3), dtype=bool)
    shear_swap_mask = np.zeros(n_points, dtype=bool)
    permutation_ambiguous_mask = np.zeros(n_points, dtype=bool)
    subspace_rotation_mask = np.zeros(n_points, dtype=bool)

    _initialize_point(
        field,
        0,
        tracked,
        branch_mode_indices,
        resolved_mask,
        segment_start_mask,
        permutation_ambiguous_mask,
    )

    for point in range(1, n_points):
        _track_shear_point(
            field,
            point,
            tracked,
            branch_mode_indices,
            sign_flip_mask,
            continuity_scores,
            resolved_mask,
            segment_start_mask,
            shear_swap_mask,
            permutation_ambiguous_mask,
            subspace_rotation_mask,
            tolerance,
        )
        _track_p_point(
            field,
            point,
            tracked,
            branch_mode_indices,
            sign_flip_mask,
            continuity_scores,
            resolved_mask,
            segment_start_mask,
        )

    arrays = (
        tracked,
        branch_mode_indices,
        sign_flip_mask,
        continuity_scores,
        resolved_mask,
        segment_start_mask,
        shear_swap_mask,
        permutation_ambiguous_mask,
        subspace_rotation_mask,
    )
    for array in arrays:
        array.setflags(write=False)

    return PolarizationTrackingResult(
        directions=field.directions,
        polarizations=tracked,
        branch_mode_indices=branch_mode_indices,
        sign_flip_mask=sign_flip_mask,
        continuity_scores=continuity_scores,
        resolved_mask=resolved_mask,
        segment_start_mask=segment_start_mask,
        shear_swap_mask=shear_swap_mask,
        shear_permutation_ambiguous_mask=permutation_ambiguous_mask,
        shear_subspace_rotation_mask=subspace_rotation_mask,
    )


def _initialize_point(
    field: PhaseFieldResult,
    point: int,
    tracked: NDArray[np.float64],
    branch_mode_indices: NDArray[np.int64],
    resolved_mask: NDArray[np.bool_],
    segment_start_mask: NDArray[np.bool_],
    permutation_ambiguous_mask: NDArray[np.bool_],
) -> None:
    """Initialize tracking at one point without a predecessor.

    Parameters
    ----------
    field : PhaseFieldResult
        Sampled phase field.
    point : int
        Point index.
    tracked, branch_mode_indices, resolved_mask, segment_start_mask : ndarray
        Mutable output arrays.
    permutation_ambiguous_mask : ndarray
        Mutable shear-permutation diagnostic.
    """
    tracked[point] = field.polarizations[point]
    branch_mode_indices[point] = [0, 1, 2]
    segment_start_mask[point] = True
    permutation_ambiguous_mask[point] = bool(field.pair_degeneracy_mask[point, 0])
    _apply_local_resolution(
        field,
        point,
        branch_mode_indices[point],
        resolved_mask[point],
    )


def _track_shear_point(
    field: PhaseFieldResult,
    point: int,
    tracked: NDArray[np.float64],
    branch_mode_indices: NDArray[np.int64],
    sign_flip_mask: NDArray[np.bool_],
    continuity_scores: NDArray[np.float64],
    resolved_mask: NDArray[np.bool_],
    segment_start_mask: NDArray[np.bool_],
    shear_swap_mask: NDArray[np.bool_],
    permutation_ambiguous_mask: NDArray[np.bool_],
    subspace_rotation_mask: NDArray[np.bool_],
    tolerance: float,
) -> None:
    """Track the two shear branches at one point.

    Parameters
    ----------
    field : PhaseFieldResult
        Sampled phase field.
    point : int
        Current point index.
    tracked, branch_mode_indices, sign_flip_mask, continuity_scores : ndarray
        Mutable tracking arrays.
    resolved_mask, segment_start_mask : ndarray
        Mutable diagnostic arrays.
    shear_swap_mask, permutation_ambiguous_mask : ndarray
        Mutable shear-permutation diagnostics.
    subspace_rotation_mask : ndarray
        Mutable degenerate-subspace diagnostic.
    tolerance : float
        Permutation-score tolerance.
    """
    current_valid = bool(np.all(field.valid_mask[point, :2]))
    previous_finite = bool(np.all(np.isfinite(tracked[point - 1, :2])))
    previous_valid = bool(np.all(field.valid_mask[point - 1, :2]))
    current_upper_degenerate = bool(field.pair_degeneracy_mask[point, 1])
    previous_upper_degenerate = bool(field.pair_degeneracy_mask[point - 1, 1])
    if (
        not current_valid
        or not previous_valid
        or not previous_finite
        or current_upper_degenerate
        or previous_upper_degenerate
    ):
        tracked[point, :2] = field.polarizations[point, :2]
        branch_mode_indices[point, :2] = [0, 1]
        segment_start_mask[point, :2] = True
        _apply_local_resolution(
            field,
            point,
            branch_mode_indices[point],
            resolved_mask[point],
        )
        return

    reference = tracked[point - 1, :2]
    current = field.polarizations[point, :2]
    shear_degenerate = bool(field.pair_degeneracy_mask[point, 0])

    if shear_degenerate:
        aligned_basis = _align_subspace(reference, current)
        tracked[point, :2] = aligned_basis
        branch_mode_indices[point, :2] = -1
        continuity_scores[point, :2] = np.abs(
            np.einsum("ij,ij->i", reference, aligned_basis)
        )
        permutation_ambiguous_mask[point] = True
        subspace_rotation_mask[point] = True
        return

    direct_score = _axis_overlap(reference[0], current[0]) + _axis_overlap(
        reference[1], current[1]
    )
    exchanged_score = _axis_overlap(reference[0], current[1]) + _axis_overlap(
        reference[1], current[0]
    )
    ambiguous = abs(direct_score - exchanged_score) <= tolerance
    permutation_ambiguous_mask[point] = ambiguous
    if exchanged_score > direct_score:
        local_indices = np.asarray([1, 0], dtype=np.int64)
        shear_swap_mask[point] = True
    else:
        local_indices = np.asarray([0, 1], dtype=np.int64)

    for branch in range(2):
        local_index = int(local_indices[branch])
        aligned, flipped = _align_axis(reference[branch], current[local_index])
        tracked[point, branch] = aligned
        branch_mode_indices[point, branch] = local_index
        sign_flip_mask[point, branch] = flipped
        continuity_scores[point, branch] = _axis_overlap(reference[branch], aligned)

    _apply_local_resolution(
        field,
        point,
        branch_mode_indices[point],
        resolved_mask[point],
    )
    if ambiguous:
        resolved_mask[point, :2] = False


def _track_p_point(
    field: PhaseFieldResult,
    point: int,
    tracked: NDArray[np.float64],
    branch_mode_indices: NDArray[np.int64],
    sign_flip_mask: NDArray[np.bool_],
    continuity_scores: NDArray[np.float64],
    resolved_mask: NDArray[np.bool_],
    segment_start_mask: NDArray[np.bool_],
) -> None:
    """Track the quasi-longitudinal branch at one point.

    Parameters
    ----------
    field : PhaseFieldResult
        Sampled phase field.
    point : int
        Current point index.
    tracked, branch_mode_indices, sign_flip_mask, continuity_scores : ndarray
        Mutable tracking arrays.
    resolved_mask, segment_start_mask : ndarray
        Mutable diagnostic arrays.
    """
    p_branch = BRANCH_INDEX[PolarizationBranch.P]
    p_mode = MODE_INDEX[WaveMode.V_P]
    current = field.polarizations[point, p_mode]
    current_valid = bool(field.valid_mask[point, p_mode])
    previous = tracked[point - 1, p_branch]
    previous_valid = bool(field.valid_mask[point - 1, p_mode])
    current_upper_degenerate = bool(field.pair_degeneracy_mask[point, 1])
    previous_upper_degenerate = bool(field.pair_degeneracy_mask[point - 1, 1])
    if (
        not current_valid
        or not previous_valid
        or not np.all(np.isfinite(previous))
        or current_upper_degenerate
        or previous_upper_degenerate
    ):
        tracked[point, p_branch] = current
        branch_mode_indices[point, p_branch] = p_mode
        segment_start_mask[point, p_branch] = True
        _apply_local_resolution(
            field,
            point,
            branch_mode_indices[point],
            resolved_mask[point],
        )
        return

    aligned, flipped = _align_axis(previous, current)
    tracked[point, p_branch] = aligned
    branch_mode_indices[point, p_branch] = p_mode
    sign_flip_mask[point, p_branch] = flipped
    continuity_scores[point, p_branch] = _axis_overlap(previous, aligned)
    _apply_local_resolution(
        field,
        point,
        branch_mode_indices[point],
        resolved_mask[point],
    )


def _apply_local_resolution(
    field: PhaseFieldResult,
    point: int,
    local_indices: NDArray[np.int64],
    output: NDArray[np.bool_],
) -> None:
    """Evaluate whether branch axes have unique local mode assignments.

    Parameters
    ----------
    field : PhaseFieldResult
        Sampled phase field.
    point : int
        Point index.
    local_indices : ndarray
        Local mode indices for each branch.
    output : ndarray
        Mutable boolean output with shape ``(3,)``.
    """
    shear_degenerate = bool(field.pair_degeneracy_mask[point, 0])
    upper_degenerate = bool(field.pair_degeneracy_mask[point, 1])
    for branch, local_index in enumerate(local_indices):
        if local_index < 0:
            output[branch] = False
            continue
        valid = bool(field.valid_mask[point, local_index])
        unique = True
        if local_index == 0:
            unique = not shear_degenerate
        elif local_index == 1:
            unique = not shear_degenerate and not upper_degenerate
        elif local_index == 2:
            unique = not upper_degenerate
        output[branch] = valid and unique
        if not unique:
            local_indices[branch] = -1


def _align_subspace(
    reference: NDArray[np.float64],
    current: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Align an orthonormal two-dimensional basis by Procrustes rotation.

    Parameters
    ----------
    reference, current : ndarray
        Row-oriented orthonormal bases with shape ``(2, 3)``.

    Returns
    -------
    ndarray
        Basis spanning the current subspace and aligned with ``reference``.
    """
    overlap = reference @ current.T
    left, _, right_transpose = np.linalg.svd(overlap)
    rotation = left @ right_transpose
    aligned = np.asarray(rotation @ current, dtype=float)
    for index in range(2):
        aligned[index], _ = _align_axis(reference[index], aligned[index])
    return aligned


def _align_axis(
    reference: NDArray[np.float64],
    vector: NDArray[np.float64],
) -> tuple[NDArray[np.float64], bool]:
    """Align one finite non-zero axis and report whether it was flipped.

    Parameters
    ----------
    reference, vector : ndarray
        Cartesian vectors with shape ``(3,)``.

    Returns
    -------
    aligned : ndarray
        Copy of ``vector`` with selected sign.
    flipped : bool
        Whether the sign was reversed.
    """
    flipped = bool(np.dot(reference, vector) < 0.0)
    aligned = np.array(-vector if flipped else vector, dtype=float, copy=True)
    return aligned, flipped


def _axis_overlap(
    first: NDArray[np.float64],
    second: NDArray[np.float64],
) -> float:
    """Return the normalized absolute overlap of two axes.

    Parameters
    ----------
    first, second : ndarray
        Non-zero Cartesian vectors.

    Returns
    -------
    float
        Absolute normalized dot product in the interval ``[0, 1]``.
    """
    denominator = float(np.linalg.norm(first) * np.linalg.norm(second))
    if denominator == 0.0 or not np.isfinite(denominator):
        return float("nan")
    value = abs(float(np.dot(first, second))) / denominator
    return float(np.clip(value, 0.0, 1.0))


def _validate_axis(axis: ArrayLike, name: str) -> NDArray[np.float64]:
    """Validate a finite non-zero Cartesian axis.

    Parameters
    ----------
    axis : array_like
        Cartesian vector.
    name : str
        Parameter name used in errors.

    Returns
    -------
    ndarray
        Floating-point vector with shape ``(3,)``.

    Raises
    ------
    ValueError
        If the axis has invalid shape, values or norm.
    """
    vector = np.asarray(axis, dtype=float)
    if vector.shape != (3,):
        raise ValueError(f"{name} must have shape (3,).")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain finite values.")
    if float(np.linalg.norm(vector)) == 0.0:
        raise ValueError(f"{name} must be non-zero.")
    return vector


def _readonly_float_view(
    array: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Return a read-only floating-point view.

    Parameters
    ----------
    array : ndarray
        Source array.

    Returns
    -------
    ndarray
        Read-only view.
    """
    view = array.view()
    view.setflags(write=False)
    return view


def _readonly_bool_view(
    array: NDArray[np.bool_],
) -> NDArray[np.bool_]:
    """Return a read-only boolean view.

    Parameters
    ----------
    array : ndarray
        Source array.

    Returns
    -------
    ndarray
        Read-only view.
    """
    view = array.view()
    view.setflags(write=False)
    return view


def _readonly_int_view(
    array: NDArray[np.int64],
) -> NDArray[np.int64]:
    """Return a read-only integer view.

    Parameters
    ----------
    array : ndarray
        Source array.

    Returns
    -------
    ndarray
        Read-only view.
    """
    view = array.view()
    view.setflags(write=False)
    return view
