# -*- coding: utf-8 -*-

"""Batched sampling of acoustic-wave fields on spherical grids."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum

import numpy as np
from numpy.typing import NDArray

from quantas.core.geometry import SphericalGrid

from .christoffel import ChristoffelSolver
from .field import (
    EnhancementFieldResult,
    GroupFieldResult,
    PhaseFieldResult,
)
from .tracking import PolarizationTrackingResult, track_polarizations


ProgressCallback = Callable[[int, int], None]


class SamplingLevel(str, Enum):
    """Highest acoustic quantity evaluated by a sampling operation."""

    PHASE = "phase"
    GROUP = "group"
    ENHANCEMENT = "enhancement"


@dataclass(frozen=True, slots=True)
class SeismicFieldResult:
    """Acoustic fields sampled on one regular spherical grid.

    Parameters
    ----------
    grid : SphericalGrid
        Angular grid and wave-normal directions.
    level : SamplingLevel
        Highest quantity evaluated by the sampler.
    phase : PhaseFieldResult
        Phase velocities, polarizations and degeneracy diagnostics.
    group : GroupFieldResult or None
        Group velocities and power-flow quantities, when requested.
    enhancement : EnhancementFieldResult or None
        Curvature and enhancement quantities, when requested.
    tracking : PolarizationTrackingResult or None
        Axial-polarization continuity data, when requested.
    batch_size : int
        Maximum number of directions evaluated in one NumPy batch.
    """

    grid: SphericalGrid
    level: SamplingLevel
    phase: PhaseFieldResult
    group: GroupFieldResult | None
    enhancement: EnhancementFieldResult | None
    tracking: PolarizationTrackingResult | None
    batch_size: int

    @property
    def n_points(self) -> int:
        """Return the number of sampled grid positions."""
        return self.grid.size

    @property
    def has_group(self) -> bool:
        """Return whether group quantities are available."""
        return self.group is not None

    @property
    def has_enhancement(self) -> bool:
        """Return whether enhancement quantities are available."""
        return self.enhancement is not None

    @property
    def has_tracking(self) -> bool:
        """Return whether polarization tracking was requested."""
        return self.tracking is not None


def sample_seismic_field(
    solver: ChristoffelSolver,
    grid: SphericalGrid,
    *,
    level: SamplingLevel | str = SamplingLevel.ENHANCEMENT,
    batch_size: int = 512,
    track_polarization_axes: bool = False,
    progress_callback: ProgressCallback | None = None,
) -> SeismicFieldResult:
    """Sample phase, group and enhancement fields in bounded NumPy batches.

    Parameters
    ----------
    solver : ChristoffelSolver
        Pointwise acoustic solver defining the material and tolerances.
    grid : SphericalGrid
        Spherical grid to sample.
    level : SamplingLevel or str, optional
        Highest quantity to calculate. Group sampling includes phase results;
        enhancement sampling includes phase and group results.
    batch_size : int, optional
        Maximum number of directions evaluated in one batch.
    track_polarization_axes : bool, optional
        Whether to align polarization axes along a deterministic serpentine
        traversal of the grid.
    progress_callback : callable or None, optional
        Callback receiving ``(current, total)`` after every completed batch.

    Returns
    -------
    SeismicFieldResult
        Read-only sampled fields and optional polarization tracking.

    Raises
    ------
    TypeError
        If ``solver`` or ``grid`` has an unsupported type.
    ValueError
        If ``level`` or ``batch_size`` is invalid.
    """
    if not isinstance(solver, ChristoffelSolver):
        raise TypeError("sample_seismic_field requires a ChristoffelSolver.")
    if not isinstance(grid, SphericalGrid):
        raise TypeError("sample_seismic_field requires a SphericalGrid.")
    calculation = SamplingLevel(level)
    size = _validate_batch_size(batch_size)

    directions = grid.flat_directions
    arrays = _allocate_arrays(grid.size, calculation)

    for start in range(0, grid.size, size):
        stop = min(start + size, grid.size)
        batch = _solve_batch(
            solver,
            directions[start:stop],
            calculation,
        )
        _store_batch(arrays, batch, slice(start, stop), calculation)
        if progress_callback is not None:
            progress_callback(stop, grid.size)

    phase = _build_phase_field(directions, arrays)
    group = (
        _build_group_field(phase, arrays)
        if calculation is not SamplingLevel.PHASE
        else None
    )
    enhancement = (
        _build_enhancement_field(group, arrays)
        if calculation is SamplingLevel.ENHANCEMENT and group is not None
        else None
    )
    tracking = (
        _track_grid_polarizations(phase, grid) if track_polarization_axes else None
    )

    return SeismicFieldResult(
        grid=grid,
        level=calculation,
        phase=phase,
        group=group,
        enhancement=enhancement,
        tracking=tracking,
        batch_size=size,
    )


@dataclass(slots=True)
class _SamplingArrays:
    """Mutable preallocated arrays used while sampling one grid."""

    eigenvalues: NDArray[np.float64]
    phase_speeds: NDArray[np.float64]
    polarizations: NDArray[np.float64]
    mode_eigenvalue_gaps: NDArray[np.float64]
    mode_relative_eigenvalue_gaps: NDArray[np.float64]
    pair_eigenvalue_gaps: NDArray[np.float64]
    pair_relative_eigenvalue_gaps: NDArray[np.float64]
    phase_valid_mask: NDArray[np.bool_]
    clamped_mask: NDArray[np.bool_]
    degeneracy_mask: NDArray[np.bool_]
    pair_degeneracy_mask: NDArray[np.bool_]
    eigenvalue_thresholds: NDArray[np.float64]
    degeneracy_thresholds: NDArray[np.float64]
    eigenvalue_gradients: NDArray[np.float64] | None = None
    group_velocities: NDArray[np.float64] | None = None
    group_speeds: NDArray[np.float64] | None = None
    ray_directions: NDArray[np.float64] | None = None
    power_flow_angles: NDArray[np.float64] | None = None
    group_valid_mask: NDArray[np.bool_] | None = None
    group_resolved_mask: NDArray[np.bool_] | None = None
    eigenvalue_hessians: NDArray[np.float64] | None = None
    ray_direction_gradients: NDArray[np.float64] | None = None
    area_factors: NDArray[np.float64] | None = None
    caustic_thresholds: NDArray[np.float64] | None = None
    enhancement: NDArray[np.float64] | None = None
    log10_enhancement: NDArray[np.float64] | None = None
    enhancement_valid_mask: NDArray[np.bool_] | None = None
    enhancement_resolved_mask: NDArray[np.bool_] | None = None
    finite_mask: NDArray[np.bool_] | None = None
    caustic_candidate_mask: NDArray[np.bool_] | None = None


@dataclass(frozen=True, slots=True)
class _BatchResult:
    """Numerical arrays produced for one independent direction batch."""

    values: dict[str, NDArray[np.generic]]


def _allocate_arrays(n_points: int, level: SamplingLevel) -> _SamplingArrays:
    """Allocate mutable arrays for one complete sampled field."""
    arrays = _SamplingArrays(
        eigenvalues=np.empty((n_points, 3), dtype=float),
        phase_speeds=np.empty((n_points, 3), dtype=float),
        polarizations=np.empty((n_points, 3, 3), dtype=float),
        mode_eigenvalue_gaps=np.empty((n_points, 3), dtype=float),
        mode_relative_eigenvalue_gaps=np.empty((n_points, 3), dtype=float),
        pair_eigenvalue_gaps=np.empty((n_points, 2), dtype=float),
        pair_relative_eigenvalue_gaps=np.empty((n_points, 2), dtype=float),
        phase_valid_mask=np.empty((n_points, 3), dtype=bool),
        clamped_mask=np.empty((n_points, 3), dtype=bool),
        degeneracy_mask=np.empty((n_points, 3), dtype=bool),
        pair_degeneracy_mask=np.empty((n_points, 2), dtype=bool),
        eigenvalue_thresholds=np.empty(n_points, dtype=float),
        degeneracy_thresholds=np.empty(n_points, dtype=float),
    )
    if level in (SamplingLevel.GROUP, SamplingLevel.ENHANCEMENT):
        arrays.eigenvalue_gradients = np.empty((n_points, 3, 3), dtype=float)
        arrays.group_velocities = np.empty((n_points, 3, 3), dtype=float)
        arrays.group_speeds = np.empty((n_points, 3), dtype=float)
        arrays.ray_directions = np.empty((n_points, 3, 3), dtype=float)
        arrays.power_flow_angles = np.empty((n_points, 3), dtype=float)
        arrays.group_valid_mask = np.empty((n_points, 3), dtype=bool)
        arrays.group_resolved_mask = np.empty((n_points, 3), dtype=bool)
    if level is SamplingLevel.ENHANCEMENT:
        arrays.eigenvalue_hessians = np.empty((n_points, 3, 3, 3), dtype=float)
        arrays.ray_direction_gradients = np.empty((n_points, 3, 3, 3), dtype=float)
        arrays.area_factors = np.empty((n_points, 3), dtype=float)
        arrays.caustic_thresholds = np.empty((n_points, 3), dtype=float)
        arrays.enhancement = np.empty((n_points, 3), dtype=float)
        arrays.log10_enhancement = np.empty((n_points, 3), dtype=float)
        arrays.enhancement_valid_mask = np.empty((n_points, 3), dtype=bool)
        arrays.enhancement_resolved_mask = np.empty((n_points, 3), dtype=bool)
        arrays.finite_mask = np.empty((n_points, 3), dtype=bool)
        arrays.caustic_candidate_mask = np.empty((n_points, 3), dtype=bool)
    return arrays


def _solve_batch(
    solver: ChristoffelSolver,
    directions: NDArray[np.float64],
    level: SamplingLevel,
) -> _BatchResult:
    """Evaluate one direction batch using vectorized NumPy kernels."""
    q = np.asarray(directions, dtype=float)
    christoffel = np.einsum(
        "nj,ijkl,nl->nik",
        q,
        solver.reduced_stiffness,
        q,
        optimize=True,
    )
    christoffel = 0.5 * (christoffel + np.swapaxes(christoffel, 1, 2))
    eigenvalues, eigenvectors = np.linalg.eigh(christoffel)
    polarizations = np.swapaxes(eigenvectors, 1, 2)

    scales = np.maximum(
        np.max(np.abs(eigenvalues), axis=1),
        np.finfo(float).tiny,
    )
    eigenvalue_thresholds = solver.eigenvalue_atol + solver.eigenvalue_rtol * scales
    phase_valid = eigenvalues >= -eigenvalue_thresholds[:, np.newaxis]
    clamped = phase_valid & (eigenvalues < 0.0)
    speed_values = np.where(clamped, 0.0, eigenvalues)
    phase_speeds = np.full_like(eigenvalues, np.nan)
    phase_speeds[phase_valid] = np.sqrt(speed_values[phase_valid])

    pair_gaps = np.diff(eigenvalues, axis=1)
    pair_relative_gaps = pair_gaps / scales[:, np.newaxis]
    mode_gaps = np.column_stack(
        (pair_gaps[:, 0], np.minimum(pair_gaps[:, 0], pair_gaps[:, 1]), pair_gaps[:, 1])
    )
    mode_relative_gaps = mode_gaps / scales[:, np.newaxis]
    degeneracy_thresholds = solver.degeneracy_atol + solver.degeneracy_rtol * scales
    pair_degeneracy = pair_gaps <= degeneracy_thresholds[:, np.newaxis]
    degeneracy = mode_gaps <= degeneracy_thresholds[:, np.newaxis]

    values: dict[str, NDArray[np.generic]] = {
        "eigenvalues": eigenvalues,
        "phase_speeds": phase_speeds,
        "polarizations": polarizations,
        "mode_eigenvalue_gaps": mode_gaps,
        "mode_relative_eigenvalue_gaps": mode_relative_gaps,
        "pair_eigenvalue_gaps": pair_gaps,
        "pair_relative_eigenvalue_gaps": pair_relative_gaps,
        "phase_valid_mask": phase_valid,
        "clamped_mask": clamped,
        "degeneracy_mask": degeneracy,
        "pair_degeneracy_mask": pair_degeneracy,
        "eigenvalue_thresholds": eigenvalue_thresholds,
        "degeneracy_thresholds": degeneracy_thresholds,
    }
    if level is SamplingLevel.PHASE:
        return _BatchResult(values)

    gradient_tensor = solver.reduced_stiffness + np.transpose(
        solver.reduced_stiffness,
        (0, 2, 1, 3),
    )
    christoffel_gradient = np.einsum(
        "nk,iakl->nail",
        q,
        gradient_tensor,
        optimize=True,
    )
    christoffel_gradient = 0.5 * (
        christoffel_gradient + np.swapaxes(christoffel_gradient, 2, 3)
    )
    eigenvalue_gradients = np.einsum(
        "nmi,naij,nmj->nma",
        polarizations,
        christoffel_gradient,
        polarizations,
        optimize=True,
    )
    group_valid = phase_valid & np.isfinite(phase_speeds) & (phase_speeds > 0.0)
    group_velocities = np.full((q.shape[0], 3, 3), np.nan, dtype=float)
    group_velocities[group_valid] = eigenvalue_gradients[group_valid] / (
        2.0 * phase_speeds[group_valid, np.newaxis]
    )
    group_speeds = np.full((q.shape[0], 3), np.nan, dtype=float)
    group_speeds[group_valid] = np.linalg.norm(group_velocities[group_valid], axis=1)
    group_valid &= np.isfinite(group_speeds) & (group_speeds > 0.0)
    ray_directions = np.full((q.shape[0], 3, 3), np.nan, dtype=float)
    ray_directions[group_valid] = (
        group_velocities[group_valid] / group_speeds[group_valid, np.newaxis]
    )
    power_flow_angles = np.full((q.shape[0], 3), np.nan, dtype=float)
    cosine = np.einsum("nmi,ni->nm", ray_directions, q, optimize=True)
    power_flow_angles[group_valid] = np.arccos(np.clip(cosine[group_valid], -1.0, 1.0))
    group_resolved = group_valid & ~degeneracy

    values.update(
        {
            "eigenvalue_gradients": eigenvalue_gradients,
            "group_velocities": group_velocities,
            "group_speeds": group_speeds,
            "ray_directions": ray_directions,
            "power_flow_angles": power_flow_angles,
            "group_valid_mask": group_valid,
            "group_resolved_mask": group_resolved,
        }
    )
    if level is SamplingLevel.GROUP:
        return _BatchResult(values)

    eigenvalue_hessians = np.full((q.shape[0], 3, 3, 3), np.nan, dtype=float)
    ray_gradients = np.full((q.shape[0], 3, 3, 3), np.nan, dtype=float)
    area_factors = np.full((q.shape[0], 3), np.nan, dtype=float)
    caustic_thresholds = np.full((q.shape[0], 3), np.nan, dtype=float)
    enhancement = np.full((q.shape[0], 3), np.nan, dtype=float)
    log10_enhancement = np.full((q.shape[0], 3), np.nan, dtype=float)
    enhancement_valid = np.zeros((q.shape[0], 3), dtype=bool)
    finite_mask = np.zeros((q.shape[0], 3), dtype=bool)
    caustic_mask = np.zeros((q.shape[0], 3), dtype=bool)

    direct = np.einsum(
        "nmi,abij,nmj->nmab",
        polarizations,
        solver.christoffel_hessian,
        polarizations,
        optimize=True,
    )
    derivative_vectors = np.einsum(
        "naij,nmj->nmai",
        christoffel_gradient,
        polarizations,
        optimize=True,
    )
    shifted = (
        eigenvalues[..., np.newaxis, np.newaxis] * np.eye(3)[np.newaxis, np.newaxis]
        - christoffel[:, np.newaxis]
    )
    inverse = np.linalg.pinv(shifted, rcond=solver.pseudoinverse_rcond)
    correction = 2.0 * np.einsum(
        "nmai,nmij,nmbj->nmab",
        derivative_vectors,
        inverse,
        derivative_vectors,
        optimize=True,
    )
    hessians = np.asarray(direct + correction, dtype=float)
    hessians = 0.5 * (hessians + np.swapaxes(hessians, 2, 3))

    denominators = 2.0 * phase_speeds * group_speeds
    curvature_valid = (
        group_resolved
        & np.all(np.isfinite(hessians), axis=(2, 3))
        & np.isfinite(denominators)
        & (denominators > 0.0)
        & np.all(np.isfinite(ray_directions), axis=2)
    )
    projections = (
        np.eye(3)[np.newaxis, np.newaxis]
        - ray_directions[..., :, np.newaxis] * ray_directions[..., np.newaxis, :]
    )
    calculated_gradients = (
        np.einsum(
            "nmij,nmjk->nmik",
            projections,
            hessians,
            optimize=True,
        )
        / denominators[..., np.newaxis, np.newaxis]
    )

    cofactors = _batch_cofactor(calculated_gradients)
    mapped = np.einsum("nmij,nj->nmi", cofactors, q, optimize=True)
    calculated_area = np.linalg.norm(mapped, axis=2)
    cofactor_scales = np.maximum(
        np.linalg.norm(cofactors, axis=(2, 3)),
        np.finfo(float).tiny,
    )
    calculated_thresholds = solver.caustic_atol + solver.caustic_rtol * cofactor_scales
    curvature_valid &= np.isfinite(calculated_area) & (calculated_area >= 0.0)

    eigenvalue_hessians[curvature_valid] = hessians[curvature_valid]
    ray_gradients[curvature_valid] = calculated_gradients[curvature_valid]
    area_factors[curvature_valid] = calculated_area[curvature_valid]
    caustic_thresholds[curvature_valid] = calculated_thresholds[curvature_valid]
    enhancement_valid[curvature_valid] = True
    caustic_mask[curvature_valid] = (
        calculated_area[curvature_valid] <= calculated_thresholds[curvature_valid]
    )
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        enhancement[curvature_valid] = 1.0 / calculated_area[curvature_valid]
        log10_enhancement[curvature_valid] = np.log10(enhancement[curvature_valid])
    finite_mask[curvature_valid] = np.isfinite(
        enhancement[curvature_valid]
    ) & np.isfinite(log10_enhancement[curvature_valid])

    values.update(
        {
            "eigenvalue_hessians": eigenvalue_hessians,
            "ray_direction_gradients": ray_gradients,
            "area_factors": area_factors,
            "caustic_thresholds": caustic_thresholds,
            "enhancement": enhancement,
            "log10_enhancement": log10_enhancement,
            "enhancement_valid_mask": enhancement_valid,
            "enhancement_resolved_mask": group_resolved,
            "finite_mask": finite_mask,
            "caustic_candidate_mask": caustic_mask,
        }
    )
    return _BatchResult(values)


def _store_batch(
    arrays: _SamplingArrays,
    batch: _BatchResult,
    target: slice,
    level: SamplingLevel,
) -> None:
    """Copy one numerical batch into the preallocated field arrays."""
    required = (
        "eigenvalues",
        "phase_speeds",
        "polarizations",
        "mode_eigenvalue_gaps",
        "mode_relative_eigenvalue_gaps",
        "pair_eigenvalue_gaps",
        "pair_relative_eigenvalue_gaps",
        "phase_valid_mask",
        "clamped_mask",
        "degeneracy_mask",
        "pair_degeneracy_mask",
        "eigenvalue_thresholds",
        "degeneracy_thresholds",
    )
    for name in required:
        getattr(arrays, name)[target] = batch.values[name]

    if level in (SamplingLevel.GROUP, SamplingLevel.ENHANCEMENT):
        for name in (
            "eigenvalue_gradients",
            "group_velocities",
            "group_speeds",
            "ray_directions",
            "power_flow_angles",
            "group_valid_mask",
            "group_resolved_mask",
        ):
            target_array = getattr(arrays, name)
            assert target_array is not None
            target_array[target] = batch.values[name]

    if level is SamplingLevel.ENHANCEMENT:
        for name in (
            "eigenvalue_hessians",
            "ray_direction_gradients",
            "area_factors",
            "caustic_thresholds",
            "enhancement",
            "log10_enhancement",
            "enhancement_valid_mask",
            "enhancement_resolved_mask",
            "finite_mask",
            "caustic_candidate_mask",
        ):
            target_array = getattr(arrays, name)
            assert target_array is not None
            target_array[target] = batch.values[name]


def _build_phase_field(
    directions: NDArray[np.float64],
    arrays: _SamplingArrays,
) -> PhaseFieldResult:
    """Create an immutable phase field from completed sampling arrays."""
    return PhaseFieldResult(
        directions=_readonly_float(directions),
        eigenvalues=_readonly_float(arrays.eigenvalues),
        phase_speeds=_readonly_float(arrays.phase_speeds),
        polarizations=_readonly_float(arrays.polarizations),
        mode_eigenvalue_gaps=_readonly_float(arrays.mode_eigenvalue_gaps),
        mode_relative_eigenvalue_gaps=_readonly_float(
            arrays.mode_relative_eigenvalue_gaps
        ),
        pair_eigenvalue_gaps=_readonly_float(arrays.pair_eigenvalue_gaps),
        pair_relative_eigenvalue_gaps=_readonly_float(
            arrays.pair_relative_eigenvalue_gaps
        ),
        valid_mask=_readonly_bool(arrays.phase_valid_mask),
        clamped_mask=_readonly_bool(arrays.clamped_mask),
        degeneracy_mask=_readonly_bool(arrays.degeneracy_mask),
        pair_degeneracy_mask=_readonly_bool(arrays.pair_degeneracy_mask),
        eigenvalue_thresholds=_readonly_float(arrays.eigenvalue_thresholds),
        degeneracy_thresholds=_readonly_float(arrays.degeneracy_thresholds),
    )


def _build_group_field(
    phase: PhaseFieldResult,
    arrays: _SamplingArrays,
) -> GroupFieldResult:
    """Create an immutable group field from completed sampling arrays."""
    assert arrays.eigenvalue_gradients is not None
    assert arrays.group_velocities is not None
    assert arrays.group_speeds is not None
    assert arrays.ray_directions is not None
    assert arrays.power_flow_angles is not None
    assert arrays.group_valid_mask is not None
    assert arrays.group_resolved_mask is not None
    return GroupFieldResult(
        phase=phase,
        eigenvalue_gradients=_readonly_float(arrays.eigenvalue_gradients),
        group_velocities=_readonly_float(arrays.group_velocities),
        group_speeds=_readonly_float(arrays.group_speeds),
        ray_directions=_readonly_float(arrays.ray_directions),
        power_flow_angles=_readonly_float(arrays.power_flow_angles),
        valid_mask=_readonly_bool(arrays.group_valid_mask),
        resolved_mask=_readonly_bool(arrays.group_resolved_mask),
    )


def _build_enhancement_field(
    group: GroupFieldResult,
    arrays: _SamplingArrays,
) -> EnhancementFieldResult:
    """Create an immutable enhancement field from completed sampling arrays."""
    assert arrays.eigenvalue_hessians is not None
    assert arrays.ray_direction_gradients is not None
    assert arrays.area_factors is not None
    assert arrays.caustic_thresholds is not None
    assert arrays.enhancement is not None
    assert arrays.log10_enhancement is not None
    assert arrays.enhancement_valid_mask is not None
    assert arrays.enhancement_resolved_mask is not None
    assert arrays.finite_mask is not None
    assert arrays.caustic_candidate_mask is not None
    return EnhancementFieldResult(
        group=group,
        eigenvalue_hessians=_readonly_float(arrays.eigenvalue_hessians),
        ray_direction_gradients=_readonly_float(arrays.ray_direction_gradients),
        area_factors=_readonly_float(arrays.area_factors),
        caustic_thresholds=_readonly_float(arrays.caustic_thresholds),
        enhancement=_readonly_float(arrays.enhancement),
        log10_enhancement=_readonly_float(arrays.log10_enhancement),
        valid_mask=_readonly_bool(arrays.enhancement_valid_mask),
        resolved_mask=_readonly_bool(arrays.enhancement_resolved_mask),
        finite_mask=_readonly_bool(arrays.finite_mask),
        caustic_candidate_mask=_readonly_bool(arrays.caustic_candidate_mask),
    )


def _track_grid_polarizations(
    phase: PhaseFieldResult,
    grid: SphericalGrid,
) -> PolarizationTrackingResult:
    """Track polarization axes along a serpentine nearest-neighbour path."""
    order = np.arange(grid.size, dtype=np.int64).reshape(grid.shape)
    order[1::2] = order[1::2, ::-1]
    traversal = order.ravel()
    inverse = np.empty_like(traversal)
    inverse[traversal] = np.arange(grid.size, dtype=np.int64)

    traversed = PhaseFieldResult(
        directions=_readonly_float(phase.directions[traversal]),
        eigenvalues=_readonly_float(phase.eigenvalues[traversal]),
        phase_speeds=_readonly_float(phase.phase_speeds[traversal]),
        polarizations=_readonly_float(phase.polarizations[traversal]),
        mode_eigenvalue_gaps=_readonly_float(phase.mode_eigenvalue_gaps[traversal]),
        mode_relative_eigenvalue_gaps=_readonly_float(
            phase.mode_relative_eigenvalue_gaps[traversal]
        ),
        pair_eigenvalue_gaps=_readonly_float(phase.pair_eigenvalue_gaps[traversal]),
        pair_relative_eigenvalue_gaps=_readonly_float(
            phase.pair_relative_eigenvalue_gaps[traversal]
        ),
        valid_mask=_readonly_bool(phase.valid_mask[traversal]),
        clamped_mask=_readonly_bool(phase.clamped_mask[traversal]),
        degeneracy_mask=_readonly_bool(phase.degeneracy_mask[traversal]),
        pair_degeneracy_mask=_readonly_bool(phase.pair_degeneracy_mask[traversal]),
        eigenvalue_thresholds=_readonly_float(phase.eigenvalue_thresholds[traversal]),
        degeneracy_thresholds=_readonly_float(phase.degeneracy_thresholds[traversal]),
    )
    tracked = track_polarizations(traversed)
    return PolarizationTrackingResult(
        directions=phase.directions,
        polarizations=_readonly_float(tracked.polarizations[inverse]),
        branch_mode_indices=_readonly_int(tracked.branch_mode_indices[inverse]),
        sign_flip_mask=_readonly_bool(tracked.sign_flip_mask[inverse]),
        continuity_scores=_readonly_float(tracked.continuity_scores[inverse]),
        resolved_mask=_readonly_bool(tracked.resolved_mask[inverse]),
        segment_start_mask=_readonly_bool(tracked.segment_start_mask[inverse]),
        shear_swap_mask=_readonly_bool(tracked.shear_swap_mask[inverse]),
        shear_permutation_ambiguous_mask=_readonly_bool(
            tracked.shear_permutation_ambiguous_mask[inverse]
        ),
        shear_subspace_rotation_mask=_readonly_bool(
            tracked.shear_subspace_rotation_mask[inverse]
        ),
    )


def _batch_cofactor(matrices: NDArray[np.float64]) -> NDArray[np.float64]:
    """Calculate cofactors for an array with final shape ``(3, 3)``."""
    value = matrices
    result = np.empty_like(value)
    result[..., 0, 0] = (
        value[..., 1, 1] * value[..., 2, 2] - value[..., 1, 2] * value[..., 2, 1]
    )
    result[..., 0, 1] = (
        value[..., 1, 2] * value[..., 2, 0] - value[..., 1, 0] * value[..., 2, 2]
    )
    result[..., 0, 2] = (
        value[..., 1, 0] * value[..., 2, 1] - value[..., 1, 1] * value[..., 2, 0]
    )
    result[..., 1, 0] = (
        value[..., 0, 2] * value[..., 2, 1] - value[..., 0, 1] * value[..., 2, 2]
    )
    result[..., 1, 1] = (
        value[..., 0, 0] * value[..., 2, 2] - value[..., 0, 2] * value[..., 2, 0]
    )
    result[..., 1, 2] = (
        value[..., 0, 1] * value[..., 2, 0] - value[..., 0, 0] * value[..., 2, 1]
    )
    result[..., 2, 0] = (
        value[..., 0, 1] * value[..., 1, 2] - value[..., 0, 2] * value[..., 1, 1]
    )
    result[..., 2, 1] = (
        value[..., 0, 2] * value[..., 1, 0] - value[..., 0, 0] * value[..., 1, 2]
    )
    result[..., 2, 2] = (
        value[..., 0, 0] * value[..., 1, 1] - value[..., 0, 1] * value[..., 1, 0]
    )
    return result


def _validate_batch_size(value: int) -> int:
    """Validate a positive integer batch size."""
    if isinstance(value, bool) or int(value) != value or int(value) < 1:
        raise ValueError("batch_size must be a positive integer.")
    return int(value)


def _readonly_float(array: NDArray[np.float64]) -> NDArray[np.float64]:
    """Return an independent read-only floating-point array."""
    result = np.array(array, dtype=float, copy=True)
    result.setflags(write=False)
    return result


def _readonly_bool(array: NDArray[np.bool_]) -> NDArray[np.bool_]:
    """Return an independent read-only boolean array."""
    result = np.array(array, dtype=bool, copy=True)
    result.setflags(write=False)
    return result


def _readonly_int(array: NDArray[np.int64]) -> NDArray[np.int64]:
    """Return an independent read-only integer array."""
    result = np.array(array, dtype=np.int64, copy=True)
    result.setflags(write=False)
    return result
