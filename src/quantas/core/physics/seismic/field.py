# -*- coding: utf-8 -*-

"""Sampled phase and group fields with propagation diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from numpy.typing import NDArray

from .enhancement import DirectionalEnhancementResult
from .group import DirectionalGroupResult
from .modes import (
    MODE_INDEX,
    MODE_PAIR_INDEX,
    DirectionalPhaseResult,
    ModePair,
    WaveMode,
)


@dataclass(frozen=True, slots=True)
class PhaseModeField:
    """Sampled phase quantities for one acoustic mode.

    Parameters
    ----------
    mode : WaveMode
        Acoustic-wave mode.
    eigenvalues : ndarray
        Christoffel eigenvalues in km^2 s^-2 with shape ``(n_points,)``.
    phase_speeds : ndarray
        Phase speeds in km s^-1 with shape ``(n_points,)``.
    polarizations : ndarray
        Polarization axes with shape ``(n_points, 3)``.
    eigenvalue_gaps : ndarray
        Smallest adjacent eigenvalue gap for the mode in km^2 s^-2.
    relative_eigenvalue_gaps : ndarray
        Eigenvalue gaps normalized by the largest directional eigenvalue.
    valid_mask : ndarray
        Boolean validity mask.
    clamped_mask : ndarray
        Boolean mask for small negative eigenvalues set to zero.
    degeneracy_mask : ndarray
        Boolean mask for modes belonging to a degenerate eigenspace.
    """

    mode: WaveMode
    eigenvalues: NDArray[np.float64]
    phase_speeds: NDArray[np.float64]
    polarizations: NDArray[np.float64]
    eigenvalue_gaps: NDArray[np.float64]
    relative_eigenvalue_gaps: NDArray[np.float64]
    valid_mask: NDArray[np.bool_]
    clamped_mask: NDArray[np.bool_]
    degeneracy_mask: NDArray[np.bool_]


@dataclass(frozen=True, slots=True)
class ModePairField:
    """Sampled separation and degeneracy state for one adjacent mode pair.

    Parameters
    ----------
    pair : ModePair
        Adjacent acoustic-mode pair.
    eigenvalue_gaps : ndarray
        Absolute eigenvalue separations in km^2 s^-2.
    relative_eigenvalue_gaps : ndarray
        Separations normalized by the largest directional eigenvalue.
    degeneracy_mask : ndarray
        Boolean mask identifying sampled degeneracy candidates.
    """

    pair: ModePair
    eigenvalue_gaps: NDArray[np.float64]
    relative_eigenvalue_gaps: NDArray[np.float64]
    degeneracy_mask: NDArray[np.bool_]


@dataclass(frozen=True, slots=True)
class GroupModeField:
    """Sampled group quantities for one acoustic mode.

    Parameters
    ----------
    mode : WaveMode
        Acoustic-wave mode.
    phase : PhaseModeField
        Sampled phase quantities for the same mode.
    eigenvalue_gradients : ndarray
        Eigenvalue gradients with shape ``(n_points, 3)``.
    group_velocities : ndarray
        Group-velocity vectors in km s^-1 with shape ``(n_points, 3)``.
    group_speeds : ndarray
        Group-speed magnitudes in km s^-1.
    ray_directions : ndarray
        Unit energy-flow directions with shape ``(n_points, 3)``.
    power_flow_angles : ndarray
        Power-flow angles in radians.
    valid_mask : ndarray
        Boolean mask identifying numerically defined group quantities.
    resolved_mask : ndarray
        Boolean mask identifying non-degenerate group solutions.
    """

    mode: WaveMode
    phase: PhaseModeField
    eigenvalue_gradients: NDArray[np.float64]
    group_velocities: NDArray[np.float64]
    group_speeds: NDArray[np.float64]
    ray_directions: NDArray[np.float64]
    power_flow_angles: NDArray[np.float64]
    valid_mask: NDArray[np.bool_]
    resolved_mask: NDArray[np.bool_]


@dataclass(frozen=True, slots=True)
class EnhancementModeField:
    """Sampled curvature and enhancement quantities for one acoustic mode.

    Parameters
    ----------
    mode : WaveMode
        Acoustic-wave mode.
    group : GroupModeField
        Sampled phase and group quantities for the same mode.
    eigenvalue_hessians : ndarray
        Eigenvalue Hessians with shape ``(n_points, 3, 3)``.
    ray_direction_gradients : ndarray
        Ray-direction gradients with shape ``(n_points, 3, 3)``.
    area_factors : ndarray
        Geometrical area factors.
    caustic_thresholds : ndarray
        Numerical thresholds for possible caustics.
    enhancement : ndarray
        Raw enhancement factors :math:`A`.
    log10_enhancement : ndarray
        Base-10 logarithms of the enhancement factors.
    valid_mask : ndarray
        Boolean mask identifying defined curvature quantities.
    resolved_mask : ndarray
        Boolean mask identifying uniquely resolved acoustic modes.
    finite_mask : ndarray
        Boolean mask identifying finite enhancement values.
    caustic_candidate_mask : ndarray
        Boolean mask identifying area factors close to zero.
    """

    mode: WaveMode
    group: GroupModeField
    eigenvalue_hessians: NDArray[np.float64]
    ray_direction_gradients: NDArray[np.float64]
    area_factors: NDArray[np.float64]
    caustic_thresholds: NDArray[np.float64]
    enhancement: NDArray[np.float64]
    log10_enhancement: NDArray[np.float64]
    valid_mask: NDArray[np.bool_]
    resolved_mask: NDArray[np.bool_]
    finite_mask: NDArray[np.bool_]
    caustic_candidate_mask: NDArray[np.bool_]


@dataclass(frozen=True, slots=True)
class PhaseFieldResult:
    """Phase solutions sampled on an ordered collection of directions.

    Parameters
    ----------
    directions : ndarray
        Unit wave-normal directions with shape ``(n_points, 3)``.
    eigenvalues : ndarray
        Christoffel eigenvalues with shape ``(n_points, 3)`` in the order
        ``V_S2``, ``V_S1`` and ``V_P``.
    phase_speeds : ndarray
        Phase speeds in km s^-1 with shape ``(n_points, 3)``.
    polarizations : ndarray
        Row-oriented polarization axes with shape ``(n_points, 3, 3)``.
    mode_eigenvalue_gaps : ndarray
        Smallest adjacent eigenvalue gap for each mode.
    mode_relative_eigenvalue_gaps : ndarray
        Mode gaps normalized by the largest directional eigenvalue.
    pair_eigenvalue_gaps : ndarray
        Adjacent pair gaps with shape ``(n_points, 2)``.
    pair_relative_eigenvalue_gaps : ndarray
        Relative adjacent pair gaps with shape ``(n_points, 2)``.
    valid_mask, clamped_mask, degeneracy_mask : ndarray
        Per-mode diagnostic masks with shape ``(n_points, 3)``.
    pair_degeneracy_mask : ndarray
        Per-pair degeneracy mask with shape ``(n_points, 2)``.
    eigenvalue_thresholds : ndarray
        Negative-eigenvalue thresholds in km^2 s^-2.
    degeneracy_thresholds : ndarray
        Eigenvalue-gap thresholds in km^2 s^-2.
    """

    directions: NDArray[np.float64]
    eigenvalues: NDArray[np.float64]
    phase_speeds: NDArray[np.float64]
    polarizations: NDArray[np.float64]
    mode_eigenvalue_gaps: NDArray[np.float64]
    mode_relative_eigenvalue_gaps: NDArray[np.float64]
    pair_eigenvalue_gaps: NDArray[np.float64]
    pair_relative_eigenvalue_gaps: NDArray[np.float64]
    valid_mask: NDArray[np.bool_]
    clamped_mask: NDArray[np.bool_]
    degeneracy_mask: NDArray[np.bool_]
    pair_degeneracy_mask: NDArray[np.bool_]
    eigenvalue_thresholds: NDArray[np.float64]
    degeneracy_thresholds: NDArray[np.float64]

    @property
    def n_points(self) -> int:
        """Return the number of sampled directions."""
        return int(self.directions.shape[0])

    def for_mode(self, mode: WaveMode | str) -> PhaseModeField:
        """Return sampled quantities for one acoustic mode.

        Parameters
        ----------
        mode : WaveMode or str
            Acoustic mode enum or its string value.

        Returns
        -------
        PhaseModeField
            Read-only views of the selected mode field.

        Raises
        ------
        ValueError
            If ``mode`` is not supported.
        """
        resolved = WaveMode(mode)
        index = MODE_INDEX[resolved]
        return PhaseModeField(
            mode=resolved,
            eigenvalues=_readonly_float_view(self.eigenvalues[:, index]),
            phase_speeds=_readonly_float_view(self.phase_speeds[:, index]),
            polarizations=_readonly_float_view(self.polarizations[:, index, :]),
            eigenvalue_gaps=_readonly_float_view(self.mode_eigenvalue_gaps[:, index]),
            relative_eigenvalue_gaps=_readonly_float_view(
                self.mode_relative_eigenvalue_gaps[:, index]
            ),
            valid_mask=_readonly_bool_view(self.valid_mask[:, index]),
            clamped_mask=_readonly_bool_view(self.clamped_mask[:, index]),
            degeneracy_mask=_readonly_bool_view(self.degeneracy_mask[:, index]),
        )

    def for_pair(self, pair: ModePair | str) -> ModePairField:
        """Return sampled diagnostics for one adjacent mode pair.

        Parameters
        ----------
        pair : ModePair or str
            Adjacent mode pair enum or its string value.

        Returns
        -------
        ModePairField
            Gap and degeneracy arrays for the selected pair.

        Raises
        ------
        ValueError
            If ``pair`` is not supported.
        """
        resolved = ModePair(pair)
        index = MODE_PAIR_INDEX[resolved]
        return ModePairField(
            pair=resolved,
            eigenvalue_gaps=_readonly_float_view(self.pair_eigenvalue_gaps[:, index]),
            relative_eigenvalue_gaps=_readonly_float_view(
                self.pair_relative_eigenvalue_gaps[:, index]
            ),
            degeneracy_mask=_readonly_bool_view(self.pair_degeneracy_mask[:, index]),
        )

    @property
    def acoustic_axis_candidate_mask(self) -> NDArray[np.bool_]:
        """Return sampled directions containing an adjacent-mode degeneracy."""
        return _readonly_bool_copy(np.any(self.pair_degeneracy_mask, axis=1))

    @property
    def shear_axis_candidate_mask(self) -> NDArray[np.bool_]:
        """Return candidate directions where ``V_S2`` and ``V_S1`` coincide."""
        index = MODE_PAIR_INDEX[ModePair.V_S2_V_S1]
        return _readonly_bool_view(self.pair_degeneracy_mask[:, index])

    @property
    def upper_axis_candidate_mask(self) -> NDArray[np.bool_]:
        """Return candidate directions where ``V_S1`` and ``V_P`` coincide."""
        index = MODE_PAIR_INDEX[ModePair.V_S1_V_P]
        return _readonly_bool_view(self.pair_degeneracy_mask[:, index])

    @property
    def triple_degeneracy_mask(self) -> NDArray[np.bool_]:
        """Return directions where both adjacent eigenvalue gaps are degenerate."""
        return _readonly_bool_copy(np.all(self.pair_degeneracy_mask, axis=1))

    def candidate_indices(
        self,
        pair: ModePair | str | None = None,
    ) -> NDArray[np.int64]:
        """Return sampled indices that are acoustic-axis candidates.

        Parameters
        ----------
        pair : ModePair, str or None, optional
            Restrict candidates to one adjacent pair. When omitted, all
            adjacent-mode degeneracies are included.

        Returns
        -------
        ndarray
            Read-only integer indices.

        Raises
        ------
        ValueError
            If ``pair`` is not supported.
        """
        if pair is None:
            mask = np.any(self.pair_degeneracy_mask, axis=1)
        else:
            resolved = ModePair(pair)
            mask = self.pair_degeneracy_mask[:, MODE_PAIR_INDEX[resolved]]
        return _readonly_int_copy(np.flatnonzero(mask).astype(np.int64, copy=False))

    def candidate_directions(
        self,
        pair: ModePair | str | None = None,
    ) -> NDArray[np.float64]:
        """Return sampled directions that are acoustic-axis candidates.

        Parameters
        ----------
        pair : ModePair, str or None, optional
            Restrict candidates to one adjacent pair.

        Returns
        -------
        ndarray
            Read-only directions with shape ``(n_candidates, 3)``.
        """
        return _readonly_float_copy(self.directions[self.candidate_indices(pair)])


@dataclass(frozen=True, slots=True)
class GroupFieldResult:
    """Group solutions sampled on an ordered collection of directions.

    Parameters
    ----------
    phase : PhaseFieldResult
        Sampled phase field associated with the group results.
    eigenvalue_gradients : ndarray
        Eigenvalue gradients with shape ``(n_points, 3, 3)``.
    group_velocities : ndarray
        Group-velocity vectors with shape ``(n_points, 3, 3)`` in km s^-1.
    group_speeds : ndarray
        Group-speed magnitudes with shape ``(n_points, 3)`` in km s^-1.
    ray_directions : ndarray
        Unit energy-flow directions with shape ``(n_points, 3, 3)``.
    power_flow_angles : ndarray
        Power-flow angles in radians with shape ``(n_points, 3)``.
    valid_mask : ndarray
        Per-mode mask for numerically defined group quantities.
    resolved_mask : ndarray
        Per-mode mask for non-degenerate group solutions.
    """

    phase: PhaseFieldResult
    eigenvalue_gradients: NDArray[np.float64]
    group_velocities: NDArray[np.float64]
    group_speeds: NDArray[np.float64]
    ray_directions: NDArray[np.float64]
    power_flow_angles: NDArray[np.float64]
    valid_mask: NDArray[np.bool_]
    resolved_mask: NDArray[np.bool_]

    @property
    def directions(self) -> NDArray[np.float64]:
        """Return the sampled unit wave-normal directions."""
        return self.phase.directions

    @property
    def n_points(self) -> int:
        """Return the number of sampled directions."""
        return self.phase.n_points

    def for_mode(self, mode: WaveMode | str) -> GroupModeField:
        """Return sampled phase and group quantities for one mode.

        Parameters
        ----------
        mode : WaveMode or str
            Acoustic mode enum or its string value.

        Returns
        -------
        GroupModeField
            Read-only views of the selected mode field.

        Raises
        ------
        ValueError
            If ``mode`` is not supported.
        """
        resolved = WaveMode(mode)
        index = MODE_INDEX[resolved]
        return GroupModeField(
            mode=resolved,
            phase=self.phase.for_mode(resolved),
            eigenvalue_gradients=_readonly_float_view(
                self.eigenvalue_gradients[:, index, :]
            ),
            group_velocities=_readonly_float_view(self.group_velocities[:, index, :]),
            group_speeds=_readonly_float_view(self.group_speeds[:, index]),
            ray_directions=_readonly_float_view(self.ray_directions[:, index, :]),
            power_flow_angles=_readonly_float_view(self.power_flow_angles[:, index]),
            valid_mask=_readonly_bool_view(self.valid_mask[:, index]),
            resolved_mask=_readonly_bool_view(self.resolved_mask[:, index]),
        )


@dataclass(frozen=True, slots=True)
class EnhancementFieldResult:
    """Enhancement solutions sampled on an ordered collection of directions.

    Parameters
    ----------
    group : GroupFieldResult
        Sampled phase and group field associated with the enhancement data.
    eigenvalue_hessians : ndarray
        Eigenvalue Hessians with shape ``(n_points, 3, 3, 3)``.
    ray_direction_gradients : ndarray
        Ray-direction gradients with shape ``(n_points, 3, 3, 3)``.
    area_factors : ndarray
        Geometrical area factors with shape ``(n_points, 3)``.
    caustic_thresholds : ndarray
        Mode-specific caustic thresholds with shape ``(n_points, 3)``.
    enhancement : ndarray
        Raw enhancement factors :math:`A` with shape ``(n_points, 3)``.
    log10_enhancement : ndarray
        Base-10 logarithms of the enhancement factors.
    valid_mask : ndarray
        Per-mode mask for defined analytical curvature quantities.
    resolved_mask : ndarray
        Per-mode mask for uniquely resolved acoustic modes.
    finite_mask : ndarray
        Per-mode mask for finite enhancement values.
    caustic_candidate_mask : ndarray
        Per-mode mask for possible caustics.
    """

    group: GroupFieldResult
    eigenvalue_hessians: NDArray[np.float64]
    ray_direction_gradients: NDArray[np.float64]
    area_factors: NDArray[np.float64]
    caustic_thresholds: NDArray[np.float64]
    enhancement: NDArray[np.float64]
    log10_enhancement: NDArray[np.float64]
    valid_mask: NDArray[np.bool_]
    resolved_mask: NDArray[np.bool_]
    finite_mask: NDArray[np.bool_]
    caustic_candidate_mask: NDArray[np.bool_]

    @property
    def directions(self) -> NDArray[np.float64]:
        """Return the sampled unit wave-normal directions."""
        return self.group.directions

    @property
    def n_points(self) -> int:
        """Return the number of sampled directions."""
        return self.group.n_points

    @property
    def has_caustic_candidate(self) -> bool:
        """Return whether any sampled mode approaches a caustic."""
        return bool(np.any(self.caustic_candidate_mask))

    def for_mode(self, mode: WaveMode | str) -> EnhancementModeField:
        """Return sampled enhancement quantities for one acoustic mode.

        Parameters
        ----------
        mode : WaveMode or str
            Acoustic mode enum or its string value.

        Returns
        -------
        EnhancementModeField
            Read-only views of the selected mode field.

        Raises
        ------
        ValueError
            If ``mode`` is not supported.
        """
        resolved = WaveMode(mode)
        index = MODE_INDEX[resolved]
        return EnhancementModeField(
            mode=resolved,
            group=self.group.for_mode(resolved),
            eigenvalue_hessians=_readonly_float_view(
                self.eigenvalue_hessians[:, index, :, :]
            ),
            ray_direction_gradients=_readonly_float_view(
                self.ray_direction_gradients[:, index, :, :]
            ),
            area_factors=_readonly_float_view(self.area_factors[:, index]),
            caustic_thresholds=_readonly_float_view(self.caustic_thresholds[:, index]),
            enhancement=_readonly_float_view(self.enhancement[:, index]),
            log10_enhancement=_readonly_float_view(self.log10_enhancement[:, index]),
            valid_mask=_readonly_bool_view(self.valid_mask[:, index]),
            resolved_mask=_readonly_bool_view(self.resolved_mask[:, index]),
            finite_mask=_readonly_bool_view(self.finite_mask[:, index]),
            caustic_candidate_mask=_readonly_bool_view(
                self.caustic_candidate_mask[:, index]
            ),
        )


def build_enhancement_field(
    results: Iterable[DirectionalEnhancementResult],
) -> EnhancementFieldResult:
    """Stack directional enhancement solutions into a sampled field.

    Parameters
    ----------
    results : iterable of DirectionalEnhancementResult
        Directional solutions in deterministic traversal order.

    Returns
    -------
    EnhancementFieldResult
        Read-only sampled phase, group and enhancement field.

    Raises
    ------
    TypeError
        If an item is not a :class:`DirectionalEnhancementResult`.
    ValueError
        If the iterable is empty.
    """
    entries = tuple(results)
    if not entries:
        raise ValueError("At least one directional enhancement result is required.")
    if not all(isinstance(item, DirectionalEnhancementResult) for item in entries):
        raise TypeError(
            "All field entries must be DirectionalEnhancementResult objects."
        )

    group = build_group_field(item.group for item in entries)
    return EnhancementFieldResult(
        group=group,
        eigenvalue_hessians=_stack_enhancement_float(entries, "eigenvalue_hessians"),
        ray_direction_gradients=_stack_enhancement_float(
            entries, "ray_direction_gradients"
        ),
        area_factors=_stack_enhancement_float(entries, "area_factors"),
        caustic_thresholds=_stack_enhancement_float(entries, "caustic_thresholds"),
        enhancement=_stack_enhancement_float(entries, "enhancement"),
        log10_enhancement=_stack_enhancement_float(entries, "log10_enhancement"),
        valid_mask=_stack_enhancement_bool(entries, "valid_mask"),
        resolved_mask=_stack_enhancement_bool(entries, "resolved_mask"),
        finite_mask=_stack_enhancement_bool(entries, "finite_mask"),
        caustic_candidate_mask=_stack_enhancement_bool(
            entries, "caustic_candidate_mask"
        ),
    )


def build_group_field(
    results: Iterable[DirectionalGroupResult],
) -> GroupFieldResult:
    """Stack directional group solutions into a sampled field.

    Parameters
    ----------
    results : iterable of DirectionalGroupResult
        Directional solutions in deterministic traversal order.

    Returns
    -------
    GroupFieldResult
        Read-only sampled phase and group field.

    Raises
    ------
    TypeError
        If an item is not a :class:`DirectionalGroupResult`.
    ValueError
        If the iterable is empty.
    """
    entries = tuple(results)
    if not entries:
        raise ValueError("At least one directional group result is required.")
    if not all(isinstance(item, DirectionalGroupResult) for item in entries):
        raise TypeError("All field entries must be DirectionalGroupResult objects.")

    phase = build_phase_field(item.phase for item in entries)
    return GroupFieldResult(
        phase=phase,
        eigenvalue_gradients=_stack_group_float(entries, "eigenvalue_gradients"),
        group_velocities=_stack_group_float(entries, "group_velocities"),
        group_speeds=_stack_group_float(entries, "group_speeds"),
        ray_directions=_stack_group_float(entries, "ray_directions"),
        power_flow_angles=_stack_group_float(entries, "power_flow_angles"),
        valid_mask=_stack_group_bool(entries, "valid_mask"),
        resolved_mask=_stack_group_bool(entries, "resolved_mask"),
    )


def build_phase_field(
    results: Iterable[DirectionalPhaseResult],
) -> PhaseFieldResult:
    """Stack directional phase solutions into a sampled field.

    Parameters
    ----------
    results : iterable of DirectionalPhaseResult
        Directional solutions in the order used by the intended path or grid.

    Returns
    -------
    PhaseFieldResult
        Read-only sampled phase field and acoustic-axis diagnostics.

    Raises
    ------
    TypeError
        If an item is not a :class:`DirectionalPhaseResult`.
    ValueError
        If the iterable is empty.
    """
    entries = tuple(results)
    if not entries:
        raise ValueError("At least one directional phase result is required.")
    if not all(isinstance(item, DirectionalPhaseResult) for item in entries):
        raise TypeError("All field entries must be DirectionalPhaseResult objects.")

    directions = _stack_float(entries, "direction")
    eigenvalues = _stack_float(entries, "eigenvalues")
    phase_speeds = _stack_float(entries, "phase_speeds")
    polarizations = _stack_float(entries, "polarizations")
    mode_gaps = _stack_float(entries, "mode_eigenvalue_gaps")
    mode_relative_gaps = _stack_float(entries, "mode_relative_eigenvalue_gaps")
    pair_gaps = _stack_float(entries, "eigenvalue_gaps")
    pair_relative_gaps = _stack_float(entries, "relative_eigenvalue_gaps")
    valid_mask = _stack_bool(entries, "valid_mask")
    clamped_mask = _stack_bool(entries, "clamped_mask")
    degeneracy_mask = _stack_bool(entries, "degeneracy_mask")
    eigenvalue_thresholds = _readonly_float_copy(
        np.asarray([item.eigenvalue_threshold for item in entries], dtype=float)
    )
    degeneracy_thresholds = _readonly_float_copy(
        np.asarray([item.degeneracy_threshold for item in entries], dtype=float)
    )
    pair_degeneracy_mask = _readonly_bool_copy(
        pair_gaps <= degeneracy_thresholds[:, np.newaxis]
    )

    return PhaseFieldResult(
        directions=directions,
        eigenvalues=eigenvalues,
        phase_speeds=phase_speeds,
        polarizations=polarizations,
        mode_eigenvalue_gaps=mode_gaps,
        mode_relative_eigenvalue_gaps=mode_relative_gaps,
        pair_eigenvalue_gaps=pair_gaps,
        pair_relative_eigenvalue_gaps=pair_relative_gaps,
        valid_mask=valid_mask,
        clamped_mask=clamped_mask,
        degeneracy_mask=degeneracy_mask,
        pair_degeneracy_mask=pair_degeneracy_mask,
        eigenvalue_thresholds=eigenvalue_thresholds,
        degeneracy_thresholds=degeneracy_thresholds,
    )


def _stack_float(
    entries: tuple[DirectionalPhaseResult, ...],
    attribute: str,
) -> NDArray[np.float64]:
    """Stack one floating-point result attribute.

    Parameters
    ----------
    entries : tuple of DirectionalPhaseResult
        Directional results.
    attribute : str
        Name of the array attribute to stack.

    Returns
    -------
    ndarray
        Read-only floating-point array.
    """
    values = [np.asarray(getattr(item, attribute), dtype=float) for item in entries]
    return _readonly_float_copy(np.stack(values, axis=0))


def _stack_bool(
    entries: tuple[DirectionalPhaseResult, ...],
    attribute: str,
) -> NDArray[np.bool_]:
    """Stack one boolean result attribute.

    Parameters
    ----------
    entries : tuple of DirectionalPhaseResult
        Directional results.
    attribute : str
        Name of the array attribute to stack.

    Returns
    -------
    ndarray
        Read-only boolean array.
    """
    values = [np.asarray(getattr(item, attribute), dtype=bool) for item in entries]
    return _readonly_bool_copy(np.stack(values, axis=0))


def _stack_group_float(
    entries: tuple[DirectionalGroupResult, ...],
    attribute: str,
) -> NDArray[np.float64]:
    """Stack one floating-point group-result attribute.

    Parameters
    ----------
    entries : tuple of DirectionalGroupResult
        Directional group results.
    attribute : str
        Name of the array attribute to stack.

    Returns
    -------
    ndarray
        Read-only floating-point array.
    """
    values = [np.asarray(getattr(item, attribute), dtype=float) for item in entries]
    return _readonly_float_copy(np.stack(values, axis=0))


def _stack_group_bool(
    entries: tuple[DirectionalGroupResult, ...],
    attribute: str,
) -> NDArray[np.bool_]:
    """Stack one boolean group-result attribute.

    Parameters
    ----------
    entries : tuple of DirectionalGroupResult
        Directional group results.
    attribute : str
        Name of the array attribute to stack.

    Returns
    -------
    ndarray
        Read-only boolean array.
    """
    values = [np.asarray(getattr(item, attribute), dtype=bool) for item in entries]
    return _readonly_bool_copy(np.stack(values, axis=0))


def _stack_enhancement_float(
    entries: tuple[DirectionalEnhancementResult, ...],
    attribute: str,
) -> NDArray[np.float64]:
    """Stack one floating-point enhancement-result attribute.

    Parameters
    ----------
    entries : tuple of DirectionalEnhancementResult
        Directional enhancement results.
    attribute : str
        Name of the array attribute to stack.

    Returns
    -------
    ndarray
        Read-only floating-point array.
    """
    values = [np.asarray(getattr(item, attribute), dtype=float) for item in entries]
    return _readonly_float_copy(np.stack(values, axis=0))


def _stack_enhancement_bool(
    entries: tuple[DirectionalEnhancementResult, ...],
    attribute: str,
) -> NDArray[np.bool_]:
    """Stack one boolean enhancement-result attribute.

    Parameters
    ----------
    entries : tuple of DirectionalEnhancementResult
        Directional enhancement results.
    attribute : str
        Name of the array attribute to stack.

    Returns
    -------
    ndarray
        Read-only boolean array.
    """
    values = [np.asarray(getattr(item, attribute), dtype=bool) for item in entries]
    return _readonly_bool_copy(np.stack(values, axis=0))


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


def _readonly_float_copy(
    array: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Return a read-only floating-point copy.

    Parameters
    ----------
    array : ndarray
        Source array.

    Returns
    -------
    ndarray
        Independent read-only copy.
    """
    copied = np.array(array, dtype=float, copy=True)
    copied.setflags(write=False)
    return copied


def _readonly_bool_copy(
    array: NDArray[np.bool_],
) -> NDArray[np.bool_]:
    """Return a read-only boolean copy.

    Parameters
    ----------
    array : ndarray
        Source array.

    Returns
    -------
    ndarray
        Independent read-only copy.
    """
    copied = np.array(array, dtype=bool, copy=True)
    copied.setflags(write=False)
    return copied


def _readonly_int_copy(
    array: NDArray[np.int64],
) -> NDArray[np.int64]:
    """Return a read-only integer copy.

    Parameters
    ----------
    array : ndarray
        Source array.

    Returns
    -------
    ndarray
        Independent read-only copy.
    """
    copied = np.array(array, dtype=np.int64, copy=True)
    copied.setflags(write=False)
    return copied
