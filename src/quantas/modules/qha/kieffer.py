# -*- coding: utf-8 -*-

"""Applicability and matching for Kieffer-enriched QHA calculations."""

from __future__ import annotations

import numpy as np

from quantas.models.kieffer import CutoffVolumeSource, KiefferVolumeSeries
from quantas.models.phonons import PhononInputData
from quantas.models.volume_matching import (
    VolumeMatch,
    VolumeMatchPolicy,
    match_sampled_volumes,
)


def validate_kieffer_qha_applicability(
    input_data: PhononInputData,
    cutoff_series: KiefferVolumeSeries,
    *,
    q_tolerance: float = 1.0e-10,
    volume_policy: VolumeMatchPolicy | None = None,
) -> tuple[VolumeMatch, ...]:
    """Validate primitive Gamma-only QHA data and match every cutoff volume.

    Every QHA volume must have one unique direct cutoff state.  Source states
    may be ordered independently; the returned matches define the explicit
    association used by the thermodynamic calculation.
    """
    if not isinstance(input_data, PhononInputData):
        raise TypeError("input_data must be a PhononInputData")
    if not isinstance(cutoff_series, KiefferVolumeSeries):
        raise TypeError("cutoff_series must be a KiefferVolumeSeries")
    if not np.isfinite(q_tolerance) or q_tolerance < 0.0:
        raise ValueError("q_tolerance must be finite and non-negative")
    if input_data.nvol < 2:
        raise ValueError("Kieffer QHA enrichment requires multiple volumes")
    if input_data.frequencies is None:
        raise ValueError("Kieffer QHA enrichment requires phonon frequencies")
    frequencies = np.asarray(input_data.frequencies)
    if frequencies.ndim != 3 or frequencies.shape[0] != 1:
        raise ValueError("Kieffer QHA enrichment requires exactly one q-point")
    if input_data.qpoints != 1:
        raise ValueError("Kieffer QHA enrichment requires qpoints == 1")
    if input_data.qcoords is None:
        raise ValueError("Kieffer QHA enrichment requires explicit Gamma coordinates")
    qcoords = np.asarray(input_data.qcoords, dtype=np.float64)
    if qcoords.shape != (1, 3) or not np.all(np.isfinite(qcoords)):
        raise ValueError("Kieffer QHA qcoords must have shape (1, 3) and be finite")
    if np.any(np.abs(qcoords - np.rint(qcoords)) > q_tolerance):
        raise ValueError("Kieffer QHA enrichment requires the Gamma q-point")
    if input_data.supercell is None:
        raise ValueError("Kieffer QHA enrichment requires an explicit 1x1x1 supercell")
    supercell = np.asarray(input_data.supercell, dtype=np.float64)
    if supercell.shape != (3, 3) or not np.array_equal(supercell, np.eye(3)):
        raise ValueError("Kieffer QHA enrichment requires the identity supercell")
    if len(cutoff_series.states) != input_data.nvol:
        raise ValueError("Kieffer QHA requires one cutoff state per sampled volume")
    if any(
        state.source is not CutoffVolumeSource.DIRECT for state in cutoff_series.states
    ):
        raise ValueError("Kieffer QHA enrichment requires direct cutoff states")

    structure = input_data.structure
    if structure is not None:
        normalization = structure.normalization
        if normalization.basis.strip().lower() != "primitive":
            raise ValueError("Kieffer QHA enrichment requires primitive normalization")
        if normalization.repetitions != 1:
            raise ValueError("Kieffer QHA enrichment requires one primitive repetition")
        if not np.array_equal(normalization.expansion_matrix, np.eye(3, dtype=int)):
            raise ValueError("Kieffer QHA enrichment requires identity cell expansion")
        if structure.nvol != input_data.nvol:
            raise ValueError("Kieffer QHA structure and phonons use different volumes")

    if input_data.volume is None:
        raise ValueError("Kieffer QHA enrichment requires explicit volumes")
    return match_sampled_volumes(
        np.asarray(input_data.volume, dtype=np.float64),
        cutoff_series.volumes,
        policy=volume_policy,
    )


def matched_kieffer_arrays(
    cutoff_series: KiefferVolumeSeries,
    matches: tuple[VolumeMatch, ...],
) -> tuple[np.ndarray, np.ndarray]:
    """Return cutoff and velocity arrays reordered onto the QHA volume axis."""
    source_indices = [match.source_index for match in matches]
    return (
        cutoff_series.frequencies_hz[:, source_indices].copy(),
        cutoff_series.effective_velocities_km_s[:, source_indices].copy(),
    )


__all__ = ["matched_kieffer_arrays", "validate_kieffer_qha_applicability"]
