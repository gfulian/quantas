# -*- coding: utf-8 -*-

"""Applicability checks for Kieffer enrichment of Gamma-only HA data."""

from __future__ import annotations

import numpy as np

from quantas.models.kieffer import CutoffVolumeSource, KiefferVolumeSeries
from quantas.models.phonons import PhononInputData
from quantas.models.volume_matching import (
    VolumeMatch,
    VolumeMatchPolicy,
    match_sampled_volumes,
)


def validate_kieffer_ha_applicability(
    input_data: PhononInputData,
    cutoff_series: KiefferVolumeSeries,
    *,
    q_tolerance: float = 1.0e-10,
    volume_policy: VolumeMatchPolicy | None = None,
) -> VolumeMatch:
    """Validate a single-volume, primitive Gamma HA enrichment.

    The acoustic branches represented by ``cutoff_series`` are additional to
    the calculated Gamma frequencies.  This function deliberately does not
    select, remove, or replace any phonon mode.

    Parameters
    ----------
    input_data : PhononInputData
        HA phonon data to enrich.
    cutoff_series : KiefferVolumeSeries
        A single direct cutoff state at the same primitive-cell volume.
    q_tolerance : float, optional
        Absolute tolerance for identifying Gamma modulo reciprocal vectors.
    volume_policy : VolumeMatchPolicy or None, optional
        Explicit tolerance policy for the independently reported volumes.

    Returns
    -------
    VolumeMatch
        Traceable association between the HA and cutoff volumes.

    Raises
    ------
    TypeError
        If either input has an unsupported type.
    ValueError
        If the calculation is not single-volume, primitive, Gamma-only, or if
        the cutoff state is not a direct volume match.
    """
    if not isinstance(input_data, PhononInputData):
        raise TypeError("input_data must be a PhononInputData")
    if not isinstance(cutoff_series, KiefferVolumeSeries):
        raise TypeError("cutoff_series must be a KiefferVolumeSeries")
    if not np.isfinite(q_tolerance) or q_tolerance < 0.0:
        raise ValueError("q_tolerance must be finite and non-negative")
    if input_data.nvol != 1:
        raise ValueError("Kieffer HA enrichment requires exactly one volume")
    if input_data.frequencies is None:
        raise ValueError("Kieffer HA enrichment requires phonon frequencies")
    frequencies = np.asarray(input_data.frequencies)
    if frequencies.ndim != 3 or frequencies.shape[0] != 1:
        raise ValueError("Kieffer HA enrichment requires exactly one q-point")
    if input_data.qpoints != 1:
        raise ValueError("Kieffer HA enrichment requires qpoints == 1")
    if input_data.qcoords is None:
        raise ValueError("Kieffer HA enrichment requires explicit Gamma coordinates")
    qcoords = np.asarray(input_data.qcoords, dtype=np.float64)
    if qcoords.shape != (1, 3) or not np.all(np.isfinite(qcoords)):
        raise ValueError("Kieffer HA qcoords must have shape (1, 3) and be finite")
    gamma_residual = qcoords - np.rint(qcoords)
    if np.any(np.abs(gamma_residual) > q_tolerance):
        raise ValueError("Kieffer HA enrichment requires the Gamma q-point")
    if input_data.supercell is None:
        raise ValueError("Kieffer HA enrichment requires an explicit 1x1x1 supercell")
    supercell = np.asarray(input_data.supercell, dtype=np.float64)
    if supercell.shape != (3, 3) or not np.array_equal(supercell, np.eye(3)):
        raise ValueError("Kieffer HA enrichment requires the identity supercell")
    if len(cutoff_series.states) != 1:
        raise ValueError("Kieffer HA enrichment requires exactly one cutoff state")
    if cutoff_series.states[0].source is not CutoffVolumeSource.DIRECT:
        raise ValueError("Kieffer HA enrichment requires a direct cutoff state")

    structure = input_data.structure
    if structure is not None:
        normalization = structure.normalization
        if normalization.basis.strip().lower() != "primitive":
            raise ValueError("Kieffer HA enrichment requires primitive normalization")
        if normalization.repetitions != 1:
            raise ValueError("Kieffer HA enrichment requires one primitive repetition")
        if not np.array_equal(normalization.expansion_matrix, np.eye(3, dtype=int)):
            raise ValueError("Kieffer HA enrichment requires identity cell expansion")
        if structure.nvol != 1:
            raise ValueError("Kieffer HA structure must contain exactly one volume")

    if input_data.volume is None:
        raise ValueError("Kieffer HA enrichment requires an explicit volume")
    matches = match_sampled_volumes(
        np.asarray(input_data.volume, dtype=np.float64),
        cutoff_series.volumes,
        policy=volume_policy,
    )
    return matches[0]


__all__ = ["validate_kieffer_ha_applicability"]
