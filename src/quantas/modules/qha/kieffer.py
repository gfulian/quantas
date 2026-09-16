# -*- coding: utf-8 -*-

"""Applicability and matching for Kieffer-enriched QHA calculations."""

from __future__ import annotations

import numpy as np

from quantas.models.kieffer import KiefferVolumeSeries
from quantas.models.kieffer_application import (
    validate_kieffer_phonon_applicability,
)
from quantas.models.phonons import PhononInputData
from quantas.models.volume_matching import (
    VolumeMatch,
    VolumeMatchPolicy,
)


def validate_kieffer_qha_applicability(
    input_data: PhononInputData,
    cutoff_series: KiefferVolumeSeries,
    *,
    q_tolerance: float = 1.0e-10,
    volume_policy: VolumeMatchPolicy | None = None,
) -> tuple[VolumeMatch, ...]:
    """Validate primitive Gamma-only QHA data and match every Kieffer cutoff volume.

    Every sampled QHA state must represent primitive-cell phonons at the single
    Gamma q-point with an identity phonon supercell. Each QHA volume must also have
    one unique *direct* cutoff state. Source cutoff states may be ordered
    independently; the returned matches define the explicit association used by the
    thermodynamic calculation and prevent accidental near-volume guessing.

    Parameters
    ----------
    input_data : PhononInputData
        Multi-volume QHA phonon data to enrich.
    cutoff_series : KiefferVolumeSeries
        Direct Kieffer cutoff states covering all sampled QHA volumes.
    q_tolerance : float, optional
        Absolute tolerance for identifying Gamma modulo reciprocal lattice vectors.
    volume_policy : VolumeMatchPolicy or None, optional
        Explicit tolerance policy for matching independently reported volumes.

    Returns
    -------
    tuple of VolumeMatch
        One traceable cutoff association for every QHA volume, ordered on the QHA
        volume axis.

    Raises
    ------
    TypeError
        If an input object has an unsupported type.
    ValueError
        If phonons are not primitive Gamma-only, the supercell is not identity,
        cutoff states are interpolated or incomplete, or a unique volume match
        cannot be established.
    """
    return validate_kieffer_phonon_applicability(
        input_data,
        cutoff_series,
        workflow="qha",
        q_tolerance=q_tolerance,
        volume_policy=volume_policy,
    )


def matched_kieffer_arrays(
    cutoff_series: KiefferVolumeSeries,
    matches: tuple[VolumeMatch, ...],
) -> tuple[np.ndarray, np.ndarray]:
    """Return Kieffer cutoff and velocity arrays reordered onto the QHA volume axis.

    Parameters
    ----------
    cutoff_series : KiefferVolumeSeries
        Source cutoff states. Their internal volume ordering need not match the QHA
        input ordering.
    matches : tuple of VolumeMatch
        Explicit associations returned by :func:`validate_kieffer_qha_applicability`.

    Returns
    -------
    tuple of ndarray
        ``(cutoff_frequencies_hz, effective_velocities_km_s)`` with shape
        ``(3, nvol)``. Acoustic branch order is ``(S1, S2, P)``.
    """
    source_indices = [match.source_index for match in matches]
    return (
        cutoff_series.frequencies_hz[:, source_indices].copy(),
        cutoff_series.effective_velocities_km_s[:, source_indices].copy(),
    )


__all__ = ["matched_kieffer_arrays", "validate_kieffer_qha_applicability"]
