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
    """Validate primitive Gamma-only QHA data and match every cutoff volume.

    Every QHA volume must have one unique direct cutoff state.  Source states
    may be ordered independently; the returned matches define the explicit
    association used by the thermodynamic calculation.
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
    """Return cutoff and velocity arrays reordered onto the QHA volume axis."""
    source_indices = [match.source_index for match in matches]
    return (
        cutoff_series.frequencies_hz[:, source_indices].copy(),
        cutoff_series.effective_velocities_km_s[:, source_indices].copy(),
    )


__all__ = ["matched_kieffer_arrays", "validate_kieffer_qha_applicability"]
