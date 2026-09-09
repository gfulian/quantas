# -*- coding: utf-8 -*-

"""Applicability checks for Kieffer enrichment of Gamma-only HA data."""

from __future__ import annotations

from quantas.models.kieffer import KiefferVolumeSeries
from quantas.models.kieffer_application import (
    validate_kieffer_phonon_applicability,
)
from quantas.models.phonons import PhononInputData
from quantas.models.volume_matching import (
    VolumeMatch,
    VolumeMatchPolicy,
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
    matches = validate_kieffer_phonon_applicability(
        input_data,
        cutoff_series,
        workflow="ha",
        q_tolerance=q_tolerance,
        volume_policy=volume_policy,
    )
    return matches[0]


__all__ = ["validate_kieffer_ha_applicability"]
