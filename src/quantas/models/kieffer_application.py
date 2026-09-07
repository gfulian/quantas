# -*- coding: utf-8 -*-

"""Shared applicability validation for Kieffer-enriched phonon inputs."""

from __future__ import annotations

from typing import Literal

import numpy as np

from quantas.models.kieffer import CutoffVolumeSource, KiefferVolumeSeries
from quantas.models.phonons import PhononInputData
from quantas.models.volume_matching import (
    VolumeMatch,
    VolumeMatchPolicy,
    match_sampled_volumes,
)


KiefferWorkflow = Literal["ha", "qha"]


def validate_kieffer_phonon_applicability(
    input_data: PhononInputData,
    cutoff_series: KiefferVolumeSeries,
    *,
    workflow: KiefferWorkflow,
    q_tolerance: float = 1.0e-10,
    volume_policy: VolumeMatchPolicy | None = None,
) -> tuple[VolumeMatch, ...]:
    """Validate a primitive Gamma-only phonon and acoustic-cutoff pairing.

    Parameters
    ----------
    input_data : PhononInputData
        HA or QHA phonon data to enrich.
    cutoff_series : KiefferVolumeSeries
        Direct acoustic cutoff state or volume series.
    workflow : {"ha", "qha"}
        Volume-cardinality policy and diagnostic context.
    q_tolerance : float, optional
        Absolute tolerance for identifying Gamma modulo reciprocal vectors.
    volume_policy : VolumeMatchPolicy or None, optional
        Explicit matching policy for independently reported volumes.

    Returns
    -------
    tuple of VolumeMatch
        One unique cutoff association for every input volume.

    Raises
    ------
    TypeError
        If the input contracts have unsupported types.
    ValueError
        If the workflow, Gamma sampling, primitive normalization, cutoff
        provenance, or volume association is invalid.
    """
    if workflow not in {"ha", "qha"}:
        raise ValueError("workflow must be 'ha' or 'qha'")
    label = workflow.upper()
    if not isinstance(input_data, PhononInputData):
        raise TypeError("input_data must be a PhononInputData")
    if not isinstance(cutoff_series, KiefferVolumeSeries):
        raise TypeError("cutoff_series must be a KiefferVolumeSeries")
    if not np.isfinite(q_tolerance) or q_tolerance < 0.0:
        raise ValueError("q_tolerance must be finite and non-negative")

    if workflow == "ha" and input_data.nvol != 1:
        raise ValueError("Kieffer HA enrichment requires exactly one volume")
    if workflow == "qha" and input_data.nvol < 2:
        raise ValueError("Kieffer QHA enrichment requires multiple volumes")
    if input_data.frequencies is None:
        raise ValueError(f"Kieffer {label} enrichment requires phonon frequencies")
    frequencies = np.asarray(input_data.frequencies)
    if frequencies.ndim != 3 or frequencies.shape[0] != 1:
        raise ValueError(
            f"Kieffer {label} enrichment requires exactly one q-point"
        )
    if input_data.qpoints != 1:
        raise ValueError(f"Kieffer {label} enrichment requires qpoints == 1")
    if input_data.qcoords is None:
        raise ValueError(
            f"Kieffer {label} enrichment requires explicit Gamma coordinates"
        )
    qcoords = np.asarray(input_data.qcoords, dtype=np.float64)
    if qcoords.shape != (1, 3) or not np.all(np.isfinite(qcoords)):
        raise ValueError(
            f"Kieffer {label} qcoords must have shape (1, 3) and be finite"
        )
    if np.any(np.abs(qcoords - np.rint(qcoords)) > q_tolerance):
        raise ValueError(f"Kieffer {label} enrichment requires the Gamma q-point")
    if input_data.supercell is None:
        raise ValueError(
            f"Kieffer {label} enrichment requires an explicit 1x1x1 supercell"
        )
    supercell = np.asarray(input_data.supercell, dtype=np.float64)
    if supercell.shape != (3, 3) or not np.array_equal(supercell, np.eye(3)):
        raise ValueError(
            f"Kieffer {label} enrichment requires the identity supercell"
        )

    if workflow == "ha" and len(cutoff_series.states) != 1:
        raise ValueError(
            "Kieffer HA enrichment requires exactly one cutoff state"
        )
    if workflow == "qha" and len(cutoff_series.states) != input_data.nvol:
        raise ValueError(
            "Kieffer QHA requires one cutoff state per sampled volume"
        )
    if any(
        state.source is not CutoffVolumeSource.DIRECT
        for state in cutoff_series.states
    ):
        state_label = (
            "a direct cutoff state"
            if workflow == "ha"
            else "direct cutoff states"
        )
        raise ValueError(f"Kieffer {label} enrichment requires {state_label}")

    structure = input_data.structure
    if structure is not None:
        normalization = structure.normalization
        if normalization.basis.strip().lower() != "primitive":
            raise ValueError(
                f"Kieffer {label} enrichment requires primitive normalization"
            )
        if normalization.repetitions != 1:
            raise ValueError(
                f"Kieffer {label} enrichment requires one primitive repetition"
            )
        if not np.array_equal(
            normalization.expansion_matrix,
            np.eye(3, dtype=int),
        ):
            raise ValueError(
                f"Kieffer {label} enrichment requires identity cell expansion"
            )
        if workflow == "ha" and structure.nvol != 1:
            raise ValueError(
                "Kieffer HA structure must contain exactly one volume"
            )
        if workflow == "qha" and structure.nvol != input_data.nvol:
            raise ValueError(
                "Kieffer QHA structure and phonons use different volumes"
            )

    if input_data.volume is None:
        volume_label = (
            "an explicit volume" if workflow == "ha" else "explicit volumes"
        )
        raise ValueError(f"Kieffer {label} enrichment requires {volume_label}")
    return match_sampled_volumes(
        np.asarray(input_data.volume, dtype=np.float64),
        cutoff_series.volumes,
        policy=volume_policy,
    )


__all__ = ["KiefferWorkflow", "validate_kieffer_phonon_applicability"]
