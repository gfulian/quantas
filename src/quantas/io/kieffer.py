# -*- coding: utf-8 -*-

"""Serialize volume-resolved Kieffer acoustic input data."""

from __future__ import annotations

from typing import Any

import numpy as np

from quantas.models.kieffer import (
    CutoffVolumeSource,
    KiefferCutoffState,
    KiefferVolumeSeries,
)


KIEFFER_INPUT_METHOD = "kieffer-sine-wave"
KIEFFER_INPUT_COMPOSITION = "additional-acoustic-branches"


def kieffer_series_to_mapping(
    series: KiefferVolumeSeries,
    *,
    provenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return the YAML mapping for one Kieffer cutoff series.

    Parameters
    ----------
    series : KiefferVolumeSeries
        Direct or explicitly interpolated acoustic cutoff states.
    provenance : dict, optional
        Enrichment and source-elastic provenance.

    Returns
    -------
    dict
        YAML-safe logical data. Presentation formatting is applied later.
    """
    if not isinstance(series, KiefferVolumeSeries):
        raise TypeError("series must be a KiefferVolumeSeries")
    states: list[dict[str, Any]] = []
    for state in series.states:
        cutoff_frequency = np.asarray(state.frequencies_hz, dtype=np.float64)
        effective_velocity = np.asarray(
            state.effective_velocities_km_s,
            dtype=np.float64,
        )
        source = CutoffVolumeSource(state.source)
        states.append(
            {
                "volume": state.volume,
                "cutoff_frequency": cutoff_frequency.tolist(),
                "effective_velocity": effective_velocity.tolist(),
                "source": source.value,
                "source_elastic_indices": list(state.source_elastic_indices),
                "metadata": dict(state.metadata),
            }
        )
    mapping: dict[str, Any] = {
        "method": KIEFFER_INPUT_METHOD,
        "composition": KIEFFER_INPUT_COMPOSITION,
        "units": {
            "volume": "angstrom^3",
            "cutoff_frequency": "Hz",
            "effective_velocity": "km/s",
            "pressure": "GPa",
        },
        "states": states,
        "metadata": dict(series.metadata),
    }
    if series.interpolation_method is not None:
        mapping["interpolation_method"] = series.interpolation_method
    if provenance:
        mapping["provenance"] = dict(provenance)
    return mapping


def kieffer_series_from_mapping(mapping: Any) -> KiefferVolumeSeries:
    """Restore and validate a Kieffer cutoff series from YAML data.

    Parameters
    ----------
    mapping : object
        Parsed value of the top-level ``kieffer`` YAML field.

    Returns
    -------
    KiefferVolumeSeries
        Validated cutoff series.

    Raises
    ------
    ValueError
        If the model identity, units, states, or numerical values are invalid.
    """
    if not isinstance(mapping, dict):
        raise ValueError("kieffer must be a YAML mapping")
    if mapping.get("method") != KIEFFER_INPUT_METHOD:
        raise ValueError(f"kieffer.method must be {KIEFFER_INPUT_METHOD!r}")
    if mapping.get("composition") != KIEFFER_INPUT_COMPOSITION:
        raise ValueError(
            f"kieffer.composition must be {KIEFFER_INPUT_COMPOSITION!r}"
        )
    units = mapping.get("units")
    expected_units = {
        "volume": "angstrom^3",
        "cutoff_frequency": "Hz",
        "effective_velocity": "km/s",
        "pressure": "GPa",
    }
    if not isinstance(units, dict) or any(
        units.get(key) != value for key, value in expected_units.items()
    ):
        raise ValueError("kieffer units do not match the canonical input units")
    raw_states = mapping.get("states")
    if not isinstance(raw_states, list) or not raw_states:
        raise ValueError("kieffer.states must be a non-empty sequence")
    states: list[KiefferCutoffState] = []
    for index, raw in enumerate(raw_states):
        if not isinstance(raw, dict):
            raise ValueError(f"kieffer state {index} must be a mapping")
        try:
            state = KiefferCutoffState(
                volume=float(raw["volume"]),
                frequencies_hz=np.asarray(raw["cutoff_frequency"], dtype=np.float64),
                effective_velocities_km_s=np.asarray(
                    raw["effective_velocity"], dtype=np.float64
                ),
                source=CutoffVolumeSource(raw.get("source", "direct")),
                source_elastic_indices=tuple(
                    int(value) for value in raw.get("source_elastic_indices", [])
                ),
                metadata=dict(raw.get("metadata", {})),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid kieffer state {index}: {exc}") from exc
        states.append(state)
    metadata = dict(mapping.get("metadata", {}))
    if "provenance" in mapping:
        metadata["input_provenance"] = dict(mapping["provenance"])
    return KiefferVolumeSeries(
        states=tuple(states),
        interpolation_method=mapping.get("interpolation_method"),
        metadata=metadata,
    )


__all__ = [
    "KIEFFER_INPUT_COMPOSITION",
    "KIEFFER_INPUT_METHOD",
    "kieffer_series_from_mapping",
    "kieffer_series_to_mapping",
]
