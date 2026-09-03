# -*- coding: utf-8 -*-

"""Passive volume-resolved contracts for Kieffer acoustic cutoffs."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray


class CutoffVolumeSource(str, Enum):
    """Identify how one volume-resolved cutoff state was obtained."""

    DIRECT = "direct"
    INTERPOLATED = "interpolated"


@dataclass(slots=True)
class KiefferCutoffState:
    """Three Kieffer acoustic cutoffs at one primitive-cell volume.

    Parameters
    ----------
    volume : float
        Primitive-cell volume in angstrom cubed.
    frequencies_hz : array_like
        Three ordinary cutoff frequencies in hertz.
    effective_velocities_km_s : array_like
        Slow-shear, fast-shear, and longitudinal effective velocities.
    source : CutoffVolumeSource or str, optional
        Whether values are direct or interpolated.
    source_elastic_indices : tuple of int, optional
        Elastic states supporting this cutoff state.
    metadata : dict, optional
        Quadrature, interpolation, and source provenance.

    Raises
    ------
    ValueError
        If values or source indices are invalid.
    """

    volume: float
    frequencies_hz: ArrayLike
    effective_velocities_km_s: ArrayLike
    source: CutoffVolumeSource | str = CutoffVolumeSource.DIRECT
    source_elastic_indices: tuple[int, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize arrays and validate one cutoff state."""
        self.volume = float(self.volume)
        if not np.isfinite(self.volume) or self.volume <= 0.0:
            raise ValueError("volume must be finite and positive")
        self.frequencies_hz = _positive_triplet(self.frequencies_hz, "frequencies_hz")
        self.effective_velocities_km_s = _positive_triplet(
            self.effective_velocities_km_s, "effective_velocities_km_s"
        )
        self.source = CutoffVolumeSource(self.source)
        self.source_elastic_indices = tuple(
            int(value) for value in self.source_elastic_indices
        )
        if any(value < 0 for value in self.source_elastic_indices):
            raise ValueError("source elastic indices must be non-negative")
        if (
            self.source is CutoffVolumeSource.DIRECT
            and len(self.source_elastic_indices) != 1
        ):
            raise ValueError("a direct cutoff state requires one source elastic index")
        if (
            self.source is CutoffVolumeSource.INTERPOLATED
            and len(self.source_elastic_indices) < 2
        ):
            raise ValueError(
                "an interpolated cutoff state requires at least two source indices"
            )
        self.metadata = dict(self.metadata)


@dataclass(slots=True)
class KiefferVolumeSeries:
    """Increasing volume series of Kieffer acoustic cutoff states.

    Parameters
    ----------
    states : tuple of KiefferCutoffState
        Cutoff states ordered by increasing volume.
    interpolation_method : str or None, optional
        Method used when at least one state is interpolated.
    metadata : dict, optional
        Series-level provenance.

    Raises
    ------
    ValueError
        If ordering or interpolation provenance is invalid.
    """

    states: tuple[KiefferCutoffState, ...]
    interpolation_method: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate ordering and interpolation metadata."""
        self.states = tuple(self.states)
        if not self.states:
            raise ValueError("Kieffer cutoff series cannot be empty")
        if np.any(np.diff(self.volumes) <= 0.0):
            raise ValueError("Kieffer states must have unique increasing volumes")
        interpolated = any(
            state.source is CutoffVolumeSource.INTERPOLATED for state in self.states
        )
        method = (
            None
            if self.interpolation_method is None
            else self.interpolation_method.strip()
        )
        if interpolated and not method:
            raise ValueError(
                "interpolated cutoff states require an interpolation method"
            )
        if not interpolated and method:
            raise ValueError(
                "interpolation method is invalid for an entirely direct series"
            )
        self.interpolation_method = method
        self.metadata = dict(self.metadata)

    @property
    def volumes(self) -> NDArray[np.float64]:
        """Return primitive-cell volumes in angstrom cubed."""
        return np.asarray([state.volume for state in self.states], dtype=np.float64)

    @property
    def frequencies_hz(self) -> NDArray[np.float64]:
        """Return cutoffs with shape ``(3, nvol)`` for thermodynamic functions."""
        return np.column_stack([state.frequencies_hz for state in self.states])

    @property
    def effective_velocities_km_s(self) -> NDArray[np.float64]:
        """Return effective velocities with shape ``(3, nvol)``."""
        return np.column_stack(
            [state.effective_velocities_km_s for state in self.states]
        )


def _positive_triplet(values: ArrayLike, name: str) -> NDArray[np.float64]:
    """Return a copied read-only positive floating-point triplet."""
    array = np.asarray(values, dtype=np.float64)
    if array.shape != (3,) or not np.all(np.isfinite(array)) or np.any(array <= 0.0):
        raise ValueError(f"{name} must contain three finite positive values")
    result = array.copy()
    result.setflags(write=False)
    return result


__all__ = ["CutoffVolumeSource", "KiefferCutoffState", "KiefferVolumeSeries"]
