# -*- coding: utf-8 -*-

"""Explicit matching of independently reported primitive-cell volumes."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike


@dataclass(frozen=True, slots=True)
class VolumeMatchPolicy:
    """Numerical policy for matching primitive-cell volumes.

    Parameters
    ----------
    relative_tolerance : float, optional
        Relative tolerance applied after unit normalization.
    absolute_tolerance : float, optional
        Absolute tolerance in cubic angstrom.
    require_unique : bool, optional
        Whether multiple source candidates constitute an error.

    Raises
    ------
    ValueError
        If a tolerance is negative or non-finite.
    """

    relative_tolerance: float = 1.0e-8
    absolute_tolerance: float = 1.0e-6
    require_unique: bool = True

    def __post_init__(self) -> None:
        """Validate tolerances."""
        for name, value in (
            ("relative_tolerance", self.relative_tolerance),
            ("absolute_tolerance", self.absolute_tolerance),
        ):
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative")


@dataclass(frozen=True, slots=True)
class VolumeMatch:
    """One traceable association between target and source volumes."""

    target_index: int
    source_index: int
    target_volume: float
    source_volume: float
    absolute_difference: float
    relative_difference: float


def match_sampled_volumes(
    target_volumes: ArrayLike,
    source_volumes: ArrayLike,
    *,
    policy: VolumeMatchPolicy | None = None,
) -> tuple[VolumeMatch, ...]:
    """Match every target volume to an explicitly unique source volume.

    Parameters
    ----------
    target_volumes, source_volumes : array_like
        One-dimensional primitive-cell volumes in cubic angstrom.
    policy : VolumeMatchPolicy or None, optional
        Explicit tolerance and uniqueness policy.

    Returns
    -------
    tuple of VolumeMatch
        Matches in target order with complete numerical diagnostics.

    Raises
    ------
    TypeError
        If ``policy`` has an unsupported type.
    ValueError
        If arrays are invalid, a target is unmatched, or a required unique
        match is ambiguous.
    """
    selected = VolumeMatchPolicy() if policy is None else policy
    if not isinstance(selected, VolumeMatchPolicy):
        raise TypeError("policy must be a VolumeMatchPolicy")
    targets = _volumes(target_volumes, "target_volumes")
    sources = _volumes(source_volumes, "source_volumes")
    matches: list[VolumeMatch] = []
    used: set[int] = set()
    for target_index, target in enumerate(targets):
        candidates = np.flatnonzero(
            np.isclose(
                sources,
                target,
                rtol=selected.relative_tolerance,
                atol=selected.absolute_tolerance,
            )
        )
        if candidates.size == 0:
            raise ValueError(f"target volume {target:.12g} has no source match")
        if selected.require_unique and candidates.size != 1:
            raise ValueError(f"target volume {target:.12g} has multiple source matches")
        source_index = int(candidates[np.argmin(np.abs(sources[candidates] - target))])
        if selected.require_unique and source_index in used:
            raise ValueError("one source volume cannot match multiple target volumes")
        used.add(source_index)
        difference = abs(float(sources[source_index] - target))
        matches.append(
            VolumeMatch(
                target_index=target_index,
                source_index=source_index,
                target_volume=float(target),
                source_volume=float(sources[source_index]),
                absolute_difference=difference,
                relative_difference=difference / abs(float(target)),
            )
        )
    return tuple(matches)


def _volumes(values: ArrayLike, name: str) -> np.ndarray:
    """Return a non-empty positive one-dimensional volume array."""
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"{name} must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(array)) or np.any(array <= 0.0):
        raise ValueError(f"{name} must contain finite positive values")
    return array


__all__ = ["VolumeMatch", "VolumeMatchPolicy", "match_sampled_volumes"]
