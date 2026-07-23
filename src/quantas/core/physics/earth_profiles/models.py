# -*- coding: utf-8 -*-

"""Passive contracts for terrestrial pressure-temperature depth profiles."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Protocol, TypeAlias, runtime_checkable

import numpy as np
from numpy.typing import NDArray


FloatArray: TypeAlias = NDArray[np.float64]


@runtime_checkable
class PressureDepthModel(Protocol):
    """Structural contract for pressure as a function of depth.

    Implementations use depth in km and return pressure in GPa.  Numerical
    models must be independent from CLI, rendering, persistence, and workflow
    objects.
    """

    @property
    def name(self) -> str:
        """Return the stable model identifier."""

    @property
    def depth_bounds(self) -> tuple[float, float]:
        """Return the supported closed depth interval in km."""

    @property
    def critical_depths(self) -> tuple[float, ...]:
        """Return depths worth preserving in generated sampling grids."""

    def pressure(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate pressure in GPa at the requested depths."""

    def metadata(self) -> dict[str, Any]:
        """Return recursively serializable scientific provenance."""


@runtime_checkable
class TemperatureDepthModel(Protocol):
    """Structural contract for temperature as a function of depth.

    Implementations use depth in km and return absolute temperature in K.
    """

    @property
    def name(self) -> str:
        """Return the stable model identifier."""

    @property
    def depth_bounds(self) -> tuple[float, float]:
        """Return the supported closed depth interval in km."""

    @property
    def critical_depths(self) -> tuple[float, ...]:
        """Return depths worth preserving in generated sampling grids."""

    def temperature(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate temperature in K at the requested depths."""

    def metadata(self) -> dict[str, Any]:
        """Return recursively serializable scientific provenance."""


@dataclass(slots=True)
class EarthDepthProfile:
    """Terrestrial pressure-temperature path parameterized by depth.

    Parameters
    ----------
    name : str
        Stable profile identifier.
    depth : ndarray
        Strictly increasing geological depths in km.
    pressure : ndarray
        Pressure values in GPa, aligned with ``depth``.
    temperature : ndarray
        Absolute temperatures in K, aligned with ``depth``.
    metadata : dict, optional
        Model composition, numerical settings, and bibliographic provenance.

    Raises
    ------
    ValueError
        If arrays are empty, misaligned, non-finite, non-monotonic in depth,
        or contain negative pressure or temperature values.
    """

    name: str
    depth: FloatArray
    pressure: FloatArray
    temperature: FloatArray
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize arrays to float64 and validate the profile contract."""
        self.name = str(self.name).strip()
        if not self.name:
            raise ValueError("profile name cannot be empty")
        self.depth = np.asarray(self.depth, dtype=np.float64).copy()
        self.pressure = np.asarray(self.pressure, dtype=np.float64).copy()
        self.temperature = np.asarray(self.temperature, dtype=np.float64).copy()
        if self.depth.ndim != 1 or self.depth.size == 0:
            raise ValueError("depth must be a non-empty one-dimensional array")
        if self.pressure.shape != self.depth.shape:
            raise ValueError("pressure and depth must have identical shapes")
        if self.temperature.shape != self.depth.shape:
            raise ValueError("temperature and depth must have identical shapes")
        if not all(
            np.all(np.isfinite(values))
            for values in (self.depth, self.pressure, self.temperature)
        ):
            raise ValueError("profile arrays must contain finite values")
        if np.any(self.depth < 0.0):
            raise ValueError("depth values must be non-negative")
        if np.any(np.diff(self.depth) <= 0.0):
            raise ValueError("depth values must be strictly increasing")
        if np.any(self.pressure < -1.0e-12):
            raise ValueError("pressure values must be non-negative")
        if np.any(self.temperature < 0.0):
            raise ValueError("temperature values must be non-negative")
        self.pressure[np.abs(self.pressure) < 1.0e-14] = 0.0
        self.metadata = dict(self.metadata)


__all__ = [
    "EarthDepthProfile",
    "FloatArray",
    "PressureDepthModel",
    "TemperatureDepthModel",
]
