# -*- coding: utf-8 -*-

"""Composition of independent terrestrial pressure and temperature models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, Sequence

import numpy as np
from numpy.typing import NDArray

from .models import (
    EarthDepthProfile,
    FloatArray,
    PressureDepthModel,
    TemperatureDepthModel,
)
from .pressure import _validated_depth


JoinMode = Literal["direct", "continuous-offset", "blend"]


@dataclass(frozen=True, slots=True)
class TemperatureSegment:
    """Passive specification of one piecewise temperature-model segment.

    Parameters
    ----------
    depth_min_km, depth_max_km : float
        Closed segment bounds in km.  Adjacent segments must share a boundary.
    model : TemperatureDepthModel
        Temperature model evaluated inside the segment.
    join_mode : {"direct", "continuous-offset", "blend"}, optional
        Policy applied when entering this segment from the previous one.
    blend_width_km : float, optional
        Width centered on the shared boundary for ``blend`` joins.

    Raises
    ------
    ValueError
        If bounds or join settings are invalid.
    """

    depth_min_km: float
    depth_max_km: float
    model: TemperatureDepthModel
    join_mode: JoinMode = "direct"
    blend_width_km: float = 0.0

    def __post_init__(self) -> None:
        """Validate the passive segment specification."""
        if (
            not np.isfinite(self.depth_min_km)
            or not np.isfinite(self.depth_max_km)
            or self.depth_min_km < 0.0
            or self.depth_max_km <= self.depth_min_km
        ):
            raise ValueError("segment depths must be finite and increasing")
        if self.join_mode not in ("direct", "continuous-offset", "blend"):
            raise ValueError("unknown temperature-segment join mode")
        if not np.isfinite(self.blend_width_km) or self.blend_width_km < 0.0:
            raise ValueError("blend_width_km must be finite and non-negative")
        if self.join_mode == "blend" and self.blend_width_km <= 0.0:
            raise ValueError("blend joins require a positive blend_width_km")
        lower, upper = self.model.depth_bounds
        tolerance = max(1.0, abs(self.depth_max_km)) * 1.0e-10
        if (
            self.depth_min_km < lower - tolerance
            or self.depth_max_km > upper + tolerance
        ):
            raise ValueError(
                f"segment [{self.depth_min_km:g}, {self.depth_max_km:g}] km "
                f"exceeds model '{self.model.name}' bounds [{lower:g}, {upper:g}] km"
            )


class PiecewiseTemperatureModel:
    """Temperature-depth model assembled from contiguous physical segments.

    ``direct`` preserves each source model exactly and may retain a temperature
    jump.  ``continuous-offset`` adds a constant to the entering segment so its
    first value matches the previous segment.  ``blend`` linearly combines the
    two source models across a finite interval centered on the join; both source
    models must support the full blend interval.  Every transformation is
    recorded in metadata.

    Parameters
    ----------
    segments : sequence of TemperatureSegment
        Contiguous segments ordered by increasing depth.
    name : str, optional
        Stable model identifier.

    Raises
    ------
    ValueError
        If segments are absent, unordered, overlapping, gapped, or cannot
        satisfy their join policy.
    """

    def __init__(
        self,
        segments: Sequence[TemperatureSegment],
        *,
        name: str = "piecewise-temperature",
    ) -> None:
        self.segments = tuple(segments)
        if not self.segments:
            raise ValueError("piecewise temperature model requires segments")
        self._name = str(name)
        for index, segment in enumerate(self.segments[1:], start=1):
            previous = self.segments[index - 1]
            tolerance = max(1.0, abs(segment.depth_min_km)) * 1.0e-10
            if not np.isclose(
                previous.depth_max_km,
                segment.depth_min_km,
                rtol=0.0,
                atol=tolerance,
            ):
                raise ValueError("piecewise temperature segments must be contiguous")
            if segment.join_mode == "blend":
                half = 0.5 * segment.blend_width_km
                blend_min = segment.depth_min_km - half
                blend_max = segment.depth_min_km + half
                if (
                    blend_min < previous.model.depth_bounds[0]
                    or blend_max > previous.model.depth_bounds[1]
                    or blend_min < segment.model.depth_bounds[0]
                    or blend_max > segment.model.depth_bounds[1]
                ):
                    raise ValueError("both models must support the full blend interval")
        self._offsets = self._calculate_offsets()

    @property
    def name(self) -> str:
        """Return the stable temperature-model identifier."""
        return self._name

    @property
    def depth_bounds(self) -> tuple[float, float]:
        """Return the full piecewise depth interval in km."""
        return self.segments[0].depth_min_km, self.segments[-1].depth_max_km

    @property
    def critical_depths(self) -> tuple[float, ...]:
        """Return segment joins and nested model critical depths."""
        values: set[float] = set()
        for segment in self.segments:
            values.add(float(segment.depth_min_km))
            values.add(float(segment.depth_max_km))
            values.update(
                value
                for value in segment.model.critical_depths
                if segment.depth_min_km <= value <= segment.depth_max_km
            )
            if segment.join_mode == "blend":
                half = 0.5 * segment.blend_width_km
                values.add(segment.depth_min_km - half)
                values.add(segment.depth_min_km + half)
        values.discard(self.depth_bounds[0])
        values.discard(self.depth_bounds[1])
        return tuple(sorted(values))

    def temperature(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate the composed temperature model in K."""
        depth = _validated_depth(depth_km, self.depth_bounds, self.name)
        flat = depth.ravel()
        result = np.empty_like(flat)
        assigned = np.zeros_like(flat, dtype=np.bool_)
        for index, segment in enumerate(self.segments):
            if index == len(self.segments) - 1:
                mask = (flat >= segment.depth_min_km) & (flat <= segment.depth_max_km)
            else:
                mask = (flat >= segment.depth_min_km) & (flat < segment.depth_max_km)
            if not np.any(mask):
                continue
            values = segment.model.temperature(flat[mask]) + self._offsets[index]
            result[mask] = values
            assigned[mask] = True
        if not np.all(assigned):
            raise RuntimeError("internal piecewise temperature assignment failure")
        for index, segment in enumerate(self.segments[1:], start=1):
            if segment.join_mode != "blend":
                continue
            previous = self.segments[index - 1]
            boundary = segment.depth_min_km
            half = 0.5 * segment.blend_width_km
            mask = (flat >= boundary - half) & (flat <= boundary + half)
            if not np.any(mask):
                continue
            blend_depth = flat[mask]
            left = previous.model.temperature(blend_depth) + self._offsets[index - 1]
            right = segment.model.temperature(blend_depth) + self._offsets[index]
            weight = (blend_depth - (boundary - half)) / (2.0 * half)
            result[mask] = (1.0 - weight) * left + weight * right
        return result.reshape(depth.shape).astype(np.float64)

    def metadata(self) -> dict[str, Any]:
        """Return segment composition and explicit join transformations."""
        return {
            "model": self.name,
            "kind": "piecewise_temperature",
            "segments": [
                {
                    "depth_min_km": segment.depth_min_km,
                    "depth_max_km": segment.depth_max_km,
                    "join_mode": segment.join_mode,
                    "blend_width_km": segment.blend_width_km,
                    "applied_offset_K": self._offsets[index],
                    "source_model": segment.model.metadata(),
                }
                for index, segment in enumerate(self.segments)
            ],
        }

    def _calculate_offsets(self) -> tuple[float, ...]:
        offsets = [0.0]
        for index, segment in enumerate(self.segments[1:], start=1):
            previous = self.segments[index - 1]
            if segment.join_mode != "continuous-offset":
                offsets.append(0.0)
                continue
            boundary = segment.depth_min_km
            shallow_depth = np.nextafter(boundary, -np.inf)
            shallow_depth = max(shallow_depth, previous.depth_min_km)
            previous_value = (
                float(previous.model.temperature(np.asarray([shallow_depth]))[0])
                + offsets[index - 1]
            )
            entering_value = float(segment.model.temperature(np.asarray([boundary]))[0])
            offsets.append(previous_value - entering_value)
        return tuple(offsets)


class EarthProfileModel:
    """Compose independent pressure and temperature models on one depth grid.

    Parameters
    ----------
    pressure_model : PressureDepthModel
        Model providing pressure in GPa as a function of depth.
    temperature_model : TemperatureDepthModel
        Model providing temperature in K as a function of depth.
    name : str, optional
        Stable composed-model identifier.

    Raises
    ------
    ValueError
        If the two model domains do not overlap.
    """

    def __init__(
        self,
        pressure_model: PressureDepthModel,
        temperature_model: TemperatureDepthModel,
        *,
        name: str = "earth-profile",
    ) -> None:
        self.pressure_model = pressure_model
        self.temperature_model = temperature_model
        self.name = str(name)
        self._bounds = (
            max(pressure_model.depth_bounds[0], temperature_model.depth_bounds[0]),
            min(pressure_model.depth_bounds[1], temperature_model.depth_bounds[1]),
        )
        if self._bounds[1] <= self._bounds[0]:
            raise ValueError("pressure and temperature model domains do not overlap")

    @property
    def depth_bounds(self) -> tuple[float, float]:
        """Return the common pressure-temperature depth domain in km."""
        return self._bounds

    @property
    def critical_depths(self) -> tuple[float, ...]:
        """Return all source-model critical depths within the common domain."""
        values = set(self.pressure_model.critical_depths)
        values.update(self.temperature_model.critical_depths)
        return tuple(
            sorted(
                value for value in values if self._bounds[0] < value < self._bounds[1]
            )
        )

    def evaluate(
        self,
        depth_km: NDArray[np.float64],
        *,
        profile_name: str | None = None,
    ) -> EarthDepthProfile:
        """Evaluate the composed model on an explicit depth array.

        Parameters
        ----------
        depth_km : ndarray
            Strictly increasing depths in km.
        profile_name : str or None, optional
            Result profile identifier.

        Returns
        -------
        EarthDepthProfile
            Aligned pressure-temperature-depth values and full provenance.

        Raises
        ------
        ValueError
            If depths are invalid or outside the common model domain.
        """
        depth = _validated_depth(depth_km, self.depth_bounds, self.name)
        if depth.ndim != 1 or depth.size == 0:
            raise ValueError("profile depth must be a non-empty one-dimensional array")
        if np.any(np.diff(depth) <= 0.0):
            raise ValueError("profile depth must be strictly increasing")
        return EarthDepthProfile(
            name=self.name if profile_name is None else profile_name,
            depth=depth,
            pressure=self.pressure_model.pressure(depth),
            temperature=self.temperature_model.temperature(depth),
            metadata=self.metadata(),
        )

    def regular_profile(
        self,
        *,
        depth_min_km: float | None = None,
        depth_max_km: float | None = None,
        depth_step_km: float = 1.0,
        include_critical_depths: bool = True,
        profile_name: str | None = None,
    ) -> EarthDepthProfile:
        """Generate and evaluate a regular depth grid.

        Parameters
        ----------
        depth_min_km, depth_max_km : float or None, optional
            Requested limits.  Common model bounds are used when omitted.
        depth_step_km : float, optional
            Positive regular spacing in km.
        include_critical_depths : bool, optional
            Insert model boundaries and transition-bracketing points even when
            they do not lie on the regular grid.
        profile_name : str or None, optional
            Result profile identifier.

        Returns
        -------
        EarthDepthProfile
            Generated pressure-temperature-depth path.

        Raises
        ------
        ValueError
            If limits or spacing are invalid.
        """
        lower = self._bounds[0] if depth_min_km is None else float(depth_min_km)
        upper = self._bounds[1] if depth_max_km is None else float(depth_max_km)
        step = float(depth_step_km)
        if not np.isfinite(step) or step <= 0.0:
            raise ValueError("depth_step_km must be finite and positive")
        if lower < self._bounds[0] or upper > self._bounds[1] or upper <= lower:
            raise ValueError(
                f"requested range must lie within [{self._bounds[0]:g}, "
                f"{self._bounds[1]:g}] km"
            )
        count = int(np.floor((upper - lower) / step + 1.0e-12)) + 1
        depth = lower + step * np.arange(count, dtype=np.float64)
        tolerance = max(1.0, abs(upper)) * 1.0e-12
        if depth[-1] < upper - tolerance:
            depth = np.append(depth, upper)
        else:
            depth[-1] = upper
        if include_critical_depths:
            critical = np.asarray(
                [value for value in self.critical_depths if lower < value < upper],
                dtype=np.float64,
            )
            if critical.size:
                depth = np.unique(np.concatenate((depth, critical)))
        return self.evaluate(depth, profile_name=profile_name)

    def metadata(self) -> dict[str, Any]:
        """Return complete pressure-temperature composition provenance."""
        return {
            "kind": "composed_earth_profile",
            "profile_schema_version": 1,
            "name": self.name,
            "depth_unit": "km",
            "pressure_unit": "GPa",
            "temperature_unit": "K",
            "depth_bounds_km": self.depth_bounds,
            "pressure": self.pressure_model.metadata(),
            "temperature": self.temperature_model.metadata(),
        }


__all__ = [
    "EarthProfileModel",
    "JoinMode",
    "PiecewiseTemperatureModel",
    "TemperatureSegment",
]
