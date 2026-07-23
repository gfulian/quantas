# -*- coding: utf-8 -*-

"""Tabulated pressure and temperature models parameterized by depth."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import PchipInterpolator, interp1d

from .models import FloatArray
from .pressure import _validated_depth


InterpolationKind = Literal["linear", "pchip"]


class _TabulatedDepthModel:
    """Shared validated interpolation for one scalar depth-dependent field."""

    def __init__(
        self,
        depth_km: NDArray[np.float64],
        values: NDArray[np.float64],
        *,
        name: str,
        field_name: str,
        field_unit: str,
        interpolation: InterpolationKind = "pchip",
        source: str | None = None,
        citation: str | None = None,
    ) -> None:
        depth = np.asarray(depth_km, dtype=np.float64)
        data = np.asarray(values, dtype=np.float64)
        if depth.ndim != 1 or depth.size < 2:
            raise ValueError("tabulated depth must contain at least two points")
        if data.shape != depth.shape:
            raise ValueError("tabulated depth and values must have identical shapes")
        if np.any(~np.isfinite(depth)) or np.any(~np.isfinite(data)):
            raise ValueError("tabulated depth and values must be finite")
        order = np.argsort(depth, kind="stable")
        depth = depth[order]
        data = data[order]
        if depth[0] < 0.0 or np.any(np.diff(depth) <= 0.0):
            raise ValueError(
                "tabulated depths must be unique, increasing, and non-negative"
            )
        if np.any(data < 0.0):
            raise ValueError(f"tabulated {field_name} values must be non-negative")
        if interpolation not in ("linear", "pchip"):
            raise ValueError("interpolation must be 'linear' or 'pchip'")
        self._depth = depth.copy()
        self._values = data.copy()
        self._name = str(name)
        self._field_name = field_name
        self._field_unit = field_unit
        self._interpolation = interpolation
        self._source = source
        self._citation = citation
        if interpolation == "pchip":
            self._interpolator = PchipInterpolator(depth, data, extrapolate=False)
        else:
            self._interpolator = interp1d(
                depth,
                data,
                kind="linear",
                bounds_error=True,
                assume_sorted=True,
            )

    @property
    def name(self) -> str:
        """Return the stable model identifier."""
        return self._name

    @property
    def depth_bounds(self) -> tuple[float, float]:
        """Return the tabulated closed depth interval in km."""
        return float(self._depth[0]), float(self._depth[-1])

    @property
    def critical_depths(self) -> tuple[float, ...]:
        """Return every original tabulated depth for exact grid preservation."""
        return tuple(float(value) for value in self._depth)

    def _evaluate(self, depth_km: NDArray[np.float64]) -> FloatArray:
        depth = _validated_depth(depth_km, self.depth_bounds, self.name)
        return np.asarray(self._interpolator(depth), dtype=np.float64)

    def metadata(self) -> dict[str, Any]:
        """Return interpolation settings and input provenance."""
        return {
            "model": self.name,
            "kind": f"tabulated_{self._field_name}",
            "depth_unit": "km",
            f"{self._field_name}_unit": self._field_unit,
            "interpolation": self._interpolation,
            "source": self._source,
            "citation": self._citation or "User-supplied tabulated profile.",
            "source_depth_km": self._depth.copy(),
            f"source_{self._field_name}": self._values.copy(),
        }


class TabulatedPressureModel(_TabulatedDepthModel):
    """Pressure-depth model interpolated from user-supplied values.

    Parameters
    ----------
    depth_km, pressure_GPa : ndarray
        Tabulated depths in km and pressures in GPa.
    name : str, optional
        Stable model identifier.
    interpolation : {"linear", "pchip"}, optional
        Interpolation scheme.  Extrapolation is never permitted.
    source, citation : str or None, optional
        File provenance and complete scientific citation.

    Raises
    ------
    ValueError
        If data or interpolation settings are invalid.
    """

    def __init__(
        self,
        depth_km: NDArray[np.float64],
        pressure_GPa: NDArray[np.float64],
        *,
        name: str = "tabulated-pressure",
        interpolation: InterpolationKind = "pchip",
        source: str | None = None,
        citation: str | None = None,
    ) -> None:
        super().__init__(
            depth_km,
            pressure_GPa,
            name=name,
            field_name="pressure",
            field_unit="GPa",
            interpolation=interpolation,
            source=source,
            citation=citation,
        )

    def pressure(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate interpolated pressure in GPa."""
        return self._evaluate(depth_km)


class TabulatedTemperatureModel(_TabulatedDepthModel):
    """Temperature-depth model interpolated from user-supplied values.

    Parameters
    ----------
    depth_km, temperature_K : ndarray
        Tabulated depths in km and absolute temperatures in K.
    name : str, optional
        Stable model identifier.
    interpolation : {"linear", "pchip"}, optional
        Interpolation scheme.  Extrapolation is never permitted.
    source, citation : str or None, optional
        File provenance and complete scientific citation.

    Raises
    ------
    ValueError
        If data or interpolation settings are invalid.
    """

    def __init__(
        self,
        depth_km: NDArray[np.float64],
        temperature_K: NDArray[np.float64],
        *,
        name: str = "tabulated-temperature",
        interpolation: InterpolationKind = "pchip",
        source: str | None = None,
        citation: str | None = None,
    ) -> None:
        super().__init__(
            depth_km,
            temperature_K,
            name=name,
            field_name="temperature",
            field_unit="K",
            interpolation=interpolation,
            source=source,
            citation=citation,
        )

    def temperature(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate interpolated temperature in K."""
        return self._evaluate(depth_km)


def read_tabulated_depth_field(
    filename: str | Path,
    *,
    field: Literal["pressure", "temperature"],
) -> tuple[FloatArray, FloatArray]:
    """Read one depth-dependent field from CSV or whitespace text.

    Parameters
    ----------
    filename : str or Path
        Input table with a header.  Depth aliases are ``depth_km``, ``depth``,
        and ``z_km``.  Pressure aliases are ``P_GPa``, ``pressure_GPa``,
        ``pressure``, and ``P``.  Temperature aliases are ``T_K``,
        ``temperature_K``, ``temperature``, and ``T``.
    field : {"pressure", "temperature"}
        Scalar field to read.

    Returns
    -------
    tuple of ndarray
        Sorted depth and field arrays.

    Raises
    ------
    ValueError
        If the table is empty, lacks columns, or contains invalid values.
    """
    path = Path(filename)
    text = path.read_text(encoding="utf-8")
    lines = [
        line
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines:
        raise ValueError(f"tabulated {field} file contains no data")
    if "," in lines[0]:
        rows = list(csv.DictReader(lines))
    else:
        headers = lines[0].split()
        rows = [dict(zip(headers, line.split(), strict=True)) for line in lines[1:]]
    if not rows:
        raise ValueError(f"tabulated {field} file contains no numerical rows")
    aliases = {
        "depth": ("depth_km", "depth", "z_km"),
        "pressure": ("P_GPa", "pressure_GPa", "pressure", "P"),
        "temperature": ("T_K", "temperature_K", "temperature", "T"),
    }

    def values_for(kind: str) -> FloatArray:
        key = next(
            (candidate for candidate in aliases[kind] if candidate in rows[0]),
            None,
        )
        if key is None:
            raise ValueError(f"missing required {kind} column")
        try:
            return np.asarray([float(row[key]) for row in rows], dtype=np.float64)
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid values in {kind} column") from exc

    depth = values_for("depth")
    data = values_for(field)
    order = np.argsort(depth, kind="stable")
    return depth[order], data[order]


__all__ = [
    "InterpolationKind",
    "TabulatedPressureModel",
    "TabulatedTemperatureModel",
    "read_tabulated_depth_field",
]
