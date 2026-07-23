# -*- coding: utf-8 -*-

"""Piecewise-linear interpolation on two-dimensional rectilinear fields.

The implementation supports non-uniform and singleton coordinate axes,
trailing field dimensions, paired point evaluation, and endpoint-slope
extrapolation.  It is independent from QHA and thermoelasticity so every
frontend and workflow uses the same numerical path.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray: TypeAlias = NDArray[np.float64]
BoolArray: TypeAlias = NDArray[np.bool_]


@dataclass(frozen=True, slots=True)
class RectilinearFieldInterpolator:
    """Interpolate fields sampled on increasing ``x`` and ``y`` axes.

    Parameters
    ----------
    x, y : array_like
        Strictly increasing source axes.
    values : array_like
        Source data with shape ``(x.size, y.size, ...)``.  Any trailing field
        dimensions are preserved.

    Notes
    -----
    Extrapolation is linear using the slope of the nearest endpoint interval.
    A singleton axis is treated as constant.  Returned arrays always use
    ``float64``.
    """

    x: FloatArray
    y: FloatArray
    values: FloatArray

    def __init__(self, x: ArrayLike, y: ArrayLike, values: ArrayLike) -> None:
        x_axis = validated_axis(x, "x")
        y_axis = validated_axis(y, "y")
        field = np.asarray(values, dtype=np.float64)
        if field.shape[:2] != (x_axis.size, y_axis.size):
            raise ValueError("values must start with shape (x.size, y.size)")
        if np.any(~np.isfinite(field)):
            raise ValueError("rectilinear source values must be finite")
        object.__setattr__(self, "x", x_axis)
        object.__setattr__(self, "y", y_axis)
        object.__setattr__(self, "values", field.copy())

    def evaluate_grid(
        self,
        target_x: ArrayLike,
        target_y: ArrayLike,
    ) -> tuple[FloatArray, BoolArray]:
        """Evaluate a Cartesian product of target coordinates.

        Parameters
        ----------
        target_x, target_y : array_like
            Strictly increasing target axes.

        Returns
        -------
        tuple of ndarray
            Interpolated field with shape ``(nx, ny, ...)`` and a Boolean
            extrapolation mask with shape ``(nx, ny)``.
        """
        tx = validated_axis(target_x, "target x")
        ty = validated_axis(target_y, "target y")
        along_y = _linear_axis(self.y, self.values, ty, axis=1)
        result = _linear_axis(self.x, along_y, tx, axis=0)
        outside = np.asarray(
            (tx[:, None] < self.x[0])
            | (tx[:, None] > self.x[-1])
            | (ty[None, :] < self.y[0])
            | (ty[None, :] > self.y[-1]),
            dtype=np.bool_,
        )
        return np.asarray(result, dtype=np.float64), outside

    def evaluate_points(
        self,
        target_x: ArrayLike,
        target_y: ArrayLike,
    ) -> tuple[FloatArray, BoolArray]:
        """Evaluate aligned coordinate pairs.

        Parameters
        ----------
        target_x, target_y : array_like
            Aligned finite one-dimensional vectors.

        Returns
        -------
        tuple of ndarray
            Interpolated values with shape ``(npoints, ...)`` and a Boolean
            extrapolation mask with shape ``(npoints,)``.
        """
        tx = np.asarray(target_x, dtype=np.float64)
        ty = np.asarray(target_y, dtype=np.float64)
        if tx.ndim != 1 or ty.shape != tx.shape:
            raise ValueError("target coordinate pairs must be aligned vectors")
        if np.any(~np.isfinite(tx)) or np.any(~np.isfinite(ty)):
            raise ValueError("target coordinate pairs must be finite")
        trailing = self.values.shape[2:]
        result = np.empty((tx.size,) + trailing, dtype=np.float64)
        for index, (x_value, y_value) in enumerate(zip(tx, ty, strict=True)):
            along_y = _linear_axis(
                self.y,
                self.values,
                np.asarray([y_value], dtype=np.float64),
                axis=1,
            )[:, 0, ...]
            value = _linear_axis(
                self.x,
                along_y,
                np.asarray([x_value], dtype=np.float64),
                axis=0,
            )
            result[index] = value[0]
        outside = np.asarray(
            (tx < self.x[0]) | (tx > self.x[-1]) | (ty < self.y[0]) | (ty > self.y[-1]),
            dtype=np.bool_,
        )
        return result, outside


def validated_axis(values: ArrayLike, name: str) -> FloatArray:
    """Return a copied finite, strictly increasing coordinate axis."""
    axis = np.asarray(values, dtype=np.float64)
    if axis.ndim != 1 or axis.size < 1 or np.any(~np.isfinite(axis)):
        raise ValueError(f"{name} must be a non-empty finite vector")
    if np.any(np.diff(axis) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    return axis.copy()


def regular_grid(start: float, stop: float, step: float) -> FloatArray:
    """Return an inclusive validated ``float64`` grid."""
    if not all(np.isfinite(value) for value in (start, stop, step)):
        raise ValueError("grid bounds and step must be finite")
    if step <= 0.0:
        raise ValueError("grid step must be positive")
    if stop < start:
        raise ValueError("grid maximum must be greater than or equal to minimum")
    count = int(np.floor((stop - start) / step + 1.0e-12)) + 1
    values = start + step * np.arange(count, dtype=np.float64)
    if values[-1] < stop - max(1.0, abs(stop)) * 1.0e-12:
        values = np.append(values, float(stop))
    else:
        values[-1] = float(stop)
    return np.asarray(values, dtype=np.float64)


def grid_step(values: ArrayLike) -> float | None:
    """Return the uniform spacing of an axis, when one exists."""
    axis = np.asarray(values, dtype=np.float64)
    if axis.ndim != 1 or axis.size < 2:
        return None
    differences = np.diff(axis)
    if np.allclose(differences, differences[0], rtol=1.0e-10, atol=1.0e-12):
        return float(differences[0])
    return None


def _linear_axis(
    x: FloatArray,
    values: FloatArray,
    target: FloatArray,
    *,
    axis: int,
) -> FloatArray:
    """Interpolate one axis while retaining all remaining dimensions."""
    moved = np.moveaxis(np.asarray(values, dtype=np.float64), axis, 0)
    if moved.shape[0] != x.size:
        raise ValueError("interpolation axis length does not match coordinates")
    if x.size == 1:
        expanded = np.broadcast_to(moved[0], (target.size,) + moved.shape[1:])
        return np.moveaxis(np.asarray(expanded, dtype=np.float64), 0, axis)
    indices = np.searchsorted(x, target, side="right")
    indices = np.clip(indices, 1, x.size - 1)
    left = indices - 1
    right = indices
    fraction = (target - x[left]) / (x[right] - x[left])
    reshape = (target.size,) + (1,) * (moved.ndim - 1)
    result = moved[left] + fraction.reshape(reshape) * (moved[right] - moved[left])
    return np.moveaxis(np.asarray(result, dtype=np.float64), 0, axis)


__all__ = [
    "RectilinearFieldInterpolator",
    "grid_step",
    "regular_grid",
    "validated_axis",
]
