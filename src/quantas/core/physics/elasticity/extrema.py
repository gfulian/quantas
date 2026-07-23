# -*- coding: utf-8 -*-

"""Global directional extrema for elastic properties."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy import optimize

from .directional import direction_vector


@dataclass(slots=True)
class DirectionalExtrema:
    """Minimum, maximum, anisotropy, and corresponding directions."""

    minimum: float
    maximum: float
    anisotropy: float
    minimum_axis: list[float]
    maximum_axis: list[float]
    minimum_measurement_axis: list[float] | None = None
    maximum_measurement_axis: list[float] | None = None


def find_directional_extrema(
    function: Callable[[np.ndarray], float],
    dimensions: int,
) -> DirectionalExtrema:
    """Calculate global extrema of a directional property.

    Parameters
    ----------
    function : callable
        Objective accepting two or three angular coordinates.
    dimensions : int
        Number of angular coordinates, either ``2`` or ``3``.

    Returns
    -------
    DirectionalExtrema
        Extrema, anisotropy, and Cartesian directions.

    Raises
    ------
    ValueError
        If ``dimensions`` is not ``2`` or ``3``.
    """
    minimum_angles, minimum = _minimize(function, dimensions)
    maximum_angles, maximum = _maximize(function, dimensions)
    anisotropy = maximum / minimum if minimum > 0.0 else np.inf

    return DirectionalExtrema(
        minimum=minimum,
        maximum=maximum,
        anisotropy=float(anisotropy),
        minimum_axis=direction_vector(minimum_angles[0], minimum_angles[1]).tolist(),
        maximum_axis=direction_vector(maximum_angles[0], maximum_angles[1]).tolist(),
        minimum_measurement_axis=(
            direction_vector(
                minimum_angles[0], minimum_angles[1], minimum_angles[2]
            ).tolist()
            if dimensions == 3
            else None
        ),
        maximum_measurement_axis=(
            direction_vector(
                maximum_angles[0], maximum_angles[1], maximum_angles[2]
            ).tolist()
            if dimensions == 3
            else None
        ),
    )


def _minimize(
    function: Callable[[np.ndarray], float],
    dimensions: int,
) -> tuple[np.ndarray, float]:
    """Minimize a directional property using the historical search grid."""
    ranges: tuple[tuple[float, float], ...]
    if dimensions == 2:
        ranges = ((0.0, np.pi), (0.0, np.pi))
        points = 25
    elif dimensions == 3:
        ranges = ((0.0, np.pi), (0.0, np.pi), (0.0, np.pi))
        points = 10
    else:
        raise ValueError("dimensions must be 2 or 3.")

    angles, value = optimize.brute(
        function,
        ranges,
        Ns=points,
        full_output=True,
        finish=optimize.fmin,
    )[0:2]
    return np.asarray(angles, dtype=float), float(value)


def _maximize(
    function: Callable[[np.ndarray], float],
    dimensions: int,
) -> tuple[np.ndarray, float]:
    """Maximize a directional property using minimization of its negative."""
    angles, value = _minimize(lambda values: -function(values), dimensions)
    return angles, -value
