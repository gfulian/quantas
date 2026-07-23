# -*- coding: utf-8 -*-

"""Shared helpers for the Matplotlib plotting backend."""

from __future__ import annotations

from typing import Any


LINE_STYLES = {
    "solid": "-",
    "dashed": "--",
    "dotted": ":",
    "dashdot": "-.",
    "none": "None",
}

SURFACE_ELEVATION = 29.0
SURFACE_AZIMUTH = 33.0


def apply_x_limits(
    axis: Any,
    limits: tuple[float | None, float | None] | None,
) -> None:
    """Apply optional horizontal axis limits.

    Parameters
    ----------
    axis : matplotlib.axes.Axes
        Axis receiving the limits.
    limits : tuple or None
        Optional lower and upper horizontal bounds.

    Returns
    -------
    None
        The axis is modified in place.
    """
    if limits is not None:
        axis.set_xlim(left=limits[0], right=limits[1])


def apply_y_limits(
    axis: Any,
    limits: tuple[float | None, float | None] | None,
) -> None:
    """Apply optional vertical axis limits.

    Parameters
    ----------
    axis : matplotlib.axes.Axes
        Axis receiving the limits.
    limits : tuple or None
        Optional lower and upper vertical bounds.

    Returns
    -------
    None
        The axis is modified in place.
    """
    if limits is not None:
        axis.set_ylim(bottom=limits[0], top=limits[1])


def pyplot() -> Any:
    """Import Matplotlib pyplot lazily.

    Returns
    -------
    module
        :mod:`matplotlib.pyplot` module.
    """
    import matplotlib.pyplot as plt

    return plt


__all__ = [
    "LINE_STYLES",
    "SURFACE_AZIMUTH",
    "SURFACE_ELEVATION",
    "apply_x_limits",
    "apply_y_limits",
    "pyplot",
]
