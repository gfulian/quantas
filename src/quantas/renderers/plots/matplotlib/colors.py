# -*- coding: utf-8 -*-

"""Colormap validation and resolution for the Matplotlib backend."""

from __future__ import annotations

from difflib import get_close_matches
from typing import Any


CUSTOM_COLORMAPS: dict[str, tuple[str, ...]] = {
    "quantas_powerflow": ("black", "cyan"),
    "quantas_enhancement": ("blue", "white", "red", "black"),
}


def available_colormap_names() -> tuple[str, ...]:
    """Return sorted Matplotlib and Quantas colormap names.

    Returns
    -------
    tuple of str
        Names accepted by :func:`resolve_colormap`.
    """
    from matplotlib import colormaps

    return tuple(sorted(set(colormaps).union(CUSTOM_COLORMAPS)))


def validate_colormap_name(name: str) -> str:
    """Validate a Matplotlib-compatible colormap name before calculation.

    Parameters
    ----------
    name : str
        Requested colormap name.

    Returns
    -------
    str
        The validated name.

    Raises
    ------
    ValueError
        If the name is unavailable. The error includes close matches when
        possible so CLI adapters can fail cleanly before expensive workflows.
    """
    available = available_colormap_names()
    if name in available:
        return name
    suggestions = get_close_matches(name, available, n=3, cutoff=0.55)
    message = f"Unknown Matplotlib colormap '{name}'."
    if suggestions:
        message += (
            " Did you mean " + ", ".join(repr(item) for item in suggestions) + "?"
        )
    raise ValueError(message)


def resolve_colormap(plt: Any, name: str) -> Any:
    """Return a Matplotlib colormap, including Quantas custom families.

    Parameters
    ----------
    plt : module
        Imported :mod:`matplotlib.pyplot` module.
    name : str
        Matplotlib colormap name or Quantas custom family.

    Returns
    -------
    object
        Matplotlib colormap instance.

    Raises
    ------
    ValueError
        If the requested colormap does not exist.
    """
    from matplotlib.colors import LinearSegmentedColormap

    validate_colormap_name(name)
    if name in CUSTOM_COLORMAPS:
        return LinearSegmentedColormap.from_list(
            name,
            CUSTOM_COLORMAPS[name],
            N=256,
        )
    return plt.get_cmap(name)


__all__ = [
    "CUSTOM_COLORMAPS",
    "available_colormap_names",
    "resolve_colormap",
    "validate_colormap_name",
]
