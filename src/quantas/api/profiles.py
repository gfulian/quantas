# -*- coding: utf-8 -*-

"""Public API for reusable terrestrial pressure-temperature profiles."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from quantas.core.physics.earth_profiles import (
    EarthDepthProfile as DepthProfile,
    EarthProfileModel as Model,
    EarthProfilePreset as Preset,
    build_earth_profile_from_mapping as _from_mapping,
    build_earth_profile_preset as _build_preset,
    earth_profile_presets as _presets,
    read_earth_profile_spec as _read_spec,
)


from .common import _public_dir

def from_mapping(specification: Mapping[str, Any]) -> DepthProfile:
    """Build a depth profile from a normalized mapping.

    Parameters
    ----------
    specification : mapping of str to Any
        Normalized pressure-temperature-depth model specification.

    Returns
    -------
    DepthProfile
        Evaluated terrestrial depth path.

    Raises
    ------
    ValueError
        If the specification is incomplete or scientifically invalid.
    """
    return _from_mapping(specification)


def build_preset(name: str) -> DepthProfile:
    """Build one named depth-profile preset.

    Parameters
    ----------
    name : str
        Stable identifier returned by :func:`preset_names`.

    Returns
    -------
    DepthProfile
        Evaluated pressure-temperature-depth path.
    """
    return _build_preset(name)


def presets() -> tuple[Preset, ...]:
    """Return available depth-profile preset contracts.

    Returns
    -------
    tuple of Preset
        Immutable preset metadata in stable display order.
    """
    return _presets()


def preset_names() -> tuple[str, ...]:
    """Return available depth-profile preset names.

    Returns
    -------
    tuple of str
        Stable identifiers accepted by :func:`build_preset`.
    """
    return tuple(item.name for item in presets())


def read_spec(source: str | Path) -> DepthProfile:
    """Read a depth-profile specification.

    Parameters
    ----------
    source : str or Path
        YAML profile specification.

    Returns
    -------
    DepthProfile
        Evaluated pressure-temperature-depth path.

    Raises
    ------
    ValueError
        If the specification cannot be parsed or validated.
    """
    return _read_spec(source)


def __dir__() -> list[str]:
    """Return only names in the supported public contract."""
    return _public_dir(__all__)


__all__ = [
    "DepthProfile",
    "Model",
    "Preset",
    "build_preset",
    "from_mapping",
    "preset_names",
    "presets",
    "read_spec",
]
