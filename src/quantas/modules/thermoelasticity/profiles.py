# -*- coding: utf-8 -*-

"""Adapters from terrestrial Earth models to thermoelastic depth profiles.

Scientific pressure and temperature models live in
:mod:`quantas.core.physics.earth_profiles`.  This module only adapts their
frontend-neutral output to the passive thermoelastic workflow contract.  Only quantitative Earth reference presets are registered publicly.  Controlled
linear paths remain available through explicit CLI/API parameters and are never
presented as geological reference models.
"""

from __future__ import annotations

from dataclasses import dataclass

from quantas.core.physics.earth_profiles import (
    build_earth_profile_preset,
    earth_profile_presets,
)
from quantas.modules.thermoelasticity.models import ThermoelasticDepthProfile


@dataclass(frozen=True, slots=True)
class ThermoelasticProfilePreset:
    """Passive specification of one built-in thermoelastic depth profile.

    Parameters
    ----------
    name : str
        Stable CLI and API identifier.
    description : str
        Human-readable scope and intended use.
    region : str
        Geological region represented by the profile.
    depth_min, depth_max, depth_step : float
        Default depth interval and spacing in km.
    pressure_model : str
        Pressure-model summary.
    temperature_model : str
        Temperature-model summary.
    references : tuple of str
        Complete bibliographic citations for quantitative models.
    core_preset : str
        Identifier delegated to the shared Earth-profile core.
    """

    name: str
    description: str
    region: str
    depth_min: float
    depth_max: float
    depth_step: float
    pressure_model: str
    temperature_model: str
    references: tuple[str, ...] = ()
    core_preset: str = ""

    def build(self, *, name: str | None = None) -> ThermoelasticDepthProfile:
        """Create a validated thermoelastic depth profile.

        Parameters
        ----------
        name : str or None, optional
            Optional result-profile name.  The preset identifier is used by
            default.

        Returns
        -------
        ThermoelasticDepthProfile
            Generated pressure-temperature-depth path.

        """
        result_name = self.name if name is None else name
        profile = build_earth_profile_preset(self.core_preset)
        metadata = dict(profile.metadata)
        metadata.update(
            {
                "thermoelastic_profile_preset": self.name,
                "description": self.description,
                "region": self.region,
                "references": self.references,
            }
        )
        return ThermoelasticDepthProfile(
            name=result_name,
            depth=profile.depth,
            pressure=profile.pressure,
            temperature=profile.temperature,
            metadata=metadata,
        )


_SCIENTIFIC_PRESETS = tuple(
    ThermoelasticProfilePreset(
        name=preset.name,
        description=preset.description,
        region=preset.region,
        depth_min=preset.depth_min_km,
        depth_max=preset.depth_max_km,
        depth_step=preset.depth_step_km,
        pressure_model=preset.pressure_model,
        temperature_model=preset.temperature_model,
        references=preset.references,
        core_preset=preset.name,
    )
    for preset in earth_profile_presets()
)

_PROFILE_PRESETS: dict[str, ThermoelasticProfilePreset] = {
    preset.name: preset for preset in _SCIENTIFIC_PRESETS
}


def thermoelastic_profile_presets() -> tuple[ThermoelasticProfilePreset, ...]:
    """Return all built-in profile specifications in stable display order.

    Returns
    -------
    tuple of ThermoelasticProfilePreset
        Quantitative scientific Earth-profile presets.
    """
    return tuple(_PROFILE_PRESETS.values())


def build_thermoelastic_profile_preset(name: str) -> ThermoelasticDepthProfile:
    """Build one named built-in geological profile.

    Parameters
    ----------
    name : str
        Preset identifier returned by :func:`thermoelastic_profile_presets`.

    Returns
    -------
    ThermoelasticDepthProfile
        Generated profile.

    Raises
    ------
    ValueError
        If ``name`` is not a known preset.
    """
    normalized = name.strip().lower()
    try:
        preset = _PROFILE_PRESETS[normalized]
    except KeyError as exc:
        available = ", ".join(_PROFILE_PRESETS)
        raise ValueError(
            f"unknown thermoelastic profile preset '{name}'; available: {available}"
        ) from exc
    return preset.build()


__all__ = [
    "ThermoelasticProfilePreset",
    "build_thermoelastic_profile_preset",
    "thermoelastic_profile_presets",
]
