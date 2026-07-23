# -*- coding: utf-8 -*-

"""Scientifically documented built-in terrestrial profile presets."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from quantas.references import render_citation_inline

from .composite import EarthProfileModel
from .models import EarthDepthProfile, TemperatureDepthModel
from .pressure import PremPressureModel
from .temperature import (
    ConductiveLayer,
    ContinentalConductiveGeotherm,
    Katsura2022MantleAdiabat,
    OceanicPlateGeotherm,
)


@dataclass(frozen=True, slots=True)
class EarthProfilePreset:
    """Passive specification of one built-in terrestrial profile.

    Parameters
    ----------
    name : str
        Stable API and CLI identifier.
    description : str
        Human-readable scientific scope.
    region : str
        Broad terrestrial region represented by the preset.
    depth_min_km, depth_max_km, depth_step_km : float
        Default sampling limits and spacing.
    pressure_model : str
        Pressure-model summary.
    temperature_model : str
        Temperature-model summary.
    references : tuple of str
        Bibliographic citations rendered from canonical Quantas citation keys.
    parameters : dict
        Explicit passive parameter mapping used to construct the preset.
    """

    name: str
    description: str
    region: str
    depth_min_km: float
    depth_max_km: float
    depth_step_km: float
    pressure_model: str
    temperature_model: str
    references: tuple[str, ...]
    parameters: dict[str, Any] = field(default_factory=dict)

    def build(self, *, name: str | None = None) -> EarthDepthProfile:
        """Construct and evaluate this preset.

        Parameters
        ----------
        name : str or None, optional
            Optional result profile name.

        Returns
        -------
        EarthDepthProfile
            Generated pressure-temperature-depth path.

        Raises
        ------
        ValueError
            If the preset parameters are internally invalid.
        """
        model = _build_model(self)
        return model.regular_profile(
            depth_min_km=self.depth_min_km,
            depth_max_km=self.depth_max_km,
            depth_step_km=self.depth_step_km,
            profile_name=self.name if name is None else name,
        )


_PREM_REFERENCE = render_citation_inline("prem_dziewonski_anderson_1981")
_HASTEROK_REFERENCE = render_citation_inline("hasterok_chapman_2011")
_PARSONS_REFERENCE = render_citation_inline("parsons_sclater_1977")
_KATSURA_REFERENCE = render_citation_inline("katsura_2022")
_KATSURA_SOFTWARE_REFERENCE = render_citation_inline("katsura_software_2022")


def _continental_parameters(
    *,
    surface_heat_flow_mW_m2: float,
    layers: list[dict[str, Any]],
) -> dict[str, Any]:
    return {
        "surface_temperature_K": 288.15,
        "surface_heat_flow_mW_m2": surface_heat_flow_mW_m2,
        "layers": layers,
    }


_PRESETS: dict[str, EarthProfilePreset] = {
    "continental-cratonic": EarthProfilePreset(
        name="continental-cratonic",
        description=(
            "Cool 200 km continental-lithosphere conductive reference with "
            "40 mW m^-2 surface heat flow."
        ),
        region="stable continental lithosphere",
        depth_min_km=0.0,
        depth_max_km=200.0,
        depth_step_km=2.0,
        pressure_model="PREM hydrostatic pressure",
        temperature_model="layered steady continental conduction",
        references=(_PREM_REFERENCE, _HASTEROK_REFERENCE),
        parameters=_continental_parameters(
            surface_heat_flow_mW_m2=40.0,
            layers=[
                {
                    "name": "upper crust",
                    "thickness_km": 15.0,
                    "conductivity_W_mK": 2.5,
                    "heat_production_uW_m3": 0.6,
                },
                {
                    "name": "lower crust",
                    "thickness_km": 25.0,
                    "conductivity_W_mK": 2.5,
                    "heat_production_uW_m3": 0.4,
                },
                {
                    "name": "lithospheric mantle",
                    "thickness_km": 160.0,
                    "conductivity_W_mK": 3.3,
                    "heat_production_uW_m3": 0.02,
                },
            ],
        ),
    ),
    "continental-reference": EarthProfilePreset(
        name="continental-reference",
        description=(
            "Intermediate 120 km continental-lithosphere conductive reference "
            "with 55 mW m^-2 surface heat flow."
        ),
        region="average continental lithosphere",
        depth_min_km=0.0,
        depth_max_km=120.0,
        depth_step_km=1.0,
        pressure_model="PREM hydrostatic pressure",
        temperature_model="layered steady continental conduction",
        references=(_PREM_REFERENCE, _HASTEROK_REFERENCE),
        parameters=_continental_parameters(
            surface_heat_flow_mW_m2=55.0,
            layers=[
                {
                    "name": "upper crust",
                    "thickness_km": 15.0,
                    "conductivity_W_mK": 2.5,
                    "heat_production_uW_m3": 0.8,
                },
                {
                    "name": "lower crust",
                    "thickness_km": 20.0,
                    "conductivity_W_mK": 2.5,
                    "heat_production_uW_m3": 0.4,
                },
                {
                    "name": "lithospheric mantle",
                    "thickness_km": 85.0,
                    "conductivity_W_mK": 3.3,
                    "heat_production_uW_m3": 0.02,
                },
            ],
        ),
    ),
    "continental-active": EarthProfilePreset(
        name="continental-active",
        description=(
            "Warm 80 km tectonically active continental conductive reference "
            "with 70 mW m^-2 surface heat flow."
        ),
        region="active continental lithosphere",
        depth_min_km=0.0,
        depth_max_km=80.0,
        depth_step_km=1.0,
        pressure_model="PREM hydrostatic pressure",
        temperature_model="layered steady continental conduction",
        references=(_PREM_REFERENCE, _HASTEROK_REFERENCE),
        parameters=_continental_parameters(
            surface_heat_flow_mW_m2=70.0,
            layers=[
                {
                    "name": "upper crust",
                    "thickness_km": 12.0,
                    "conductivity_W_mK": 2.5,
                    "heat_production_uW_m3": 1.2,
                },
                {
                    "name": "lower crust",
                    "thickness_km": 18.0,
                    "conductivity_W_mK": 2.5,
                    "heat_production_uW_m3": 0.4,
                },
                {
                    "name": "lithospheric mantle",
                    "thickness_km": 50.0,
                    "conductivity_W_mK": 3.3,
                    "heat_production_uW_m3": 0.02,
                },
            ],
        ),
    ),
    **{
        f"oceanic-{age}ma": EarthProfilePreset(
            name=f"oceanic-{age}ma",
            description=(f"Finite 125 km oceanic plate cooling profile at {age} Ma."),
            region="oceanic lithosphere",
            depth_min_km=0.0,
            depth_max_km=125.0,
            depth_step_km=1.0,
            pressure_model="PREM hydrostatic pressure",
            temperature_model="Parsons-Sclater finite plate cooling",
            references=(_PREM_REFERENCE, _PARSONS_REFERENCE),
            parameters={
                "age_Ma": float(age),
                "plate_thickness_km": 125.0,
                "surface_temperature_K": 273.15,
                "mantle_temperature_K": 1623.15,
                "diffusivity_m2_s": 1.0e-6,
                "series_terms": 200,
            },
        )
        for age in (10, 50, 100)
    },
    "mantle-katsura-2022": EarthProfilePreset(
        name="mantle-katsura-2022",
        description=(
            "Dry-pyrolite mantle adiabat from 50 to 2800 km, paired with PREM "
            "pressure and explicit phase-boundary sampling."
        ),
        region="convecting mantle",
        depth_min_km=50.0,
        depth_max_km=2800.0,
        depth_step_km=10.0,
        pressure_model="PREM hydrostatic pressure",
        temperature_model="Katsura (2022) published-constraint mantle adiabat",
        references=(
            _PREM_REFERENCE,
            _KATSURA_REFERENCE,
            _KATSURA_SOFTWARE_REFERENCE,
        ),
        parameters={"transition_epsilon_km": 1.0e-3},
    ),
}


def earth_profile_presets() -> tuple[EarthProfilePreset, ...]:
    """Return scientific terrestrial presets in stable display order.

    Returns
    -------
    tuple of EarthProfilePreset
        Immutable preset records.
    """
    return tuple(_PRESETS.values())


def build_earth_profile_preset(name: str) -> EarthDepthProfile:
    """Build one named scientific terrestrial profile.

    Parameters
    ----------
    name : str
        Identifier returned by :func:`earth_profile_presets`.

    Returns
    -------
    EarthDepthProfile
        Generated terrestrial profile.

    Raises
    ------
    ValueError
        If ``name`` is unknown.
    """
    normalized = name.strip().lower()
    try:
        preset = _PRESETS[normalized]
    except KeyError as exc:
        available = ", ".join(_PRESETS)
        raise ValueError(
            f"unknown Earth profile preset '{name}'; available: {available}"
        ) from exc
    return preset.build()


def _build_model(preset: EarthProfilePreset) -> EarthProfileModel:
    pressure = PremPressureModel(max_depth_km=max(2891.0, preset.depth_max_km))
    temperature: TemperatureDepthModel
    if preset.name.startswith("continental-"):
        layers = tuple(
            ConductiveLayer(
                thickness_km=float(item["thickness_km"]),
                conductivity_W_mK=float(item["conductivity_W_mK"]),
                heat_production_uW_m3=float(item["heat_production_uW_m3"]),
                name=str(item["name"]),
            )
            for item in preset.parameters["layers"]
        )
        temperature = ContinentalConductiveGeotherm(
            layers,
            surface_temperature_K=float(preset.parameters["surface_temperature_K"]),
            surface_heat_flow_mW_m2=float(preset.parameters["surface_heat_flow_mW_m2"]),
            name=preset.name,
        )
    elif preset.name.startswith("oceanic-"):
        temperature = OceanicPlateGeotherm(
            age_Ma=float(preset.parameters["age_Ma"]),
            plate_thickness_km=float(preset.parameters["plate_thickness_km"]),
            surface_temperature_K=float(preset.parameters["surface_temperature_K"]),
            mantle_temperature_K=float(preset.parameters["mantle_temperature_K"]),
            diffusivity_m2_s=float(preset.parameters["diffusivity_m2_s"]),
            series_terms=int(preset.parameters["series_terms"]),
        )
    elif preset.name == "mantle-katsura-2022":
        temperature = Katsura2022MantleAdiabat(
            transition_epsilon_km=float(preset.parameters["transition_epsilon_km"])
        )
    else:  # pragma: no cover - protected by the private preset table
        raise ValueError(f"unsupported Earth profile preset '{preset.name}'")
    return EarthProfileModel(pressure, temperature, name=preset.name)


__all__ = [
    "EarthProfilePreset",
    "build_earth_profile_preset",
    "earth_profile_presets",
]
