# -*- coding: utf-8 -*-

"""YAML specifications for composed terrestrial depth profiles."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, cast

import yaml

from .composite import (
    EarthProfileModel,
    PiecewiseTemperatureModel,
    TemperatureSegment,
)
from .models import EarthDepthProfile, PressureDepthModel, TemperatureDepthModel
from .pressure import (
    LayeredLithostaticPressureModel,
    LithostaticLayer,
    PremPressureModel,
)
from .tabulated import (
    TabulatedPressureModel,
    TabulatedTemperatureModel,
    read_tabulated_depth_field,
)
from .temperature import (
    ConductiveLayer,
    ContinentalConductiveGeotherm,
    Katsura2022MantleAdiabat,
    LinearThermalBoundaryLayer,
    OceanicHalfSpaceGeotherm,
    OceanicPlateGeotherm,
)


PROFILE_SPEC_SCHEMA_VERSION = 1


def read_earth_profile_spec(filename: str | Path) -> EarthDepthProfile:
    """Read and evaluate a composed terrestrial profile specification.

    Parameters
    ----------
    filename : str or Path
        YAML file with top-level ``name``, ``depth``, ``pressure``, and
        ``temperature`` mappings.  Relative tabulated paths are resolved from
        the YAML file directory.

    Returns
    -------
    EarthDepthProfile
        Evaluated pressure-temperature-depth path with full model provenance.

    Raises
    ------
    ValueError
        If the YAML structure, schema version, model settings, or paths are
        invalid.
    """
    path = Path(filename)
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, Mapping):
        raise ValueError("Earth profile specification must be a YAML mapping")
    profile = build_earth_profile_from_mapping(raw, base_path=path.parent)
    profile.metadata["specification_source"] = str(path)
    return profile


def build_earth_profile_from_mapping(
    specification: Mapping[str, Any],
    *,
    base_path: str | Path = ".",
) -> EarthDepthProfile:
    """Build a terrestrial profile from a passive mapping.

    Parameters
    ----------
    specification : mapping
        Parsed profile specification.
    base_path : str or Path, optional
        Directory used to resolve relative tabulated files.

    Returns
    -------
    EarthDepthProfile
        Evaluated composed profile.

    Raises
    ------
    ValueError
        If required fields or model options are invalid.
    """
    version = int(specification.get("schema_version", PROFILE_SPEC_SCHEMA_VERSION))
    if version != PROFILE_SPEC_SCHEMA_VERSION:
        raise ValueError(
            f"unsupported Earth profile schema_version {version}; "
            f"expected {PROFILE_SPEC_SCHEMA_VERSION}"
        )
    name = str(specification.get("name", "custom-earth-profile")).strip()
    if not name:
        raise ValueError("Earth profile name cannot be empty")
    depth_spec = _mapping(specification.get("depth"), "depth")
    pressure_spec = _mapping(specification.get("pressure"), "pressure")
    temperature_spec = _mapping(specification.get("temperature"), "temperature")
    base = Path(base_path)
    pressure = _build_pressure_model(pressure_spec, base_path=base)
    temperature = _build_temperature_model(temperature_spec, base_path=base)
    model = EarthProfileModel(pressure, temperature, name=name)
    lower = float(depth_spec.get("min_km", model.depth_bounds[0]))
    upper = float(depth_spec.get("max_km", model.depth_bounds[1]))
    step = float(depth_spec.get("step_km", 1.0))
    include_critical = bool(depth_spec.get("include_critical_depths", True))
    profile = model.regular_profile(
        depth_min_km=lower,
        depth_max_km=upper,
        depth_step_km=step,
        include_critical_depths=include_critical,
        profile_name=name,
    )
    profile.metadata["profile_schema_version"] = version
    return profile


def _build_pressure_model(
    specification: Mapping[str, Any],
    *,
    base_path: Path,
) -> PressureDepthModel:
    model = _model_name(specification)
    if model == "prem":
        return PremPressureModel(
            max_depth_km=float(specification.get("max_depth_km", 2891.0)),
            integration_step_km=float(specification.get("integration_step_km", 0.25)),
        )
    if model in {"layered-lithostatic", "lithostatic"}:
        items = _sequence(specification.get("layers"), "pressure.layers")
        layers = tuple(
            LithostaticLayer(
                thickness_km=float(_mapping(item, "pressure layer")["thickness_km"]),
                density_kg_m3=float(_mapping(item, "pressure layer")["density_kg_m3"]),
                name=str(_mapping(item, "pressure layer").get("name", "layer")),
            )
            for item in items
        )
        return LayeredLithostaticPressureModel(
            layers,
            gravity_m_s2=float(specification.get("gravity_m_s2", 9.80665)),
            name=str(specification.get("name", "layered-lithostatic")),
            citation=_optional_string(specification.get("citation")),
        )
    if model in {"table", "tabulated"}:
        path = _resolve_path(specification, base_path)
        depth, values = read_tabulated_depth_field(path, field="pressure")
        return TabulatedPressureModel(
            depth,
            values,
            name=str(specification.get("name", path.stem)),
            interpolation=cast(Any, specification.get("interpolation", "pchip")),
            source=str(path),
            citation=_optional_string(specification.get("citation")),
        )
    raise ValueError(f"unknown pressure model '{model}'")


def _build_temperature_model(
    specification: Mapping[str, Any],
    *,
    base_path: Path,
) -> TemperatureDepthModel:
    model = _model_name(specification)
    if model in {"continental-conductive", "conductive"}:
        items = _sequence(specification.get("layers"), "temperature.layers")
        layers = tuple(
            ConductiveLayer(
                thickness_km=float(_mapping(item, "temperature layer")["thickness_km"]),
                conductivity_W_mK=float(
                    _mapping(item, "temperature layer")["conductivity_W_mK"]
                ),
                heat_production_uW_m3=float(
                    _mapping(item, "temperature layer")["heat_production_uW_m3"]
                ),
                name=str(_mapping(item, "temperature layer").get("name", "layer")),
            )
            for item in items
        )
        return ContinentalConductiveGeotherm(
            layers,
            surface_temperature_K=float(
                specification.get("surface_temperature_K", 288.15)
            ),
            surface_heat_flow_mW_m2=float(
                specification.get("surface_heat_flow_mW_m2", 40.0)
            ),
            name=str(specification.get("name", "continental-conductive")),
            citation=_optional_string(specification.get("citation")),
        )
    if model == "oceanic-half-space":
        return OceanicHalfSpaceGeotherm(
            age_Ma=float(specification["age_Ma"]),
            mantle_temperature_K=float(
                specification.get("mantle_temperature_K", 1623.15)
            ),
            surface_temperature_K=float(
                specification.get("surface_temperature_K", 273.15)
            ),
            diffusivity_m2_s=float(specification.get("diffusivity_m2_s", 1.0e-6)),
            max_depth_km=float(specification.get("max_depth_km", 300.0)),
        )
    if model == "oceanic-plate":
        return OceanicPlateGeotherm(
            age_Ma=float(specification["age_Ma"]),
            plate_thickness_km=float(specification.get("plate_thickness_km", 125.0)),
            mantle_temperature_K=float(
                specification.get("mantle_temperature_K", 1623.15)
            ),
            surface_temperature_K=float(
                specification.get("surface_temperature_K", 273.15)
            ),
            diffusivity_m2_s=float(specification.get("diffusivity_m2_s", 1.0e-6)),
            series_terms=int(specification.get("series_terms", 200)),
            max_depth_km=(
                None
                if specification.get("max_depth_km") is None
                else float(specification["max_depth_km"])
            ),
        )
    if model in {"katsura-2022", "mantle-katsura-2022"}:
        return Katsura2022MantleAdiabat(
            transition_epsilon_km=float(
                specification.get("transition_epsilon_km", 1.0e-3)
            )
        )
    if model in {"linear-boundary-layer", "thermal-boundary-layer"}:
        return LinearThermalBoundaryLayer(
            depth_top_km=float(specification["depth_top_km"]),
            depth_bottom_km=float(specification["depth_bottom_km"]),
            temperature_top_K=float(specification["temperature_top_K"]),
            temperature_bottom_K=float(specification["temperature_bottom_K"]),
            exponent=float(specification.get("exponent", 1.0)),
            name=str(specification.get("name", "linear-thermal-boundary-layer")),
            citation=_optional_string(specification.get("citation")),
        )
    if model in {"table", "tabulated"}:
        path = _resolve_path(specification, base_path)
        depth, values = read_tabulated_depth_field(path, field="temperature")
        return TabulatedTemperatureModel(
            depth,
            values,
            name=str(specification.get("name", path.stem)),
            interpolation=cast(Any, specification.get("interpolation", "pchip")),
            source=str(path),
            citation=_optional_string(specification.get("citation")),
        )
    if model == "piecewise":
        items = _sequence(specification.get("segments"), "temperature.segments")
        segments: list[TemperatureSegment] = []
        for item in items:
            segment_spec = dict(_mapping(item, "temperature segment"))
            depth_min = float(segment_spec.pop("depth_min_km"))
            depth_max = float(segment_spec.pop("depth_max_km"))
            join = segment_spec.pop("join", None)
            join_mapping = {} if join is None else _mapping(join, "segment.join")
            join_mode = str(
                join_mapping.get(
                    "mode",
                    segment_spec.pop("join_mode", "direct"),
                )
            )
            blend_width = float(
                join_mapping.get(
                    "width_km",
                    segment_spec.pop("blend_width_km", 0.0),
                )
            )
            nested = _build_temperature_model(segment_spec, base_path=base_path)
            segments.append(
                TemperatureSegment(
                    depth_min_km=depth_min,
                    depth_max_km=depth_max,
                    model=nested,
                    join_mode=cast(Any, join_mode),
                    blend_width_km=blend_width,
                )
            )
        return PiecewiseTemperatureModel(
            segments,
            name=str(specification.get("name", "piecewise-temperature")),
        )
    raise ValueError(f"unknown temperature model '{model}'")


def _model_name(specification: Mapping[str, Any]) -> str:
    source = specification.get("source")
    model = specification.get("model", source)
    if model is None:
        raise ValueError("profile component requires a 'model' or 'source' field")
    return str(model).strip().lower()


def _resolve_path(specification: Mapping[str, Any], base_path: Path) -> Path:
    value = specification.get("file")
    if value is None:
        raise ValueError("tabulated profile component requires a 'file' path")
    path = Path(str(value))
    if not path.is_absolute():
        path = base_path / path
    if not path.is_file():
        raise ValueError(f"tabulated profile file does not exist: {path}")
    return path


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return cast(Mapping[str, Any], value)


def _sequence(value: Any, name: str) -> Sequence[Any]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a sequence")
    if len(value) == 0:
        raise ValueError(f"{name} cannot be empty")
    return cast(Sequence[Any], value)


def _optional_string(value: Any) -> str | None:
    return None if value is None else str(value)


__all__ = [
    "PROFILE_SPEC_SCHEMA_VERSION",
    "build_earth_profile_from_mapping",
    "read_earth_profile_spec",
]
