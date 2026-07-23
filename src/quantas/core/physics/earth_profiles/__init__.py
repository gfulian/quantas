# -*- coding: utf-8 -*-

"""Reusable terrestrial pressure-temperature profiles parameterized by depth."""

from .composite import (
    EarthProfileModel,
    JoinMode,
    PiecewiseTemperatureModel,
    TemperatureSegment,
)
from .models import (
    EarthDepthProfile,
    PressureDepthModel,
    TemperatureDepthModel,
)
from .presets import (
    EarthProfilePreset,
    build_earth_profile_preset,
    earth_profile_presets,
)
from .pressure import (
    LayeredLithostaticPressureModel,
    LithostaticLayer,
    PremPressureModel,
)
from .spec import (
    PROFILE_SPEC_SCHEMA_VERSION,
    build_earth_profile_from_mapping,
    read_earth_profile_spec,
)
from .tabulated import (
    InterpolationKind,
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

__all__ = [
    "ConductiveLayer",
    "ContinentalConductiveGeotherm",
    "EarthDepthProfile",
    "EarthProfileModel",
    "EarthProfilePreset",
    "InterpolationKind",
    "JoinMode",
    "Katsura2022MantleAdiabat",
    "LayeredLithostaticPressureModel",
    "LinearThermalBoundaryLayer",
    "LithostaticLayer",
    "OceanicHalfSpaceGeotherm",
    "OceanicPlateGeotherm",
    "PROFILE_SPEC_SCHEMA_VERSION",
    "PiecewiseTemperatureModel",
    "PremPressureModel",
    "PressureDepthModel",
    "TabulatedPressureModel",
    "TabulatedTemperatureModel",
    "TemperatureDepthModel",
    "TemperatureSegment",
    "build_earth_profile_from_mapping",
    "build_earth_profile_preset",
    "earth_profile_presets",
    "read_earth_profile_spec",
    "read_tabulated_depth_field",
]
