# -*- coding: utf-8 -*-

"""Domain-specific EOS fitting adapters.

The package groups the fitting layer by scientific domain while preserving the
historical import modules ``quantas.modules.eos.fitting``,
``quantas.modules.eos.temperature_fitting``, and
``quantas.modules.eos.pvt_fitting`` as compatibility facades.
"""

from .pv import (
    AxialEOSFitModel,
    PressureEOSFitModel,
    axial_to_volume_parameters,
    build_axial_parameter_map,
    build_pressure_parameter_map,
    estimate_axial_parameters,
    estimate_pressure_parameters,
)
from .pvt import (
    PVTEOSFitModel,
    build_pvt_parameter_map,
    estimate_pvt_parameters,
    pvt_parameter_names,
)
from .vt import (
    TemperatureEOSFitModel,
    build_temperature_parameter_map,
    estimate_temperature_parameters,
    temperature_parameter_names,
)

__all__ = [
    "AxialEOSFitModel",
    "PressureEOSFitModel",
    "PVTEOSFitModel",
    "TemperatureEOSFitModel",
    "axial_to_volume_parameters",
    "build_axial_parameter_map",
    "build_pressure_parameter_map",
    "build_pvt_parameter_map",
    "build_temperature_parameter_map",
    "estimate_axial_parameters",
    "estimate_pressure_parameters",
    "estimate_pvt_parameters",
    "estimate_temperature_parameters",
    "pvt_parameter_names",
    "temperature_parameter_names",
]
