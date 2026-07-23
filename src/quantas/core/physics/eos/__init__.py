# -*- coding: utf-8 -*-

"""Equation-of-state models shared by Quantas scientific workflows."""

from .diagnostics import EOSStrainFamily, EOSStrainTransform, PressureEOSDiagnostics
from .energy import EnergyEOS, EnergyEOSFitModel
from .fitted import FittedEnergyEOS
from .parameters import (
    EOSParameters,
    free_energy_parameters,
    free_pressure_parameters,
    implied_kp,
    implied_kpp,
    resolve_energy_parameters,
    resolve_pressure_parameters,
    resolved_energy_parameter_covariance,
    resolved_energy_parameter_jacobian,
)
from .pressure import PressureEOS
from .pvt import (
    PVTCouplingFamily,
    PVTEOS,
    PVTModel,
    available_pvt_couplings,
    parse_pvt_coupling,
)
from .spec import (
    EOSFamily,
    EOSModel,
    available_eos_models,
    available_eos_tags,
    parse_eos_model,
)
from .state import EOSState, EOSStateUncertainty
from .thermal_pressure import (
    MGDNormalization,
    MGDThermalPressure,
    MGDParameters,
    MGDVariant,
    MGDVolumeBasis,
    ThermalPressureFamily,
    ThermalPressureModel,
    debye_function_3,
    parse_thermal_pressure_model,
    thermal_pressure_model_contracts,
    thermal_pressure_parameter_names,
)
from .temperature import (
    TemperatureEOS,
    TemperatureEOSFamily,
    TemperatureEOSModel,
    TemperatureEOSParameters,
    TemperatureEOSVariant,
    available_temperature_eos_models,
    parse_temperature_eos_model,
)

__all__ = [
    "EOSFamily",
    "EOSStrainFamily",
    "EOSStrainTransform",
    "EOSModel",
    "EOSParameters",
    "EOSState",
    "EOSStateUncertainty",
    "EnergyEOS",
    "EnergyEOSFitModel",
    "FittedEnergyEOS",
    "PressureEOS",
    "PressureEOSDiagnostics",
    "MGDNormalization",
    "MGDThermalPressure",
    "MGDParameters",
    "MGDVariant",
    "MGDVolumeBasis",
    "ThermalPressureFamily",
    "ThermalPressureModel",
    "debye_function_3",
    "PVTCouplingFamily",
    "PVTEOS",
    "PVTModel",
    "TemperatureEOS",
    "TemperatureEOSFamily",
    "TemperatureEOSModel",
    "TemperatureEOSParameters",
    "TemperatureEOSVariant",
    "available_eos_models",
    "available_eos_tags",
    "available_temperature_eos_models",
    "available_pvt_couplings",
    "thermal_pressure_model_contracts",
    "thermal_pressure_parameter_names",
    "free_energy_parameters",
    "free_pressure_parameters",
    "implied_kp",
    "implied_kpp",
    "parse_eos_model",
    "parse_temperature_eos_model",
    "parse_thermal_pressure_model",
    "parse_pvt_coupling",
    "resolve_energy_parameters",
    "resolve_pressure_parameters",
    "resolved_energy_parameter_covariance",
    "resolved_energy_parameter_jacobian",
]
