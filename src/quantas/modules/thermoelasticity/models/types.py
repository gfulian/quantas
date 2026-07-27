# -*- coding: utf-8 -*-

"""Shared type aliases and QHA protocol for thermoelastic contracts."""

from __future__ import annotations
from typing import Literal, Protocol, TypeAlias
import numpy as np
from numpy.typing import NDArray

FloatArray: TypeAlias = NDArray[np.float64]
ThermoelasticMethod = Literal["quasistatic"]
ThermoelasticQualityLevel = Literal[
    "supported", "caution", "unsupported", "not_applicable"
]
ThermoelasticQualityPolicy = Literal["fail", "warn", "allow"]
ThermoelasticStabilityPolicy = Literal["fail", "warn", "allow"]
ThermoelasticAdiabaticMode = Literal["auto", "off", "require"]
ThermoelasticTensorCondition = Literal["isothermal", "adiabatic"]
ThermoelasticReportLevel = Literal["standard", "extended", "debug"]
ThermoelasticExtrapolationPolicy = Literal["fail", "warn", "allow"]
ThermoelasticFitFailurePolicy = Literal["stop", "continue", "raise"]


class QHAThermoelasticPayload(Protocol):
    """Structural protocol for QHA fields consumed by thermoelasticity."""

    temperature: FloatArray | None
    pressure: FloatArray | None
    volume: FloatArray | None
    static_energy: FloatArray | None
    equilibrium_volume: FloatArray | None
    isochoric_heat_capacity: FloatArray | None
    isothermal_bulk_modulus: FloatArray | None
    bulk_modulus_derivative: FloatArray | None
    thermal_expansion: FloatArray | None
    axial_thermal_expansion: FloatArray | None
    thermal_expansion_tensor: FloatArray | None
    equilibrium_lattice: FloatArray | None
    lattice_parameters: FloatArray | None
    uncertainties: dict[str, FloatArray]


# Internal compatibility alias retained while the public protocol name is
# adopted throughout the refactored codebase.
_QHAThermoelasticPayload = QHAThermoelasticPayload


__all__ = [
    "FloatArray",
    "ThermoelasticAdiabaticMode",
    "ThermoelasticExtrapolationPolicy",
    "ThermoelasticFitFailurePolicy",
    "ThermoelasticMethod",
    "ThermoelasticQualityLevel",
    "ThermoelasticQualityPolicy",
    "ThermoelasticReportLevel",
    "ThermoelasticStabilityPolicy",
    "ThermoelasticTensorCondition",
    "QHAThermoelasticPayload",
    "_QHAThermoelasticPayload",
]
