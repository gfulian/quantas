# -*- coding: utf-8 -*-

"""Shared property schema and formatting for QHA table export."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from quantas.modules.qha.formatting import (
    QHATableFormat,
    format_number,
    format_spec_for_property,
    property_unit_label,
)

EXPORT_PROPERTY_ORDER: tuple[str, ...] = (
    "equilibrium_volume",
    "isothermal_bulk_modulus",
    "adiabatic_bulk_modulus",
    "bulk_modulus_derivative",
    "thermal_expansion",
    "zero_point_energy",
    "thermal_energy",
    "internal_energy",
    "vibrational_free_energy",
    "free_energy",
    "enthalpy",
    "gibbs_free_energy",
    "entropy",
    "isochoric_heat_capacity",
    "isobaric_heat_capacity",
    "heat_capacity_difference",
    "gruneisen",
    "mode_weighted_gruneisen",
)

THERMOELASTIC_PROPERTIES = {
    "equilibrium_volume",
    "isothermal_bulk_modulus",
    "adiabatic_bulk_modulus",
    "bulk_modulus_derivative",
    "thermal_expansion",
}

GRUNEISEN_PROPERTIES = {
    "gruneisen",
    "mode_weighted_gruneisen",
}

ENERGY_PROPERTIES = {
    "zero_point_energy",
    "thermal_energy",
    "internal_energy",
    "vibrational_free_energy",
    "free_energy",
    "enthalpy",
    "gibbs_free_energy",
}

ENERGY_CONTRIBUTION_PROPERTIES = {
    "zero_point_energy",
    "thermal_energy",
    "internal_energy",
    "vibrational_free_energy",
}

THERMODYNAMIC_POTENTIAL_PROPERTIES = {
    "free_energy",
    "enthalpy",
    "gibbs_free_energy",
}

THERMAL_PROPERTIES = {
    "entropy",
    "isochoric_heat_capacity",
    "isobaric_heat_capacity",
    "heat_capacity_difference",
}

STRUCTURAL_PROPERTY_ALIASES = {
    "structure",
    "structural",
    "cell",
    "cell_parameters",
    "lattice",
    "lattice_parameters",
    "alphaabc",
    "axial_thermal_expansion",
}

@dataclass(slots=True)
class _StructuralExportData:
    """Validated structural arrays prepared for QHA table export."""

    lattice_parameters: np.ndarray
    lattice_parameter_uncertainties: np.ndarray | None
    axial_thermal_expansion: np.ndarray | None
    axial_thermal_expansion_uncertainties: np.ndarray | None
    thermal_expansion_tensor: np.ndarray | None
    volumetric_thermal_expansion: np.ndarray | None
    tensor_trace: np.ndarray | None
    trace_residual: np.ndarray | None
    extrapolation_mask: np.ndarray | None

def _format_for_property(attr: str, table_format: QHATableFormat) -> str:
    """Return the numerical format associated with a property attribute."""
    return format_spec_for_property(attr, table_format)

def _format_number(value: float, spec: str) -> str:
    """Format a floating-point value using a table format specification."""
    return format_number(value, spec)

def _unit_for_property(attr: str, units: dict[str, Any]) -> str | None:
    """Return a unit label for a QHA symbol property.

    Parameters
    ----------
    attr : str
        QHAResult attribute name.
    units : dict
        Unit metadata stored in the result.

    Returns
    -------
    str or None
        Unit label suitable for reports and HDF5 attributes.
    """
    return property_unit_label(attr, units)

__all__ = [
    "ENERGY_CONTRIBUTION_PROPERTIES",
    "ENERGY_PROPERTIES",
    "EXPORT_PROPERTY_ORDER",
    "GRUNEISEN_PROPERTIES",
    "STRUCTURAL_PROPERTY_ALIASES",
    "THERMAL_PROPERTIES",
    "THERMODYNAMIC_POTENTIAL_PROPERTIES",
    "THERMOELASTIC_PROPERTIES",
]
