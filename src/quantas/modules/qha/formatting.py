# -*- coding: utf-8 -*-

"""Formatting helpers shared by QHA reports and text exporters.

The utilities in this module define the presentation conventions used for
pressure-temperature tables. Numerical results remain stored in their native
units; only labels and string formatting are handled here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(slots=True)
class QHATableFormat:
    """Numerical formats used by QHA reports and text tables.

    Parameters
    ----------
    temperature : str, optional
        Format specification for temperatures.
    pressure : str, optional
        Format specification for pressures.
    volume : str, optional
        Format specification for volumes.
    bulk_modulus : str, optional
        Format specification for isothermal and adiabatic bulk moduli.
    bulk_derivative : str, optional
        Format specification for pressure derivatives of the bulk modulus.
    thermal_expansion : str, optional
        Format specification for unscaled volumetric thermal expansion.
    thermal_expansion_scaled : str, optional
        Format specification for volumetric thermal expansion multiplied by
        ``1e5``.
    lattice_length : str, optional
        Format specification for crystallographic lengths and their
        uncertainties.
    lattice_angle : str, optional
        Format specification for crystallographic angles and their
        uncertainties.
    energy : str, optional
        Format specification for energies and thermodynamic quantities.
    default : str, optional
        Format specification used for other numerical values.
    signed_scientific : str, optional
        Scientific format used in fixed-width human-readable tables. A leading
        blank is reserved for positive values so that positive and negative
        numbers occupy the same width.
    energy_width : int, optional
        Fixed width reserved for energy, entropy and heat-capacity columns in
        human-readable tables.
    """

    temperature: str = ".2f"
    pressure: str = ".2f"
    volume: str = ".6f"
    bulk_modulus: str = ".5f"
    bulk_derivative: str = ".5f"
    thermal_expansion: str = ".6E"
    thermal_expansion_scaled: str = ".6f"
    lattice_length: str = ".8f"
    lattice_angle: str = ".8f"
    energy: str = ".12E"
    default: str = ".12E"
    signed_scientific: str = " .12E"
    temperature_width: int = 9
    pressure_width: int = 10
    volume_width: int = 12
    bulk_modulus_width: int = 12
    bulk_derivative_width: int = 11
    thermal_expansion_width: int = 14
    lattice_length_width: int = 13
    lattice_angle_width: int = 13
    energy_width: int = 20
    default_width: int = 20


def volume_unit_label(unit: str) -> str:
    """Return a volume unit from a length or volume label.

    Parameters
    ----------
    unit : str
        Length unit such as ``"A"`` or an already expanded volume unit such as
        ``"A^3"``.

    Returns
    -------
    str
        Canonical volume-unit label.
    """
    label = canonical_unit_label(unit)
    lowered = label.lower().replace(" ", "")
    if lowered in {"", "unknown", "-"}:
        return label or "unknown"
    if "^3" in lowered or "³" in lowered:
        return label
    return f"{label}^3"


def native_energy_unit_label(unit: str) -> str:
    """Return a native energy label with explicit cell normalization.

    Parameters
    ----------
    unit : str
        Energy-unit label.

    Returns
    -------
    str
        Canonical energy label.
    """
    label = canonical_unit_label(unit)
    lowered = label.lower()
    if lowered in {"", "unknown", "-"}:
        return label or "unknown"
    if "mol" in lowered or "cell" in lowered or "f.u." in lowered:
        return label
    return f"{label} cell^-1"


def canonical_unit_label(unit: str) -> str:
    """Normalize a unit label to the Quantas text-report convention.

    Parameters
    ----------
    unit : str
        Unit label to normalize.

    Returns
    -------
    str
        Unit label using multiplicative notation and negative powers.
    """
    label = str(unit).strip()
    replacements = {
        "kJ/mol/K": "kJ mol^-1 K^-1",
        "J/mol/K": "J mol^-1 K^-1",
        "kJ/mol": "kJ mol^-1",
        "J/mol": "J mol^-1",
        "Ha/K": "Ha K^-1",
        "eV/K": "eV K^-1",
        "Ry/K": "Ry K^-1",
    }
    if label in replacements:
        return replacements[label]
    label = label.replace(" per ", " ")
    label = label.replace("/cell", " cell^-1")
    label = label.replace("/mol", " mol^-1")
    label = label.replace("/K", " K^-1")
    label = " ".join(label.split())
    return label


def property_unit_label(attr: str, units: dict[str, Any]) -> str:
    """Return the display unit associated with a QHA property.

    Parameters
    ----------
    attr : str
        :class:`~quantas.modules.qha.models.QHAResult` attribute name.
    units : dict
        Unit metadata stored in the result.

    Returns
    -------
    str
        Canonical display-unit label. Dimensionless quantities are returned as
        ``"-"``.
    """
    energy = str(units.get("energy", "unknown"))
    pressure = canonical_unit_label(str(units.get("pressure", "unknown")))
    volume = str(units.get("volume", "unknown"))
    temperature = canonical_unit_label(str(units.get("temperature", "K")))

    if attr in {"volume", "equilibrium_volume"}:
        return volume_unit_label(volume)
    if attr == "pressure":
        return pressure
    if attr == "temperature":
        return temperature
    if attr == "entropy":
        return canonical_unit_label(
            str(units.get("entropy", f"{energy} cell^-1 {temperature}^-1"))
        )
    if attr in {
        "isochoric_heat_capacity",
        "isobaric_heat_capacity",
        "heat_capacity_difference",
    }:
        return canonical_unit_label(
            str(
                units.get(
                    "heat_capacity",
                    f"{energy} cell^-1 {temperature}^-1",
                )
            )
        )
    if attr in {"isothermal_bulk_modulus", "adiabatic_bulk_modulus"}:
        return pressure
    if attr in {
        "thermal_expansion",
        "thermal_expansion_mixed",
        "thermal_expansion_mode",
        "thermal_expansion_numerical",
    }:
        return f"{temperature}^-1"
    if attr in {
        "bulk_modulus_derivative",
        "gruneisen",
        "mode_weighted_gruneisen",
        "mode_gruneisen",
    }:
        return "-"
    if (
        attr == "static_energy"
        or attr.endswith("energy")
        or attr in {"enthalpy", "gibbs_free_energy", "free_energy"}
    ):
        return native_energy_unit_label(energy)
    return "-"


def format_spec_for_property(attr: str, table_format: QHATableFormat) -> str:
    """Return the numerical format associated with a QHA property.

    Parameters
    ----------
    attr : str
        QHA result attribute name.
    table_format : QHATableFormat
        Formatting configuration.

    Returns
    -------
    str
        Python format specification.
    """
    if attr == "equilibrium_volume":
        return table_format.volume
    if attr in {"isothermal_bulk_modulus", "adiabatic_bulk_modulus"}:
        return table_format.bulk_modulus
    if attr in {
        "bulk_modulus_derivative",
        "gruneisen",
        "mode_weighted_gruneisen",
    }:
        return table_format.bulk_derivative
    if attr in {
        "thermal_expansion",
        "thermal_expansion_mixed",
        "thermal_expansion_mode",
        "thermal_expansion_numerical",
    }:
        return table_format.thermal_expansion
    if attr in {
        "zero_point_energy",
        "thermal_energy",
        "internal_energy",
        "entropy",
        "vibrational_free_energy",
        "free_energy",
        "isochoric_heat_capacity",
        "isobaric_heat_capacity",
        "heat_capacity_difference",
        "enthalpy",
        "gibbs_free_energy",
    }:
        return table_format.energy
    return table_format.default


def format_number(value: float, spec: str) -> str:
    """Format a finite floating-point value.

    Parameters
    ----------
    value : float
        Value to format.
    spec : str
        Python format specification.

    Returns
    -------
    str
        Formatted value, or ``"nan"`` for non-finite input.
    """
    import numpy as np

    if not np.isfinite(value):
        return "nan"
    return format(value, spec)
