# -*- coding: utf-8 -*-

"""Labels and property metadata for QHA plots.

The definitions in this module provide a single source for plot symbols, unit
selection, scaling factors, and filename-safe property names used by API and CLI
plotting frontends.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from quantas.modules.qha.report import QHA_PROPERTY_LABELS, resolve_property_name


@dataclass(frozen=True, slots=True)
class QHAPlotProperty:
    """Description of a plottable QHA property.

    Parameters
    ----------
    attribute : str
        Name of the :class:`~quantas.modules.qha.models.QHAResult` attribute.
    key : str
        Historical short key used in filenames and CLI options.
    symbol : str
        Matplotlib mathtext symbol used in axis labels and titles.
    symbol_plain : str
        Unicode or plain-text symbol independent of a concrete renderer.
    description : str
        Human-readable description of the property.
    unit_key : str
        Metadata key used to retrieve the physical unit.
    scale : float, optional
        Multiplicative factor applied only for plotting.
    scale_label : str, optional
        Text appended to labels when a scaled quantity is shown.
    """

    attribute: str
    key: str
    symbol: str
    symbol_plain: str
    description: str
    unit_key: str
    scale: float = 1.0
    scale_label: str = ""


_PROPERTY_OVERRIDES: dict[str, tuple[str, str, str, float, str]] = {
    "equilibrium_volume": (r"$V$", "V", "volume", 1.0, ""),
    "zero_point_energy": (r"$U_{zp}$", "U_ZP", "energy", 1.0, ""),
    "thermal_energy": (r"$U_{th}$", "U_th", "energy", 1.0, ""),
    "internal_energy": (r"$U_{tot}$", "U", "energy", 1.0, ""),
    "entropy": (r"$S$", "S", "entropy", 1.0, ""),
    "vibrational_free_energy": (r"$F_{vib}$", "F_vib", "energy", 1.0, ""),
    "free_energy": (r"$F$", "F", "energy", 1.0, ""),
    "isochoric_heat_capacity": (r"$C_V$", "Cᵥ", "heat_capacity", 1.0, ""),
    "isobaric_heat_capacity": (r"$C_P$", "Cₚ", "heat_capacity", 1.0, ""),
    "heat_capacity_difference": (
        r"$C_P-C_V$",
        "Cₚ−Cᵥ",
        "heat_capacity",
        1.0,
        "",
    ),
    "isothermal_bulk_modulus": (r"$K_T$", "Kₜ", "pressure", 1.0, ""),
    "adiabatic_bulk_modulus": (r"$K_S$", "Kₛ", "pressure", 1.0, ""),
    "bulk_modulus_derivative": (
        r"$K^{\prime}$",
        "K′",
        "dimensionless",
        1.0,
        "",
    ),
    "thermal_expansion": (
        r"$\alpha_V$",
        "αᵥ",
        "temperature_inverse",
        1.0e5,
        r"$\times 10^5$",
    ),
    "enthalpy": (r"$H$", "H", "energy", 1.0, ""),
    "gibbs_free_energy": (r"$G$", "G", "energy", 1.0, ""),
    "gruneisen": (r"$\gamma_{th}$", "γ_th", "dimensionless", 1.0, ""),
    "mode_weighted_gruneisen": (
        r"$\gamma_{mode}$",
        "γ_mode",
        "dimensionless",
        1.0,
        "",
    ),
}


def available_plot_properties() -> dict[str, QHAPlotProperty]:
    """Return metadata for all standard QHA plot properties.

    Returns
    -------
    dict
        Mapping from QHA result attribute name to plotting metadata.
    """
    properties: dict[str, QHAPlotProperty] = {}
    for attribute, (key, description) in QHA_PROPERTY_LABELS.items():
        if attribute not in _PROPERTY_OVERRIDES:
            continue
        symbol, symbol_plain, unit_key, scale, scale_label = _PROPERTY_OVERRIDES[
            attribute
        ]
        properties[attribute] = QHAPlotProperty(
            attribute=attribute,
            key=key,
            symbol=symbol,
            symbol_plain=symbol_plain,
            description=description,
            unit_key=unit_key,
            scale=scale,
            scale_label=scale_label,
        )
    return properties


def resolve_plot_property(name: str) -> QHAPlotProperty:
    """Resolve a user property name into plotting metadata.

    Parameters
    ----------
    name : str
        Historical key, attribute name, or supported alias.

    Returns
    -------
    QHAPlotProperty
        Resolved plotting metadata.

    Raises
    ------
    KeyError
        If the requested property cannot be plotted.
    """
    normalized = name.strip()
    aliases = {
        "v": "VT",
        "volume": "VT",
        "alpha": "alphaV",
        "alpha_v": "alphaV",
        "thermal_expansion": "alphaV",
        "cv": "Cv",
        "cp": "Cp",
        "heat_capacity": "Cp",
        "heat_capacities": "heat_capacities",
        "cp_cv": "heat_capacities",
    }
    lowered = normalized.lower()
    if lowered in {"heat_capacities", "cp_cv"}:
        raise KeyError("heat_capacities is a combined plot, not a single property")
    attr = resolve_property_name(aliases.get(lowered, normalized))
    properties = available_plot_properties()
    if attr not in properties:
        raise KeyError(f"QHA property '{name}' cannot be plotted")
    return properties[attr]


def unit_for_property(
    property_info: QHAPlotProperty, units: Mapping[str, object]
) -> str:
    """Return the display unit for a plot property.

    Parameters
    ----------
    property_info : QHAPlotProperty
        Plot property metadata.
    units : mapping
        Unit metadata stored in a QHA result or plotting display context.

    Returns
    -------
    str
        Unit label suitable for plot axes.
    """
    key = property_info.unit_key
    if key == "volume":
        volume_unit = str(units.get("volume", "A"))

        if volume_unit in {"A", "angstrom", "Angstrom", "Å"}:
            return r"$\AA^3$"
        if volume_unit.endswith("^3") or volume_unit.endswith("³"):
            return volume_unit
        return f"{volume_unit}^3"
    if key == "pressure":
        return str(units.get("pressure", "GPa"))
    if key == "energy":
        return str(units.get("energy", "energy unit"))
    if key == "entropy":
        return str(units.get("entropy", units.get("heat_capacity", "energy/K")))
    if key == "heat_capacity":
        return str(
            units.get("heat_capacity", str(units.get("energy", "energy unit")) + "/K")
        )
    if key == "temperature_inverse":
        temperature_unit = str(units.get("temperature", "K"))
        return f"10$^{{-5}}$ {temperature_unit}$^{{-1}}$"
    if key == "dimensionless":
        return "dimensionless"
    return str(units.get(key, "unknown"))


def ylabel_for_property(
    property_info: QHAPlotProperty, units: Mapping[str, object]
) -> str:
    """Build a complete y-axis label for a QHA property.

    Parameters
    ----------
    property_info : QHAPlotProperty
        Plot property metadata.
    units : mapping
        Unit metadata stored in a QHA result or plotting display context.

    Returns
    -------
    str
        Axis label containing symbol and unit.
    """
    unit = unit_for_property(property_info, units)
    if property_info.scale_label and property_info.unit_key != "temperature_inverse":
        return f"{property_info.symbol} {property_info.scale_label} ({unit})"
    return f"{property_info.symbol} ({unit})"
