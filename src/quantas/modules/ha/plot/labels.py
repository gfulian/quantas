# -*- coding: utf-8 -*-

"""Scientific property metadata and label helpers for HA plots.

This module is the single authoritative catalogue for harmonic plot property
keys, names, renderer-neutral symbols, established MathText labels, and
scientific categories.  It also formats the legacy Matplotlib-compatible axis
labels used by the current static renderer.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from quantas.modules.ha.report import THERMODYNAMIC_LABELS


@dataclass(frozen=True, slots=True)
class HAPlotProperty:
    """Describe one standard harmonic plot property.

    Parameters
    ----------
    attribute : str
        Name of the :class:`~quantas.modules.ha.models.HAResult` attribute.
    key : str
        Stable short key accepted by the public HA builder.
    symbol : str
        Established MathText symbol used by the static renderer.
    symbol_math : str
        Renderer-neutral mathematical symbol without delimiters.
    symbol_plain : str
        Unicode or plain-text symbol.
    name : str
        Extended human-readable property name.
    description : str
        Short scientific description.
    category : str
        Stable scientific grouping key.
    """

    attribute: str
    key: str
    symbol: str
    symbol_math: str
    symbol_plain: str
    name: str
    description: str
    category: str


_PROPERTY_DETAILS: dict[str, tuple[str, str, str, str, str]] = {
    "static_energy": (
        r"$U_0$",
        "U_0",
        "U₀",
        "Electronic static energy at each sampled volume.",
        "energy",
    ),
    "zero_point_energy": (
        r"$U_\mathrm{ZP}$",
        r"U_{\mathrm{ZP}}",
        "U_ZP",
        "Temperature-independent vibrational zero-point energy.",
        "energy",
    ),
    "thermal_energy": (
        r"$U_\mathrm{th}$",
        r"U_{\mathrm{th}}",
        "U_th",
        "Thermal vibrational contribution to the internal energy.",
        "energy",
    ),
    "internal_energy": (
        r"$U$",
        "U",
        "U",
        "Total static, zero-point, and thermal internal energy.",
        "energy",
    ),
    "entropy": (
        r"$S$",
        "S",
        "S",
        "Harmonic vibrational entropy.",
        "entropy",
    ),
    "vibrational_free_energy": (
        r"$F_\mathrm{vib}$",
        r"F_{\mathrm{vib}}",
        "F_vib",
        "Vibrational contribution to the Helmholtz free energy.",
        "energy",
    ),
    "free_energy": (
        r"$F$",
        "F",
        "F",
        "Total harmonic Helmholtz free energy.",
        "energy",
    ),
    "isochoric_heat_capacity": (
        r"$C_V$",
        "C_V",
        "Cᵥ",
        "Heat capacity at constant volume.",
        "heat_capacity",
    ),
}


def available_plot_properties() -> dict[str, HAPlotProperty]:
    """Return the authoritative standard HA plot-property catalogue."""
    labels = {"static_energy": ("U0", "Static energy"), **THERMODYNAMIC_LABELS}
    properties: dict[str, HAPlotProperty] = {}
    for attribute, (key, name) in labels.items():
        symbol, symbol_math, symbol_plain, description, category = (
            _PROPERTY_DETAILS[attribute]
        )
        properties[attribute] = HAPlotProperty(
            attribute=attribute,
            key=key,
            symbol=symbol,
            symbol_math=symbol_math,
            symbol_plain=symbol_plain,
            name=name,
            description=description,
            category=category,
        )
    return properties


def resolve_plot_property(name: str) -> HAPlotProperty:
    """Resolve one HA result attribute or stable short plotting key."""
    properties = available_plot_properties()
    if name in properties:
        return properties[name]
    for property_info in properties.values():
        if name == property_info.key:
            return property_info
    raise KeyError(f"unknown HA plot property: {name}")


THERMODYNAMIC_SYMBOLS = {
    attribute: property_info.symbol
    for attribute, property_info in available_plot_properties().items()
}


TEMPERATURE_UNITS = {
    "K": "K",
    "C": "°C",
    "c": "°C",
    "Celsius": "°C",
    "celsius": "°C",
    "F": "°F",
    "f": "°F",
    "Fahrenheit": "°F",
    "fahrenheit": "°F",
}


def format_temperature_unit(unit: str | None) -> str:
    """
    Return a human-readable temperature unit label.

    Parameters
    ----------
    unit : str or None
        Stored or requested temperature unit.

    Returns
    -------
    str
        Temperature unit label suitable for axis labels.
    """
    if unit is None:
        return "K"
    return TEMPERATURE_UNITS.get(unit, unit)


def format_temperature_axis_label(unit: str | None) -> str:
    """
    Return the x-axis label for HA temperature plots.

    Parameters
    ----------
    unit : str or None
        Temperature unit.

    Returns
    -------
    str
        Axis label in the form ``"Temperature (unit)"``.
    """
    return f"Temperature ({format_temperature_unit(unit)})"


def format_property_axis_label(attr: str, unit: str | None) -> str:
    """
    Return the y-axis label for a HA thermodynamic property.

    Parameters
    ----------
    attr : str
        HAResult attribute name.
    unit : str or None
        Unit label associated with the plotted values.

    Returns
    -------
    str
        Axis label combining the property symbol and the unit in parentheses.

    Raises
    ------
    KeyError
        If ``attr`` is not a known HA thermodynamic property.
    """
    symbol = THERMODYNAMIC_SYMBOLS[attr]
    if unit is None or str(unit).strip() == "":
        return symbol
    return f"{symbol} ({format_unit_mathtext(str(unit))})"


def format_volume_legend_label(unit: str | None) -> str:
    """
    Return the volume label for HA plots.

    Parameters
    ----------
    unit : str or None
        Volume unit.

    Returns
    -------
    str
        Legend label in the form ``"Volume (unit)"``.
    """
    return (
        "Volume"
        if unit is None or not str(unit).strip()
        else f"Volume ({format_unit_mathtext(str(unit))})"
    )


def format_unit_mathtext(unit: str) -> str:
    r"""
    Return a unit label with Matplotlib MathText exponents.

    Parameters
    ----------
    unit : str
        Unit label, such as ``"Ha"``, ``"kJ mol^-1"`` or
        ``"kJ mol^-1 K^-1"``.

    Returns
    -------
    str
        Unit label with negative powers formatted as MathText, for example
        ``"kJ mol$^{-1}$ K$^{-1}$"``.
    """
    label = unit.replace("−", "-").replace("**", "^")
    label = label.replace("/mol", " mol^-1")
    label = label.replace("/K", " K^-1")

    def repl(match: re.Match[str]) -> str:
        return f"$^{{{match.group(1)}}}$"

    label = re.sub(r"\^\{?(-?\d+)\}?", repl, label)
    return label
