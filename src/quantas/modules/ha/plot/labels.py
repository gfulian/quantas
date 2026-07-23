# -*- coding: utf-8 -*-

"""
Axis-label utilities for harmonic-approximation plots.

This module centralizes scientific symbols and unit formatting used by neutral
plot specifications. It does not perform plotting or unit conversion. The
returned labels use portable mathematical text understood by the supported
static renderer and reusable by future frontend adapters.
"""

from __future__ import annotations

import re


THERMODYNAMIC_SYMBOLS = {
    "static_energy": r"$U_0$",
    "zero_point_energy": r"$U_\mathrm{ZP}$",
    "thermal_energy": r"$U_\mathrm{th}$",
    "internal_energy": r"$U$",
    "entropy": r"$S$",
    "vibrational_free_energy": r"$F_\mathrm{vib}$",
    "free_energy": r"$F$",
    "isochoric_heat_capacity": r"$C_V$",
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
