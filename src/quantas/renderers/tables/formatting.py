# -*- coding: utf-8 -*-

"""Central numeric formatting profiles for tabular renderers.

Formatting is intentionally independent from in-memory and HDF5 precision.
Renderers may display a compact or publication-oriented representation without
altering the underlying numerical values.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np


NUMERIC_FORMAT_PROFILES: dict[str, str] = {
    "general": ".8g",
    "integer": "d",
    "energy": ".12E",
    "energy_ha": ".12E",
    "energy_ev": ".10E",
    "volume": ".8f",
    "pressure": ".8f",
    "modulus": ".8f",
    "thermoelastic_modulus": ".4f",
    "thermoelastic_uncertainty": ".4f",
    "temperature": ".2f",
    "frequency": ".6f",
    "dimensionless": ".8f",
    "uncertainty": ".6E",
    "angle": ".6f",
    "eos_pressure": ".4f",
    "eos_temperature": ".2f",
    "eos_structural": ".6f",
    "eos_correlation": ".6f",
    "eos_covariance": ".6e",
}


NumericFormatter = Callable[[Any], str]


def _adaptive_six(value: Any) -> str:
    """Format EOS values with six decimals and scientific fallback.

    Ordinary values use fixed-point notation for visually aligned tables.
    Very small or very large values use scientific notation so significant
    information is not hidden by a fixed decimal width.
    """
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    number = float(value)
    magnitude = abs(number)
    if number != 0.0 and (magnitude < 1.0e-4 or magnitude >= 1.0e6):
        return f"{number:.6e}"
    return f"{number:.6f}"


def _pressure_uncertainty(value: Any) -> str:
    """Format pressure uncertainties without rounding finite values to zero."""
    number = float(value)
    if number != 0.0 and abs(number) < 5.0e-5:
        return f"{number:.6e}"
    return f"{number:.4f}"


def _temperature_uncertainty(value: Any) -> str:
    """Format temperature uncertainties with a two-decimal default."""
    number = float(value)
    if number != 0.0 and abs(number) < 5.0e-3:
        return f"{number:.6e}"
    return f"{number:.2f}"


ADAPTIVE_NUMERIC_FORMATTERS: dict[str, NumericFormatter] = {
    "eos_parameter": _adaptive_six,
    "eos_uncertainty": _adaptive_six,
    "eos_residual": _adaptive_six,
    "eos_statistic": _adaptive_six,
    "eos_pressure_uncertainty": _pressure_uncertainty,
    "eos_temperature_uncertainty": _temperature_uncertainty,
}


def resolve_numeric_format(format_name: str | None) -> str:
    """Resolve a profile name or explicit Python format specification.

    Parameters
    ----------
    format_name : str or None
        Registered profile name, explicit format specification such as
        ``".12e"``, or ``None`` for the general profile.

    Returns
    -------
    str
        Python format specification without surrounding braces.

    Raises
    ------
    ValueError
        If a named profile is unknown or an explicit specification is invalid.
    """
    if format_name is None:
        return NUMERIC_FORMAT_PROFILES["general"]
    if format_name in NUMERIC_FORMAT_PROFILES:
        return NUMERIC_FORMAT_PROFILES[format_name]
    if format_name in ADAPTIVE_NUMERIC_FORMATTERS:
        return format_name
    try:
        format(1.2345, format_name)
    except (ValueError, TypeError) as exc:
        available = ", ".join(sorted(NUMERIC_FORMAT_PROFILES))
        raise ValueError(
            f"Unknown numeric format profile or invalid format specification "
            f"'{format_name}'. Available profiles: {available}."
        ) from exc
    return format_name


def format_numeric(value: Any, format_name: str | None = None) -> str:
    """Format one real numerical value without modifying the source value."""
    specification = resolve_numeric_format(format_name)
    if specification in ADAPTIVE_NUMERIC_FORMATTERS:
        return ADAPTIVE_NUMERIC_FORMATTERS[specification](value)
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value))
    if isinstance(value, (int, np.integer)):
        if specification == "d":
            return format(int(value), specification)
        return format(float(value), specification)
    if isinstance(value, (float, np.floating)):
        return format(float(value), specification)
    raise TypeError(f"value of type {type(value).__name__!r} is not a real number")


def format_cell(value: Any, format_name: str | None = None) -> str:
    """Convert one neutral table cell to display text.

    Parameters
    ----------
    value : Any
        Cell value to render.
    format_name : str or None, optional
        Numeric profile name or explicit Python format specification.

    Returns
    -------
    str
        Display representation that leaves the source value unchanged.
    """
    if value is None:
        return ""
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value))
    if isinstance(value, (int, float, np.integer, np.floating)):
        return format_numeric(value, format_name)
    return str(value)


__all__ = [
    "ADAPTIVE_NUMERIC_FORMATTERS",
    "NUMERIC_FORMAT_PROFILES",
    "format_cell",
    "format_numeric",
    "resolve_numeric_format",
]
