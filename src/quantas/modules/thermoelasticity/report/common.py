# -*- coding: utf-8 -*-

"""Shared report-formatting helpers."""

from __future__ import annotations
from typing import Any
import numpy as np
from quantas.modules.thermoelasticity.models import ThermoelasticResult


def _outside_distance(value: float, lower: Any, upper: Any) -> float | None:
    try:
        lower_value = float(lower)
        upper_value = float(upper)
    except (TypeError, ValueError):
        return None
    if value < lower_value:
        return value - lower_value
    if value > upper_value:
        return value - upper_value
    return 0.0


def _metadata_float(result: ThermoelasticResult, key: str) -> float | None:
    try:
        value = float(result.metadata[key])
    except (KeyError, TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _metadata_range_text(
    result: ThermoelasticResult,
    lower_key: str,
    upper_key: str,
) -> str:
    lower = _metadata_float(result, lower_key)
    upper = _metadata_float(result, upper_key)
    if lower is None or upper is None:
        return "unavailable"
    return f"{lower:.8g} – {upper:.8g}"


def _range_text(values: Any) -> str:
    array = np.asarray(values, dtype=np.float64)
    if array.size == 0 or not np.any(np.isfinite(array)):
        return "unavailable"
    return f"{float(np.nanmin(array)):.8g} – {float(np.nanmax(array)):.8g}"


def _maximum_relative_uncertainty(values: Any, sigma: Any) -> float | None:
    if values is None or sigma is None:
        return None
    value_array = np.asarray(values, dtype=np.float64)
    sigma_array = np.asarray(sigma, dtype=np.float64)
    denominator = np.abs(value_array)
    valid = (
        np.isfinite(value_array) & np.isfinite(sigma_array) & (denominator > 1.0e-12)
    )
    if not np.any(valid):
        return None
    return 100.0 * float(np.nanmax(sigma_array[valid] / denominator[valid]))


def _minimum_stability_eigenvalue(stability: Any) -> float | None:
    """Return the minimum finite stiffness eigenvalue from a field."""
    if stability is None:
        return None
    values = np.asarray(stability.minimum_eigenvalue, dtype=np.float64)
    finite = values[np.isfinite(values)]
    return None if finite.size == 0 else float(np.min(finite))


__all__ = [
    "_maximum_relative_uncertainty",
    "_metadata_float",
    "_metadata_range_text",
    "_minimum_stability_eigenvalue",
    "_outside_distance",
    "_range_text",
]
