# -*- coding: utf-8 -*-

"""Shared labels and array helpers for neutral QHA reports."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from quantas.modules.qha.models import QHAResult

QHA_PROPERTY_LABELS: dict[str, tuple[str, str]] = {
    "equilibrium_volume": ("VT", "Equilibrium volume"),
    "zero_point_energy": ("Uzp", "Zero-point energy"),
    "thermal_energy": ("Uth", "Thermal energy"),
    "internal_energy": ("Utot", "Internal energy"),
    "entropy": ("S", "Entropy"),
    "vibrational_free_energy": ("Fvib", "Vibrational Helmholtz free energy"),
    "free_energy": ("F", "Helmholtz free energy"),
    "isochoric_heat_capacity": ("Cv", "Isochoric heat capacity"),
    "isobaric_heat_capacity": ("Cp", "Isobaric heat capacity"),
    "heat_capacity_difference": ("Cp-Cv", "Heat-capacity correction"),
    "isothermal_bulk_modulus": ("KT", "Isothermal bulk modulus"),
    "adiabatic_bulk_modulus": ("KS", "Adiabatic bulk modulus"),
    "bulk_modulus_derivative": ("Kp", "Bulk-modulus pressure derivative"),
    "thermal_expansion": ("alphaV", "Volumetric thermal expansion"),
    "thermal_expansion_mixed": (
        "alphaV_mixed",
        "Volumetric thermal expansion from the mixed derivative",
    ),
    "thermal_expansion_mode": (
        "alphaV_mode",
        "Volumetric thermal expansion from mode Gruneisen parameters",
    ),
    "thermal_expansion_numerical": (
        "alphaV_numerical",
        "Volumetric thermal expansion from the equilibrium-volume derivative",
    ),
    "enthalpy": ("H", "Enthalpy"),
    "gibbs_free_energy": ("G", "Gibbs free energy"),
    "gruneisen": ("gamma_th", "Thermodynamic Gruneisen parameter"),
    "mode_weighted_gruneisen": (
        "gamma_mode",
        "Heat-capacity-weighted mode Gruneisen parameter",
    ),
}

def resolve_property_name(name: str) -> str:
    """Resolve a result attribute or historical key to a QHAResult attribute.

    Parameters
    ----------
    name : str
        Attribute name or historical key.

    Returns
    -------
    str
        QHAResult attribute name.

    Raises
    ------
    KeyError
        If the property is unknown.
    """
    if name in QHA_PROPERTY_LABELS:
        return name
    for attr, (symbol, _) in QHA_PROPERTY_LABELS.items():
        if name == symbol:
            return attr
    raise KeyError(f"unknown QHA property: {name}")

def _property_rows(
    result: QHAResult,
    data: np.ndarray,
    sigma: np.ndarray | None,
    *,
    max_rows: int | None,
) -> list[list[Any]]:
    """Return rows for a QHA property array.

    Parameters
    ----------
    result : QHAResult
        Result containing pressure and temperature axes.
    data : ndarray
        Property data.
    sigma : ndarray or None
        Optional property uncertainty with the same shape as ``data``.
    max_rows : int or None
        Maximum number of rows.

    Returns
    -------
    list of list
        Table rows.

    Raises
    ------
    ValueError
        If the property shape cannot be mapped to pressure-temperature axes.
    """
    temperature = (
        np.asarray(result.temperature, dtype=np.float64)
        if result.temperature is not None
        else None
    )
    pressure = (
        np.asarray(result.pressure, dtype=np.float64)
        if result.pressure is not None
        else None
    )

    rows: list[list[Any]] = []
    if data.ndim == 0:
        rows.append([None, None, float(data), None if sigma is None else float(sigma)])
    elif data.ndim == 1:
        axis = (
            temperature
            if temperature is not None and temperature.shape[0] == data.shape[0]
            else None
        )
        if axis is None:
            axis = (
                pressure
                if pressure is not None and pressure.shape[0] == data.shape[0]
                else None
            )
        for index, value in enumerate(data):
            t_value = None if axis is None or axis is pressure else float(axis[index])
            p_value = (
                None if axis is None or axis is temperature else float(axis[index])
            )
            row = [t_value, p_value, float(value)]
            if sigma is not None:
                row.append(float(sigma[index]))
            rows.append(row)
    elif data.ndim == 2:
        if temperature is None or pressure is None:
            raise ValueError(
                "two-dimensional QHA properties require temperature and pressure axes"
            )
        if data.shape == (temperature.shape[0], pressure.shape[0]):
            for i, temp in enumerate(temperature):
                for j, press in enumerate(pressure):
                    row = [float(temp), float(press), float(data[i, j])]
                    if sigma is not None:
                        row.append(float(sigma[i, j]))
                    rows.append(row)
        elif data.shape == (pressure.shape[0], temperature.shape[0]):
            for i, press in enumerate(pressure):
                for j, temp in enumerate(temperature):
                    row = [float(temp), float(press), float(data[i, j])]
                    if sigma is not None:
                        row.append(float(sigma[i, j]))
                    rows.append(row)
        else:
            raise ValueError(
                "two-dimensional QHA property shape is incompatible with pressure-temperature axes"
            )
    else:
        raise ValueError(
            "QHA properties must be scalar, one-dimensional, or two-dimensional"
        )

    if max_rows is not None:
        return rows[:max_rows]
    return rows

def _uncertainty_for(
    result: QHAResult,
    attr: str,
    symbol: str,
) -> np.ndarray | None:
    """Return the uncertainty array associated with a property.

    Parameters
    ----------
    result : QHAResult
        Result containing uncertainty mappings.
    attr : str
        Result attribute name.
    symbol : str
        Historical property key.

    Returns
    -------
    ndarray or None
        Matching uncertainty array when available.
    """
    candidates = (
        f"sigma_{attr}",
        f"sigma_{symbol}",
        f"sigma_{symbol.replace('-', '_')}",
        attr,
        symbol,
    )
    for key in candidates:
        value = result.uncertainties.get(key)
        if value is not None:
            return np.asarray(value, dtype=np.float64)
    return None

def _normalized_indices(
    size: int,
    indices: Sequence[int] | None,
    label: str,
) -> list[int]:
    """Return unique validated indices in requested order.

    Parameters
    ----------
    size : int
        Number of entries in the source grid.
    indices : sequence of int or None
        Requested indices. ``None`` selects the complete grid.
    label : str
        Human-readable grid label used in error messages.

    Returns
    -------
    list of int
        Normalized indices.

    Raises
    ------
    IndexError
        If a requested index lies outside the source grid.
    """
    if indices is None:
        return list(range(size))
    normalized: list[int] = []
    for raw_index in indices:
        index = raw_index if raw_index >= 0 else size + raw_index
        if index < 0 or index >= size:
            raise IndexError(f"QHA {label} index out of range: {raw_index}")
        if index not in normalized:
            normalized.append(index)
    return normalized

__all__ = ["QHA_PROPERTY_LABELS", "resolve_property_name"]
