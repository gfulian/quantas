# -*- coding: utf-8 -*-

"""Density calculations from chemical formulas and cell volumes."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np
from scipy import constants as cs

from quantas.core.physics.units import convert_volume

from .formula import molar_mass

_DENSITY_TO_KG_M3 = {
    "kg m^-3": 1.0,
    "kg/m^3": 1.0,
    "kg m-3": 1.0,
    "kgm^-3": 1.0,
    "g cm^-3": 1000.0,
    "g/cm^3": 1000.0,
    "g cm-3": 1000.0,
    "gcm^-3": 1000.0,
}


def _normalize_density_unit(unit: str) -> str:
    if not isinstance(unit, str) or not unit.strip():
        raise ValueError("density unit must be a non-empty string")
    return " ".join(unit.strip().lower().replace("³", "^3").split())


def convert_density(
    value: Any,
    from_unit: str,
    to_unit: str,
) -> float | np.ndarray:
    """Convert mass densities between supported units.

    Parameters
    ----------
    value : scalar or array_like
        Density value or values.
    from_unit, to_unit : str
        Supported units are ``"kg m^-3"`` and ``"g cm^-3"`` together with
        common slash variants.

    Returns
    -------
    float or ndarray
        Converted density value or values.

    Raises
    ------
    NotImplementedError
        If one of the density units is unsupported.
    """
    from_key = _normalize_density_unit(from_unit)
    to_key = _normalize_density_unit(to_unit)
    if from_key not in _DENSITY_TO_KG_M3:
        raise NotImplementedError(f"unsupported density unit: {from_unit!r}")
    if to_key not in _DENSITY_TO_KG_M3:
        raise NotImplementedError(f"unsupported density unit: {to_unit!r}")
    array = np.asarray(value, dtype=np.float64)
    converted = array * _DENSITY_TO_KG_M3[from_key] / _DENSITY_TO_KG_M3[to_key]
    return float(converted) if np.isscalar(value) else converted


def density_from_formula(
    formula: str | Mapping[str, float],
    volume: Any,
    *,
    z: float = 1.0,
    volume_unit: str = "angstrom",
    density_unit: str = "kg m^-3",
) -> float | np.ndarray:
    """Return crystal density from formula, cell volume and formula units.

    Parameters
    ----------
    formula : str or mapping
        Chemical formula of the formula unit, or an element-count mapping. If
        the formula already represents the full unit-cell composition, use
        ``z=1``.
    volume : scalar or array_like
        Unit-cell volume.
    z : float, optional
        Number of formula units in the unit cell.
    volume_unit : str, optional
        Length unit defining the volume unit. For example, ``"angstrom"``
        means Angstrom cubed, matching :func:`quantas.core.physics.units.convert_volume`.
    density_unit : str, optional
        Output density unit. Supported values include ``"kg m^-3"`` and
        ``"g cm^-3"``.

    Returns
    -------
    float or ndarray
        Density in the requested unit.

    Raises
    ------
    ValueError
        If ``z`` or any volume is not finite and positive.
    NotImplementedError
        If a requested unit is unsupported.
    """
    z_value = float(z)
    if not np.isfinite(z_value) or z_value <= 0.0:
        raise ValueError("z must be finite and positive")
    volume_m3 = np.asarray(convert_volume(volume, volume_unit, "m"), dtype=np.float64)
    if not np.all(np.isfinite(volume_m3)) or np.any(volume_m3 <= 0.0):
        raise ValueError("cell volumes must be finite and positive")
    molar_mass_kg_mol = molar_mass(formula) / 1000.0
    mass_cell = z_value * molar_mass_kg_mol / cs.Avogadro
    density_si = mass_cell / volume_m3
    converted = convert_density(density_si, "kg m^-3", density_unit)
    return float(converted) if np.isscalar(volume) else converted
