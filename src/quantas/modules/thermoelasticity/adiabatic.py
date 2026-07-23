# -*- coding: utf-8 -*-

"""Thermoelastic adapters for anisotropic adiabatic stiffness conversion."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from quantas.core.physics.elasticity import (
    AdiabaticStiffnessFieldResult,
    adiabatic_stiffness_field,
)
from quantas.core.physics.units import (
    convert_energy_per_temperature,
    convert_volume,
)
from quantas.modules.thermoelasticity.models import ThermoelasticAdiabaticMode


@dataclass(slots=True)
class StandardAdiabaticInputs:
    r"""QHA thermoelastic inputs normalized to per-cell SI units.

    Parameters
    ----------
    heat_capacity : ndarray
        Isochoric heat capacity in J K\ :sup:`-1` per normalized cell.
    thermal_expansion_tensor : ndarray
        Cartesian thermal-expansion tensor in K\ :sup:`-1`.
    sigma_heat_capacity : ndarray or None
        One-standard-deviation heat-capacity uncertainty in matching units.
    sigma_thermal_expansion_tensor : ndarray or None
        One-standard-deviation tensor uncertainty in K\ :sup:`-1`.
    metadata : dict
        Unit and provenance information.
    """

    heat_capacity: np.ndarray
    thermal_expansion_tensor: np.ndarray
    sigma_heat_capacity: np.ndarray | None = None
    sigma_thermal_expansion_tensor: np.ndarray | None = None
    metadata: dict[str, Any] | None = None


def standard_adiabatic_inputs_from_qha(
    qha: Any,
    qha_options: Mapping[str, Any],
    target_shape: tuple[int, ...],
) -> StandardAdiabaticInputs | None:
    """Return QHA ``C_V`` and expansion tensors in per-cell SI units.

    Parameters
    ----------
    qha : object
        QHA result payload.
    qha_options : mapping
        QHA unit options. Heat capacity is assumed to use the documented native
        normalization ``energy_unit cell^-1 K^-1``.
    target_shape : tuple of int
        Required pressure-temperature field shape.

    Returns
    -------
    StandardAdiabaticInputs or None
        Normalized inputs, or ``None`` when either mandatory field is absent or
        shape-incompatible.
    """
    cv_value = getattr(qha, "isochoric_heat_capacity", None)
    alpha_value = getattr(qha, "thermal_expansion_tensor", None)
    if cv_value is None or alpha_value is None:
        return None
    cv = np.asarray(cv_value, dtype=np.float64)
    alpha = np.asarray(alpha_value, dtype=np.float64)
    if cv.shape != target_shape or alpha.shape != target_shape + (3, 3):
        return None
    energy_unit = str(qha_options.get("energy_unit", "Ha"))
    cv_j = np.asarray(
        convert_energy_per_temperature(
            cv,
            f"{energy_unit} cell^-1 K^-1",
            "J cell^-1 K^-1",
        ),
        dtype=np.float64,
    )
    uncertainties = getattr(qha, "uncertainties", {})
    sigma_cv = _uncertainty_array(
        uncertainties,
        ("sigma_Cv", "sigma_isochoric_heat_capacity", "Cv"),
        target_shape,
    )
    if sigma_cv is not None:
        sigma_cv = np.asarray(
            convert_energy_per_temperature(
                sigma_cv,
                f"{energy_unit} cell^-1 K^-1",
                "J cell^-1 K^-1",
            ),
            dtype=np.float64,
        )
    sigma_alpha = _uncertainty_array(
        uncertainties,
        (
            "sigma_alpha_tensor",
            "sigma_thermal_expansion_tensor",
            "alpha_tensor",
        ),
        target_shape + (3, 3),
    )
    return StandardAdiabaticInputs(
        heat_capacity=cv_j,
        thermal_expansion_tensor=alpha.copy(),
        sigma_heat_capacity=sigma_cv,
        sigma_thermal_expansion_tensor=sigma_alpha,
        metadata={
            "source": "QHA equilibrium thermodynamic and structural fields",
            "native_heat_capacity_unit": f"{energy_unit} cell^-1 K^-1",
            "stored_heat_capacity_unit": "J cell^-1 K^-1",
            "thermal_expansion_unit": "K^-1",
            "normalization": "same normalized cell as equilibrium_volume",
        },
    )


def calculate_adiabatic_field(
    *,
    stiffness_isothermal: np.ndarray,
    sigma_stiffness_isothermal: np.ndarray | None,
    temperature_k: np.ndarray,
    volume_a3: np.ndarray,
    sigma_volume_a3: np.ndarray | None,
    inputs: StandardAdiabaticInputs | None,
    mode: ThermoelasticAdiabaticMode,
    propagate_uncertainty: bool,
) -> AdiabaticStiffnessFieldResult | None:
    """Calculate an adiabatic field according to the configured availability mode.

    Parameters
    ----------
    stiffness_isothermal, sigma_stiffness_isothermal : ndarray or None
        Isothermal stiffness and uncertainty in GPa.
    temperature_k : ndarray
        Absolute temperature field in K.
    volume_a3, sigma_volume_a3 : ndarray or None
        Cell volume and uncertainty in angstrom cubed.
    inputs : StandardAdiabaticInputs or None
        Standardized QHA heat capacity and expansion data.
    mode : {"auto", "off", "require"}
        Availability policy.
    propagate_uncertainty : bool
        Include available first-order uncertainties.

    Returns
    -------
    AdiabaticStiffnessFieldResult or None
        Converted field, or ``None`` when disabled/unavailable in ``auto`` mode.

    Raises
    ------
    ValueError
        If ``require`` is selected and complete inputs are unavailable or any
        nonzero-temperature state is invalid.
    """
    if mode == "off":
        return None
    if inputs is None:
        if mode == "require":
            raise ValueError(
                "adiabatic conversion requires QHA isochoric heat capacity and "
                "Cartesian thermal-expansion tensor"
            )
        return None
    volume_m3 = np.asarray(convert_volume(volume_a3, "A", "m"), dtype=np.float64)
    sigma_volume_m3 = None
    if sigma_volume_a3 is not None and propagate_uncertainty:
        sigma_volume_m3 = np.asarray(
            convert_volume(sigma_volume_a3, "A", "m"), dtype=np.float64
        )
    result = adiabatic_stiffness_field(
        stiffness_isothermal,
        temperature_k,
        volume_m3,
        inputs.heat_capacity,
        inputs.thermal_expansion_tensor,
        sigma_stiffness_isothermal=(
            sigma_stiffness_isothermal if propagate_uncertainty else None
        ),
        sigma_volume_m3=sigma_volume_m3,
        sigma_heat_capacity_j_per_k=(
            inputs.sigma_heat_capacity if propagate_uncertainty else None
        ),
        sigma_thermal_expansion_tensor=(
            inputs.sigma_thermal_expansion_tensor if propagate_uncertainty else None
        ),
    )
    if mode == "require" and np.any(result.invalid_mask):
        count = int(np.count_nonzero(result.invalid_mask))
        raise ValueError(f"adiabatic conversion is invalid at {count} requested states")
    result.metadata.update(inputs.metadata or {})
    return result


def _uncertainty_array(
    uncertainties: object,
    keys: tuple[str, ...],
    shape: tuple[int, ...],
) -> np.ndarray | None:
    if not isinstance(uncertainties, Mapping):
        return None
    for key in keys:
        if key not in uncertainties:
            continue
        value = np.asarray(uncertainties[key], dtype=np.float64)
        if value.shape != shape or np.any(value < 0.0):
            return None
        return value.copy()
    return None


__all__ = [
    "StandardAdiabaticInputs",
    "calculate_adiabatic_field",
    "standard_adiabatic_inputs_from_qha",
]
