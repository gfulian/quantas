# -*- coding: utf-8 -*-

"""Aggregated calibration and analysis result contracts."""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any
import numpy as np
from numpy.typing import NDArray
from quantas.core.physics.elasticity import StabilityFieldResult
from .fields import ThermoelasticProfileResult
from .fitting import ElasticComponentFit, ReferenceEOSFit
from .types import FloatArray, ThermoelasticTensorCondition


@dataclass(slots=True)
class ThermoelasticResult:
    """Complete quasi-static thermoelastic fit and reconstruction result.

    Parameters
    ----------
    jobname : str
        Workflow description.
    reference_eos : ReferenceEOSFit
        Fixed static EOS used by every elastic-component fit.
    component_fits : dict
        Independent component fit records keyed by label.
    independent_labels : tuple of str
        Component order used in reconstructed arrays.
    temperature, pressure : ndarray
        QHA grid coordinates in K and GPa.
    equilibrium_volume : ndarray
        QHA equilibrium volumes in angstrom cubed.
    density : ndarray
        Density in kg m^-3.
    independent_stiffness, sigma_independent_stiffness : ndarray
        Independent component values and one-sigma uncertainties.
    independent_stiffness_covariance : ndarray
        Full covariance among independent components at every grid point.
    stiffness_isothermal, sigma_stiffness_isothermal : ndarray
        Full symmetry-constrained Voigt tensors and one-sigma uncertainties.
    extrapolation_mask : ndarray
        Elastic-volume extrapolation mask.
    sigma_equilibrium_volume : ndarray or None, optional
        One-sigma QHA equilibrium-volume uncertainty in angstrom cubed.
    qha_extrapolation_mask : ndarray or None, optional
        Mask marking requested states outside the archived QHA coordinate grid.
    profiles : dict, optional
        Named geological depth-profile results.
    stability : StabilityFieldResult or None, optional
        Positive-definiteness diagnostics for reconstructed Wallace stiffness
        matrices on the pressure-temperature grid.
    completed : bool, optional
        Whether all active fits and grid reconstruction completed.
    metadata : dict, optional
        Scientific and uncertainty-propagation metadata.
    """

    jobname: str
    reference_eos: ReferenceEOSFit
    component_fits: dict[str, ElasticComponentFit]
    independent_labels: tuple[str, ...]
    temperature: FloatArray
    pressure: FloatArray
    equilibrium_volume: FloatArray
    density: FloatArray
    independent_stiffness: FloatArray | None
    sigma_independent_stiffness: FloatArray | None
    independent_stiffness_covariance: FloatArray | None
    stiffness_isothermal: FloatArray | None
    sigma_stiffness_isothermal: FloatArray | None
    extrapolation_mask: NDArray[np.bool_]
    sigma_equilibrium_volume: FloatArray | None = None
    qha_extrapolation_mask: NDArray[np.bool_] | None = None
    profiles: dict[str, ThermoelasticProfileResult] = field(default_factory=dict)
    isochoric_heat_capacity_cell: FloatArray | None = None
    sigma_isochoric_heat_capacity_cell: FloatArray | None = None
    thermal_expansion_tensor: FloatArray | None = None
    sigma_thermal_expansion_tensor: FloatArray | None = None
    stiffness_adiabatic: FloatArray | None = None
    sigma_stiffness_adiabatic: FloatArray | None = None
    adiabatic_correction: FloatArray | None = None
    adiabatic_thermal_stress: FloatArray | None = None
    adiabatic_valid_mask: NDArray[np.bool_] | None = None
    stability: StabilityFieldResult | None = None
    completed: bool = True
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize principal arrays and validate grid dimensions."""
        self.temperature = np.asarray(self.temperature, dtype=np.float64).copy()
        self.pressure = np.asarray(self.pressure, dtype=np.float64).copy()
        self.equilibrium_volume = np.asarray(
            self.equilibrium_volume, dtype=np.float64
        ).copy()
        self.density = np.asarray(self.density, dtype=np.float64).copy()
        shape = (self.temperature.size, self.pressure.size)
        if self.equilibrium_volume.shape != shape or self.density.shape != shape:
            raise ValueError("thermoelastic volume and density must match the T-P grid")
        self.extrapolation_mask = np.asarray(
            self.extrapolation_mask, dtype=np.bool_
        ).copy()
        if self.extrapolation_mask.shape != shape:
            raise ValueError("extrapolation_mask must match the T-P grid")
        if self.sigma_equilibrium_volume is not None:
            sigma_volume = np.asarray(self.sigma_equilibrium_volume, dtype=np.float64)
            if sigma_volume.shape != shape:
                raise ValueError("sigma_equilibrium_volume must match the T-P grid")
            if np.any(~np.isfinite(sigma_volume)) or np.any(sigma_volume < 0.0):
                raise ValueError(
                    "sigma_equilibrium_volume must be finite and non-negative"
                )
            self.sigma_equilibrium_volume = sigma_volume.copy()
        if self.qha_extrapolation_mask is None:
            self.qha_extrapolation_mask = np.zeros(shape, dtype=np.bool_)
        else:
            qha_extrapolation = np.asarray(self.qha_extrapolation_mask, dtype=np.bool_)
            if qha_extrapolation.shape != shape:
                raise ValueError("qha_extrapolation_mask must match the T-P grid")
            self.qha_extrapolation_mask = qha_extrapolation.copy()
        ncomponents = len(self.independent_labels)
        if self.independent_stiffness is not None:
            self.independent_stiffness = np.asarray(
                self.independent_stiffness, dtype=np.float64
            ).copy()
            if self.independent_stiffness.shape != shape + (ncomponents,):
                raise ValueError("independent_stiffness has an invalid shape")
        if self.sigma_independent_stiffness is not None:
            self.sigma_independent_stiffness = np.asarray(
                self.sigma_independent_stiffness, dtype=np.float64
            ).copy()
            if self.sigma_independent_stiffness.shape != shape + (ncomponents,):
                raise ValueError("sigma_independent_stiffness has an invalid shape")
        if self.independent_stiffness_covariance is not None:
            self.independent_stiffness_covariance = np.asarray(
                self.independent_stiffness_covariance, dtype=np.float64
            ).copy()
            if self.independent_stiffness_covariance.shape != shape + (
                ncomponents,
                ncomponents,
            ):
                raise ValueError(
                    "independent_stiffness_covariance has an invalid shape"
                )
        for name in ("stiffness_isothermal", "sigma_stiffness_isothermal"):
            value = getattr(self, name)
            if value is not None:
                array = np.asarray(value, dtype=np.float64)
                if array.shape != shape + (6, 6):
                    raise ValueError(f"{name} has an invalid shape")
                setattr(self, name, array.copy())
        for name in (
            "isochoric_heat_capacity_cell",
            "sigma_isochoric_heat_capacity_cell",
        ):
            value = getattr(self, name)
            if value is not None:
                array = np.asarray(value, dtype=np.float64)
                if array.shape != shape:
                    raise ValueError(f"{name} must match the T-P grid")
                if name.startswith("sigma_") and np.any(array < 0.0):
                    raise ValueError(f"{name} must be non-negative")
                setattr(self, name, array.copy())
        for name in (
            "thermal_expansion_tensor",
            "sigma_thermal_expansion_tensor",
        ):
            value = getattr(self, name)
            if value is not None:
                array = np.asarray(value, dtype=np.float64)
                if array.shape != shape + (3, 3):
                    raise ValueError(f"{name} must have shape T-P + (3, 3)")
                if name.startswith("sigma_") and np.any(array < 0.0):
                    raise ValueError(f"{name} must be non-negative")
                setattr(self, name, array.copy())
        for name in (
            "stiffness_adiabatic",
            "sigma_stiffness_adiabatic",
            "adiabatic_correction",
        ):
            value = getattr(self, name)
            if value is not None:
                array = np.asarray(value, dtype=np.float64)
                if array.shape != shape + (6, 6):
                    raise ValueError(f"{name} has an invalid shape")
                setattr(self, name, array.copy())
        if self.adiabatic_thermal_stress is not None:
            thermal_stress = np.asarray(self.adiabatic_thermal_stress, dtype=np.float64)
            if thermal_stress.shape != shape + (6,):
                raise ValueError("adiabatic_thermal_stress has an invalid shape")
            self.adiabatic_thermal_stress = thermal_stress.copy()
        if self.adiabatic_valid_mask is not None:
            valid = np.asarray(self.adiabatic_valid_mask, dtype=np.bool_)
            if valid.shape != shape:
                raise ValueError("adiabatic_valid_mask must match the T-P grid")
            self.adiabatic_valid_mask = valid.copy()
        if self.stability is not None:
            if self.stability.minimum_eigenvalue.shape != shape:
                raise ValueError("stability field must match the T-P grid")
        self.component_fits = dict(self.component_fits)
        self.independent_labels = tuple(str(value) for value in self.independent_labels)
        self.profiles = dict(self.profiles)
        self.metadata = dict(self.metadata)


def select_stiffness_tensor(
    result: ThermoelasticResult | ThermoelasticProfileResult,
    condition: ThermoelasticTensorCondition = "isothermal",
) -> tuple[FloatArray, FloatArray | None]:
    """Return stiffness and uncertainty for one thermodynamic condition.

    Parameters
    ----------
    result : ThermoelasticResult or ThermoelasticProfileResult
        Reconstructed grid or depth-profile result.
    condition : {"isothermal", "adiabatic"}, optional
        Requested thermodynamic condition.

    Returns
    -------
    tuple of ndarray and ndarray or None
        Stiffness and optional one-sigma uncertainty in GPa.

    Raises
    ------
    ValueError
        If the requested tensor is unavailable or contains invalid states.
    """
    if condition == "isothermal":
        stiffness = result.stiffness_isothermal
        sigma = result.sigma_stiffness_isothermal
    elif condition == "adiabatic":
        stiffness = result.stiffness_adiabatic
        sigma = result.sigma_stiffness_adiabatic
        valid = result.adiabatic_valid_mask
        if valid is not None and not np.all(valid):
            raise ValueError(
                "adiabatic stiffness contains invalid states; select a valid "
                "pressure-temperature point explicitly"
            )
    else:
        raise ValueError("condition must be 'isothermal' or 'adiabatic'")
    if stiffness is None:
        raise ValueError(f"{condition} stiffness tensor is unavailable")
    return np.asarray(stiffness, dtype=np.float64), (
        None if sigma is None else np.asarray(sigma, dtype=np.float64)
    )


__all__ = ["ThermoelasticResult", "select_stiffness_tensor"]
