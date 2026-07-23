# -*- coding: utf-8 -*-

"""Passive depth-path and evaluated-profile contracts."""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any
import numpy as np
from numpy.typing import NDArray
from quantas.core.physics.elasticity import StabilityFieldResult
from .types import FloatArray


@dataclass(slots=True)
class ThermoelasticDepthProfile:
    """Pressure-temperature path parameterized by geological depth.

    Parameters
    ----------
    name : str
        Profile identifier.
    depth : ndarray
        Depth values in km, strictly increasing.
    pressure : ndarray
        Pressure in GPa.
    temperature : ndarray
        Temperature in K.
    metadata : dict, optional
        Profile provenance and gradient descriptions.
    """

    name: str
    depth: FloatArray
    pressure: FloatArray
    temperature: FloatArray
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate aligned and monotonic path arrays."""
        self.name = str(self.name).strip()
        if not self.name or "/" in self.name or self.name in {".", ".."}:
            raise ValueError("profile name must be non-empty and cannot contain '/'")
        self.depth = np.asarray(self.depth, dtype=np.float64).copy()
        self.pressure = np.asarray(self.pressure, dtype=np.float64).copy()
        self.temperature = np.asarray(self.temperature, dtype=np.float64).copy()
        if self.depth.ndim != 1 or self.depth.size < 1:
            raise ValueError("depth must be a non-empty one-dimensional array")
        if (
            self.pressure.shape != self.depth.shape
            or self.temperature.shape != self.depth.shape
        ):
            raise ValueError(
                "depth, pressure, and temperature must have identical shapes"
            )
        if not all(
            np.all(np.isfinite(value))
            for value in (self.depth, self.pressure, self.temperature)
        ):
            raise ValueError("profile arrays must contain finite values")
        if np.any(np.diff(self.depth) <= 0.0):
            raise ValueError("depth values must be strictly increasing")
        if np.any(self.temperature < 0.0):
            raise ValueError("profile temperatures must be non-negative")
        self.metadata = dict(self.metadata)

    @classmethod
    def linear(
        cls,
        *,
        name: str,
        depth_min: float,
        depth_max: float,
        npoints: int,
        pressure_at_depth_min: float = 0.0,
        pressure_gradient: float = 0.03,
        temperature_at_depth_min: float = 298.15,
        temperature_gradient: float = 0.5,
    ) -> "ThermoelasticDepthProfile":
        """Build a profile from constant pressure and temperature gradients.

        Gradients are expressed in GPa km^-1 and K km^-1.
        """
        if npoints < 2:
            raise ValueError("npoints must be at least 2")
        depth = np.linspace(
            float(depth_min), float(depth_max), int(npoints), dtype=np.float64
        )
        offset = depth - depth[0]
        pressure = float(pressure_at_depth_min) + float(pressure_gradient) * offset
        temperature = (
            float(temperature_at_depth_min) + float(temperature_gradient) * offset
        )
        return cls(
            name=name,
            depth=depth,
            pressure=pressure,
            temperature=temperature,
            metadata={
                "kind": "linear_gradients",
                "pressure_gradient_GPa_per_km": float(pressure_gradient),
                "temperature_gradient_K_per_km": float(temperature_gradient),
            },
        )


@dataclass(slots=True)
class ThermoelasticProfileResult:
    """Reconstructed quasi-static tensors along one depth profile."""

    name: str
    depth: FloatArray
    pressure: FloatArray
    temperature: FloatArray
    volume: FloatArray
    density: FloatArray
    independent_stiffness: FloatArray
    sigma_independent_stiffness: FloatArray
    independent_stiffness_covariance: FloatArray
    stiffness_isothermal: FloatArray
    sigma_stiffness_isothermal: FloatArray
    qha_extrapolation_mask: NDArray[np.bool_]
    elastic_extrapolation_mask: NDArray[np.bool_]
    stiffness_adiabatic: FloatArray | None = None
    sigma_stiffness_adiabatic: FloatArray | None = None
    adiabatic_correction: FloatArray | None = None
    adiabatic_thermal_stress: FloatArray | None = None
    adiabatic_valid_mask: NDArray[np.bool_] | None = None
    stability: StabilityFieldResult | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize profile result arrays."""
        self.depth = np.asarray(self.depth, dtype=np.float64).copy()
        npoints = self.depth.size
        for name in ("pressure", "temperature", "volume", "density"):
            value = np.asarray(getattr(self, name), dtype=np.float64)
            if value.shape != (npoints,):
                raise ValueError(f"{name} must have shape (npoints,)")
            setattr(self, name, value.copy())
        independent = np.asarray(self.independent_stiffness, dtype=np.float64)
        sigma_independent = np.asarray(
            self.sigma_independent_stiffness, dtype=np.float64
        )
        if independent.ndim != 2 or independent.shape[0] != npoints:
            raise ValueError(
                "independent_stiffness must have shape (npoints, ncomponents)"
            )
        if sigma_independent.shape != independent.shape:
            raise ValueError(
                "sigma_independent_stiffness must match independent_stiffness"
            )
        covariance = np.asarray(self.independent_stiffness_covariance, dtype=np.float64)
        ncomponents = independent.shape[1]
        if covariance.shape != (npoints, ncomponents, ncomponents):
            raise ValueError(
                "independent_stiffness_covariance must have shape "
                "(npoints, ncomponents, ncomponents)"
            )
        self.independent_stiffness = independent.copy()
        self.sigma_independent_stiffness = sigma_independent.copy()
        self.independent_stiffness_covariance = covariance.copy()
        for name in ("stiffness_isothermal", "sigma_stiffness_isothermal"):
            value = np.asarray(getattr(self, name), dtype=np.float64)
            if value.shape != (npoints, 6, 6):
                raise ValueError(f"{name} must have shape (npoints, 6, 6)")
            setattr(self, name, value.copy())
        for name in (
            "stiffness_adiabatic",
            "sigma_stiffness_adiabatic",
            "adiabatic_correction",
        ):
            value = getattr(self, name)
            if value is not None:
                array = np.asarray(value, dtype=np.float64)
                if array.shape != (npoints, 6, 6):
                    raise ValueError(f"{name} must have shape (npoints, 6, 6)")
                setattr(self, name, array.copy())
        if self.adiabatic_thermal_stress is not None:
            thermal_stress = np.asarray(self.adiabatic_thermal_stress, dtype=np.float64)
            if thermal_stress.shape != (npoints, 6):
                raise ValueError(
                    "adiabatic_thermal_stress must have shape (npoints, 6)"
                )
            self.adiabatic_thermal_stress = thermal_stress.copy()
        if self.adiabatic_valid_mask is not None:
            valid = np.asarray(self.adiabatic_valid_mask, dtype=np.bool_)
            if valid.shape != (npoints,):
                raise ValueError("adiabatic_valid_mask must have shape (npoints,)")
            self.adiabatic_valid_mask = valid.copy()
        self.qha_extrapolation_mask = np.asarray(
            self.qha_extrapolation_mask, dtype=np.bool_
        ).copy()
        self.elastic_extrapolation_mask = np.asarray(
            self.elastic_extrapolation_mask, dtype=np.bool_
        ).copy()
        if self.qha_extrapolation_mask.shape != (
            npoints,
        ) or self.elastic_extrapolation_mask.shape != (npoints,):
            raise ValueError("profile masks must have shape (npoints,)")
        if self.stability is not None:
            if self.stability.minimum_eigenvalue.shape != (npoints,):
                raise ValueError("profile stability must have shape (npoints,)")
        self.metadata = dict(self.metadata)


__all__ = ["ThermoelasticDepthProfile", "ThermoelasticProfileResult"]
