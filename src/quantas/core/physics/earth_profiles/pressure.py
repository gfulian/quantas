# -*- coding: utf-8 -*-

"""Pressure-depth models for terrestrial thermoelastic paths."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

from quantas.references import get_citation, render_citation_inline

from .models import FloatArray


_PREM_KEY = "prem_dziewonski_anderson_1981"
_PREM_REFERENCE = get_citation(_PREM_KEY)
_PREM_CITATION = render_citation_inline(_PREM_KEY)


@dataclass(frozen=True, slots=True)
class LithostaticLayer:
    """One constant-density layer in a lithostatic pressure model.

    Parameters
    ----------
    thickness_km : float
        Layer thickness in km.
    density_kg_m3 : float
        Constant mass density in kg m^-3.
    name : str, optional
        Human-readable layer identifier.

    Raises
    ------
    ValueError
        If thickness or density is not finite and positive.
    """

    thickness_km: float
    density_kg_m3: float
    name: str = "layer"

    def __post_init__(self) -> None:
        """Validate the passive layer specification."""
        if not np.isfinite(self.thickness_km) or self.thickness_km <= 0.0:
            raise ValueError("layer thickness_km must be finite and positive")
        if not np.isfinite(self.density_kg_m3) or self.density_kg_m3 <= 0.0:
            raise ValueError("layer density_kg_m3 must be finite and positive")


class PremPressureModel:
    """Hydrostatic pressure reconstructed from the PREM density profile.

    The density polynomials are those in Table I of the Preliminary Reference
    Earth Model.  Enclosed mass is integrated from the centre, gravity is
    evaluated as ``G m(r) / r**2``, and pressure is integrated inward from the
    surface according to ``dP/dr = -rho(r) g(r)``.  The complete PREM density
    distribution is used internally so that mantle gravity includes the core
    mass, even when public evaluation is restricted to the surface--CMB range.

    References
    ----------
    Canonical citation key: ``prem_dziewonski_anderson_1981``.

    Parameters
    ----------
    max_depth_km : float, optional
        Deepest public evaluation depth.  The default, 2891 km, approximately
        reaches the core-mantle boundary while excluding core profiles from
        mineral thermoelastic workflows.
    integration_step_km : float, optional
        Maximum radial spacing used for mass and pressure integration.

    Raises
    ------
    ValueError
        If the public depth range or integration spacing is invalid.
    """

    earth_radius_km = 6371.0
    core_mantle_boundary_depth_km = 2891.0
    gravitational_constant = 6.67430e-11

    def __init__(
        self,
        *,
        max_depth_km: float = core_mantle_boundary_depth_km,
        integration_step_km: float = 0.25,
    ) -> None:
        self._max_depth_km = float(max_depth_km)
        self._integration_step_km = float(integration_step_km)
        if (
            not np.isfinite(self._max_depth_km)
            or self._max_depth_km <= 0.0
            or self._max_depth_km > self.earth_radius_km
        ):
            raise ValueError("max_depth_km must lie in (0, 6371]")
        if (
            not np.isfinite(self._integration_step_km)
            or self._integration_step_km <= 0.0
            or self._integration_step_km > 10.0
        ):
            raise ValueError("integration_step_km must lie in (0, 10]")
        self._prepare_reference_profile()

    @property
    def name(self) -> str:
        """Return the stable pressure-model identifier."""
        return "prem"

    @property
    def depth_bounds(self) -> tuple[float, float]:
        """Return the supported public depth interval in km."""
        return 0.0, self._max_depth_km

    @property
    def critical_depths(self) -> tuple[float, ...]:
        """Return PREM layer boundaries within the public depth domain."""
        radius_boundaries = (
            6368.0,
            6356.0,
            6346.6,
            6291.0,
            6151.0,
            5971.0,
            5771.0,
            5701.0,
            5600.0,
            3630.0,
            3480.0,
        )
        depths = tuple(self.earth_radius_km - value for value in radius_boundaries)
        return tuple(value for value in depths if value <= self._max_depth_km)

    @property
    def earth_mass_kg(self) -> float:
        """Return the mass implied by the integrated PREM density profile."""
        return self._earth_mass_kg

    @property
    def surface_gravity_m_s2(self) -> float:
        """Return surface gravity implied by the integrated PREM density."""
        return self._surface_gravity_m_s2

    def pressure(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate PREM hydrostatic pressure.

        Parameters
        ----------
        depth_km : ndarray
            Geological depth in km.

        Returns
        -------
        ndarray
            Hydrostatic pressure in GPa with the input shape.

        Raises
        ------
        ValueError
            If depths are non-finite or outside ``depth_bounds``.
        """
        depth = _validated_depth(depth_km, self.depth_bounds, self.name)
        result = np.asarray(self._pressure_interpolator(depth), dtype=np.float64)
        result[np.abs(result) < 1.0e-13] = 0.0
        return result

    def density(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate PREM density in kg m^-3 within the public domain."""
        depth = _validated_depth(depth_km, self.depth_bounds, self.name)
        radius = self.earth_radius_km - depth
        return 1000.0 * self._density_g_cm3(radius)

    def gravity(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate gravity in m s^-2 from the integrated PREM mass profile."""
        depth = _validated_depth(depth_km, self.depth_bounds, self.name)
        return np.asarray(self._gravity_interpolator(depth), dtype=np.float64)

    def metadata(self) -> dict[str, Any]:
        """Return model parameters and complete bibliographic provenance."""
        return {
            "model": self.name,
            "kind": "hydrostatic_reference_earth_model",
            "depth_unit": "km",
            "pressure_unit": "GPa",
            "density_unit": "kg m^-3",
            "max_depth_km": self._max_depth_km,
            "integration_step_km": self._integration_step_km,
            "earth_radius_km": self.earth_radius_km,
            "earth_mass_kg": self.earth_mass_kg,
            "surface_gravity_m_s2": self.surface_gravity_m_s2,
            "citation": _PREM_CITATION,
            "citation_key": _PREM_KEY,
            "doi": _PREM_REFERENCE.doi,
            "scientific_scope": (
                "One-dimensional global-average radial pressure reference; "
                "not a site-specific shallow-crust lithostatic model."
            ),
        }

    def _prepare_reference_profile(self) -> None:
        boundaries = np.asarray(
            [
                0.0,
                1221.5,
                3480.0,
                3630.0,
                5600.0,
                5701.0,
                5771.0,
                5971.0,
                6151.0,
                6291.0,
                6346.6,
                6356.0,
                6368.0,
                6371.0,
            ],
            dtype=np.float64,
        )
        regular = np.arange(
            0.0,
            self.earth_radius_km,
            self._integration_step_km,
            dtype=np.float64,
        )
        radius_km = np.unique(np.concatenate((regular, boundaries)))
        if radius_km[-1] != self.earth_radius_km:
            radius_km = np.append(radius_km, self.earth_radius_km)
        radius_m = radius_km * 1000.0
        density = 1000.0 * self._density_g_cm3(radius_km)
        shell_integrand = 4.0 * np.pi * radius_m**2 * density
        enclosed_mass = cumulative_trapezoid(
            shell_integrand,
            radius_m,
            initial=0.0,
        )
        gravity = np.zeros_like(radius_m)
        nonzero = radius_m > 0.0
        gravity[nonzero] = (
            self.gravitational_constant
            * enclosed_mass[nonzero]
            / radius_m[nonzero] ** 2
        )
        pressure_integrand = density * gravity
        reversed_integral = cumulative_trapezoid(
            pressure_integrand[::-1],
            radius_m[::-1],
            initial=0.0,
        )
        pressure_pa = -reversed_integral[::-1]
        depth_km = self.earth_radius_km - radius_km
        order = np.argsort(depth_km)
        depth_sorted = depth_km[order]
        pressure_gpa = pressure_pa[order] / 1.0e9
        gravity_sorted = gravity[order]
        self._pressure_interpolator = PchipInterpolator(
            depth_sorted,
            pressure_gpa,
            extrapolate=False,
        )
        self._gravity_interpolator = PchipInterpolator(
            depth_sorted,
            gravity_sorted,
            extrapolate=False,
        )
        self._earth_mass_kg = float(enclosed_mass[-1])
        self._surface_gravity_m_s2 = float(gravity[-1])

    @classmethod
    def _density_g_cm3(cls, radius_km: NDArray[np.float64]) -> FloatArray:
        radius = np.asarray(radius_km, dtype=np.float64)
        x = radius / cls.earth_radius_km
        rho = np.empty_like(x)
        regions = (
            (radius <= 1221.5, 13.0885 - 8.8381 * x**2),
            (
                (radius > 1221.5) & (radius <= 3480.0),
                12.5815 - 1.2638 * x - 3.6426 * x**2 - 5.5281 * x**3,
            ),
            (
                (radius > 3480.0) & (radius <= 5701.0),
                7.9565 - 6.4761 * x + 5.5283 * x**2 - 3.0807 * x**3,
            ),
            (
                (radius > 5701.0) & (radius <= 5771.0),
                5.3197 - 1.4836 * x,
            ),
            (
                (radius > 5771.0) & (radius <= 5971.0),
                11.2494 - 8.0298 * x,
            ),
            (
                (radius > 5971.0) & (radius <= 6151.0),
                7.1089 - 3.8045 * x,
            ),
            (
                (radius > 6151.0) & (radius <= 6346.6),
                2.6910 + 0.6924 * x,
            ),
            ((radius > 6346.6) & (radius <= 6356.0), np.full_like(x, 2.900)),
            ((radius > 6356.0) & (radius <= 6368.0), np.full_like(x, 2.600)),
            ((radius > 6368.0) & (radius <= 6371.0), np.full_like(x, 1.020)),
        )
        assigned = np.zeros_like(x, dtype=np.bool_)
        for mask, values in regions:
            rho[mask] = values[mask] if values.shape == x.shape else values
            assigned |= mask
        if np.any(~assigned):
            raise ValueError("radius lies outside the PREM Earth radius")
        return rho


class LayeredLithostaticPressureModel:
    """Piecewise-constant-density lithostatic pressure profile.

    Pressure is accumulated downward as ``dP/dz = rho g`` with a constant
    gravitational acceleration.  This model is intended for site-specific
    crustal or lithospheric profiles where PREM's global-average shallow
    layering is inappropriate.

    Parameters
    ----------
    layers : sequence of LithostaticLayer
        Ordered layers from the surface downward.
    gravity_m_s2 : float, optional
        Constant gravitational acceleration.
    name : str, optional
        Stable model identifier.
    citation : str or None, optional
        Complete citation or provenance for user-supplied layer parameters.

    Raises
    ------
    ValueError
        If no layers are provided or gravity is invalid.
    """

    def __init__(
        self,
        layers: Sequence[LithostaticLayer],
        *,
        gravity_m_s2: float = 9.80665,
        name: str = "layered-lithostatic",
        citation: str | None = None,
    ) -> None:
        self.layers = tuple(layers)
        if not self.layers:
            raise ValueError("at least one lithostatic layer is required")
        self.gravity_m_s2 = float(gravity_m_s2)
        if not np.isfinite(self.gravity_m_s2) or self.gravity_m_s2 <= 0.0:
            raise ValueError("gravity_m_s2 must be finite and positive")
        self._name = str(name)
        self._citation = citation
        self._boundaries = np.cumsum(
            [layer.thickness_km for layer in self.layers], dtype=np.float64
        )

    @property
    def name(self) -> str:
        """Return the stable pressure-model identifier."""
        return self._name

    @property
    def depth_bounds(self) -> tuple[float, float]:
        """Return the supported depth interval in km."""
        return 0.0, float(self._boundaries[-1])

    @property
    def critical_depths(self) -> tuple[float, ...]:
        """Return internal layer boundaries in km."""
        return tuple(float(value) for value in self._boundaries[:-1])

    def pressure(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate lithostatic pressure in GPa."""
        depth = _validated_depth(depth_km, self.depth_bounds, self.name)
        flat = depth.ravel()
        pressure_pa = np.zeros_like(flat)
        top = 0.0
        for layer, bottom in zip(self.layers, self._boundaries, strict=True):
            within = np.clip(flat - top, 0.0, bottom - top)
            pressure_pa += layer.density_kg_m3 * self.gravity_m_s2 * within * 1000.0
            top = float(bottom)
        return (pressure_pa.reshape(depth.shape) / 1.0e9).astype(np.float64)

    def metadata(self) -> dict[str, Any]:
        """Return layer parameters and optional user provenance."""
        return {
            "model": self.name,
            "kind": "piecewise_constant_density_lithostatic",
            "depth_unit": "km",
            "pressure_unit": "GPa",
            "gravity_m_s2": self.gravity_m_s2,
            "layers": [
                {
                    "name": layer.name,
                    "thickness_km": layer.thickness_km,
                    "density_kg_m3": layer.density_kg_m3,
                }
                for layer in self.layers
            ],
            "citation": self._citation or "User-specified layered lithostatic model.",
        }


def _validated_depth(
    values: NDArray[np.float64],
    bounds: tuple[float, float],
    model_name: str,
) -> FloatArray:
    depth = np.asarray(values, dtype=np.float64)
    if np.any(~np.isfinite(depth)):
        raise ValueError(f"{model_name} depths must be finite")
    tolerance = max(1.0, abs(bounds[1])) * 1.0e-10
    if np.any(depth < bounds[0] - tolerance) or np.any(depth > bounds[1] + tolerance):
        raise ValueError(
            f"{model_name} supports depths from {bounds[0]:g} to {bounds[1]:g} km"
        )
    return np.clip(depth, bounds[0], bounds[1]).astype(np.float64, copy=False)


__all__ = [
    "LayeredLithostaticPressureModel",
    "LithostaticLayer",
    "PremPressureModel",
]
