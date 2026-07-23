# -*- coding: utf-8 -*-

"""Temperature-depth models for terrestrial thermoelastic paths."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import CubicHermiteSpline
from scipy.special import erf

from quantas.references import get_citation, render_citation_inline

from .models import FloatArray
from .pressure import _validated_depth


_HASTEROK_CHAPMAN_KEY = "hasterok_chapman_2011"
_HASTEROK_CHAPMAN_REFERENCE = get_citation(_HASTEROK_CHAPMAN_KEY)
_HASTEROK_CHAPMAN_CITATION = render_citation_inline(_HASTEROK_CHAPMAN_KEY)
_PARSONS_SCLATER_KEY = "parsons_sclater_1977"
_PARSONS_SCLATER_REFERENCE = get_citation(_PARSONS_SCLATER_KEY)
_PARSONS_SCLATER_CITATION = render_citation_inline(_PARSONS_SCLATER_KEY)
_KATSURA_KEY = "katsura_2022"
_KATSURA_SOFTWARE_KEY = "katsura_software_2022"
_KATSURA_REFERENCE = get_citation(_KATSURA_KEY)
_KATSURA_SOFTWARE_REFERENCE = get_citation(_KATSURA_SOFTWARE_KEY)
_KATSURA_CITATION = render_citation_inline(_KATSURA_KEY)
_KATSURA_SOFTWARE_CITATION = render_citation_inline(_KATSURA_SOFTWARE_KEY)


@dataclass(frozen=True, slots=True)
class ConductiveLayer:
    """One homogeneous layer in a steady conductive geotherm.

    Parameters
    ----------
    thickness_km : float
        Layer thickness in km.
    conductivity_W_mK : float
        Thermal conductivity in W m^-1 K^-1.
    heat_production_uW_m3 : float
        Volumetric radiogenic heat production in microW m^-3.
    name : str, optional
        Human-readable layer identifier.

    Raises
    ------
    ValueError
        If thickness or conductivity is not finite and positive, or heat
        production is non-finite or negative.
    """

    thickness_km: float
    conductivity_W_mK: float
    heat_production_uW_m3: float
    name: str = "layer"

    def __post_init__(self) -> None:
        """Validate the passive conductive-layer specification."""
        if not np.isfinite(self.thickness_km) or self.thickness_km <= 0.0:
            raise ValueError("layer thickness_km must be finite and positive")
        if not np.isfinite(self.conductivity_W_mK) or self.conductivity_W_mK <= 0.0:
            raise ValueError("layer conductivity_W_mK must be finite and positive")
        if (
            not np.isfinite(self.heat_production_uW_m3)
            or self.heat_production_uW_m3 < 0.0
        ):
            raise ValueError(
                "layer heat_production_uW_m3 must be finite and non-negative"
            )


class ContinentalConductiveGeotherm:
    """Steady one-dimensional conductive continental geotherm.

    Each layer solves

    ``d/dz(k dT/dz) + A = 0``

    for constant conductivity ``k`` and volumetric heat production ``A``.
    Temperature and upward heat flow are propagated continuously from the
    surface through the ordered layers.  The implementation is a general
    layered conductive model; parameter sets supplied as Quantas presets are
    representative members of the surface-heat-flow geotherm family discussed
    by Hasterok and Chapman rather than digitized site-specific curves.

    References
    ----------
    Canonical citation key: ``hasterok_chapman_2011``.

    Parameters
    ----------
    layers : sequence of ConductiveLayer
        Layers ordered from the surface downward.
    surface_temperature_K : float, optional
        Temperature at zero depth in K.
    surface_heat_flow_mW_m2 : float, optional
        Upward surface heat flow in mW m^-2.
    name : str, optional
        Stable model identifier.
    citation : str or None, optional
        Alternative complete provenance for a custom parameterization.

    Raises
    ------
    ValueError
        If layers, surface temperature, or heat flow are invalid, or if the
        resulting model predicts negative absolute temperatures.
    """

    def __init__(
        self,
        layers: Sequence[ConductiveLayer],
        *,
        surface_temperature_K: float = 288.15,
        surface_heat_flow_mW_m2: float = 40.0,
        name: str = "continental-conductive",
        citation: str | None = None,
    ) -> None:
        self.layers = tuple(layers)
        if not self.layers:
            raise ValueError("at least one conductive layer is required")
        self.surface_temperature_K = float(surface_temperature_K)
        self.surface_heat_flow_mW_m2 = float(surface_heat_flow_mW_m2)
        if (
            not np.isfinite(self.surface_temperature_K)
            or self.surface_temperature_K < 0.0
        ):
            raise ValueError("surface_temperature_K must be finite and non-negative")
        if (
            not np.isfinite(self.surface_heat_flow_mW_m2)
            or self.surface_heat_flow_mW_m2 <= 0.0
        ):
            raise ValueError("surface_heat_flow_mW_m2 must be finite and positive")
        self._name = str(name)
        self._citation = citation or _HASTEROK_CHAPMAN_CITATION
        self._boundaries = np.cumsum(
            [layer.thickness_km for layer in self.layers], dtype=np.float64
        )
        bottom_temperature = self.temperature(
            np.asarray([self._boundaries[-1]], dtype=np.float64)
        )[0]
        if bottom_temperature < 0.0:
            raise ValueError("conductive model predicts a negative temperature")

    @property
    def name(self) -> str:
        """Return the stable temperature-model identifier."""
        return self._name

    @property
    def depth_bounds(self) -> tuple[float, float]:
        """Return the supported depth interval in km."""
        return 0.0, float(self._boundaries[-1])

    @property
    def critical_depths(self) -> tuple[float, ...]:
        """Return internal conductive-layer boundaries in km."""
        return tuple(float(value) for value in self._boundaries[:-1])

    def temperature(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate the layered conductive temperature in K."""
        depth = _validated_depth(depth_km, self.depth_bounds, self.name)
        flat = depth.ravel()
        result = np.full_like(flat, self.surface_temperature_K)
        heat_flow = self.surface_heat_flow_mW_m2 * 1.0e-3
        top_km = 0.0
        active = np.ones_like(flat, dtype=np.bool_)
        for layer, bottom_km in zip(self.layers, self._boundaries, strict=True):
            local_km = np.clip(flat - top_km, 0.0, bottom_km - top_km)
            local_m = local_km * 1000.0
            production = layer.heat_production_uW_m3 * 1.0e-6
            increment = (
                heat_flow * local_m / layer.conductivity_W_mK
                - production * local_m** 2 / (2.0 * layer.conductivity_W_mK)
            )
            result[active] += increment[active]
            fully_below = flat >= bottom_km
            partial = active & ~fully_below
            active[partial] = False
            heat_flow -= production * layer.thickness_km * 1000.0
            top_km = float(bottom_km)
        return result.reshape(depth.shape).astype(np.float64)

    def metadata(self) -> dict[str, Any]:
        """Return conductive parameters and complete bibliographic provenance."""
        return {
            "model": self.name,
            "kind": "steady_layered_conduction",
            "depth_unit": "km",
            "temperature_unit": "K",
            "surface_temperature_K": self.surface_temperature_K,
            "surface_heat_flow_mW_m2": self.surface_heat_flow_mW_m2,
            "layers": [
                {
                    "name": layer.name,
                    "thickness_km": layer.thickness_km,
                    "conductivity_W_mK": layer.conductivity_W_mK,
                    "heat_production_uW_m3": layer.heat_production_uW_m3,
                }
                for layer in self.layers
            ],
            "citation": self._citation,
            "citation_key": _HASTEROK_CHAPMAN_KEY,
            "doi": _HASTEROK_CHAPMAN_REFERENCE.doi,
            "scientific_scope": (
                "Steady one-dimensional conductive approximation; advection, "
                "transient tectonics, and temperature-dependent conductivity "
                "are not included."
            ),
        }


class OceanicHalfSpaceGeotherm:
    """Oceanic half-space cooling temperature profile.

    The model evaluates

    ``T = Ts + (Tm - Ts) erf[z / (2 sqrt(kappa t))]``.

    It describes conductive cooling of an initially hot semi-infinite medium.
    Parsons and Sclater discuss the simple cooling solution and its relation to
    the finite-thickness plate model when interpreting ocean-floor heat flow
    and bathymetry.

    References
    ----------
    Canonical citation key: ``parsons_sclater_1977``.

    Parameters
    ----------
    age_Ma : float
        Oceanic-lithosphere age in million years.
    mantle_temperature_K : float, optional
        Initial and asymptotic mantle temperature in K.
    surface_temperature_K : float, optional
        Surface boundary temperature in K.
    diffusivity_m2_s : float, optional
        Thermal diffusivity in m^2 s^-1.
    max_depth_km : float, optional
        Numerical depth domain exposed by the model.

    Raises
    ------
    ValueError
        If parameters are non-finite or outside their physical domains.
    """

    def __init__(
        self,
        *,
        age_Ma: float,
        mantle_temperature_K: float = 1623.15,
        surface_temperature_K: float = 273.15,
        diffusivity_m2_s: float = 1.0e-6,
        max_depth_km: float = 300.0,
    ) -> None:
        self.age_Ma = float(age_Ma)
        self.mantle_temperature_K = float(mantle_temperature_K)
        self.surface_temperature_K = float(surface_temperature_K)
        self.diffusivity_m2_s = float(diffusivity_m2_s)
        self.max_depth_km = float(max_depth_km)
        if not np.isfinite(self.age_Ma) or self.age_Ma <= 0.0:
            raise ValueError("age_Ma must be finite and positive")
        if (
            not np.isfinite(self.mantle_temperature_K)
            or self.mantle_temperature_K <= 0.0
        ):
            raise ValueError("mantle_temperature_K must be finite and positive")
        if (
            not np.isfinite(self.surface_temperature_K)
            or self.surface_temperature_K < 0.0
            or self.surface_temperature_K >= self.mantle_temperature_K
        ):
            raise ValueError(
                "surface_temperature_K must be non-negative and below mantle temperature"
            )
        if not np.isfinite(self.diffusivity_m2_s) or self.diffusivity_m2_s <= 0.0:
            raise ValueError("diffusivity_m2_s must be finite and positive")
        if not np.isfinite(self.max_depth_km) or self.max_depth_km <= 0.0:
            raise ValueError("max_depth_km must be finite and positive")

    @property
    def name(self) -> str:
        """Return the stable temperature-model identifier."""
        return "oceanic-half-space"

    @property
    def depth_bounds(self) -> tuple[float, float]:
        """Return the exposed numerical depth interval in km."""
        return 0.0, self.max_depth_km

    @property
    def critical_depths(self) -> tuple[float, ...]:
        """Return an empty tuple because the solution is continuous."""
        return ()

    def temperature(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate half-space cooling temperature in K."""
        depth = _validated_depth(depth_km, self.depth_bounds, self.name)
        age_seconds = self.age_Ma * 1.0e6 * 365.25 * 24.0 * 3600.0
        argument = depth * 1000.0 / (2.0 * np.sqrt(self.diffusivity_m2_s * age_seconds))
        return (
            self.surface_temperature_K
            + (self.mantle_temperature_K - self.surface_temperature_K) * erf(argument)
        ).astype(np.float64)

    def metadata(self) -> dict[str, Any]:
        """Return model parameters and complete bibliographic provenance."""
        return {
            "model": self.name,
            "kind": "oceanic_half_space_cooling",
            "depth_unit": "km",
            "temperature_unit": "K",
            "age_Ma": self.age_Ma,
            "surface_temperature_K": self.surface_temperature_K,
            "mantle_temperature_K": self.mantle_temperature_K,
            "diffusivity_m2_s": self.diffusivity_m2_s,
            "citation": _PARSONS_SCLATER_CITATION,
            "citation_key": _PARSONS_SCLATER_KEY,
            "doi": _PARSONS_SCLATER_REFERENCE.doi,
            "scientific_scope": (
                "Semi-infinite conductive cooling approximation; no finite "
                "plate thickness or hydrothermal circulation."
            ),
        }


class OceanicPlateGeotherm:
    """Finite-thickness oceanic plate cooling temperature profile.

    The conductive solution satisfies fixed surface and basal temperatures for
    a plate of thickness ``L``:

    ``(T-Ts)/(Tm-Ts) = z/L + sum[2 sin(n*pi*z/L)/(n*pi) * exp(-kappa*n^2*pi^2*t/L^2)]``.

    At and below the plate base the temperature is set to ``Tm``.  The default
    plate thickness and basal temperature are close to the classical values
    inferred by Parsons and Sclater, while every parameter remains explicit.

    References
    ----------
    Canonical citation key: ``parsons_sclater_1977``.

    Parameters
    ----------
    age_Ma : float
        Oceanic-lithosphere age in million years.
    plate_thickness_km : float, optional
        Fixed thermal plate thickness in km.
    mantle_temperature_K : float, optional
        Basal plate temperature in K.
    surface_temperature_K : float, optional
        Surface boundary temperature in K.
    diffusivity_m2_s : float, optional
        Thermal diffusivity in m^2 s^-1.
    series_terms : int, optional
        Number of Fourier terms retained.
    max_depth_km : float or None, optional
        Exposed domain.  Defaults to the plate thickness.

    Raises
    ------
    ValueError
        If parameters are invalid.
    """

    def __init__(
        self,
        *,
        age_Ma: float,
        plate_thickness_km: float = 125.0,
        mantle_temperature_K: float = 1623.15,
        surface_temperature_K: float = 273.15,
        diffusivity_m2_s: float = 1.0e-6,
        series_terms: int = 200,
        max_depth_km: float | None = None,
    ) -> None:
        self.age_Ma = float(age_Ma)
        self.plate_thickness_km = float(plate_thickness_km)
        self.mantle_temperature_K = float(mantle_temperature_K)
        self.surface_temperature_K = float(surface_temperature_K)
        self.diffusivity_m2_s = float(diffusivity_m2_s)
        self.series_terms = int(series_terms)
        self.max_depth_km = (
            self.plate_thickness_km if max_depth_km is None else float(max_depth_km)
        )
        if not np.isfinite(self.age_Ma) or self.age_Ma <= 0.0:
            raise ValueError("age_Ma must be finite and positive")
        if not np.isfinite(self.plate_thickness_km) or self.plate_thickness_km <= 0.0:
            raise ValueError("plate_thickness_km must be finite and positive")
        if (
            not np.isfinite(self.mantle_temperature_K)
            or self.mantle_temperature_K <= 0.0
        ):
            raise ValueError("mantle_temperature_K must be finite and positive")
        if (
            not np.isfinite(self.surface_temperature_K)
            or self.surface_temperature_K < 0.0
            or self.surface_temperature_K >= self.mantle_temperature_K
        ):
            raise ValueError(
                "surface_temperature_K must be non-negative and below mantle temperature"
            )
        if not np.isfinite(self.diffusivity_m2_s) or self.diffusivity_m2_s <= 0.0:
            raise ValueError("diffusivity_m2_s must be finite and positive")
        if self.series_terms < 1:
            raise ValueError("series_terms must be positive")
        if (
            not np.isfinite(self.max_depth_km)
            or self.max_depth_km < self.plate_thickness_km
        ):
            raise ValueError("max_depth_km must be at least plate_thickness_km")

    @property
    def name(self) -> str:
        """Return the stable temperature-model identifier."""
        return "oceanic-plate"

    @property
    def depth_bounds(self) -> tuple[float, float]:
        """Return the exposed depth interval in km."""
        return 0.0, self.max_depth_km

    @property
    def critical_depths(self) -> tuple[float, ...]:
        """Return the finite thermal-plate base."""
        return (self.plate_thickness_km,)

    def temperature(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate finite-plate cooling temperature in K."""
        depth = _validated_depth(depth_km, self.depth_bounds, self.name)
        flat = depth.ravel()
        normalized = np.clip(flat / self.plate_thickness_km, 0.0, 1.0)
        age_seconds = self.age_Ma * 1.0e6 * 365.25 * 24.0 * 3600.0
        plate_m = self.plate_thickness_km * 1000.0
        n = np.arange(1, self.series_terms + 1, dtype=np.float64)[:, None]
        z = normalized[None, :]
        decay = np.exp(
            -self.diffusivity_m2_s * n**2 * np.pi**2 * age_seconds / plate_m**2
        )
        series = np.sum(
            2.0 * np.sin(n * np.pi * z) * decay / (n * np.pi),
            axis=0,
        )
        theta = normalized + series
        theta[flat >= self.plate_thickness_km] = 1.0
        theta[flat == 0.0] = 0.0
        result = (
            self.surface_temperature_K
            + (self.mantle_temperature_K - self.surface_temperature_K) * theta
        )
        return result.reshape(depth.shape).astype(np.float64)

    def metadata(self) -> dict[str, Any]:
        """Return model parameters and complete bibliographic provenance."""
        return {
            "model": self.name,
            "kind": "oceanic_finite_plate_cooling",
            "depth_unit": "km",
            "temperature_unit": "K",
            "age_Ma": self.age_Ma,
            "plate_thickness_km": self.plate_thickness_km,
            "surface_temperature_K": self.surface_temperature_K,
            "mantle_temperature_K": self.mantle_temperature_K,
            "diffusivity_m2_s": self.diffusivity_m2_s,
            "series_terms": self.series_terms,
            "citation": _PARSONS_SCLATER_CITATION,
            "citation_key": _PARSONS_SCLATER_KEY,
            "doi": _PARSONS_SCLATER_REFERENCE.doi,
            "scientific_scope": (
                "One-dimensional fixed-thickness conductive plate; hydrothermal "
                "circulation and lateral variations are not included."
            ),
        }


class Katsura2022MantleAdiabat:
    """Published-constraint reconstruction of the Katsura (2022) mantle adiabat.

    Quantas reconstructs the most-probable dry-pyrolite profile with monotonic
    cubic Hermite segments constrained by the temperatures and local gradients
    reported in the article.  The phase-boundary temperatures are represented
    on their shallow and deep sides at 410, 520, and 660 km.  At a depth exactly
    equal to a transition, the deeper assemblage is selected.  This class is a
    deterministic reconstruction of the published profile, not a re-execution
    of the author's full Monte Carlo MATLAB workflow.

    References
    ----------
    Canonical citation keys: ``katsura_2022`` and
    ``katsura_software_2022``.

    Parameters
    ----------
    transition_epsilon_km : float, optional
        Half-width used only to insert paired sampling points around abrupt
        phase-boundary temperature changes in generated profiles.

    Raises
    ------
    ValueError
        If ``transition_epsilon_km`` is not finite and positive.
    """

    _transitions = (410.0, 520.0, 660.0)

    def __init__(self, *, transition_epsilon_km: float = 1.0e-3) -> None:
        self.transition_epsilon_km = float(transition_epsilon_km)
        if (
            not np.isfinite(self.transition_epsilon_km)
            or self.transition_epsilon_km <= 0.0
            or self.transition_epsilon_km >= 1.0
        ):
            raise ValueError("transition_epsilon_km must lie in (0, 1) km")
        self._segments = (
            CubicHermiteSpline([50.0, 410.0], [1646.0, 1799.0], [0.54, 0.36]),
            CubicHermiteSpline([410.0, 520.0], [1860.0, 1899.0], [0.36, 0.35]),
            CubicHermiteSpline([520.0, 660.0], [1942.0, 1994.0], [0.39, 0.37]),
            CubicHermiteSpline([660.0, 2400.0], [1960.0, 2490.0], [0.41, 0.26]),
            CubicHermiteSpline([2400.0, 2800.0], [2490.0, 2587.0], [0.26, 0.23]),
        )

    @property
    def name(self) -> str:
        """Return the stable temperature-model identifier."""
        return "katsura-2022"

    @property
    def depth_bounds(self) -> tuple[float, float]:
        """Return the published-profile reconstruction domain in km."""
        return 50.0, 2800.0

    @property
    def critical_depths(self) -> tuple[float, ...]:
        """Return paired sampling points around phase boundaries and 2400 km."""
        values: list[float] = []
        for depth in self._transitions:
            values.extend(
                [depth - self.transition_epsilon_km, depth + self.transition_epsilon_km]
            )
        values.append(2400.0)
        return tuple(values)

    def temperature(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate the reconstructed mantle adiabat in K."""
        depth = _validated_depth(depth_km, self.depth_bounds, self.name)
        flat = depth.ravel()
        result = np.empty_like(flat)
        masks = (
            flat < 410.0,
            (flat >= 410.0) & (flat < 520.0),
            (flat >= 520.0) & (flat < 660.0),
            (flat >= 660.0) & (flat < 2400.0),
            flat >= 2400.0,
        )
        for mask, segment in zip(masks, self._segments, strict=True):
            result[mask] = segment(flat[mask])
        return result.reshape(depth.shape).astype(np.float64)

    def metadata(self) -> dict[str, Any]:
        """Return reconstruction constraints and complete provenance."""
        return {
            "model": self.name,
            "kind": "published_constraint_mantle_adiabat",
            "depth_unit": "km",
            "temperature_unit": "K",
            "transition_side_at_exact_depth": "deep",
            "transition_epsilon_km": self.transition_epsilon_km,
            "temperature_constraints_K": {
                "50_km": 1646.0,
                "410_km_shallow": 1799.0,
                "410_km_deep": 1860.0,
                "520_km_shallow": 1899.0,
                "520_km_deep": 1942.0,
                "660_km_shallow": 1994.0,
                "660_km_deep": 1960.0,
                "2400_km": 2490.0,
                "2800_km": 2587.0,
            },
            "citation": _KATSURA_CITATION,
            "software_citation": _KATSURA_SOFTWARE_CITATION,
            "citation_key": _KATSURA_KEY,
            "software_citation_key": _KATSURA_SOFTWARE_KEY,
            "doi": _KATSURA_REFERENCE.doi,
            "software_doi": _KATSURA_SOFTWARE_REFERENCE.doi,
            "scientific_scope": (
                "Most-probable dry-pyrolite adiabatic mantle profile; not a "
                "lithospheric conductive geotherm and not a core profile."
            ),
            "reconstruction_note": (
                "Cubic Hermite segments reproduce published temperatures and "
                "reported endpoint gradients; uncertainty envelopes are not "
                "propagated by this deterministic preset."
            ),
        }


class LinearThermalBoundaryLayer:
    """User-parameterized thermal boundary layer between two depths.

    Parameters
    ----------
    depth_top_km, depth_bottom_km : float
        Boundary-layer depth limits in km.
    temperature_top_K, temperature_bottom_K : float
        Endpoint temperatures in K.
    exponent : float, optional
        Shape exponent applied to normalized depth.  One gives a linear ramp.
    name : str, optional
        Stable model identifier.
    citation : str or None, optional
        Complete provenance for the chosen boundary conditions.

    Raises
    ------
    ValueError
        If bounds, temperatures, or exponent are invalid.
    """

    def __init__(
        self,
        *,
        depth_top_km: float,
        depth_bottom_km: float,
        temperature_top_K: float,
        temperature_bottom_K: float,
        exponent: float = 1.0,
        name: str = "linear-thermal-boundary-layer",
        citation: str | None = None,
    ) -> None:
        self.depth_top_km = float(depth_top_km)
        self.depth_bottom_km = float(depth_bottom_km)
        self.temperature_top_K = float(temperature_top_K)
        self.temperature_bottom_K = float(temperature_bottom_K)
        self.exponent = float(exponent)
        self._name = str(name)
        self._citation = citation
        if (
            not np.isfinite(self.depth_top_km)
            or not np.isfinite(self.depth_bottom_km)
            or self.depth_top_km < 0.0
            or self.depth_bottom_km <= self.depth_top_km
        ):
            raise ValueError("boundary-layer depths must be finite and increasing")
        if (
            not np.isfinite(self.temperature_top_K)
            or not np.isfinite(self.temperature_bottom_K)
            or self.temperature_top_K < 0.0
            or self.temperature_bottom_K < 0.0
        ):
            raise ValueError("boundary-layer temperatures must be non-negative")
        if not np.isfinite(self.exponent) or self.exponent <= 0.0:
            raise ValueError("boundary-layer exponent must be finite and positive")

    @property
    def name(self) -> str:
        """Return the stable temperature-model identifier."""
        return self._name

    @property
    def depth_bounds(self) -> tuple[float, float]:
        """Return the supported boundary-layer depth interval."""
        return self.depth_top_km, self.depth_bottom_km

    @property
    def critical_depths(self) -> tuple[float, ...]:
        """Return the boundary-layer endpoints."""
        return (self.depth_top_km, self.depth_bottom_km)

    def temperature(self, depth_km: NDArray[np.float64]) -> FloatArray:
        """Evaluate the parameterized boundary-layer temperature in K."""
        depth = _validated_depth(depth_km, self.depth_bounds, self.name)
        fraction = (depth - self.depth_top_km) / (
            self.depth_bottom_km - self.depth_top_km
        )
        return (
            self.temperature_top_K
            + (self.temperature_bottom_K - self.temperature_top_K)
            * fraction**self.exponent
        ).astype(np.float64)

    def metadata(self) -> dict[str, Any]:
        """Return endpoint parameters and user-supplied provenance."""
        return {
            "model": self.name,
            "kind": "parameterized_thermal_boundary_layer",
            "depth_top_km": self.depth_top_km,
            "depth_bottom_km": self.depth_bottom_km,
            "temperature_top_K": self.temperature_top_K,
            "temperature_bottom_K": self.temperature_bottom_K,
            "exponent": self.exponent,
            "citation": self._citation or "User-specified boundary-layer model.",
        }


__all__ = [
    "ConductiveLayer",
    "ContinentalConductiveGeotherm",
    "Katsura2022MantleAdiabat",
    "LinearThermalBoundaryLayer",
    "OceanicHalfSpaceGeotherm",
    "OceanicPlateGeotherm",
]
