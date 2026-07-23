# -*- coding: utf-8 -*-

"""Volume minimization utilities for quasi-harmonic calculations.

The functions in this module evaluate equilibrium volumes from
volume-dependent free-energy data using either polynomial or equation-of-state
fits.  The routines return structured results containing the optimized volume,
fit diagnostics and flags describing whether the minimum lies inside the
sampled volume interval.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping, Sequence

import numpy as np

from quantas.core.math.polynomials import (
    FittedPolynomial,
    fit_polynomial,
    fit_polynomial_result,
)
from quantas.core.physics.units import energy_to_pressure, pressure_to_energy
from quantas.core.physics.eos import FittedEnergyEOS
from quantas.core.physics.eos import EnergyEOS
from quantas.core.math.fitting import FitQuality, FitResult, FitStatus, validate_xy

ArrayLike = np.ndarray | Sequence[float]


class MinimumStatus(str, Enum):
    """Execution status of a volume minimization.

    Attributes
    ----------
    SUCCESS
        The minimization returned a usable equilibrium volume.
    FAILED
        The minimization did not return a usable equilibrium volume.
    INVALID_INPUT
        The supplied data are not suitable for minimization.
    """

    SUCCESS = "success"
    FAILED = "failed"
    INVALID_INPUT = "invalid_input"


class VolumeRangeStatus(str, Enum):
    """Location of a minimum relative to the sampled volume interval.

    Attributes
    ----------
    INSIDE
        The minimum lies inside the sampled interval.
    NEAR_BOUNDARY
        The minimum lies inside the interval but close to one boundary.
    OUTSIDE
        The minimum lies outside the sampled interval.
    """

    INSIDE = "inside"
    NEAR_BOUNDARY = "near_boundary"
    OUTSIDE = "outside"


@dataclass(slots=True)
class VolumeMinimumResult:
    """Result of a free-energy volume minimization.

    Parameters
    ----------
    success : bool
        Whether the minimization returned a usable result.
    status : MinimumStatus
        Execution status of the minimization.
    method : str
        Minimization method, such as ``"polynomial"`` or ``"eos"``.
    volume : float or None
        Equilibrium volume at the requested pressure and temperature.
    bulk_modulus : float or None
        Isothermal bulk modulus, when available.
    bulk_modulus_derivative : float or None
        First pressure derivative of the bulk modulus, when available.
    sigma_volume : float or None
        One-sigma uncertainty of the equilibrium volume, when available.
    sigma_bulk_modulus : float or None
        One-sigma uncertainty of the isothermal bulk modulus, when available.
    sigma_bulk_modulus_derivative : float or None
        One-sigma uncertainty of the pressure derivative of the bulk modulus,
        when available.
    message : str
        Human-readable status message.
    range_status : VolumeRangeStatus or None
        Location of the minimum relative to the sampled volumes.
    fit : FitResult or None
        Fit diagnostics associated with the minimization.
    warnings : list of str
        Non-fatal diagnostic messages.
    metadata : dict
        Additional caller-defined diagnostic information.
    """

    success: bool
    status: MinimumStatus
    method: str
    volume: float | None = None
    bulk_modulus: float | None = None
    bulk_modulus_derivative: float | None = None
    sigma_volume: float | None = None
    sigma_bulk_modulus: float | None = None
    sigma_bulk_modulus_derivative: float | None = None
    message: str = ""
    range_status: VolumeRangeStatus | None = None
    fit: FitResult | None = None
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def failed(
        cls,
        message: str,
        *,
        method: str,
        status: MinimumStatus = MinimumStatus.FAILED,
        fit: FitResult | None = None,
        metadata: Mapping[str, Any] | None = None,
    ) -> VolumeMinimumResult:
        """Create a failed minimization result.

        Parameters
        ----------
        message : str
            Explanation of the failure.
        method : str
            Minimization method.
        status : MinimumStatus, optional
            Failure category.
        fit : FitResult, optional
            Fit diagnostics, when a fit was attempted.
        metadata : mapping, optional
            Additional diagnostic information.

        Returns
        -------
        VolumeMinimumResult
            Failed minimization result.
        """
        return cls(
            success=False,
            status=status,
            method=method,
            message=message,
            fit=fit,
            metadata=dict(metadata or {}),
        )

    def as_dict(self) -> dict[str, Any]:
        """Return the minimization result as a serializable dictionary.

        Returns
        -------
        dict
            Dictionary representation of the minimization result.
        """
        return {
            "success": self.success,
            "status": self.status.value,
            "method": self.method,
            "volume": self.volume,
            "bulk_modulus": self.bulk_modulus,
            "bulk_modulus_derivative": self.bulk_modulus_derivative,
            "sigma_volume": self.sigma_volume,
            "sigma_bulk_modulus": self.sigma_bulk_modulus,
            "sigma_bulk_modulus_derivative": self.sigma_bulk_modulus_derivative,
            "message": self.message,
            "range_status": None
            if self.range_status is None
            else self.range_status.value,
            "fit": None if self.fit is None else self.fit.as_dict(),
            "warnings": list(self.warnings),
            "metadata": dict(self.metadata),
        }


@dataclass(slots=True)
class EOSModelFitResult:
    """Result of fitting one energy EOS to a free-energy curve.

    Parameters
    ----------
    success : bool
        Whether a usable fitted EOS model was created.
    eos : str
        Canonical equation-of-state name.
    fit : FitResult
        Structured diagnostics from the energy-volume fit.
    model : FittedEnergyEOS or None
        Pressure-volume model built from the fitted parameters.
    message : str
        Human-readable status message.
    warnings : list of str
        Non-fatal diagnostics associated with model construction.
    metadata : dict
        Unit conversions and model details.
    """

    success: bool
    eos: str
    fit: FitResult
    model: FittedEnergyEOS | None = None
    message: str = ""
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        """Return a serializable representation of the fitted EOS model.

        Returns
        -------
        dict
            Fit status, diagnostics and model metadata.
        """
        return {
            "success": self.success,
            "eos": self.eos,
            "fit": self.fit.as_dict(),
            "message": self.message,
            "warnings": list(self.warnings),
            "metadata": dict(self.metadata),
        }


@dataclass(slots=True)
class PolynomialModelFitResult:
    """Result of fitting one polynomial free-energy model.

    Parameters
    ----------
    success : bool
        Whether a usable polynomial model was created.
    fit : FitResult
        Structured diagnostics from the scaled polynomial fit.
    model : FittedPolynomialFreeEnergy or None
        Reusable polynomial free-energy model.
    message : str
        Human-readable status message.
    warnings : list of str
        Non-fatal diagnostics associated with model construction.
    metadata : dict
        Scaling and model details.
    """

    success: bool
    fit: FitResult
    model: "FittedPolynomialFreeEnergy | None" = None
    message: str = ""
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        """Return a serializable representation of the fitted model.

        Returns
        -------
        dict
            Fit status, diagnostics and scaling metadata.
        """
        return {
            "success": self.success,
            "fit": self.fit.as_dict(),
            "message": self.message,
            "warnings": list(self.warnings),
            "metadata": dict(self.metadata),
        }


class FittedPolynomialFreeEnergy(FittedPolynomial):
    """Represent a Helmholtz free-energy curve with a scaled polynomial.

    The fitted coordinate is

    .. math::

        x = \frac{V - V_c}{V_s},

    where ``V_c`` and ``V_s`` are the volume center and scale. Derivatives are
    evaluated with respect to the physical volume ``V``.

    Parameters
    ----------
    parameters : array-like
        Polynomial coefficients in ascending order in the scaled coordinate.
    center : float
        Volume used as the coordinate origin.
    scale : float
        Positive volume scale.
    sampled_volumes : array-like
        Volume values used to determine the fit.
    """

    def __init__(
        self,
        parameters: ArrayLike,
        *,
        center: float,
        scale: float,
        sampled_volumes: ArrayLike,
    ) -> None:
        super().__init__(
            parameters,
            center=center,
            scale=scale,
            sampled_coordinates=sampled_volumes,
        )
        assert self.sampled_coordinates is not None
        self.sampled_volumes = self.sampled_coordinates.copy()

    def scaled_volume(self, volume: ArrayLike | float) -> np.ndarray:
        """Transform volumes to the scaled polynomial coordinate.

        Parameters
        ----------
        volume : scalar or array-like
            Physical volume values.

        Returns
        -------
        ndarray
            Scaled volume coordinates.
        """
        return self.scaled_coordinate(volume)

    def free_energy(self, volume: ArrayLike | float) -> np.ndarray:
        """Evaluate the fitted Helmholtz free energy.

        Parameters
        ----------
        volume : scalar or array-like
            Physical volume values.

        Returns
        -------
        ndarray
            Fitted free-energy values.
        """
        return self.evaluate(volume)

    def stationary_minima(
        self,
        pressure_energy_density: float = 0.0,
    ) -> np.ndarray:
        """Return positive-volume minima of ``F(V) + P V``.

        Parameters
        ----------
        pressure_energy_density : float, optional
            Pressure expressed in the native energy-density scale.

        Returns
        -------
        ndarray
            Positive-volume local minima.
        """
        minima = self.local_minima(pressure_energy_density)
        return np.asarray(
            minima[np.isfinite(minima) & (minima > 0.0)], dtype=np.float64
        )

    def equilibrium_volume(self, pressure_energy_density: float = 0.0) -> float:
        """Return the minimum connected to the sampled free-energy basin.

        The sampled volume with the smallest value of ``F(V) + P V`` defines
        the reference basin. If several stationary minima are present, the
        closest minimum is selected and objective values break exact distance
        ties.

        Parameters
        ----------
        pressure_energy_density : float, optional
            Pressure expressed in the native energy-density scale.

        Returns
        -------
        float
            Equilibrium volume.

        Raises
        ------
        RuntimeError
            If no physical minimum is available.
        """
        minima = self.stationary_minima(pressure_energy_density)
        if minima.size == 0:
            raise RuntimeError(
                "the polynomial curve did not produce a physical minimum at the requested pressure"
            )
        sampled_objective = self.free_energy(self.sampled_volumes) + (
            float(pressure_energy_density) * self.sampled_volumes
        )
        reference = float(self.sampled_volumes[int(np.argmin(sampled_objective))])
        distance = np.abs(minima - reference)
        closest = np.flatnonzero(
            np.isclose(distance, np.min(distance), rtol=1.0e-12, atol=1.0e-12)
        )
        if closest.size == 1:
            return float(minima[int(closest[0])])
        objective = self.free_energy(minima[closest]) + (
            float(pressure_energy_density) * minima[closest]
        )
        return float(minima[int(closest[int(np.argmin(objective))])])

    def analytic_bulk_properties(self, volume: float) -> tuple[float, float]:
        r"""Return ``K_T`` and ``K'_T`` from polynomial derivatives.

        The properties are calculated from

        .. math::

            K_T = V \frac{\partial^2 F}{\partial V^2},

        and

        .. math::

            K'_T = -1 - V
            \frac{\partial^3 F / \partial V^3}
                  {\partial^2 F / \partial V^2}.

        Parameters
        ----------
        volume : float
            Volume at which the properties are evaluated.

        Returns
        -------
        tuple of float
            Bulk modulus in native energy-density units and dimensionless
            pressure derivative.

        Raises
        ------
        RuntimeError
            If the fitted curvature is not positive and finite.
        """
        second = float(self.derivative(volume, 2))
        third = float(self.derivative(volume, 3))
        if not np.isfinite(second) or second <= 0.0:
            raise RuntimeError(
                "the polynomial curvature is not positive at the minimum"
            )
        bulk_modulus = float(volume * second)
        derivative = float(-1.0 - volume * third / second)
        return bulk_modulus, derivative


def fine_grid(npoints: int = 5, separation: float = 0.05) -> tuple[np.ndarray, int]:
    """Return relative volume factors around a central point.

    Parameters
    ----------
    npoints : int, optional
        Number of grid points.  The value must be an odd positive integer.
    separation : float, optional
        Separation between adjacent points, expressed as a percentage.

    Returns
    -------
    tuple of ndarray and int
        Relative factors and index of the central point.

    Raises
    ------
    ValueError
        If ``npoints`` is not a positive odd integer or ``separation`` is not
        positive.
    """
    if npoints <= 0 or npoints % 2 == 0:
        raise ValueError("npoints must be a positive odd integer")
    if separation <= 0.0:
        raise ValueError("separation must be positive")
    spacing = separation / 100.0
    half = npoints // 2
    expansion = np.arange(-half, half + 1, dtype=np.float64) * spacing
    return expansion, half


def target_pressure_energy_density(
    pressure: float,
    *,
    energy_unit: str,
    volume_unit: str,
    pressure_unit: str,
) -> float:
    """Convert a pressure to the energy-density scale used by free energies.

    Parameters
    ----------
    pressure : float
        Pressure value.
    energy_unit : str
        Energy unit of the free-energy data.
    volume_unit : str
        Length unit defining the volume unit.
    pressure_unit : str
        Unit of the pressure value.

    Returns
    -------
    float
        Pressure converted to energy per volume.
    """
    return float(pressure_to_energy(pressure, energy_unit, volume_unit, pressure_unit))


def pressure_shifted_free_energy(
    free_energy: ArrayLike,
    volume: ArrayLike,
    pressure_energy_density: float,
) -> np.ndarray:
    """Add the pressure-volume term to a free-energy curve.

    Parameters
    ----------
    free_energy : array-like
        Helmholtz free-energy values.
    volume : array-like
        Volume values.
    pressure_energy_density : float
        Pressure in the energy-density scale of the free-energy data.

    Returns
    -------
    ndarray
        Pressure-shifted free energy ``F + P V``.

    Raises
    ------
    ValueError
        If the arrays are not suitable for element-wise operations.
    """
    volume_array, free_energy_array = validate_xy(volume, free_energy)
    return free_energy_array + float(pressure_energy_density) * volume_array


def classify_volume(
    volume: float,
    sampled_volumes: ArrayLike,
    *,
    boundary_fraction: float = 0.05,
) -> VolumeRangeStatus:
    """Classify a volume relative to a sampled volume interval.

    Parameters
    ----------
    volume : float
        Candidate equilibrium volume.
    sampled_volumes : array-like
        Volumes used by the fit.
    boundary_fraction : float, optional
        Fraction of the sampled interval treated as a boundary region.

    Returns
    -------
    VolumeRangeStatus
        Location of ``volume`` relative to the sampled interval.

    Raises
    ------
    ValueError
        If the sampled volumes are not finite or the interval is degenerate.
    """
    values = np.asarray(sampled_volumes, dtype=np.float64)
    if values.ndim != 1 or values.size < 2 or not np.all(np.isfinite(values)):
        raise ValueError("sampled volumes must be a finite one-dimensional array")
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    width = vmax - vmin
    if width <= 0.0:
        raise ValueError("sampled volume interval must be non-degenerate")
    if volume < vmin or volume > vmax:
        return VolumeRangeStatus.OUTSIDE
    boundary = boundary_fraction * width
    if (volume - vmin) <= boundary or (vmax - volume) <= boundary:
        return VolumeRangeStatus.NEAR_BOUNDARY
    return VolumeRangeStatus.INSIDE


def polynomial_fit(
    volume: ArrayLike,
    free_energy: ArrayLike,
    degree: int,
) -> FitResult:
    """Fit a polynomial to a free-energy curve.

    Parameters
    ----------
    volume : array-like
        Volume values.
    free_energy : array-like
        Free-energy values.
    degree : int
        Polynomial degree.

    Returns
    -------
    FitResult
        Polynomial coefficients in ascending order and diagnostics.
    """
    return fit_polynomial_result(
        volume,
        free_energy,
        degree,
        scale_coordinate=False,
    )


def fit_polynomial_free_energy_model(
    volume: ArrayLike,
    free_energy: ArrayLike,
    *,
    degree: int = 4,
) -> PolynomialModelFitResult:
    """Fit one scaled polynomial model to a free-energy curve.

    The coordinate transformation

    .. math::

        x = \frac{V - V_c}{V_s}

    improves the conditioning of the polynomial fit over narrow volume
    intervals.

    Parameters
    ----------
    volume : array-like
        Sampled volume values.
    free_energy : array-like
        Helmholtz free-energy values.
    degree : int, optional
        Polynomial degree.

    Returns
    -------
    PolynomialModelFitResult
        Fit diagnostics and reusable free-energy model.
    """
    fit, generic_model = fit_polynomial(
        volume,
        free_energy,
        degree,
        scale_coordinate=True,
    )
    if not fit.success or generic_model is None or fit.parameters is None:
        return PolynomialModelFitResult(False, fit, message=fit.message)

    assert generic_model.sampled_coordinates is not None
    try:
        model = FittedPolynomialFreeEnergy(
            fit.parameters,
            center=generic_model.center,
            scale=generic_model.scale,
            sampled_volumes=generic_model.sampled_coordinates,
        )
    except ValueError as exc:
        return PolynomialModelFitResult(
            False,
            fit,
            message=f"could not construct the polynomial free-energy model: {exc}",
        )

    vmin = float(np.min(model.sampled_volumes))
    vmax = float(np.max(model.sampled_volumes))
    metadata = {
        "degree": int(degree),
        "volume_center": model.center,
        "volume_scale": model.scale,
        "sampled_volume_range": [vmin, vmax],
    }
    fit.metadata.update(
        {
            "volume_center": model.center,
            "volume_scale": model.scale,
            "sampled_volume_min": vmin,
            "sampled_volume_max": vmax,
            "coordinate": "x=(V-center)/scale",
        }
    )
    return PolynomialModelFitResult(
        True,
        fit,
        model=model,
        message=fit.message,
        warnings=list(fit.warnings),
        metadata=metadata,
    )


def _local_polynomial_bulk_properties(
    volume: ArrayLike,
    free_energy: ArrayLike,
    *,
    index: int,
    degree: int,
) -> tuple[float, float, PolynomialModelFitResult, FitResult]:
    """Return local ``K_T`` and ``K'_T`` from a small volume grid.

    Parameters
    ----------
    volume : array-like
        Local volume grid.
    free_energy : array-like
        Helmholtz free energies evaluated on the local grid.
    index : int
        Index of the central volume.
    degree : int
        Polynomial degree used for the local free-energy fit.

    Returns
    -------
    tuple
        Bulk modulus, pressure derivative, local free-energy fit and
        quadratic ``K(P)`` fit.

    Raises
    ------
    RuntimeError
        If either local fit fails or the central curvature is not positive.
    IndexError
        If ``index`` is outside the local grid.
    """
    volume_array, free_energy_array = validate_xy(volume, free_energy)
    if index < 0 or index >= volume_array.size:
        raise IndexError("index is outside the local volume grid")
    effective_degree = min(int(degree), int(volume_array.size) - 1)
    local_fit = fit_polynomial_free_energy_model(
        volume_array,
        free_energy_array,
        degree=effective_degree,
    )
    if not local_fit.success or local_fit.model is None:
        raise RuntimeError(local_fit.message or "local polynomial fit failed")

    model = local_fit.model
    pressure_values = -model.derivative(volume_array, 1)
    bulk_values = volume_array * model.derivative(volume_array, 2)
    if not np.all(np.isfinite(pressure_values)) or not np.all(np.isfinite(bulk_values)):
        raise RuntimeError("local polynomial derivatives are not finite")
    if float(bulk_values[index]) <= 0.0:
        raise RuntimeError("local polynomial curvature is not positive")

    pressure_degree = min(2, int(volume_array.size) - 1)
    pressure_fit = polynomial_fit(pressure_values, bulk_values, pressure_degree)
    if not pressure_fit.success or pressure_fit.parameters is None:
        raise RuntimeError(pressure_fit.message or "local K(P) fit failed")
    derivative_parameters = np.polynomial.polynomial.polyder(
        pressure_fit.parameters,
        1,
    )
    kp = float(
        np.polynomial.polynomial.polyval(
            pressure_values[index],
            derivative_parameters,
        )
    )
    return float(bulk_values[index]), kp, local_fit, pressure_fit


def evaluate_fitted_polynomial_at_pressure(
    fitted: PolynomialModelFitResult,
    pressure_energy_density: float,
    *,
    derivative_method: str = "local_grid",
    local_free_energy: Callable[[np.ndarray], np.ndarray] | None = None,
    local_grid_points: int = 5,
    local_grid_separation: float = 0.05,
    local_degree: int = 3,
) -> VolumeMinimumResult:
    """Evaluate one fitted polynomial model at a target pressure.

    Parameters
    ----------
    fitted : PolynomialModelFitResult
        Reusable polynomial free-energy model.
    pressure_energy_density : float
        Target pressure in the energy-density scale of the fit.
    derivative_method : {"local_grid", "analytic"}, optional
        Method used to obtain ``K_T`` and ``K'_T``.
    local_free_energy : callable or None, optional
        Function receiving a local volume array and returning Helmholtz free
        energies. When omitted, the fitted global polynomial is evaluated.
    local_grid_points : int, optional
        Number of points in the local volume grid.
    local_grid_separation : float, optional
        Adjacent local-grid spacing as a percentage of the central volume.
    local_degree : int, optional
        Polynomial degree used for the local free-energy fit.

    Returns
    -------
    VolumeMinimumResult
        Equilibrium volume and thermoelastic properties.
    """
    if not fitted.success or fitted.model is None:
        return VolumeMinimumResult.failed(
            fitted.message or "the polynomial free-energy fit failed",
            method="polynomial",
            fit=fitted.fit,
            metadata={"failure_stage": "polynomial_fit"},
        )
    if derivative_method not in {"local_grid", "analytic"}:
        return VolumeMinimumResult.failed(
            f"unsupported polynomial derivative method: {derivative_method}",
            method="polynomial",
            status=MinimumStatus.INVALID_INPUT,
            fit=fitted.fit,
        )

    model = fitted.model
    try:
        volume_minimum = model.equilibrium_volume(pressure_energy_density)
    except RuntimeError as exc:
        return VolumeMinimumResult.failed(
            str(exc),
            method="polynomial",
            fit=fitted.fit,
            metadata={
                "pressure_energy_density": float(pressure_energy_density),
                "failure_stage": "polynomial_minimum",
            },
        )

    try:
        range_status = classify_volume(volume_minimum, model.sampled_volumes)
    except ValueError:
        range_status = None
    warnings_: list[str] = []
    if range_status is VolumeRangeStatus.NEAR_BOUNDARY:
        warnings_.append(
            "the equilibrium volume is close to the sampled volume boundary"
        )
    elif range_status is VolumeRangeStatus.OUTSIDE:
        warnings_.append(
            "the equilibrium volume is outside the sampled volume interval"
        )

    metadata: dict[str, Any] = {
        "pressure_energy_density": float(pressure_energy_density),
        "derivative_method": derivative_method,
        "global_fit": {
            "success": bool(fitted.success),
            "quality": fitted.fit.quality.value,
            "degree": fitted.metadata.get("degree"),
            "volume_center": fitted.metadata.get("volume_center"),
            "volume_scale": fitted.metadata.get("volume_scale"),
            "sampled_volume_range": fitted.metadata.get("sampled_volume_range"),
        },
    }
    try:
        analytic_kt, analytic_kp = model.analytic_bulk_properties(volume_minimum)
        metadata["analytic_bulk_modulus"] = analytic_kt
        metadata["analytic_bulk_modulus_derivative"] = analytic_kp
        if derivative_method == "analytic":
            bulk_modulus = analytic_kt
            bulk_modulus_derivative = analytic_kp
        else:
            expansion, center_index = fine_grid(
                npoints=local_grid_points,
                separation=local_grid_separation,
            )
            local_volumes = volume_minimum * (1.0 + expansion)
            if local_free_energy is None:
                local_energies = model.free_energy(local_volumes)
            else:
                local_energies = np.asarray(
                    local_free_energy(local_volumes),
                    dtype=np.float64,
                )
            if local_energies.shape != local_volumes.shape:
                raise RuntimeError(
                    "the local free-energy evaluator returned an incompatible shape"
                )
            bulk_modulus, bulk_modulus_derivative, local_fit, pressure_fit = (
                _local_polynomial_bulk_properties(
                    local_volumes,
                    local_energies,
                    index=center_index,
                    degree=local_degree,
                )
            )
            if local_fit.model is None:
                raise RuntimeError("local polynomial model is unavailable")
            local_pressure = float(-local_fit.model.derivative(volume_minimum, 1))
            metadata["local_grid"] = {
                "npoints": int(local_grid_points),
                "separation_percent": float(local_grid_separation),
                "degree_requested": int(local_degree),
                "degree_effective": int(min(local_degree, local_grid_points - 1)),
                "volumes": local_volumes.tolist(),
                "free_energy": local_energies.tolist(),
                "central_pressure_energy_density": local_pressure,
                "target_pressure_energy_density": float(pressure_energy_density),
                "pressure_residual_energy_density": float(
                    local_pressure - pressure_energy_density
                ),
                "free_energy_fit": local_fit.as_dict(),
                "bulk_pressure_fit": pressure_fit.as_dict(),
            }
    except (ValueError, RuntimeError, FloatingPointError, ZeroDivisionError) as exc:
        return VolumeMinimumResult.failed(
            f"polynomial thermoelastic derivatives failed: {exc}",
            method="polynomial",
            fit=fitted.fit,
            metadata={**metadata, "failure_stage": "polynomial_derivatives"},
        )

    return VolumeMinimumResult(
        success=True,
        status=MinimumStatus.SUCCESS,
        method="polynomial",
        volume=volume_minimum,
        bulk_modulus=float(bulk_modulus),
        bulk_modulus_derivative=float(bulk_modulus_derivative),
        message="polynomial state evaluated",
        range_status=range_status,
        fit=None,
        warnings=warnings_,
        metadata=metadata,
    )


def polynomial_stationary_points(
    parameters: ArrayLike, *, addendum: float = 0.0
) -> np.ndarray:
    """Return minima of a polynomial after adding a constant derivative term.

    Parameters
    ----------
    parameters : array-like
        Polynomial coefficients in ascending order.
    addendum : float, optional
        Constant added to the first derivative.  For QHA minimization this is
        the pressure expressed as energy per volume.

    Returns
    -------
    ndarray
        Real stationary points with positive second derivative.
    """
    polynomial = np.polynomial.polynomial.Polynomial(
        np.asarray(parameters, dtype=np.float64)
    )
    derivative = polynomial.deriv(1)
    derivative = derivative + float(addendum)
    roots = derivative.roots()
    real_roots = np.real(roots[np.isclose(np.imag(roots), 0.0)])
    if real_roots.size == 0:
        return np.asarray([], dtype=np.float64)
    curvature = polynomial.deriv(2)(real_roots)
    return np.asarray(real_roots[curvature > 0.0], dtype=np.float64)


def minimize_polynomial(
    volume: ArrayLike,
    free_energy: ArrayLike,
    *,
    pressure_energy_density: float = 0.0,
    degree: int = 4,
    parameters: ArrayLike | None = None,
    derivative_method: str = "analytic",
    local_grid_points: int = 5,
    local_grid_separation: float = 0.05,
    local_degree: int | None = None,
) -> VolumeMinimumResult:
    """Minimize a pressure-shifted polynomial free-energy curve.

    Parameters
    ----------
    volume : array-like
        Volume values.
    free_energy : array-like
        Free-energy values.
    pressure_energy_density : float, optional
        Pressure in the energy-density scale of the free-energy data.
    degree : int, optional
        Polynomial degree used when coefficients are not provided.
    parameters : array-like, optional
        Polynomial coefficients in ascending order for the physical volume
        coordinate. This argument is retained for direct API compatibility.
    derivative_method : {"analytic", "local_grid"}, optional
        Method used to calculate ``K_T`` and ``K'_T``.
    local_grid_points : int, optional
        Number of points used by the local-grid derivative method.
    local_grid_separation : float, optional
        Adjacent local-grid spacing as a percentage of the central volume.
    local_degree : int or None, optional
        Degree of the local free-energy polynomial. Defaults to ``degree``.

    Returns
    -------
    VolumeMinimumResult
        Equilibrium volume and fit diagnostics.
    """
    try:
        volume_array, free_energy_array = validate_xy(volume, free_energy)
    except ValueError as exc:
        return VolumeMinimumResult.failed(
            str(exc),
            method="polynomial",
            status=MinimumStatus.INVALID_INPUT,
        )

    if parameters is None:
        fitted = fit_polynomial_free_energy_model(
            volume_array,
            free_energy_array,
            degree=degree,
        )
    else:
        coefficients = np.asarray(parameters, dtype=np.float64)
        fit = FitResult(
            success=True,
            status=FitStatus.SUCCESS,
            quality=FitQuality.GOOD,
            message="polynomial coefficients supplied by caller",
            parameters=coefficients,
            n_points=int(volume_array.size),
            n_parameters=int(coefficients.size),
            metadata={
                "model": "polynomial",
                "coordinate": "physical_volume",
                "provided_parameters": True,
            },
        )
        try:
            model = FittedPolynomialFreeEnergy(
                coefficients,
                center=0.0,
                scale=1.0,
                sampled_volumes=volume_array,
            )
        except ValueError as exc:
            return VolumeMinimumResult.failed(
                str(exc),
                method="polynomial",
                status=MinimumStatus.INVALID_INPUT,
                fit=fit,
            )
        fitted = PolynomialModelFitResult(
            True,
            fit,
            model=model,
            message=fit.message,
            metadata={"degree": int(coefficients.size - 1)},
        )

    result = evaluate_fitted_polynomial_at_pressure(
        fitted,
        pressure_energy_density,
        derivative_method=derivative_method,
        local_grid_points=local_grid_points,
        local_grid_separation=local_grid_separation,
        local_degree=degree if local_degree is None else int(local_degree),
    )
    result.fit = fitted.fit
    return result


def fit_eos_free_energy_model(
    volume: ArrayLike,
    free_energy: ArrayLike,
    *,
    eos: str = "BM3",
    energy_unit: str | None = None,
    volume_unit: str | None = None,
    pressure_unit: str | None = None,
    maxfev: int | None = None,
) -> EOSModelFitResult:
    """Fit one energy EOS and build the matching pressure-volume model.

    The EOS family and order are resolved once and shared by the integrated
    energy fit and the pressure evaluator.  Free parameters are converted from
    energy-density units to the requested pressure scale by name: ``K0`` is
    multiplied by the conversion factor, while a fitted ``KPP`` is divided by
    the same factor.
    """
    units = (energy_unit, volume_unit, pressure_unit)
    if any(value is not None for value in units) and not all(
        value is not None for value in units
    ):
        raise ValueError(
            "energy_unit, volume_unit, and pressure_unit must be supplied together"
        )

    energy_eos = EnergyEOS()
    try:
        model_spec = energy_eos.model(eos)
    except ValueError as exc:
        failed = FitResult.failed(str(exc), status=FitStatus.INVALID_INPUT)
        return EOSModelFitResult(False, str(eos), failed, message=str(exc))

    fit = energy_eos.fit(model_spec, volume, free_energy, maxfev=maxfev)
    if not fit.success or fit.parameters is None:
        return EOSModelFitResult(
            success=False,
            eos=model_spec.tag,
            fit=fit,
            message=fit.message,
        )

    parameter_names = tuple(
        str(name)
        for name in fit.metadata.get(
            "parameter_order", model_spec.energy_parameter_names
        )
    )
    parameters = np.asarray(fit.parameters, dtype=np.float64).copy()
    pressure_factor = 1.0
    if energy_unit is not None:
        if volume_unit is None or pressure_unit is None:
            raise ValueError(
                "volume_unit and pressure_unit are required with energy_unit"
            )
        pressure_factor = float(
            energy_to_pressure(1.0, energy_unit, volume_unit, pressure_unit)
        )

    scales = np.ones(parameters.size, dtype=np.float64)
    for index, name in enumerate(parameter_names):
        if name == "K0":
            scales[index] = pressure_factor
        elif name == "KPP":
            scales[index] = 1.0 / pressure_factor
    parameters *= scales

    covariance = None
    warnings_: list[str] = []
    if fit.covariance is not None:
        candidate = np.asarray(fit.covariance, dtype=np.float64)
        expected = (parameters.size, parameters.size)
        if candidate.shape == expected and np.all(np.isfinite(candidate)):
            transform = np.diag(scales)
            covariance = transform @ candidate @ transform.T
        else:
            warnings_.append(
                "the EOS covariance is unavailable or non-finite; "
                "pressure-state uncertainties cannot be propagated"
            )

    try:
        model = FittedEnergyEOS(
            model_spec,
            parameters,
            sampled_volumes=volume,
            covariance=covariance,
        )
    except ValueError as exc:
        if covariance is None:
            return EOSModelFitResult(
                success=False,
                eos=model_spec.tag,
                fit=fit,
                message=f"could not construct the fitted pressure EOS: {exc}",
                warnings=warnings_,
            )
        warnings_.append(
            f"the EOS covariance could not be stabilized and was discarded: {exc}"
        )
        try:
            model = FittedEnergyEOS(
                model_spec,
                parameters,
                sampled_volumes=volume,
                covariance=None,
            )
        except ValueError as model_exc:
            return EOSModelFitResult(
                success=False,
                eos=model_spec.tag,
                fit=fit,
                message=f"could not construct the fitted pressure EOS: {model_exc}",
                warnings=warnings_,
            )

    if fit.quality is FitQuality.POOR:
        warnings_.append("the EOS fit converged with diagnostic warnings")

    resolved = model.resolved_parameters
    return EOSModelFitResult(
        success=True,
        eos=model_spec.tag,
        fit=fit,
        model=model,
        message="EOS model fitted",
        warnings=warnings_,
        metadata={
            "eos_family": model_spec.family.value,
            "eos_tag": model_spec.tag,
            "eos_order": model_spec.order,
            "parameter_order": list(parameter_names),
            "parameter_sources": {"E0": "fitted", **model_spec.parameter_sources},
            "pressure_factor": pressure_factor,
            "energy_unit": energy_unit,
            "volume_unit": volume_unit,
            "pressure_unit": pressure_unit,
            "model_parameters": parameters.tolist(),
            "resolved_model_parameters": resolved.as_dict(),
            "sampled_volume_range": model.sampled_volume_range,
        },
    )


def evaluate_fitted_eos_at_pressure(
    fitted: EOSModelFitResult,
    pressure: float,
    *,
    uncertainty_method: str = "none",
    relative_step: float = 1.0e-5,
    confidence_level: float = 0.95,
    monte_carlo_samples: int = 10000,
    random_state: int | np.random.Generator | None = None,
    minimum_success_fraction: float = 0.8,
) -> VolumeMinimumResult:
    """Evaluate one fitted EOS model at a target pressure.

    Parameters
    ----------
    fitted : EOSModelFitResult
        Successful free-energy EOS fit.
    pressure : float
        Target pressure in the scale used by the fitted model.
    uncertainty_method : str, optional
        ``"none"``, ``"covariance"`` or ``"montecarlo"``.
    relative_step : float, optional
        Relative parameter step used for covariance propagation.
    confidence_level : float, optional
        Confidence probability used for propagated intervals.
    monte_carlo_samples : int, optional
        Number of correlated EOS parameter samples.
    random_state : int, Generator or None, optional
        Random-number source used by Monte Carlo propagation.
    minimum_success_fraction : float, optional
        Minimum accepted fraction of physical Monte Carlo states.

    Returns
    -------
    VolumeMinimumResult
        Pressure-dependent volume, bulk modulus, derivative and uncertainties.
    """
    if not fitted.success or fitted.model is None:
        return VolumeMinimumResult.failed(
            fitted.message or "the EOS model is not available",
            method="eos",
            fit=fitted.fit,
            metadata=fitted.as_dict(),
        )

    warnings_: list[str] = []
    method = str(uncertainty_method).strip().lower()
    if method == "bootstrap":
        warnings_.append(
            "bootstrap uncertainty propagation is not available for EOS states; "
            "state uncertainties were not calculated"
        )
        method = "none"
    if method != "none" and fitted.model.covariance is None:
        warnings_.append(
            "the EOS parameter covariance is unavailable; state uncertainties "
            "were not calculated"
        )
        method = "none"

    try:
        state = fitted.model.state_at_pressure(
            float(pressure),
            uncertainty_method=method,
            relative_step=relative_step,
            confidence_level=confidence_level,
            monte_carlo_samples=monte_carlo_samples,
            random_state=random_state,
            minimum_success_fraction=minimum_success_fraction,
        )
    except ValueError as exc:
        return VolumeMinimumResult.failed(
            str(exc),
            method="eos",
            fit=fitted.fit,
            metadata={
                "pressure": float(pressure),
                "eos_fit": fitted.as_dict(),
            },
        )

    range_status = VolumeRangeStatus.OUTSIDE if state.extrapolated else None
    if range_status is None and fitted.model.sampled_volume_range is not None:
        try:
            range_status = classify_volume(
                state.volume,
                fitted.model.sampled_volume_range,
            )
        except (TypeError, ValueError):
            range_status = None
    if range_status is VolumeRangeStatus.NEAR_BOUNDARY:
        warnings_.append(
            "the equilibrium volume is close to the sampled volume boundary"
        )
    elif range_status is VolumeRangeStatus.OUTSIDE:
        warnings_.append(
            "the equilibrium volume is outside the sampled volume interval"
        )

    return VolumeMinimumResult(
        success=True,
        status=MinimumStatus.SUCCESS,
        method="eos",
        volume=state.volume,
        bulk_modulus=state.bulk_modulus,
        bulk_modulus_derivative=state.bulk_modulus_derivative,
        sigma_volume=state.sigma_volume,
        sigma_bulk_modulus=state.sigma_bulk_modulus,
        sigma_bulk_modulus_derivative=state.sigma_bulk_modulus_derivative,
        message="EOS state evaluated",
        range_status=range_status,
        fit=None,
        warnings=warnings_,
        metadata={
            "pressure": float(pressure),
            "eos": fitted.eos,
            "state": state.as_dict(),
        },
    )


def minimize_eos(
    volume: ArrayLike,
    free_energy: ArrayLike,
    *,
    pressure_energy_density: float = 0.0,
    eos: str = "BM3",
    energy_unit: str | None = None,
    volume_unit: str | None = None,
    pressure_unit: str | None = None,
    maxfev: int | None = None,
) -> VolumeMinimumResult:
    """Evaluate a free-energy EOS at the requested pressure.

    Parameters
    ----------
    volume : array-like
        Volume values.
    free_energy : array-like
        Free-energy values.
    pressure_energy_density : float, optional
        Pressure in the energy-density scale of the free-energy data.
    eos : str, optional
        Equation-of-state model.
    energy_unit, volume_unit, pressure_unit : str, optional
        Units used to convert the EOS bulk modulus from energy density to a
        pressure unit.  If omitted, ``bulk_modulus`` is returned in the energy
        density scale of the free-energy data.
    maxfev : int, optional
        Maximum number of optimizer evaluations.

    Returns
    -------
    VolumeMinimumResult
        Equilibrium volume, EOS parameters and fit diagnostics.
    """
    try:
        fitted = fit_eos_free_energy_model(
            volume,
            free_energy,
            eos=eos,
            energy_unit=energy_unit,
            volume_unit=volume_unit,
            pressure_unit=pressure_unit,
            maxfev=maxfev,
        )
    except ValueError as exc:
        return VolumeMinimumResult.failed(
            str(exc),
            method="eos",
            status=MinimumStatus.INVALID_INPUT,
        )

    target_pressure = float(pressure_energy_density)
    if (
        energy_unit is not None
        and volume_unit is not None
        and pressure_unit is not None
    ):
        target_pressure = float(
            energy_to_pressure(
                pressure_energy_density,
                energy_unit,
                volume_unit,
                pressure_unit,
            )
        )
    minimum = evaluate_fitted_eos_at_pressure(fitted, target_pressure)
    minimum.fit = fitted.fit
    minimum.metadata["pressure_energy_density"] = float(pressure_energy_density)
    return minimum


def numerical_bulk_modulus(
    volume: ArrayLike,
    free_energy: ArrayLike,
    *,
    index: int,
    degree: int = 4,
    energy_unit: str | None = None,
    volume_unit: str | None = None,
    pressure_unit: str | None = None,
) -> tuple[float, float, FitResult]:
    """Calculate ``K_T`` and ``K'_T`` from a local polynomial grid.

    Parameters
    ----------
    volume : array-like
        Local volume values.
    free_energy : array-like
        Helmholtz free energies evaluated on ``volume``.
    index : int
        Index of the volume at which the quantities are evaluated.
    degree : int, optional
        Polynomial degree used for the local free-energy fit.
    energy_unit, volume_unit, pressure_unit : str, optional
        Units used to convert the bulk modulus to pressure. All three values
        must be supplied together.

    Returns
    -------
    tuple
        Isothermal bulk modulus, pressure derivative of the bulk modulus and
        local polynomial fit diagnostics.

    Raises
    ------
    RuntimeError
        If the local polynomial analysis fails.
    IndexError
        If ``index`` is outside the volume array.
    ValueError
        If only part of the unit triplet is supplied.
    """
    kt, kp, local_fit, _ = _local_polynomial_bulk_properties(
        volume,
        free_energy,
        index=index,
        degree=degree,
    )
    units = (energy_unit, volume_unit, pressure_unit)
    if any(value is not None for value in units) and not all(
        value is not None for value in units
    ):
        raise ValueError(
            "energy_unit, volume_unit, and pressure_unit must be supplied together"
        )
    if energy_unit is not None:
        if volume_unit is None or pressure_unit is None:
            raise ValueError(
                "volume_unit and pressure_unit are required with energy_unit"
            )
        kt = float(
            energy_to_pressure(
                kt,
                energy_unit,
                volume_unit,
                pressure_unit,
            )
        )
    return float(kt), float(kp), local_fit.fit
