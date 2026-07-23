# -*- coding: utf-8 -*-

"""Passive fit-result contracts for thermoelastic calibration."""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any
import numpy as np
from quantas.core.math.fitting import FitResult
from .types import FloatArray, ThermoelasticQualityLevel


@dataclass(slots=True)
class ThermoelasticFitQuality:
    """Scientific support diagnostics for one elastic-component fit.

    Parameters
    ----------
    level : {"supported", "caution", "unsupported", "not_applicable"}
        Overall diagnostic classification.  The classification never modifies
        fitted parameters or reconstructed tensors.
    issues : tuple of str
        Stable machine-readable diagnostic identifiers.
    n_observations, n_parameters, degrees_of_freedom : int
        Data and model dimensions.
    eulerian_strain_min, eulerian_strain_max, eulerian_strain_span : float
        Sampled finite-strain support.
    reference_volume_bracketed : bool
        Whether the fitted static-EOS ``V0`` lies inside the elastic-volume
        interval.
    reference_volume_distance_fraction : float
        Distance of ``V0`` outside the interval, normalized by interval width.
    design_rank : int
        Rank of the component design matrix.
    design_condition_number : float
        Two-norm condition number after independent column normalization.
    leverage : ndarray
        Diagonal of the ordinary-least-squares hat matrix.
    maximum_leverage : float
        Maximum point leverage.
    maximum_relative_symmetry_spread : float
        Largest normalized disagreement among symmetry-equivalent entries.
    maximum_leave_one_out_parameter_change : float or None
        Largest scale-aware parameter change after omitting one point.
    maximum_order_parameter_change : float or None
        Scale-aware parameter change between finite-strain orders.
    thresholds : dict
        User-visible thresholds used for the classification.
    """

    level: ThermoelasticQualityLevel
    issues: tuple[str, ...]
    n_observations: int
    n_parameters: int
    degrees_of_freedom: int
    eulerian_strain_min: float
    eulerian_strain_max: float
    eulerian_strain_span: float
    reference_volume_bracketed: bool
    reference_volume_distance_fraction: float
    design_rank: int
    design_condition_number: float
    leverage: FloatArray
    maximum_leverage: float
    maximum_relative_symmetry_spread: float
    maximum_leave_one_out_parameter_change: float | None
    maximum_order_parameter_change: float | None
    thresholds: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize diagnostic values and validate dimensions."""
        if self.level not in {
            "supported",
            "caution",
            "unsupported",
            "not_applicable",
        }:
            raise ValueError("invalid thermoelastic fit-quality level")
        self.issues = tuple(str(value) for value in self.issues)
        for name in (
            "n_observations",
            "n_parameters",
            "degrees_of_freedom",
            "design_rank",
        ):
            value = int(getattr(self, name))
            if value < 0:
                raise ValueError(f"{name} must be non-negative")
            setattr(self, name, value)
        self.leverage = np.asarray(self.leverage, dtype=np.float64).copy()
        if self.leverage.shape != (self.n_observations,):
            raise ValueError("leverage must have shape (n_observations,)")
        if np.any(~np.isfinite(self.leverage)) or np.any(self.leverage < 0.0):
            raise ValueError("leverage must be finite and non-negative")
        for name in (
            "eulerian_strain_min",
            "eulerian_strain_max",
            "eulerian_strain_span",
            "reference_volume_distance_fraction",
            "maximum_leverage",
            "maximum_relative_symmetry_spread",
        ):
            float_value = float(getattr(self, name))
            if not np.isfinite(float_value):
                raise ValueError(f"{name} must be finite")
            setattr(self, name, float_value)
        self.design_condition_number = float(self.design_condition_number)
        if np.isnan(self.design_condition_number) or self.design_condition_number < 0.0:
            raise ValueError(
                "design_condition_number must be non-negative or positive infinity"
            )
        for name in (
            "maximum_leave_one_out_parameter_change",
            "maximum_order_parameter_change",
        ):
            value = getattr(self, name)
            if value is not None:
                normalized = float(value)
                if not np.isfinite(normalized) or normalized < 0.0:
                    raise ValueError(f"{name} must be finite and non-negative")
                setattr(self, name, normalized)
        self.reference_volume_bracketed = bool(self.reference_volume_bracketed)
        self.thresholds = dict(self.thresholds)

    def as_dict(self) -> dict[str, Any]:
        """Return a recursively serializable diagnostic mapping."""
        return {
            "level": self.level,
            "issues": list(self.issues),
            "n_observations": self.n_observations,
            "n_parameters": self.n_parameters,
            "degrees_of_freedom": self.degrees_of_freedom,
            "eulerian_strain_min": self.eulerian_strain_min,
            "eulerian_strain_max": self.eulerian_strain_max,
            "eulerian_strain_span": self.eulerian_strain_span,
            "reference_volume_bracketed": self.reference_volume_bracketed,
            "reference_volume_distance_fraction": (
                self.reference_volume_distance_fraction
            ),
            "design_rank": self.design_rank,
            "design_condition_number": self.design_condition_number,
            "leverage": self.leverage.copy(),
            "maximum_leverage": self.maximum_leverage,
            "maximum_relative_symmetry_spread": (self.maximum_relative_symmetry_spread),
            "maximum_leave_one_out_parameter_change": (
                self.maximum_leave_one_out_parameter_change
            ),
            "maximum_order_parameter_change": self.maximum_order_parameter_change,
            "thresholds": dict(self.thresholds),
        }


@dataclass(slots=True)
class ReferenceEOSFit:
    """Static reference EOS parameters and complete fit diagnostics.

    Parameters
    ----------
    eos : str
        Canonical EOS tag.
    reference_volume : float
        ``V0`` in angstrom cubed per normalized cell.
    bulk_modulus : float
        ``K0`` in GPa.
    bulk_modulus_derivative : float
        ``Kprime`` at the reference state.
    bulk_modulus_second_derivative : float
        ``Ksecond`` in inverse GPa.
    covariance : ndarray or None
        Covariance of ``V0, K0, Kprime`` in that order.
    fit : FitResult
        Complete integrated energy-EOS fit result.
    metadata : dict, optional
        Unit conversions and parameter provenance.
    """

    eos: str
    reference_volume: float
    bulk_modulus: float
    bulk_modulus_derivative: float
    bulk_modulus_second_derivative: float
    covariance: FloatArray | None
    fit: FitResult
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate reference parameters and normalize covariance."""
        values = np.asarray(
            [
                self.reference_volume,
                self.bulk_modulus,
                self.bulk_modulus_derivative,
                self.bulk_modulus_second_derivative,
            ],
            dtype=np.float64,
        )
        if (
            not np.all(np.isfinite(values))
            or self.reference_volume <= 0.0
            or self.bulk_modulus <= 0.0
        ):
            raise ValueError(
                "reference EOS parameters must be finite with positive V0 and K0"
            )
        if self.covariance is not None:
            covariance = np.asarray(self.covariance, dtype=np.float64)
            if covariance.shape != (3, 3) or not np.all(np.isfinite(covariance)):
                raise ValueError(
                    "reference EOS covariance must be finite with shape (3, 3)"
                )
            self.covariance = covariance.copy()
        self.metadata = dict(self.metadata)


@dataclass(slots=True)
class ElasticComponentFit:
    """Fit result and point-level diagnostics for one independent component.

    Parameters
    ----------
    label : str
        Canonical Voigt label.
    entries : tuple of tuple
        Symmetry-equivalent ``(row, column, multiplier)`` entries.
    wallace_delta : float
        Wallace hydrostatic-delta entry used by the model.
    volumes, pressures : ndarray
        Sampled primitive volumes and CRYSTAL pressures.
    observed, fitted, residuals, relative_residuals : ndarray
        Point-level fit diagnostics in GPa.
    symmetry_spread : ndarray
        Maximum disagreement among symmetry-equivalent observations at each
        volume after sign normalization.
    fit : FitResult or None
        General fitting result. ``None`` for numerically zero components.
    active : bool
        Whether the component was fitted.
    zero_by_tolerance : bool
        Whether it was retained as an exact numerical zero.
    metadata : dict, optional
        Component-specific diagnostics.
    quality : ThermoelasticFitQuality or None, optional
        Scientific support diagnostics for the fitted component.
    """

    label: str
    entries: tuple[tuple[int, int, float], ...]
    wallace_delta: float
    volumes: FloatArray
    pressures: FloatArray
    observed: FloatArray
    fitted: FloatArray
    residuals: FloatArray
    relative_residuals: FloatArray
    symmetry_spread: FloatArray
    fit: FitResult | None
    active: bool = True
    zero_by_tolerance: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)
    quality: ThermoelasticFitQuality | None = None

    def __post_init__(self) -> None:
        """Normalize aligned diagnostic arrays."""
        arrays = []
        finite_required = {"volumes", "pressures", "observed", "symmetry_spread"}
        for name in (
            "volumes",
            "pressures",
            "observed",
            "fitted",
            "residuals",
            "relative_residuals",
            "symmetry_spread",
        ):
            value = np.asarray(getattr(self, name), dtype=np.float64)
            if value.ndim != 1:
                raise ValueError(f"{name} must be a one-dimensional array")
            if name in finite_required and not np.all(np.isfinite(value)):
                raise ValueError(f"{name} must contain finite values")
            setattr(self, name, value.copy())
            arrays.append(value)
        if len({array.size for array in arrays}) != 1:
            raise ValueError("component diagnostic arrays must have equal length")
        self.entries = tuple(
            (int(i), int(j), float(multiplier)) for i, j, multiplier in self.entries
        )
        self.wallace_delta = float(self.wallace_delta)
        self.metadata = dict(self.metadata)

    @property
    def parameters(self) -> FloatArray | None:
        """Return ``C0, Cprime`` when available."""
        if self.fit is None or self.fit.parameters is None:
            return None
        return np.asarray(self.fit.parameters, dtype=np.float64).copy()

    @property
    def covariance(self) -> FloatArray | None:
        """Return the fitted ``2 x 2`` parameter covariance."""
        if self.fit is None or self.fit.covariance is None:
            return None
        return np.asarray(self.fit.covariance, dtype=np.float64).copy()


__all__ = ["ElasticComponentFit", "ReferenceEOSFit", "ThermoelasticFitQuality"]
