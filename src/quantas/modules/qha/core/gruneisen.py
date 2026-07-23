# -*- coding: utf-8 -*-

"""Grüneisen-parameter utilities for quasi-harmonic workflows.

The functions in this module evaluate mode-resolved and heat-capacity-weighted
Grüneisen parameters from volume-dependent phonon frequencies.  They operate on
plain arrays and return structured diagnostic objects so that callers can decide
whether a dataset is suitable for a mode-continuous QHA analysis.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping, Sequence

import numpy as np
from scipy import constants as cs

from quantas.core.math.fitting import FitResult
from quantas.core.math.derivative import polynomial_derivative_from_coefficients
from quantas.core.math.polynomials import fit_polynomial_result as polynomial_fit

ArrayLike = np.ndarray | Sequence[float]


class GruneisenStatus(str, Enum):
    """Execution status of a Grüneisen-parameter calculation.

    Attributes
    ----------
    SUCCESS
        The requested quantities were calculated successfully.
    PARTIAL
        At least one mode failed while at least one result is usable.
    FAILED
        No usable result was produced.
    INVALID_INPUT
        The input arrays are not suitable for the calculation.
    """

    SUCCESS = "success"
    PARTIAL = "partial"
    FAILED = "failed"
    INVALID_INPUT = "invalid_input"


@dataclass(slots=True)
class ModeGruneisenResult:
    """Mode-resolved Grüneisen parameters and diagnostics.

    Parameters
    ----------
    success : bool
        Whether all requested modes were evaluated successfully.
    status : GruneisenStatus
        Execution status of the calculation.
    volume : ndarray or None
        Volumes at which the Grüneisen parameters were evaluated.
    gamma : ndarray or None
        Mode Grüneisen parameters.  The leading dimensions match the input
        phonon-mode dimensions and the last dimension corresponds to volume.
    degree : int
        Polynomial degree used to fit ``ln(frequency)`` as a function of
        ``ln(volume)``.
    fits : dict
        Fit diagnostics keyed by mode index.
    failed_modes : list of tuple
        Indices of modes that could not be evaluated.
    warnings : list of str
        Non-fatal diagnostic messages.
    metadata : dict
        Additional caller-defined diagnostic information.
    """

    success: bool
    status: GruneisenStatus
    volume: np.ndarray | None = None
    gamma: np.ndarray | None = None
    degree: int = 2
    fits: dict[tuple[int, ...], FitResult] = field(default_factory=dict)
    failed_modes: list[tuple[int, ...]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def failed(
        cls,
        message: str,
        *,
        degree: int = 2,
        status: GruneisenStatus = GruneisenStatus.FAILED,
        metadata: Mapping[str, Any] | None = None,
    ) -> ModeGruneisenResult:
        """Create a failed mode-Grüneisen result.

        Parameters
        ----------
        message : str
            Explanation of the failure.
        degree : int, optional
            Requested polynomial degree.
        status : GruneisenStatus, optional
            Failure category.
        metadata : mapping, optional
            Additional diagnostic information.

        Returns
        -------
        ModeGruneisenResult
            Failed result with no numerical arrays.
        """
        return cls(
            success=False,
            status=status,
            degree=int(degree),
            warnings=[message],
            metadata=dict(metadata or {}),
        )

    def as_dict(self) -> dict[str, Any]:
        """Return the result as a serializable dictionary.

        Returns
        -------
        dict
            Dictionary representation of the result.
        """
        return {
            "success": self.success,
            "status": self.status.value,
            "volume": None if self.volume is None else self.volume.tolist(),
            "gamma": None if self.gamma is None else self.gamma.tolist(),
            "degree": self.degree,
            "fits": {str(key): value.as_dict() for key, value in self.fits.items()},
            "failed_modes": [list(index) for index in self.failed_modes],
            "warnings": list(self.warnings),
            "metadata": dict(self.metadata),
        }


@dataclass(slots=True)
class ThermalGruneisenResult:
    """Heat-capacity-weighted Grüneisen parameter.

    Parameters
    ----------
    success : bool
        Whether the calculation produced finite values.
    status : GruneisenStatus
        Execution status of the calculation.
    gamma : ndarray or None
        Heat-capacity-weighted Grüneisen parameter.
    denominator : ndarray or None
        Sum of heat-capacity weights used in the average.
    message : str
        Human-readable status message.
    warnings : list of str
        Non-fatal diagnostic messages.
    metadata : dict
        Additional caller-defined diagnostic information.
    """

    success: bool
    status: GruneisenStatus
    gamma: np.ndarray | None = None
    denominator: np.ndarray | None = None
    message: str = ""
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def failed(
        cls,
        message: str,
        *,
        status: GruneisenStatus = GruneisenStatus.FAILED,
        metadata: Mapping[str, Any] | None = None,
    ) -> ThermalGruneisenResult:
        """Create a failed thermal-Grüneisen result.

        Parameters
        ----------
        message : str
            Explanation of the failure.
        status : GruneisenStatus, optional
            Failure category.
        metadata : mapping, optional
            Additional diagnostic information.

        Returns
        -------
        ThermalGruneisenResult
            Failed result with no numerical arrays.
        """
        return cls(
            success=False,
            status=status,
            message=message,
            warnings=[message],
            metadata=dict(metadata or {}),
        )

    def as_dict(self) -> dict[str, Any]:
        """Return the result as a serializable dictionary.

        Returns
        -------
        dict
            Dictionary representation of the result.
        """
        return {
            "success": self.success,
            "status": self.status.value,
            "gamma": None if self.gamma is None else self.gamma.tolist(),
            "denominator": None
            if self.denominator is None
            else self.denominator.tolist(),
            "message": self.message,
            "warnings": list(self.warnings),
            "metadata": dict(self.metadata),
        }


def validate_mode_continuity(
    frequencies: ArrayLike, volumes: ArrayLike
) -> tuple[np.ndarray, np.ndarray]:
    """Validate volume-dependent phonon frequencies for mode analysis.

    Parameters
    ----------
    frequencies : array-like
        Phonon frequencies with volume as the last dimension.
    volumes : array-like
        Volumes corresponding to the last dimension of ``frequencies``.

    Returns
    -------
    tuple of ndarray
        Frequency and volume arrays converted to ``float64``.

    Raises
    ------
    ValueError
        If the arrays are not finite, have incompatible shapes, or contain
        non-positive volumes.
    """
    volume_array = np.asarray(volumes, dtype=np.float64)
    frequency_array = np.asarray(frequencies, dtype=np.float64)
    if volume_array.ndim != 1:
        raise ValueError("volumes must be a one-dimensional array")
    if volume_array.size < 2:
        raise ValueError("at least two volumes are required")
    if frequency_array.ndim < 1:
        raise ValueError("frequencies must have at least one dimension")
    if frequency_array.shape[-1] != volume_array.size:
        raise ValueError(
            "the last frequency dimension must match the number of volumes"
        )
    if not np.all(np.isfinite(volume_array)) or np.any(volume_array <= 0.0):
        raise ValueError("volumes must be finite and positive")
    if not np.all(np.isfinite(frequency_array)):
        raise ValueError("frequencies must be finite")
    return frequency_array, volume_array


class ModeGruneisenEvaluator:
    """Fit and evaluate mode-resolved Grüneisen parameters.

    Parameters
    ----------
    frequencies : array-like
        Mode-continuous frequencies with volume as the last dimension.
    volumes : array-like
        Positive sampled volumes.
    degree : int, optional
        Polynomial degree used to fit ``ln(frequency)`` against ``ln(volume)``.
    allow_nonpositive : bool, optional
        If ``True``, modes containing non-positive values are excluded.

    Raises
    ------
    ValueError
        If the arrays are inconsistent or no usable mode can be fitted.
    """

    def __init__(
        self,
        frequencies: ArrayLike,
        volumes: ArrayLike,
        *,
        degree: int = 2,
        allow_nonpositive: bool = False,
    ) -> None:
        frequency_array, volume_array = validate_mode_continuity(frequencies, volumes)
        degree = int(degree)
        if degree < 1:
            raise ValueError("degree must be at least one")
        if volume_array.size <= degree:
            raise ValueError(
                "number of volumes must be larger than the polynomial degree"
            )
        if np.any(frequency_array <= 0.0) and not allow_nonpositive:
            raise ValueError("frequencies must be positive for logarithmic derivatives")

        self.frequencies = frequency_array
        self.volumes = volume_array
        self.degree = degree
        self.allow_nonpositive = bool(allow_nonpositive)
        self.mode_shape = frequency_array.shape[:-1]
        self.fits: dict[tuple[int, ...], FitResult] = {}
        self.failed_modes: list[tuple[int, ...]] = []
        self.warnings: list[str] = []
        self._coefficients = np.full(
            (int(np.prod(self.mode_shape)), degree + 1),
            np.nan,
            dtype=np.float64,
        )
        self._fit_modes()
        if not np.any(np.all(np.isfinite(self._coefficients), axis=1)):
            raise ValueError("no usable mode Grüneisen fit was produced")

    def _fit_modes(self) -> None:
        """Fit logarithmic frequency-volume curves for every mode."""
        log_volume = np.log(self.volumes)
        flat = self.frequencies.reshape((-1, self.volumes.size))
        for flat_index, values in enumerate(flat):
            mode_index = tuple(
                int(item) for item in np.unravel_index(flat_index, self.mode_shape)
            )
            if np.any(values <= 0.0):
                self.failed_modes.append(mode_index)
                self.warnings.append(
                    f"mode {mode_index} contains non-positive frequencies"
                )
                continue
            fit = polynomial_fit(log_volume, np.log(values), degree=self.degree)
            self.fits[mode_index] = fit
            if not fit.success or fit.parameters is None:
                self.failed_modes.append(mode_index)
                self.warnings.append(f"mode {mode_index} logarithmic fit failed")
                continue
            self._coefficients[flat_index] = fit.parameters

    @property
    def usable_modes(self) -> int:
        """Return the number of successfully fitted modes.

        Returns
        -------
        int
            Number of modes with finite logarithmic-fit coefficients.
        """
        return int(np.count_nonzero(np.all(np.isfinite(self._coefficients), axis=1)))

    def gamma_at(self, volume: ArrayLike | float) -> np.ndarray:
        """Evaluate mode Grüneisen parameters at selected volumes.

        Parameters
        ----------
        volume : array-like or float
            Positive volume values.

        Returns
        -------
        ndarray
            Mode Grüneisen parameters with shape ``mode_shape + volume_shape``.

        Raises
        ------
        ValueError
            If a requested volume is non-positive or non-finite.
        """
        volume_array = np.asarray(volume, dtype=np.float64)
        if not np.all(np.isfinite(volume_array)) or np.any(volume_array <= 0.0):
            raise ValueError("volumes must be finite and positive")
        log_volume = np.log(volume_array)
        output_shape = self.mode_shape + volume_array.shape
        gamma = np.full(output_shape, np.nan, dtype=np.float64)
        flat_gamma = gamma.reshape((self._coefficients.shape[0],) + volume_array.shape)
        for index, coefficients in enumerate(self._coefficients):
            if not np.all(np.isfinite(coefficients)):
                continue
            flat_gamma[index] = -polynomial_derivative_from_coefficients(
                coefficients,
                log_volume,
            )
        return gamma

    def result(self) -> ModeGruneisenResult:
        """Return sampled mode parameters and fit diagnostics.

        Returns
        -------
        ModeGruneisenResult
            Structured mode-Grüneisen result on the sampled volume grid.
        """
        gamma = self.gamma_at(self.volumes)
        status = (
            GruneisenStatus.SUCCESS
            if not self.failed_modes
            else GruneisenStatus.PARTIAL
        )
        r_squared = np.asarray(
            [
                fit.r_squared
                for fit in self.fits.values()
                if fit.r_squared is not None and np.isfinite(fit.r_squared)
            ],
            dtype=np.float64,
        )
        metadata: dict[str, Any] = {
            "n_modes": int(np.prod(self.mode_shape)),
            "n_usable_modes": self.usable_modes,
            "n_failed_modes": len(self.failed_modes),
        }
        if r_squared.size:
            metadata["fit_r_squared"] = {
                "minimum": float(np.min(r_squared)),
                "mean": float(np.mean(r_squared)),
                "median": float(np.median(r_squared)),
                "maximum": float(np.max(r_squared)),
            }
        return ModeGruneisenResult(
            success=not self.failed_modes,
            status=status,
            volume=self.volumes.copy(),
            gamma=gamma,
            degree=self.degree,
            fits=dict(self.fits),
            failed_modes=list(self.failed_modes),
            warnings=list(self.warnings),
            metadata=metadata,
        )


def thermal_gruneisen_from_modes(
    gamma: ArrayLike,
    frequencies_hz: ArrayLike,
    temperature_k: ArrayLike | float,
    qpoint_weights: ArrayLike,
) -> ThermalGruneisenResult:
    """Calculate heat-capacity-weighted mode Grüneisen parameters.

    Parameters
    ----------
    gamma : array-like
        Mode Grüneisen parameters with shape ``(nq, nmodes, nstate)``.
    frequencies_hz : array-like
        Frequencies in hertz with the same shape as ``gamma``.
    temperature_k : array-like or float
        Temperatures in kelvin. A scalar applies to all states; otherwise the
        last dimension must match ``nstate``.
    qpoint_weights : array-like
        Q-point weights with shape ``(nq,)``.

    Returns
    -------
    ThermalGruneisenResult
        Heat-capacity-weighted mode average for each state.

    Raises
    ------
    ValueError
        Never raised directly; invalid inputs are reported in the result.
    """
    try:
        gamma_array = np.asarray(gamma, dtype=np.float64)
        frequency_array = np.asarray(frequencies_hz, dtype=np.float64)
        weights = np.asarray(qpoint_weights, dtype=np.float64)
        temperature = np.asarray(temperature_k, dtype=np.float64)
        if gamma_array.shape != frequency_array.shape or gamma_array.ndim != 3:
            raise ValueError(
                "gamma and frequencies must have shape (qpoints, modes, states)"
            )
        if weights.ndim != 1 or weights.size != gamma_array.shape[0]:
            raise ValueError("q-point weights must match the q-point axis")
        if not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
            raise ValueError("q-point weights must be finite and non-negative")
        if weights.sum() <= 0.0:
            raise ValueError("q-point weights must have a positive sum")
        if temperature.ndim == 0:
            temperature = np.full(gamma_array.shape[-1], float(temperature))
        if temperature.ndim != 1 or temperature.size != gamma_array.shape[-1]:
            raise ValueError("temperature must be scalar or match the state axis")
    except ValueError as exc:
        return ThermalGruneisenResult.failed(
            str(exc), status=GruneisenStatus.INVALID_INPUT
        )

    normalized_weights = weights / weights.sum()
    temperature_grid = temperature[None, None, :]
    valid = (
        np.isfinite(gamma_array)
        & np.isfinite(frequency_array)
        & (frequency_array > 0.0)
        & (temperature_grid > 0.0)
    )
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        x = cs.Planck * frequency_array / (cs.Boltzmann * temperature_grid)
        exponential = np.exp(x)
        contribution = exponential * (x / np.expm1(x)) ** 2
    contribution = np.where(valid & np.isfinite(contribution), contribution, 0.0)
    weighted_contribution = contribution * normalized_weights[:, None, None]
    denominator = np.sum(weighted_contribution, axis=(0, 1))
    numerator = np.sum(
        np.where(valid, gamma_array, 0.0) * weighted_contribution,
        axis=(0, 1),
    )
    average = np.zeros_like(denominator, dtype=np.float64)
    nonzero = denominator > 0.0
    average[nonzero] = numerator[nonzero] / denominator[nonzero]
    if not np.all(np.isfinite(average)):
        return ThermalGruneisenResult.failed(
            "weighted mode Grüneisen parameter is not finite"
        )
    return ThermalGruneisenResult(
        success=True,
        status=GruneisenStatus.SUCCESS,
        gamma=average,
        denominator=denominator,
        message="weighted mode Grüneisen parameter calculated successfully",
        metadata={
            "zero_weight_states": int(np.count_nonzero(~nonzero)),
            "weighting": "harmonic_mode_Cv",
        },
    )


def mode_gruneisen(
    frequencies: ArrayLike,
    volumes: ArrayLike,
    *,
    degree: int = 2,
    allow_nonpositive: bool = False,
) -> ModeGruneisenResult:
    """Calculate mode Grüneisen parameters from frequency-volume data.

    Parameters
    ----------
    frequencies : array-like
        Mode-continuous frequencies with volume as the last dimension.
    volumes : array-like
        Volumes corresponding to the frequency samples.
    degree : int, optional
        Polynomial degree used to fit ``ln(frequency)`` against ``ln(volume)``.
    allow_nonpositive : bool, optional
        If ``True``, non-positive modes are excluded and reported.

    Returns
    -------
    ModeGruneisenResult
        Mode-resolved Grüneisen parameters and fit diagnostics.
    """
    try:
        evaluator = ModeGruneisenEvaluator(
            frequencies,
            volumes,
            degree=degree,
            allow_nonpositive=allow_nonpositive,
        )
    except ValueError as exc:
        return ModeGruneisenResult.failed(
            str(exc),
            degree=int(degree),
            status=GruneisenStatus.INVALID_INPUT,
        )
    return evaluator.result()


def thermal_gruneisen(
    gamma: ArrayLike, heat_capacity: ArrayLike, *, axis: int = 0
) -> ThermalGruneisenResult:
    """Calculate a heat-capacity-weighted Grüneisen parameter.

    Parameters
    ----------
    gamma : array-like
        Mode Grüneisen parameters.
    heat_capacity : array-like
        Mode heat capacities or equivalent positive weights.  The shape must be
        broadcast-compatible with ``gamma``.
    axis : int, optional
        Axis corresponding to the mode summation.

    Returns
    -------
    ThermalGruneisenResult
        Weighted Grüneisen parameter and denominator of the average.
    """
    try:
        gamma_array = np.asarray(gamma, dtype=np.float64)
        weight_array = np.asarray(heat_capacity, dtype=np.float64)
        gamma_b, weight_b = np.broadcast_arrays(gamma_array, weight_array)
    except ValueError as exc:
        return ThermalGruneisenResult.failed(
            str(exc), status=GruneisenStatus.INVALID_INPUT
        )

    if not np.all(np.isfinite(gamma_b)) or not np.all(np.isfinite(weight_b)):
        return ThermalGruneisenResult.failed(
            "gamma and heat-capacity arrays must be finite",
            status=GruneisenStatus.INVALID_INPUT,
        )
    denominator = np.sum(weight_b, axis=axis)
    if np.any(denominator == 0.0):
        return ThermalGruneisenResult.failed(
            "heat-capacity weights must not sum to zero",
            status=GruneisenStatus.INVALID_INPUT,
        )
    weighted = np.sum(gamma_b * weight_b, axis=axis) / denominator
    if not np.all(np.isfinite(weighted)):
        return ThermalGruneisenResult.failed(
            "weighted Grüneisen parameter is not finite"
        )

    return ThermalGruneisenResult(
        success=True,
        status=GruneisenStatus.SUCCESS,
        gamma=np.asarray(weighted, dtype=np.float64),
        denominator=np.asarray(denominator, dtype=np.float64),
        message="weighted Grüneisen parameter calculated successfully",
    )


def gruneisen_from_power_law(
    volumes: ArrayLike,
    gamma: float,
    reference_frequency: float,
    reference_volume: float,
) -> np.ndarray:
    """Generate frequencies following a constant-Grüneisen power law.

    Parameters
    ----------
    volumes : array-like
        Volumes at which frequencies are evaluated.
    gamma : float
        Constant mode Grüneisen parameter.
    reference_frequency : float
        Frequency at ``reference_volume``.
    reference_volume : float
        Reference volume.

    Returns
    -------
    ndarray
        Frequencies satisfying ``nu(V) = nu0 * (V / V0)**(-gamma)``.

    Raises
    ------
    ValueError
        If the input values are not physically meaningful.
    """
    volume_array = np.asarray(volumes, dtype=np.float64)
    if not np.all(np.isfinite(volume_array)) or np.any(volume_array <= 0.0):
        raise ValueError("volumes must be finite and positive")
    if reference_frequency <= 0.0 or reference_volume <= 0.0:
        raise ValueError("reference frequency and volume must be positive")
    return float(reference_frequency) * (volume_array / float(reference_volume)) ** (
        -float(gamma)
    )
