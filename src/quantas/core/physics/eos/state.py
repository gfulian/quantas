# -*- coding: utf-8 -*-

"""Passive state containers for fitted equations of state."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

_STATE_ORDER = ("volume", "bulk_modulus", "bulk_modulus_derivative")


@dataclass(slots=True, frozen=True)
class EOSStateUncertainty:
    """Uncertainty estimate for a pressure-dependent EOS state.

    Parameters
    ----------
    method : str
        Propagation method, currently ``"covariance"`` or ``"montecarlo"``.
    mean : ndarray
        Mean values ordered as volume, bulk modulus and bulk-modulus
        derivative.
    covariance : ndarray
        State covariance matrix with shape ``(3, 3)``.
    correlation : ndarray
        State correlation matrix with shape ``(3, 3)``.
    confidence_intervals : ndarray, optional
        Lower and upper confidence bounds with shape ``(3, 2)``.
    n_samples : int, optional
        Number of Monte Carlo parameter samples requested.
    n_successful : int, optional
        Number of Monte Carlo samples producing a physical EOS state.
    metadata : dict, optional
        Additional propagation diagnostics.
    """

    method: str
    mean: np.ndarray
    covariance: np.ndarray
    correlation: np.ndarray
    confidence_intervals: np.ndarray | None = None
    n_samples: int | None = None
    n_successful: int | None = None
    metadata: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        """Validate and detach uncertainty arrays from caller-owned data.

        Raises
        ------
        ValueError
            If an uncertainty array has an incompatible shape or contains
            non-finite values.
        """
        mean = np.asarray(self.mean, dtype=np.float64)
        covariance = np.asarray(self.covariance, dtype=np.float64)
        correlation = np.asarray(self.correlation, dtype=np.float64)
        intervals = None
        if self.confidence_intervals is not None:
            intervals = np.asarray(self.confidence_intervals, dtype=np.float64)

        if mean.shape != (3,):
            raise ValueError("uncertainty mean must have shape (3,)")
        if covariance.shape != (3, 3):
            raise ValueError("state covariance must have shape (3, 3)")
        if correlation.shape != (3, 3):
            raise ValueError("state correlation must have shape (3, 3)")
        if intervals is not None and intervals.shape != (3, 2):
            raise ValueError("confidence intervals must have shape (3, 2)")
        arrays = [mean, covariance, correlation]
        if intervals is not None:
            arrays.append(intervals)
        if not all(np.all(np.isfinite(values)) for values in arrays):
            raise ValueError("uncertainty arrays must contain only finite values")

        object.__setattr__(self, "mean", mean.copy())
        object.__setattr__(self, "covariance", covariance.copy())
        object.__setattr__(self, "correlation", correlation.copy())
        object.__setattr__(
            self,
            "confidence_intervals",
            None if intervals is None else intervals.copy(),
        )
        object.__setattr__(self, "metadata", dict(self.metadata or {}))

    @property
    def standard_deviations(self) -> np.ndarray:
        """Return one-sigma uncertainties for the three state quantities.

        Returns
        -------
        ndarray
            Standard deviations ordered as volume, bulk modulus and
            bulk-modulus derivative.
        """
        diagonal = np.clip(np.diag(self.covariance), 0.0, None)
        return np.sqrt(diagonal)

    @property
    def sigma_volume(self) -> float:
        """Return the one-sigma volume uncertainty.

        Returns
        -------
        float
            Volume standard deviation.
        """
        return float(self.standard_deviations[0])

    @property
    def sigma_bulk_modulus(self) -> float:
        """Return the one-sigma bulk-modulus uncertainty.

        Returns
        -------
        float
            Bulk-modulus standard deviation.
        """
        return float(self.standard_deviations[1])

    @property
    def sigma_bulk_modulus_derivative(self) -> float:
        """Return the one-sigma uncertainty of the bulk-modulus derivative.

        Returns
        -------
        float
            Standard deviation of ``dK/dP``.
        """
        return float(self.standard_deviations[2])

    def as_dict(self) -> dict[str, Any]:
        """Return a serializable representation of the uncertainty result.

        Returns
        -------
        dict
            Dictionary containing covariance, correlation and diagnostics.
        """
        return {
            "method": self.method,
            "state_order": list(_STATE_ORDER),
            "mean": self.mean.tolist(),
            "standard_deviations": self.standard_deviations.tolist(),
            "covariance": self.covariance.tolist(),
            "correlation": self.correlation.tolist(),
            "confidence_intervals": (
                None
                if self.confidence_intervals is None
                else self.confidence_intervals.tolist()
            ),
            "n_samples": self.n_samples,
            "n_successful": self.n_successful,
            "metadata": dict(self.metadata or {}),
        }


@dataclass(slots=True, frozen=True)
class EOSState:
    """Thermoelastic state evaluated from one fitted equation of state.

    Parameters
    ----------
    pressure : float
        Requested pressure in the scale defined by the EOS bulk modulus.
    volume : float
        Volume satisfying ``P(V) = pressure``.
    bulk_modulus : float
        Isothermal bulk modulus at ``volume``.
    bulk_modulus_derivative : float
        First pressure derivative of the bulk modulus at ``volume``.
    extrapolated : bool
        Whether the volume lies outside the sampled fit interval.
    uncertainty : EOSStateUncertainty, optional
        Propagated uncertainty of the pressure-dependent state.
    metadata : dict, optional
        Additional information about the root solution.
    """

    pressure: float
    volume: float
    bulk_modulus: float
    bulk_modulus_derivative: float
    extrapolated: bool = False
    uncertainty: EOSStateUncertainty | None = None
    metadata: dict[str, Any] | None = None

    @property
    def sigma_volume(self) -> float | None:
        """Return the propagated one-sigma volume uncertainty.

        Returns
        -------
        float or None
            Volume uncertainty when propagation was requested.
        """
        return None if self.uncertainty is None else self.uncertainty.sigma_volume

    @property
    def sigma_bulk_modulus(self) -> float | None:
        """Return the propagated one-sigma bulk-modulus uncertainty.

        Returns
        -------
        float or None
            Bulk-modulus uncertainty when propagation was requested.
        """
        return None if self.uncertainty is None else self.uncertainty.sigma_bulk_modulus

    @property
    def sigma_bulk_modulus_derivative(self) -> float | None:
        """Return the propagated one-sigma uncertainty of ``dK/dP``.

        Returns
        -------
        float or None
            Uncertainty of the bulk-modulus derivative when available.
        """
        if self.uncertainty is None:
            return None
        return self.uncertainty.sigma_bulk_modulus_derivative

    def as_dict(self) -> dict[str, Any]:
        """Return the state as a serializable dictionary.

        Returns
        -------
        dict
            Dictionary representation of the state.
        """
        return {
            "pressure": self.pressure,
            "volume": self.volume,
            "bulk_modulus": self.bulk_modulus,
            "bulk_modulus_derivative": self.bulk_modulus_derivative,
            "extrapolated": self.extrapolated,
            "uncertainty": None
            if self.uncertainty is None
            else self.uncertainty.as_dict(),
            "metadata": dict(self.metadata or {}),
        }
