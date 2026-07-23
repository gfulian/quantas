# -*- coding: utf-8 -*-

"""Validation and diagnostics shared by numerical fitting algorithms."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TypeAlias

import numpy as np

from .options import CovarianceScaling
from .result import FitQuality

ArrayLike: TypeAlias = np.ndarray | Sequence[float]


def validate_xy(x: ArrayLike, y: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    """Validate one- or multi-coordinate fitting data and observations.

    Parameters
    ----------
    x, y : array-like
        Input coordinates and observed values.

    Returns
    -------
    tuple of ndarray
        Converted arrays with ``float64`` dtype.

    Raises
    ------
    ValueError
        If ``y`` is not a finite vector, ``x`` is not a finite vector or
        coordinate matrix, their observation counts differ, or fewer than two
        observations are provided.
    """
    x_array = np.asarray(x, dtype=np.float64)
    y_array = np.asarray(y, dtype=np.float64)
    if x_array.ndim not in {1, 2} or y_array.ndim != 1:
        raise ValueError(
            "fit x must be a vector or coordinate matrix and y must be a vector"
        )
    if x_array.shape[-1] != y_array.size:
        raise ValueError("fit x and y must contain the same number of observations")
    if y_array.size < 2:
        raise ValueError("at least two points are required for a fit")
    if not np.all(np.isfinite(x_array)) or not np.all(np.isfinite(y_array)):
        raise ValueError("fit input arrays must contain only finite values")
    return x_array, y_array


def residual_metrics(observed: ArrayLike, fitted: ArrayLike) -> dict[str, float]:
    """Calculate residual statistics for fitted data.

    Parameters
    ----------
    observed, fitted : array-like
        Observed and calculated values.

    Returns
    -------
    dict
        ``rmse``, ``mae``, ``max_abs_error`` and ``r_squared``.

    Raises
    ------
    ValueError
        If the arrays cannot be compared element by element.
    """
    y_obs, y_fit = validate_xy(observed, fitted)
    residuals = y_obs - y_fit
    ss_res = float(np.sum(residuals**2))
    centred = y_obs - np.mean(y_obs)
    ss_tot = float(np.sum(centred**2))
    r_squared = float(1.0 - ss_res / ss_tot) if ss_tot > 0.0 else np.nan
    return {
        "rmse": float(np.sqrt(np.mean(residuals**2))),
        "mae": float(np.mean(np.abs(residuals))),
        "max_abs_error": float(np.max(np.abs(residuals))),
        "r_squared": r_squared,
    }


def covariance_errors(
    covariance: np.ndarray | None,
) -> tuple[np.ndarray | None, float | None, bool]:
    """Calculate parameter errors and covariance diagnostics.

    Parameters
    ----------
    covariance : ndarray or None
        Parameter covariance matrix.

    Returns
    -------
    tuple
        Parameter errors, condition number and a flag indicating whether the
        covariance matrix is numerically usable.
    """
    if covariance is None:
        return None, None, False
    matrix = np.asarray(covariance, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        return None, None, False
    if not np.all(np.isfinite(matrix)):
        return np.full(matrix.shape[0], np.inf), np.inf, False
    errors = np.sqrt(np.maximum(np.diag(matrix), 0.0))
    condition_number = float(np.linalg.cond(matrix))
    usable = bool(np.all(np.isfinite(errors)) and np.isfinite(condition_number))
    return errors, condition_number, usable


def covariance_correlation(covariance: np.ndarray | None) -> np.ndarray | None:
    """Convert a covariance matrix to a parameter correlation matrix.

    Parameters
    ----------
    covariance : ndarray or None
        Square covariance matrix.

    Returns
    -------
    ndarray or None
        Correlation matrix. Rows or columns with zero variance, such as fixed
        parameters, are represented by zeros; positive-variance diagonal
        elements are exactly one.

    Raises
    ------
    ValueError
        If the covariance matrix is not finite and square.
    """
    if covariance is None:
        return None
    matrix = np.asarray(covariance, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("covariance matrix must be square")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("covariance matrix must contain finite values")
    errors = np.sqrt(np.maximum(np.diag(matrix), 0.0))
    denominator = np.outer(errors, errors)
    correlation = np.zeros_like(matrix)
    np.divide(matrix, denominator, out=correlation, where=denominator > 0.0)
    positive = errors > 0.0
    diagonal = np.diag_indices_from(correlation)
    correlation[diagonal] = positive.astype(np.float64)
    return np.clip(correlation, -1.0, 1.0)


def standardized_residuals(
    residuals: ArrayLike,
    sigma: ArrayLike,
) -> np.ndarray:
    """Return residuals normalized by strictly positive standard deviations."""
    residual_array, sigma_array = validate_xy(residuals, sigma)
    if np.any(sigma_array <= 0.0):
        raise ValueError("standard deviations must be strictly positive")
    return residual_array / sigma_array


def assess_fit_quality(
    *,
    success: bool,
    covariance_usable: bool,
    warnings_: Sequence[str],
    r_squared: float | None,
    max_condition_number: float = 1.0e12,
    condition_number: float | None = None,
) -> FitQuality:
    """Assign a qualitative label to a fit result.

    Parameters
    ----------
    success : bool
        Whether the optimizer converged.
    covariance_usable : bool
        Whether the covariance matrix is usable.
    warnings_ : sequence of str
        Warning messages raised by the optimizer.
    r_squared : float or None
        Coefficient of determination.
    max_condition_number : float, optional
        Maximum accepted condition number.
    condition_number : float or None, optional
        Condition number of the fitted problem.

    Returns
    -------
    FitQuality
        Qualitative fit assessment.
    """
    if not success:
        return FitQuality.FAILED
    if warnings_ or not covariance_usable:
        return FitQuality.POOR
    if condition_number is not None and condition_number > max_condition_number:
        return FitQuality.POOR
    if r_squared is not None and np.isfinite(r_squared) and r_squared < 0.95:
        return FitQuality.POOR
    return FitQuality.GOOD


def covariance_scale_factor(
    policy: CovarianceScaling,
    reduced_chi_square: float | None,
    *,
    weighted: bool,
) -> float:
    """Return the multiplicative scale for a backend covariance matrix.

    Parameters
    ----------
    policy : CovarianceScaling
        Requested covariance convention.
    reduced_chi_square : float or None
        Weighted objective divided by residual degrees of freedom.
    weighted : bool
        Whether supplied observation uncertainties entered the objective.

    Returns
    -------
    float
        Scale applied to an absolute covariance matrix. Unweighted OLS
        covariance is assumed to have been scaled by its backend already.

    Raises
    ------
    ValueError
        If the policy is unsupported or reduced chi-square is non-finite.
    """
    if not weighted or reduced_chi_square is None:
        return 1.0
    reduced = float(reduced_chi_square)
    if not np.isfinite(reduced) or reduced < 0.0:
        raise ValueError("reduced chi-square must be finite and non-negative")
    policy = CovarianceScaling(policy)
    if policy is CovarianceScaling.ABSOLUTE:
        return 1.0
    if policy is CovarianceScaling.REDUCED_CHI_SQUARE:
        return reduced
    if policy is CovarianceScaling.INFLATE_ONLY:
        return max(1.0, reduced)
    raise ValueError(f"unsupported covariance-scaling policy: {policy}")


def parameters_at_bounds(
    parameters: ArrayLike,
    bounds: tuple[ArrayLike, ArrayLike],
    *,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-12,
) -> tuple[bool, ...]:
    """Return flags identifying fitted parameters at either declared bound.

    Parameters
    ----------
    parameters : array-like
        Optimized parameter vector.
    bounds : tuple of array-like
        Lower and upper bounds in the same order.
    rtol, atol : float, optional
        Relative and absolute comparison tolerances.

    Returns
    -------
    tuple of bool
        One flag per parameter.

    Raises
    ------
    ValueError
        If parameter and bound vectors are incompatible.
    """
    values = np.asarray(parameters, dtype=np.float64)
    lower = np.asarray(bounds[0], dtype=np.float64)
    upper = np.asarray(bounds[1], dtype=np.float64)
    if values.ndim != 1 or lower.shape != values.shape or upper.shape != values.shape:
        raise ValueError("parameters and bounds must be aligned vectors")
    return tuple(
        bool(
            np.isclose(value, lower_value, rtol=rtol, atol=atol)
            or np.isclose(value, upper_value, rtol=rtol, atol=atol)
        )
        for value, lower_value, upper_value in zip(values, lower, upper, strict=True)
    )
