# -*- coding: utf-8 -*-

"""Polynomial fitting, evaluation, and derivative utilities.

The helpers support scalar and batched datasets used to represent physical
properties as functions of volume, pressure, temperature, and related variables.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.polynomial.polynomial import polyfit as _polyfit
from numpy.polynomial.polynomial import polyval as _polyval
from numpy.polynomial import Polynomial
from numpy.typing import ArrayLike, NDArray


def polyfit(
    x: ArrayLike,
    y: ArrayLike,
    degree: int,
) -> NDArray[np.float64]:
    """
    Fit one or more datasets with a polynomial.

    Coefficients follow the convention used by
    :mod:`numpy.polynomial.polynomial`, where the coefficient of degree zero is
    stored first.

    Parameters
    ----------
    x : array_like
        Values of the independent variable.
    y : array_like
        Values of the dependent variable. If two-dimensional, each row is
        treated as an independent dataset.
    degree : int
        Degree of the polynomial.

    Returns
    -------
    ndarray
        Polynomial coefficients. For a single dataset, the shape is
        ``(degree + 1,)``. For multiple datasets, the shape is
        ``(n_datasets, degree + 1)``.
    """
    x_array = np.asarray(x, dtype=float)
    y_array = np.asarray(y, dtype=float)

    if y_array.ndim == 1:
        return _polyfit(x_array, y_array, degree)

    if y_array.ndim == 2:
        return _polyfit(x_array, y_array.T, degree).T

    raise ValueError("y must be a one- or two-dimensional array.")


def interpolate(
    x: ArrayLike | float,
    parameters: ArrayLike,
) -> float | NDArray[np.float64]:
    """
    Evaluate one or more polynomials.

    Parameters
    ----------
    x : scalar or array_like
        Point or points where the polynomial has to be evaluated.
    parameters : array_like
        Polynomial coefficients. Accepted shapes are ``(degree + 1,)`` for a
        single polynomial, ``(n_polynomials, degree + 1)`` for several
        polynomials, or ``(n_blocks, n_polynomials, degree + 1)`` for grouped
        polynomial datasets.

    Returns
    -------
    float or ndarray
        Interpolated value or values.
    """
    x_array = np.asarray(x, dtype=float)
    pars = np.asarray(parameters, dtype=float)

    if pars.ndim == 1:
        values = _polyval(x_array, pars)
    elif pars.ndim == 2:
        values = np.array([_polyval(x_array, row) for row in pars])
    elif pars.ndim == 3:
        values = np.array([[_polyval(x_array, row) for row in block] for block in pars])
    else:
        raise ValueError("parameters must have one, two, or three dimensions.")

    if np.isscalar(x) and pars.ndim == 1:
        return float(values)

    return values


def single_R_squared(
    x: ArrayLike,
    y: ArrayLike,
    parameters: ArrayLike,
) -> float:
    """
    Calculate the coefficient of determination for a polynomial fit.

    Parameters
    ----------
    x : array_like
        Values of the independent variable.
    y : array_like
        Observed values of the dependent variable.
    parameters : array_like
        Polynomial coefficients.

    Returns
    -------
    float
        Coefficient of determination.
    """
    x_array = np.asarray(x, dtype=float)
    y_array = np.asarray(y, dtype=float)
    y_hat = interpolate(x_array, parameters)
    y_bar = np.mean(y_array)

    ss_res = np.sum((y_array - y_hat) ** 2)
    ss_tot = np.sum((y_array - y_bar) ** 2)

    if ss_tot == 0.0:
        return 0.0

    return float(1.0 - ss_res / ss_tot)


def R_squared(
    x: ArrayLike,
    y: ArrayLike,
    parameters: ArrayLike,
) -> float | NDArray[np.float64]:
    """
    Calculate R-squared values for one or more polynomial fits.

    Parameters
    ----------
    x : array_like
        Values of the independent variable.
    y : array_like
        Observed values. If two-dimensional, each row is treated as an
        independent dataset.
    parameters : array_like
        Polynomial coefficients.

    Returns
    -------
    float or ndarray
        R-squared value or values.
    """
    y_array = np.asarray(y, dtype=float)
    pars = np.asarray(parameters, dtype=float)

    if y_array.ndim == 1:
        return single_R_squared(x, y_array, pars)

    if y_array.ndim == 2:
        return np.array(
            [single_R_squared(x, y_array[i], pars[i]) for i in range(y_array.shape[0])]
        )

    raise ValueError("y must be a one- or two-dimensional array.")


def find_polynomial_minimum(
    parameters: ArrayLike,
    addendum: Polynomial | None = None,
) -> NDArray[np.float64]:
    """
    Find the real minima of a polynomial.

    Parameters
    ----------
    parameters : array_like
        Polynomial coefficients.
    addendum : Polynomial or None, optional
        Additional polynomial contribution added to the first derivative before
        searching for critical points.

    Returns
    -------
    ndarray
        Real points where the polynomial has a local minimum.
    """
    polynomial = Polynomial(np.asarray(parameters, dtype=float))
    first_derivative = polynomial.deriv()

    if addendum is not None:
        first_derivative = first_derivative + addendum

    critical_points = first_derivative.roots()
    real_points = critical_points[np.isclose(critical_points.imag, 0.0)].real

    second_derivative = polynomial.deriv(2)
    return real_points[second_derivative(real_points) > 0.0]


class FittedPolynomial:
    """Represent a polynomial fitted in a scaled physical coordinate.

    The independent coordinate is transformed according to

    .. math::

        s = \frac{x - x_c}{x_s},

    where ``x_c`` is the coordinate center and ``x_s`` is a positive scale.
    Derivatives returned by :meth:`derivative` are evaluated with respect to
    the original physical coordinate ``x``.

    Parameters
    ----------
    parameters : array-like
        Polynomial coefficients in ascending order in the scaled coordinate.
    center : float
        Origin of the scaled coordinate.
    scale : float
        Positive coordinate scale.
    sampled_coordinates : array-like, optional
        Coordinates used to determine the fit.

    Raises
    ------
    ValueError
        If parameters, scaling values or sampled coordinates are invalid.
    """

    def __init__(
        self,
        parameters: ArrayLike,
        *,
        center: float = 0.0,
        scale: float = 1.0,
        sampled_coordinates: ArrayLike | None = None,
    ) -> None:
        coefficients = np.asarray(parameters, dtype=np.float64)
        if coefficients.ndim != 1 or coefficients.size < 1:
            raise ValueError("polynomial parameters must be a one-dimensional array")
        if not np.all(np.isfinite(coefficients)):
            raise ValueError("polynomial parameters must be finite")
        if not np.isfinite(center):
            raise ValueError("polynomial center must be finite")
        if not np.isfinite(scale) or scale <= 0.0:
            raise ValueError("polynomial scale must be positive")
        if sampled_coordinates is None:
            samples = None
        else:
            samples = np.asarray(sampled_coordinates, dtype=np.float64)
            if samples.ndim != 1 or samples.size < 2:
                raise ValueError("sampled_coordinates must be one-dimensional")
            if not np.all(np.isfinite(samples)):
                raise ValueError("sampled_coordinates must be finite")
        self.parameters = coefficients.copy()
        self.center = float(center)
        self.scale = float(scale)
        self.sampled_coordinates = None if samples is None else samples.copy()
        self._polynomial = Polynomial(self.parameters)

    def scaled_coordinate(self, x: ArrayLike | float) -> NDArray[np.float64]:
        """Transform physical coordinates to the scaled fit coordinate.

        Parameters
        ----------
        x : scalar or array-like
            Physical coordinate values.

        Returns
        -------
        ndarray
            Scaled coordinate values.
        """
        return (np.asarray(x, dtype=np.float64) - self.center) / self.scale

    def evaluate(self, x: ArrayLike | float) -> NDArray[np.float64]:
        """Evaluate the polynomial in the physical coordinate.

        Parameters
        ----------
        x : scalar or array-like
            Physical coordinate values.

        Returns
        -------
        ndarray
            Polynomial values.
        """
        return np.asarray(self._polynomial(self.scaled_coordinate(x)), dtype=np.float64)

    def derivative(
        self,
        x: ArrayLike | float,
        order: int = 1,
    ) -> NDArray[np.float64]:
        """Evaluate a derivative with respect to the physical coordinate.

        Parameters
        ----------
        x : scalar or array-like
            Physical coordinate values.
        order : int, optional
            Derivative order.

        Returns
        -------
        ndarray
            Derivative values.

        Raises
        ------
        ValueError
            If ``order`` is negative.
        """
        if order < 0:
            raise ValueError("derivative order must be non-negative")
        derivative = self._polynomial.deriv(order)
        values = derivative(self.scaled_coordinate(x)) / self.scale**order
        return np.asarray(values, dtype=np.float64)

    def stationary_points(
        self, derivative_addendum: float = 0.0
    ) -> NDArray[np.float64]:
        """Return real stationary points after adding a constant derivative.

        The roots satisfy

        .. math::

            \frac{df}{dx} + a = 0,

        where ``a`` is ``derivative_addendum``.

        Parameters
        ----------
        derivative_addendum : float, optional
            Constant added to the derivative in physical-coordinate units.

        Returns
        -------
        ndarray
            Real stationary coordinates.
        """
        derivative = self._polynomial.deriv(1) + (
            float(derivative_addendum) * self.scale
        )
        roots = derivative.roots()
        real = np.real(roots[np.isclose(np.imag(roots), 0.0)])
        return np.asarray(self.center + self.scale * real, dtype=np.float64)

    def local_minima(self, derivative_addendum: float = 0.0) -> NDArray[np.float64]:
        """Return stationary points with positive curvature.

        Parameters
        ----------
        derivative_addendum : float, optional
            Constant added to the physical derivative.

        Returns
        -------
        ndarray
            Coordinates of local minima.
        """
        points = self.stationary_points(derivative_addendum)
        if points.size == 0:
            return points
        curvature = self.derivative(points, 2)
        mask = np.isfinite(points) & np.isfinite(curvature) & (curvature > 0.0)
        return np.asarray(points[mask], dtype=np.float64)


def fit_polynomial(
    x: ArrayLike,
    y: ArrayLike,
    degree: int,
    *,
    scale_coordinate: bool = False,
):
    """Fit a polynomial and return diagnostics with a reusable model.

    Parameters
    ----------
    x, y : array-like
        Coordinates and observed values.
    degree : int
        Polynomial degree.
    scale_coordinate : bool, optional
        Transform coordinates to ``(x - center) / scale`` before fitting.

    Returns
    -------
    tuple
        ``(FitResult, FittedPolynomial or None)``.

    Raises
    ------
    ValueError
        This function reports invalid inputs through ``FitResult`` rather than
        raising, except for internal programming errors.
    """
    from .fitting import (
        FitQuality,
        FitResult,
        FitStatus,
        residual_metrics,
        validate_xy,
    )

    metadata: dict[str, Any] = {
        "model": "scaled_polynomial" if scale_coordinate else "polynomial",
        "degree": int(degree),
        "parameter_order": [f"c{i}" for i in range(degree + 1)],
    }
    try:
        x_array, y_array = validate_xy(
            np.asarray(x, dtype=np.float64),
            np.asarray(y, dtype=np.float64),
        )
    except ValueError as exc:
        return (
            FitResult.failed(
                str(exc), status=FitStatus.INVALID_INPUT, metadata=metadata
            ),
            None,
        )
    if degree < 1:
        return (
            FitResult.failed(
                "polynomial degree must be at least one",
                status=FitStatus.INVALID_INPUT,
                metadata=metadata,
            ),
            None,
        )
    if x_array.size < degree + 1:
        return (
            FitResult.failed(
                "the number of fit points is smaller than the number of polynomial coefficients",
                status=FitStatus.INVALID_INPUT,
                n_points=int(x_array.size),
                n_parameters=degree + 1,
                metadata=metadata,
            ),
            None,
        )

    if scale_coordinate:
        xmin = float(np.min(x_array))
        xmax = float(np.max(x_array))
        center = 0.5 * (xmin + xmax)
        scale = 0.5 * (xmax - xmin)
        if not np.isfinite(scale) or scale <= 0.0:
            return (
                FitResult.failed(
                    "sampled coordinate interval must be non-degenerate",
                    status=FitStatus.INVALID_INPUT,
                    metadata=metadata,
                ),
                None,
            )
        fit_x = (x_array - center) / scale
        metadata.update(
            {
                "coordinate_center": center,
                "coordinate_scale": scale,
                "sampled_coordinate_min": xmin,
                "sampled_coordinate_max": xmax,
                "coordinate": "s=(x-center)/scale",
            }
        )
    else:
        center = 0.0
        scale = 1.0
        fit_x = x_array

    with np.errstate(all="ignore"):
        coefficients, full_output = np.polynomial.polynomial.polyfit(
            fit_x, y_array, degree, full=True
        )
    residual_sum, rank, singular_values, rcond = full_output
    coefficients = np.asarray(coefficients, dtype=np.float64)
    if not np.all(np.isfinite(coefficients)):
        return (
            FitResult.failed(
                "the polynomial fit returned non-finite coefficients",
                n_points=int(x_array.size),
                n_parameters=degree + 1,
                metadata=metadata,
            ),
            None,
        )

    model = FittedPolynomial(
        coefficients,
        center=center,
        scale=scale,
        sampled_coordinates=x_array,
    )
    fitted = model.evaluate(x_array)
    metrics = residual_metrics(y_array, fitted)
    dof = max(int(x_array.size) - (degree + 1), 0)
    variance = (
        float(residual_sum[0] / dof)
        if residual_sum.size and dof > 0
        else float(metrics["rmse"] ** 2)
    )
    try:
        vandermonde = np.polynomial.polynomial.polyvander(fit_x, degree)
        covariance = variance * np.linalg.pinv(vandermonde.T @ vandermonde)
        errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
        condition_number = float(np.linalg.cond(vandermonde))
    except np.linalg.LinAlgError:
        covariance = None
        errors = None
        condition_number = np.inf

    warnings_: list[str] = []
    quality = FitQuality.GOOD
    if rank < degree + 1:
        warnings_.append("polynomial design matrix is rank deficient")
        quality = FitQuality.POOR
    if np.isfinite(condition_number) and condition_number > 1.0e12:
        warnings_.append("polynomial design matrix is ill-conditioned")
        quality = FitQuality.POOR
    if np.isfinite(metrics["r_squared"]) and metrics["r_squared"] < 0.95:
        quality = FitQuality.POOR

    result = FitResult(
        success=True,
        status=FitStatus.SUCCESS,
        quality=quality,
        message=(
            "fit converged"
            if quality is FitQuality.GOOD
            else "fit converged with diagnostic warnings"
        ),
        parameters=coefficients,
        covariance=covariance,
        errors=errors,
        fitted=fitted,
        residuals=y_array - fitted,
        rmse=metrics["rmse"],
        mae=metrics["mae"],
        max_abs_error=metrics["max_abs_error"],
        r_squared=metrics["r_squared"],
        n_points=int(x_array.size),
        n_parameters=degree + 1,
        dof=dof,
        condition_number=condition_number,
        warnings=warnings_,
        metadata={
            **metadata,
            "rank": int(rank),
            "rcond": float(rcond),
            "singular_values": singular_values.tolist(),
        },
    )
    return result, model


def fit_polynomial_result(
    x: ArrayLike,
    y: ArrayLike,
    degree: int,
    *,
    scale_coordinate: bool = False,
):
    """Fit a polynomial and return only the structured fit diagnostics.

    Parameters
    ----------
    x, y : array_like
        Coordinates and observed values.
    degree : int
        Polynomial degree.
    scale_coordinate : bool, optional
        Transform coordinates to ``(x - center) / scale`` before fitting.

    Returns
    -------
    FitResult
        Polynomial coefficients and fit diagnostics.

    Notes
    -----
    This convenience wrapper is used by workflows that need polynomial
    coefficients and diagnostics but do not need the reusable
    :class:`FittedPolynomial` model returned by :func:`fit_polynomial`.
    """
    fit, _ = fit_polynomial(
        x,
        y,
        degree,
        scale_coordinate=scale_coordinate,
    )
    return fit
