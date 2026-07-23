# -*- coding: utf-8 -*-

"""Numerical and polynomial derivative utilities.

The functions defined here centralize common derivative operations used by
scientific workflows.  They intentionally operate on plain arrays and do not
know about thermodynamic, elastic, or seismic semantics; callers are
responsible for unit conversions and physical interpretation.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from typing import Literal, cast


class LenError(ValueError):
    """Exception raised when too few points are available for a derivative."""


def _as_matching_1d_arrays(
    x: ArrayLike,
    y: ArrayLike,
    *,
    min_points: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return validated one-dimensional coordinate and value arrays.

    Parameters
    ----------
    x, y : array_like
        Coordinate and sampled values.
    min_points : int
        Minimum required number of points.

    Returns
    -------
    tuple of ndarray
        Validated ``x`` and ``y`` arrays.

    Raises
    ------
    LenError
        If too few points are provided.
    ValueError
        If the arrays are not one-dimensional or do not have the same shape.
    """
    x_array = np.asarray(x, dtype=np.float64)
    y_array = np.asarray(y, dtype=np.float64)

    if x_array.ndim != 1 or y_array.ndim != 1:
        raise ValueError("x and y must be one-dimensional arrays.")
    if x_array.shape != y_array.shape:
        raise ValueError("x and y must have the same shape.")
    if x_array.size < min_points:
        raise LenError(
            f"At least {min_points} points are required to calculate a derivative."
        )
    return x_array, y_array


def derivative(
    x: ArrayLike,
    y: ArrayLike,
    edge_order: int = 1,
) -> NDArray[np.float64]:
    """Calculate the numerical derivative of a one-dimensional dataset.

    Parameters
    ----------
    x : array_like
        Values of the independent variable.
    y : array_like
        Values of the dependent variable.
    edge_order : int, optional
        Gradient accuracy at the boundaries. This value is passed to
        :func:`numpy.gradient`.

    Returns
    -------
    ndarray
        Numerical derivative ``dy/dx``.

    Raises
    ------
    LenError
        If fewer than two points are provided.
    ValueError
        If ``x`` and ``y`` do not have the same shape, or if they are not
        one-dimensional arrays.
    """
    x_array, y_array = _as_matching_1d_arrays(x, y, min_points=2)
    if edge_order not in (1, 2):
        raise ValueError("edge_order must be 1 or 2")
    return np.asarray(
        np.gradient(
            y_array,
            x_array,
            edge_order=cast(Literal[1, 2], edge_order),
        )
    )


def relative_derivative(
    x: ArrayLike,
    y: ArrayLike,
    edge_order: int = 1,
) -> NDArray[np.float64]:
    """Return the logarithmic-style derivative ``(1/y) dy/dx``.

    Parameters
    ----------
    x : array_like
        Values of the independent variable.
    y : array_like
        Values of the dependent variable.
    edge_order : int, optional
        Gradient accuracy at the boundaries. This value is passed to
        :func:`numpy.gradient`.

    Returns
    -------
    ndarray
        Relative derivative. Points with non-finite or zero ``y`` are returned
        as ``NaN``.

    Raises
    ------
    LenError
        If fewer than two points are provided.
    ValueError
        If ``x`` and ``y`` do not have the same shape, or if they are not
        one-dimensional arrays.
    """
    y_array = np.asarray(y, dtype=np.float64)
    values = derivative(x, y_array, edge_order=edge_order)
    out = np.full_like(values, np.nan, dtype=np.float64)
    mask = np.isfinite(values) & np.isfinite(y_array) & (y_array != 0.0)
    out[mask] = values[mask] / y_array[mask]
    return out


def polynomial_derivative(
    x: ArrayLike,
    y: ArrayLike,
    x_eval: ArrayLike | float,
    degree: int,
    *,
    order: int = 1,
) -> NDArray[np.float64]:
    """Fit a polynomial and evaluate its derivative.

    Parameters
    ----------
    x, y : array_like
        One-dimensional sampled coordinates and values. Non-finite points are
        discarded before fitting.
    x_eval : scalar or array_like
        Coordinates at which the derivative is evaluated.
    degree : int
        Polynomial degree.
    order : int, optional
        Derivative order.

    Returns
    -------
    ndarray
        Derivative values at ``x_eval``.

    Raises
    ------
    LenError
        If too few finite points are available for the requested polynomial.
    ValueError
        If inputs are malformed or ``degree``/``order`` is invalid.
    """
    if degree < 0:
        raise ValueError("polynomial degree must be non-negative")
    if order < 0:
        raise ValueError("derivative order must be non-negative")
    x_array, y_array = _as_matching_1d_arrays(x, y, min_points=degree + 1)
    finite = np.isfinite(x_array) & np.isfinite(y_array)
    if np.count_nonzero(finite) < degree + 1:
        raise LenError(
            "too few finite points are available for the requested polynomial"
        )
    coefficients = np.polynomial.polynomial.polyfit(
        x_array[finite],
        y_array[finite],
        deg=int(degree),
    )
    derivative_coefficients = np.polynomial.polynomial.polyder(
        coefficients,
        m=int(order),
    )
    return np.asarray(
        np.polynomial.polynomial.polyval(
            np.asarray(x_eval, dtype=np.float64),
            derivative_coefficients,
        ),
        dtype=np.float64,
    )


def polynomial_derivative_from_coefficients(
    coefficients: ArrayLike,
    x_eval: ArrayLike | float,
    *,
    order: int = 1,
    axis: int = 0,
) -> NDArray[np.float64]:
    """Evaluate a derivative from existing polynomial coefficients.

    Parameters
    ----------
    coefficients : array_like
        Polynomial coefficients in ascending order along ``axis``.
    x_eval : scalar or array_like
        Coordinates at which the derivative is evaluated.
    order : int, optional
        Derivative order.
    axis : int, optional
        Axis storing the polynomial coefficients.

    Returns
    -------
    ndarray
        Derivative values evaluated at ``x_eval``.

    Raises
    ------
    ValueError
        If ``order`` is negative or coefficients are malformed.
    """
    if order < 0:
        raise ValueError("derivative order must be non-negative")
    coefficient_array = np.asarray(coefficients, dtype=np.float64)
    if coefficient_array.ndim == 0:
        raise ValueError("polynomial coefficients must be at least one-dimensional")
    derivative_coefficients = np.polynomial.polynomial.polyder(
        coefficient_array,
        m=int(order),
        axis=axis,
    )
    return np.asarray(
        np.polynomial.polynomial.polyval(
            np.asarray(x_eval, dtype=np.float64),
            derivative_coefficients,
        ),
        dtype=np.float64,
    )
