# -*- coding: utf-8 -*-

"""Polynomial interpolation utilities for QHA volume-dependent data.

The functions in this module fit and evaluate scalar datasets sampled on a
common volume grid.  They are designed for quantities such as phonon
frequencies, static energies and harmonic thermodynamic properties whose QHA
workflow requires repeated interpolation along the volume coordinate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from itertools import product
from typing import Any, Mapping, Sequence

import numpy as np

from quantas.core.math.fitting import FitQuality, FitResult
from quantas.core.math.polynomials import fit_polynomial_result as polynomial_fit

ArrayLike = np.ndarray | Sequence[float]


class InterpolationStatus(str, Enum):
    """Execution status of an interpolation operation.

    Attributes
    ----------
    SUCCESS
        All requested datasets were fitted or evaluated successfully.
    PARTIAL
        At least one dataset failed, while at least one result is usable.
    FAILED
        No usable interpolation result was produced.
    INVALID_INPUT
        The input arrays are not suitable for interpolation.
    """

    SUCCESS = "success"
    PARTIAL = "partial"
    FAILED = "failed"
    INVALID_INPUT = "invalid_input"


@dataclass(slots=True)
class PolynomialSeriesFit:
    """Polynomial fit of one or more datasets sampled on a common grid.

    Parameters
    ----------
    success : bool
        Whether all datasets were fitted successfully.
    status : InterpolationStatus
        Execution status of the series fit.
    degree : int
        Polynomial degree used for the fit.
    coefficients : ndarray or None
        Polynomial coefficients in ascending order.  The leading dimensions
        match the input datasets, while the last dimension stores the
        coefficients.
    fits : dict
        Individual fit diagnostics keyed by dataset index.
    failed_indices : list of tuple
        Dataset indices whose fit failed.
    warnings : list of str
        Non-fatal diagnostic messages.
    metadata : dict
        Additional caller-defined diagnostic information.
    """

    success: bool
    status: InterpolationStatus
    degree: int
    coefficients: np.ndarray | None = None
    fits: dict[tuple[int, ...], FitResult] = field(default_factory=dict)
    failed_indices: list[tuple[int, ...]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def failed(
        cls,
        message: str,
        *,
        degree: int,
        status: InterpolationStatus = InterpolationStatus.FAILED,
        metadata: Mapping[str, Any] | None = None,
    ) -> PolynomialSeriesFit:
        """Create a failed series-fit result.

        Parameters
        ----------
        message : str
            Explanation of the failure.
        degree : int
            Requested polynomial degree.
        status : InterpolationStatus, optional
            Failure category.
        metadata : mapping, optional
            Additional diagnostic information.

        Returns
        -------
        PolynomialSeriesFit
            Failed series-fit result.
        """
        return cls(
            success=False,
            status=status,
            degree=int(degree),
            warnings=[message],
            metadata=dict(metadata or {}),
        )

    def as_dict(self) -> dict[str, Any]:
        """Return the series fit as a serializable dictionary.

        Returns
        -------
        dict
            Dictionary representation of the series fit.
        """
        return {
            "success": self.success,
            "status": self.status.value,
            "degree": self.degree,
            "coefficients": None
            if self.coefficients is None
            else self.coefficients.tolist(),
            "fits": {str(key): value.as_dict() for key, value in self.fits.items()},
            "failed_indices": [list(index) for index in self.failed_indices],
            "warnings": list(self.warnings),
            "metadata": dict(self.metadata),
        }


@dataclass(slots=True)
class InterpolationResult:
    """Values obtained by evaluating fitted interpolation functions.

    Parameters
    ----------
    success : bool
        Whether the interpolation produced finite values.
    status : InterpolationStatus
        Execution status of the interpolation.
    values : ndarray or None
        Interpolated values.
    message : str
        Human-readable status message.
    warnings : list of str
        Non-fatal diagnostic messages.
    metadata : dict
        Additional caller-defined diagnostic information.
    """

    success: bool
    status: InterpolationStatus
    values: np.ndarray | None = None
    message: str = ""
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def failed(
        cls,
        message: str,
        *,
        status: InterpolationStatus = InterpolationStatus.FAILED,
        metadata: Mapping[str, Any] | None = None,
    ) -> InterpolationResult:
        """Create a failed interpolation result.

        Parameters
        ----------
        message : str
            Explanation of the failure.
        status : InterpolationStatus, optional
            Failure category.
        metadata : mapping, optional
            Additional diagnostic information.

        Returns
        -------
        InterpolationResult
            Failed interpolation result.
        """
        return cls(
            success=False,
            status=status,
            message=message,
            warnings=[message],
            metadata=dict(metadata or {}),
        )

    def as_dict(self) -> dict[str, Any]:
        """Return the interpolation result as a serializable dictionary.

        Returns
        -------
        dict
            Dictionary representation of the interpolation result.
        """
        return {
            "success": self.success,
            "status": self.status.value,
            "values": None if self.values is None else self.values.tolist(),
            "message": self.message,
            "warnings": list(self.warnings),
            "metadata": dict(self.metadata),
        }


def validate_volume_grid(volume: ArrayLike) -> np.ndarray:
    """Validate a one-dimensional volume grid.

    Parameters
    ----------
    volume : array-like
        Volume values.

    Returns
    -------
    ndarray
        Finite volume values as a one-dimensional ``float64`` array.

    Raises
    ------
    ValueError
        If the volume grid is not one-dimensional, contains fewer than two
        points, contains non-finite values or contains duplicates.
    """
    volume_array = np.asarray(volume, dtype=np.float64)
    if volume_array.ndim != 1:
        raise ValueError("volume grid must be one-dimensional")
    if volume_array.size < 2:
        raise ValueError("at least two volume points are required")
    if not np.all(np.isfinite(volume_array)):
        raise ValueError("volume grid must contain only finite values")
    if np.unique(volume_array).size != volume_array.size:
        raise ValueError("volume grid must not contain duplicate values")
    return volume_array


def move_axis_to_last(values: ArrayLike, axis: int) -> np.ndarray:
    """Move an interpolation axis to the last position.

    Parameters
    ----------
    values : array-like
        Input data array.
    axis : int
        Axis corresponding to the sampled volume coordinate.

    Returns
    -------
    ndarray
        Input values with the interpolation axis moved to the last position.

    Raises
    ------
    ValueError
        If the input values are empty or the requested axis is invalid.
    """
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 0:
        raise ValueError("values must have at least one dimension")
    if not np.all(np.isfinite(array)):
        raise ValueError("values must contain only finite numbers")
    try:
        return np.moveaxis(array, axis, -1)
    except np.exceptions.AxisError as exc:
        raise ValueError("interpolation axis is out of bounds") from exc


def fit_polynomial_series(
    volume: ArrayLike,
    values: ArrayLike,
    degree: int,
    *,
    axis: int = -1,
    metadata: Mapping[str, Any] | None = None,
) -> PolynomialSeriesFit:
    """Fit polynomial functions to datasets sampled on a volume grid.

    Parameters
    ----------
    volume : array-like
        Volume values used as the independent variable.
    values : array-like
        Values sampled at each volume.  The volume axis is selected with
        ``axis``; all remaining axes define independent datasets.
    degree : int
        Polynomial degree.
    axis : int, optional
        Axis of ``values`` corresponding to ``volume``.
    metadata : mapping, optional
        Additional diagnostic information stored in the returned result.

    Returns
    -------
    PolynomialSeriesFit
        Polynomial coefficients and per-dataset fit diagnostics.
    """
    if degree < 0:
        return PolynomialSeriesFit.failed(
            "polynomial degree must be non-negative",
            degree=degree,
            status=InterpolationStatus.INVALID_INPUT,
            metadata=metadata,
        )
    try:
        volume_array = validate_volume_grid(volume)
        values_array = move_axis_to_last(values, axis)
    except ValueError as exc:
        return PolynomialSeriesFit.failed(
            str(exc),
            degree=degree,
            status=InterpolationStatus.INVALID_INPUT,
            metadata=metadata,
        )
    if values_array.shape[-1] != volume_array.size:
        return PolynomialSeriesFit.failed(
            "the interpolation axis length must match the number of volumes",
            degree=degree,
            status=InterpolationStatus.INVALID_INPUT,
            metadata=metadata,
        )
    if volume_array.size < degree + 1:
        return PolynomialSeriesFit.failed(
            "the number of volume points is smaller than the number of polynomial coefficients",
            degree=degree,
            status=InterpolationStatus.INVALID_INPUT,
            metadata=metadata,
        )

    dataset_shape = values_array.shape[:-1]
    coefficients = np.full(dataset_shape + (degree + 1,), np.nan, dtype=np.float64)
    fits: dict[tuple[int, ...], FitResult] = {}
    failed_indices: list[tuple[int, ...]] = []
    warnings_: list[str] = []

    if dataset_shape == ():
        iterator: list[tuple[int, ...]] = [()]
    else:
        iterator = list(
            product(*(range(int(length)) for length in dataset_shape))
        )

    for index in iterator:
        y = values_array[index] if index else values_array
        fit = polynomial_fit(volume_array, y, degree)
        fits[index] = fit
        if fit.success and fit.parameters is not None:
            coefficients[index] = fit.parameters
            if fit.quality is not FitQuality.GOOD:
                warnings_.append(
                    f"polynomial fit at index {index} has quality {fit.quality.value}"
                )
        else:
            failed_indices.append(index)
            warnings_.append(f"polynomial fit at index {index} failed: {fit.message}")

    if len(failed_indices) == 0:
        status = InterpolationStatus.SUCCESS
    elif len(failed_indices) == len(iterator):
        status = InterpolationStatus.FAILED
    else:
        status = InterpolationStatus.PARTIAL

    return PolynomialSeriesFit(
        success=status is InterpolationStatus.SUCCESS,
        status=status,
        degree=int(degree),
        coefficients=coefficients,
        fits=fits,
        failed_indices=failed_indices,
        warnings=warnings_,
        metadata={
            "axis": int(axis),
            "dataset_shape": dataset_shape,
            **dict(metadata or {}),
        },
    )


def evaluate_polynomial_series(
    points: ArrayLike | float, coefficients: ArrayLike
) -> InterpolationResult:
    """Evaluate one or more polynomial interpolation functions.

    Parameters
    ----------
    points : scalar or array-like
        Points where the polynomials are evaluated.
    coefficients : array-like
        Polynomial coefficients in ascending order.  The last axis stores the
        coefficients; all preceding axes define independent datasets.

    Returns
    -------
    InterpolationResult
        Interpolated values.  If ``points`` is an array, its shape is appended
        to the dataset shape.
    """
    coefficient_array = np.asarray(coefficients, dtype=np.float64)
    if coefficient_array.ndim == 0:
        return InterpolationResult.failed(
            "coefficient array must have at least one dimension",
            status=InterpolationStatus.INVALID_INPUT,
        )
    if not np.all(np.isfinite(coefficient_array)):
        return InterpolationResult.failed(
            "coefficient array contains non-finite values",
            status=InterpolationStatus.INVALID_INPUT,
        )
    point_array = np.asarray(points, dtype=np.float64)
    if not np.all(np.isfinite(point_array)):
        return InterpolationResult.failed(
            "interpolation points must contain only finite values",
            status=InterpolationStatus.INVALID_INPUT,
        )

    values = np.polynomial.polynomial.polyval(point_array, coefficient_array.T)
    if coefficient_array.ndim == 1 and np.isscalar(points):
        values = np.asarray(float(values), dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(values)):
        return InterpolationResult.failed("interpolation produced non-finite values")
    return InterpolationResult(
        success=True,
        status=InterpolationStatus.SUCCESS,
        values=values,
        message="interpolation completed",
        metadata={
            "points_shape": point_array.shape,
            "coefficient_shape": coefficient_array.shape,
        },
    )


def interpolate_polynomial_series(
    volume: ArrayLike,
    values: ArrayLike,
    points: ArrayLike | float,
    degree: int,
    *,
    axis: int = -1,
) -> tuple[PolynomialSeriesFit, InterpolationResult]:
    """Fit and evaluate polynomial interpolation functions.

    Parameters
    ----------
    volume : array-like
        Volume values used as the independent variable.
    values : array-like
        Values sampled at each volume.
    points : scalar or array-like
        Points where the fitted polynomials are evaluated.
    degree : int
        Polynomial degree.
    axis : int, optional
        Axis of ``values`` corresponding to ``volume``.

    Returns
    -------
    tuple
        Series fit result and interpolation result.
    """
    fit = fit_polynomial_series(volume, values, degree, axis=axis)
    if fit.coefficients is None or fit.status is InterpolationStatus.FAILED:
        return fit, InterpolationResult.failed("polynomial interpolation fit failed")
    result = evaluate_polynomial_series(points, fit.coefficients)
    return fit, result


def interpolate_static_energy(
    volume: ArrayLike, static_energy: ArrayLike, points: ArrayLike | float, degree: int
) -> tuple[FitResult, InterpolationResult]:
    """Interpolate static energy sampled as a function of volume.

    Parameters
    ----------
    volume : array-like
        Volume values.
    static_energy : array-like
        Static-energy values.
    points : scalar or array-like
        Volumes where static energy is evaluated.
    degree : int
        Polynomial degree.

    Returns
    -------
    tuple
        Polynomial fit diagnostics and interpolation result.
    """
    fit = polynomial_fit(volume, static_energy, degree)
    if not fit.success or fit.parameters is None:
        return fit, InterpolationResult.failed(fit.message)
    return fit, evaluate_polynomial_series(points, fit.parameters)


def interpolate_frequency_grid(
    volume: ArrayLike,
    frequencies: ArrayLike,
    points: ArrayLike | float,
    degree: int,
    *,
    axis: int = -1,
) -> tuple[PolynomialSeriesFit, InterpolationResult]:
    """Interpolate phonon frequencies along the volume coordinate.

    Parameters
    ----------
    volume : array-like
        Volume values.
    frequencies : array-like
        Frequency values.  Typical QHA data have shape
        ``(n_qpoints, n_bands, n_volumes)``.
    points : scalar or array-like
        Volumes where frequencies are evaluated.
    degree : int
        Polynomial degree.
    axis : int, optional
        Axis of ``frequencies`` corresponding to ``volume``.

    Returns
    -------
    tuple
        Series fit result and interpolation result.
    """
    fit, interpolation = interpolate_polynomial_series(
        volume, frequencies, points, degree, axis=axis
    )
    if (
        interpolation.success
        and interpolation.values is not None
        and np.any(interpolation.values < 0.0)
    ):
        interpolation.warnings.append(
            "interpolated frequency grid contains negative values"
        )
        interpolation.metadata["contains_negative_values"] = True
    return fit, interpolation


def interpolate_property_grid(
    volume: ArrayLike,
    properties: ArrayLike,
    points: ArrayLike | float,
    degree: int,
    *,
    axis: int = -1,
) -> tuple[PolynomialSeriesFit, InterpolationResult]:
    """Interpolate thermodynamic properties along the volume coordinate.

    Parameters
    ----------
    volume : array-like
        Volume values.
    properties : array-like
        Property values sampled on the same volume grid.
    points : scalar or array-like
        Volumes where the properties are evaluated.
    degree : int
        Polynomial degree.
    axis : int, optional
        Axis of ``properties`` corresponding to ``volume``.

    Returns
    -------
    tuple
        Series fit result and interpolation result.
    """
    return interpolate_polynomial_series(volume, properties, points, degree, axis=axis)
