# -*- coding: utf-8 -*-

"""Validation utilities for elastic stiffness and compliance tensors."""

from __future__ import annotations

from enum import Enum

import numpy as np
from numpy.typing import ArrayLike, NDArray


DEFAULT_STIFFNESS_SYMMETRY_TOLERANCE = 1.0e-3


class StiffnessSymmetryCriterion(str, Enum):
    """Numerical criteria used to validate stiffness-matrix symmetry."""

    FROBENIUS = "frobenius"
    ELEMENTWISE = "elementwise"


def stiffness_symmetry_residual(matrix: ArrayLike) -> float:
    """Return the Frobenius norm of the stiffness antisymmetric part.

    Parameters
    ----------
    matrix : array_like
        Candidate stiffness matrix with shape ``(6, 6)``.

    Returns
    -------
    float
        ``||C - C.T||`` evaluated with the Euclidean/Frobenius norm.

    Raises
    ------
    ValueError
        If the input does not have shape ``(6, 6)`` or contains non-finite
        values.
    """
    array = np.asarray(matrix, dtype=float)
    if array.shape != (6, 6):
        raise ValueError("The elastic stiffness matrix must have shape (6, 6).")
    if not np.all(np.isfinite(array)):
        raise ValueError("The elastic stiffness matrix must contain finite values.")
    return float(np.linalg.norm(array - array.T))


def stiffness_maximum_asymmetry(matrix: ArrayLike) -> float:
    """Return the largest absolute antisymmetric stiffness component.

    Parameters
    ----------
    matrix : array_like
        Candidate stiffness matrix with shape ``(6, 6)``.

    Returns
    -------
    float
        Maximum value of ``abs(C_ij - C_ji)``.

    Raises
    ------
    ValueError
        If the input does not have shape ``(6, 6)`` or contains non-finite
        values.
    """
    array = np.asarray(matrix, dtype=float)
    if array.shape != (6, 6):
        raise ValueError("The elastic stiffness matrix must have shape (6, 6).")
    if not np.all(np.isfinite(array)):
        raise ValueError("The elastic stiffness matrix must contain finite values.")
    return float(np.max(np.abs(array - array.T)))


def validate_stiffness_matrix(
    stiffness: ArrayLike,
    *,
    symmetry_tolerance: float = DEFAULT_STIFFNESS_SYMMETRY_TOLERANCE,
    symmetry_criterion: StiffnessSymmetryCriterion | str = (
        StiffnessSymmetryCriterion.FROBENIUS
    ),
    copy: bool = True,
) -> NDArray[np.float64]:
    """Validate and normalize an elastic stiffness matrix.

    The default symmetry criterion intentionally preserves the historical
    Quantas behavior, based on the Frobenius norm ``||C - C.T||`` rather than
    an element-wise comparison.

    Parameters
    ----------
    stiffness : array_like
        Elastic stiffness matrix in Voigt notation.
    symmetry_tolerance : float, optional
        Maximum accepted symmetry residual.
    symmetry_criterion : StiffnessSymmetryCriterion or str, optional
        ``"frobenius"`` preserves the historical Elasticity criterion based
        on ``||C - C.T||``. ``"elementwise"`` preserves the historical
        SEISMIC criterion based on the largest individual mismatch.
    copy : bool, optional
        Whether to return an independent copy of the input matrix.

    Returns
    -------
    numpy.ndarray
        Validated floating-point matrix with shape ``(6, 6)``.

    Raises
    ------
    ValueError
        If the matrix has an invalid shape, contains non-finite values, is not
        symmetric within the requested tolerance, or if the tolerance is
        negative or non-finite.
    """
    if not np.isfinite(symmetry_tolerance) or symmetry_tolerance < 0.0:
        raise ValueError("symmetry_tolerance must be finite and non-negative.")
    criterion = StiffnessSymmetryCriterion(symmetry_criterion)

    if copy:
        matrix = np.array(stiffness, dtype=float, copy=True)
    else:
        matrix = np.asarray(stiffness, dtype=float)
    if criterion is StiffnessSymmetryCriterion.FROBENIUS:
        residual = stiffness_symmetry_residual(matrix)
    else:
        residual = stiffness_maximum_asymmetry(matrix)
    if residual > symmetry_tolerance:
        raise ValueError("The elastic stiffness matrix must be symmetric.")
    return matrix


def invert_stiffness_matrix(
    stiffness: ArrayLike,
) -> NDArray[np.float64]:
    """Return the compliance matrix associated with a stiffness matrix.

    Parameters
    ----------
    stiffness : array_like
        Valid stiffness matrix with shape ``(6, 6)``.

    Returns
    -------
    numpy.ndarray
        Inverse matrix in engineering Voigt notation.

    Raises
    ------
    ValueError
        If the matrix is singular.
    """
    matrix = np.asarray(stiffness, dtype=float)
    try:
        return np.linalg.inv(matrix)
    except np.linalg.LinAlgError as exc:
        raise ValueError("matrix is singular") from exc
