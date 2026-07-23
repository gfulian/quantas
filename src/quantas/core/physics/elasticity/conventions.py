# -*- coding: utf-8 -*-

"""Voigt notation and Cartesian elastic-tensor conversions.

Quantas uses the conventional ordering ``11, 22, 33, 23, 13, 12``. Stiffness
matrices are paired with engineering shear strains, while compliance matrices
therefore require factors of two for each shear index when converted to or
from a full Cartesian fourth-rank tensor.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from quantas.core.geometry.rotations import rotate_cartesian_rank4


VOIGT_LABELS: tuple[str, ...] = ("11", "22", "33", "23", "13", "12")
VOIGT_PAIRS: tuple[tuple[int, int], ...] = (
    (0, 0),
    (1, 1),
    (2, 2),
    (1, 2),
    (0, 2),
    (0, 1),
)
ENGINEERING_SHEAR_FACTORS: tuple[float, ...] = (1.0, 1.0, 1.0, 2.0, 2.0, 2.0)

VOIGT_INDEX_MAP = np.asarray(
    [
        [0, 5, 4],
        [5, 1, 3],
        [4, 3, 2],
    ],
    dtype=int,
)
VOIGT_INDEX_MAP.setflags(write=False)


def _validate_voigt_matrix(matrix: ArrayLike, name: str) -> NDArray[np.float64]:
    """Return a finite ``6 x 6`` floating-point matrix."""
    array = np.asarray(matrix, dtype=float)
    if array.shape != (6, 6):
        raise ValueError(f"The {name} matrix must have shape (6, 6).")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"The {name} matrix must contain finite values.")
    return array


def _validate_cartesian_tensor(tensor: ArrayLike, name: str) -> NDArray[np.float64]:
    """Return a finite Cartesian fourth-rank tensor."""
    array = np.asarray(tensor, dtype=float)
    if array.shape != (3, 3, 3, 3):
        raise ValueError(f"The Cartesian {name} tensor must have shape (3, 3, 3, 3).")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"The Cartesian {name} tensor must contain finite values.")
    return array


def voigt_stiffness_to_cartesian(
    stiffness: ArrayLike,
) -> NDArray[np.float64]:
    """Convert a Voigt stiffness matrix to ``C_ijkl``.

    Parameters
    ----------
    stiffness : array_like
        Stiffness matrix in the Quantas Voigt ordering.

    Returns
    -------
    numpy.ndarray
        Cartesian fourth-rank stiffness tensor with shape ``(3, 3, 3, 3)``.

    Raises
    ------
    ValueError
        If the matrix has an invalid shape or contains non-finite values.
    """
    matrix = _validate_voigt_matrix(stiffness, "stiffness")
    tensor = matrix[
        VOIGT_INDEX_MAP[:, :, None, None], VOIGT_INDEX_MAP[None, None, :, :]
    ]
    return np.asarray(tensor, dtype=float)


def cartesian_stiffness_to_voigt(
    tensor: ArrayLike,
) -> NDArray[np.float64]:
    """Convert ``C_ijkl`` to a Voigt stiffness matrix.

    Parameters
    ----------
    tensor : array_like
        Cartesian stiffness tensor with shape ``(3, 3, 3, 3)``.

    Returns
    -------
    numpy.ndarray
        Stiffness matrix in the ordering ``11, 22, 33, 23, 13, 12``.

    Raises
    ------
    ValueError
        If the tensor has an invalid shape or contains non-finite values.
    """
    cartesian = _validate_cartesian_tensor(tensor, "stiffness")
    matrix = np.empty((6, 6), dtype=float)
    for p, (i, j) in enumerate(VOIGT_PAIRS):
        for q, (k, ell) in enumerate(VOIGT_PAIRS):
            matrix[p, q] = cartesian[i, j, k, ell]
    return matrix


def voigt_compliance_to_cartesian(
    compliance: ArrayLike,
) -> NDArray[np.float64]:
    """Convert an engineering Voigt compliance matrix to ``S_ijkl``.

    Parameters
    ----------
    compliance : array_like
        Compliance matrix in the Quantas Voigt ordering.

    Returns
    -------
    numpy.ndarray
        Cartesian fourth-rank compliance tensor.

    Raises
    ------
    ValueError
        If the matrix has an invalid shape or contains non-finite values.
    """
    matrix = _validate_voigt_matrix(compliance, "compliance")
    factors = np.asarray(ENGINEERING_SHEAR_FACTORS, dtype=float)
    scaled = matrix / (factors[:, None] * factors[None, :])
    tensor = scaled[
        VOIGT_INDEX_MAP[:, :, None, None],
        VOIGT_INDEX_MAP[None, None, :, :],
    ]
    return np.asarray(tensor, dtype=float)


def cartesian_compliance_to_voigt(
    tensor: ArrayLike,
) -> NDArray[np.float64]:
    """Convert ``S_ijkl`` to an engineering Voigt compliance matrix.

    Parameters
    ----------
    tensor : array_like
        Cartesian compliance tensor with shape ``(3, 3, 3, 3)``.

    Returns
    -------
    numpy.ndarray
        Engineering compliance matrix in the Quantas Voigt ordering.

    Raises
    ------
    ValueError
        If the tensor has an invalid shape or contains non-finite values.
    """
    cartesian = _validate_cartesian_tensor(tensor, "compliance")
    matrix = np.empty((6, 6), dtype=float)
    factors = np.asarray(ENGINEERING_SHEAR_FACTORS, dtype=float)
    for p, (i, j) in enumerate(VOIGT_PAIRS):
        for q, (k, ell) in enumerate(VOIGT_PAIRS):
            matrix[p, q] = factors[p] * factors[q] * cartesian[i, j, k, ell]
    return matrix


def rotate_voigt_stiffness(
    stiffness: ArrayLike,
    rotation: ArrayLike,
) -> NDArray[np.float64]:
    """Transform stiffness components to a rotated Cartesian frame.

    The supplied matrix ``R`` follows the component transformation

    ``C'_{ijkl} = R_ia R_jb R_kc R_ld C_abcd``.

    When the rows of ``R`` are the target Cartesian basis vectors expressed in
    the source frame, the result contains the components of the same physical
    tensor in the target frame. The same equation can also be interpreted as
    an active right-handed rotation when ``R`` is an active rotation matrix.

    Parameters
    ----------
    stiffness : array_like
        Stiffness matrix in Voigt notation.
    rotation : array_like
        Proper orthogonal ``3 x 3`` transformation matrix.

    Returns
    -------
    numpy.ndarray
        Transformed stiffness matrix in the same Voigt ordering.

    Raises
    ------
    ValueError
        If the stiffness matrix or rotation is invalid.
    """
    cartesian = voigt_stiffness_to_cartesian(stiffness)
    rotated = rotate_cartesian_rank4(cartesian, rotation)
    return cartesian_stiffness_to_voigt(rotated)
