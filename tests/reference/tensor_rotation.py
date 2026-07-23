"""Independent fourth-rank tensor transformations for regression tests.

The helpers in this module define the passive/active transformation convention
that future production rotation APIs must state explicitly.  They remain in
the test tree until that public API is implemented.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

_VOIGT_PAIRS = (
    (0, 0),
    (1, 1),
    (2, 2),
    (1, 2),
    (0, 2),
    (0, 1),
)


def voigt_stiffness_to_cartesian(
    stiffness: ArrayLike,
) -> NDArray[np.float64]:
    """Convert a Voigt stiffness matrix to ``C_ijkl`` independently.

    Parameters
    ----------
    stiffness : array_like
        Stiffness matrix with shape ``(6, 6)``.

    Returns
    -------
    numpy.ndarray
        Cartesian fourth-rank tensor with shape ``(3, 3, 3, 3)``.

    Raises
    ------
    ValueError
        If the input does not have shape ``(6, 6)``.
    """
    matrix = np.asarray(stiffness, dtype=float)
    if matrix.shape != (6, 6):
        raise ValueError("The stiffness matrix must have shape (6, 6).")
    tensor = np.empty((3, 3, 3, 3), dtype=float)
    for p, (i, j) in enumerate(_VOIGT_PAIRS):
        for q, (k, ell) in enumerate(_VOIGT_PAIRS):
            for first in ((i, j), (j, i)):
                for second in ((k, ell), (ell, k)):
                    tensor[first[0], first[1], second[0], second[1]] = matrix[p, q]
    return tensor


def cartesian_stiffness_to_voigt(
    tensor: ArrayLike,
) -> NDArray[np.float64]:
    """Convert ``C_ijkl`` to the Quantas Voigt stiffness convention.

    Parameters
    ----------
    tensor : array_like
        Cartesian stiffness tensor with shape ``(3, 3, 3, 3)``.

    Returns
    -------
    numpy.ndarray
        Stiffness matrix with shape ``(6, 6)``.

    Raises
    ------
    ValueError
        If the input does not have shape ``(3, 3, 3, 3)``.
    """
    cartesian = np.asarray(tensor, dtype=float)
    if cartesian.shape != (3, 3, 3, 3):
        raise ValueError("The Cartesian tensor must have shape (3, 3, 3, 3).")
    matrix = np.empty((6, 6), dtype=float)
    for p, (i, j) in enumerate(_VOIGT_PAIRS):
        for q, (k, ell) in enumerate(_VOIGT_PAIRS):
            matrix[p, q] = cartesian[i, j, k, ell]
    return matrix


def rotation_matrix(axis: ArrayLike, angle: float) -> NDArray[np.float64]:
    """Return a right-handed active rotation matrix.

    Parameters
    ----------
    axis : array_like
        Cartesian rotation axis with shape ``(3,)``.
    angle : float
        Right-handed rotation angle in radians.

    Returns
    -------
    numpy.ndarray
        Proper orthogonal matrix with shape ``(3, 3)``.

    Raises
    ------
    ValueError
        If the axis is invalid or has zero length.
    """
    vector = np.asarray(axis, dtype=float)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError("The rotation axis must be a finite vector of shape (3,).")
    norm = np.linalg.norm(vector)
    if norm == 0.0:
        raise ValueError("The rotation axis must have non-zero length.")
    x, y, z = vector / norm
    cosine = float(np.cos(angle))
    sine = float(np.sin(angle))
    complement = 1.0 - cosine
    return np.asarray(
        [
            [
                cosine + x * x * complement,
                x * y * complement - z * sine,
                x * z * complement + y * sine,
            ],
            [
                y * x * complement + z * sine,
                cosine + y * y * complement,
                y * z * complement - x * sine,
            ],
            [
                z * x * complement - y * sine,
                z * y * complement + x * sine,
                cosine + z * z * complement,
            ],
        ],
        dtype=float,
    )


def rotate_cartesian_stiffness(
    tensor: ArrayLike,
    rotation: ArrayLike,
) -> NDArray[np.float64]:
    """Actively rotate a fourth-rank Cartesian stiffness tensor.

    The convention is

    ``C'_{ijkl} = R_{ia} R_{jb} R_{kc} R_{ld} C_{abcd}``.

    Parameters
    ----------
    tensor : array_like
        Cartesian stiffness tensor with shape ``(3, 3, 3, 3)``.
    rotation : array_like
        Proper orthogonal rotation matrix with shape ``(3, 3)``.

    Returns
    -------
    numpy.ndarray
        Rotated Cartesian stiffness tensor.

    Raises
    ------
    ValueError
        If the tensor or rotation has an invalid shape, or if the matrix is not
        a proper orthogonal rotation.
    """
    cartesian = np.asarray(tensor, dtype=float)
    matrix = np.asarray(rotation, dtype=float)
    if cartesian.shape != (3, 3, 3, 3):
        raise ValueError("The Cartesian tensor must have shape (3, 3, 3, 3).")
    if matrix.shape != (3, 3):
        raise ValueError("The rotation matrix must have shape (3, 3).")
    if not np.allclose(matrix @ matrix.T, np.eye(3), atol=1.0e-12):
        raise ValueError("The rotation matrix must be orthogonal.")
    if not np.isclose(np.linalg.det(matrix), 1.0, atol=1.0e-12):
        raise ValueError("The rotation matrix must be right-handed.")
    return np.einsum(
        "ia,jb,kc,ld,abcd->ijkl",
        matrix,
        matrix,
        matrix,
        matrix,
        cartesian,
        optimize=True,
    )


def rotate_voigt_stiffness(
    stiffness: ArrayLike,
    rotation: ArrayLike,
) -> NDArray[np.float64]:
    """Rotate a Voigt stiffness matrix through its Cartesian tensor.

    Parameters
    ----------
    stiffness : array_like
        Stiffness matrix with shape ``(6, 6)``.
    rotation : array_like
        Proper orthogonal rotation matrix with shape ``(3, 3)``.

    Returns
    -------
    numpy.ndarray
        Rotated stiffness matrix in the same Voigt ordering.
    """
    cartesian = voigt_stiffness_to_cartesian(stiffness)
    rotated = rotate_cartesian_stiffness(cartesian, rotation)
    return cartesian_stiffness_to_voigt(rotated)
