# -*- coding: utf-8 -*-

"""Cell-matrix construction, conversion, and periodic-coordinate utilities."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


def cell_volume(lattice: ArrayLike) -> float:
    """Return the positive volume of a direct lattice.

    Parameters
    ----------
    lattice : array_like
        Direct lattice vectors stored by rows with shape ``(3, 3)``.

    Returns
    -------
    float
        Absolute determinant of the lattice matrix.

    Raises
    ------
    ValueError
        If the lattice does not have shape ``(3, 3)`` or is singular.
    """
    matrix = np.asarray(lattice, dtype=np.float64)
    if matrix.shape != (3, 3):
        raise ValueError("lattice must have shape (3, 3)")
    volume = abs(float(np.linalg.det(matrix)))
    if volume <= np.finfo(np.float64).eps:
        raise ValueError("lattice must be non-singular")
    return volume


def lattice_from_parameters(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> FloatArray:
    """Construct row-wise direct lattice vectors from cell parameters.

    The first vector is placed along Cartesian ``x`` and the second in the
    ``xy`` plane. Angles are interpreted in degrees.

    Parameters
    ----------
    a, b, c : float
        Lattice-vector lengths in angstrom.
    alpha, beta, gamma : float
        Inter-vector angles in degrees. ``alpha`` is the angle between
        ``b`` and ``c``, ``beta`` between ``a`` and ``c``, and ``gamma``
        between ``a`` and ``b``.

    Returns
    -------
    ndarray
        Direct lattice matrix with shape ``(3, 3)`` and vectors stored by rows.

    Raises
    ------
    ValueError
        If lengths are non-positive or the parameters define a singular cell.
    """
    lengths = np.asarray([a, b, c], dtype=np.float64)
    if np.any(lengths <= 0.0):
        raise ValueError("lattice lengths must be positive")
    ar, br, gr = np.deg2rad([alpha, beta, gamma])
    sin_gamma = float(np.sin(gr))
    if abs(sin_gamma) <= np.finfo(np.float64).eps:
        raise ValueError("gamma produces a singular lattice")

    lattice = np.zeros((3, 3), dtype=np.float64)
    lattice[0] = [a, 0.0, 0.0]
    lattice[1] = [b * np.cos(gr), b * sin_gamma, 0.0]
    lattice[2, 0] = c * np.cos(br)
    lattice[2, 1] = c * (np.cos(ar) - np.cos(br) * np.cos(gr)) / sin_gamma
    z2 = c * c - lattice[2, 0] ** 2 - lattice[2, 1] ** 2
    if z2 < -1.0e-10:
        raise ValueError("cell parameters do not define a real lattice")
    lattice[2, 2] = np.sqrt(max(z2, 0.0))
    cell_volume(lattice)
    return lattice


def lattice_parameters(lattice: ArrayLike) -> FloatArray:
    """Return lengths and angles for a row-wise direct lattice.

    Parameters
    ----------
    lattice : array_like
        Direct lattice matrix with shape ``(3, 3)``.

    Returns
    -------
    ndarray
        ``[a, b, c, alpha, beta, gamma]`` with lengths in angstrom and angles
        in degrees.
    """
    matrix = np.asarray(lattice, dtype=np.float64)
    cell_volume(matrix)
    lengths = np.linalg.norm(matrix, axis=1)

    def _angle(first: FloatArray, second: FloatArray) -> float:
        cosine = float(
            np.dot(first, second) / (np.linalg.norm(first) * np.linalg.norm(second))
        )
        return float(np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0))))

    alpha = _angle(matrix[1], matrix[2])
    beta = _angle(matrix[0], matrix[2])
    gamma = _angle(matrix[0], matrix[1])
    return np.asarray([*lengths, alpha, beta, gamma], dtype=np.float64)


def wrap_fractional(positions: ArrayLike) -> FloatArray:
    """Wrap fractional coordinates into the half-open interval ``[0, 1)``.

    Parameters
    ----------
    positions : array_like
        Fractional coordinates with final dimension three.

    Returns
    -------
    ndarray
        Wrapped fractional coordinates.
    """
    values = np.asarray(positions, dtype=np.float64)
    return values - np.floor(values)


def minimum_image_delta(delta: ArrayLike) -> FloatArray:
    """Map fractional displacements to the nearest periodic image.

    Parameters
    ----------
    delta : array_like
        Fractional displacement vectors with final dimension three.

    Returns
    -------
    ndarray
        Displacements in the interval ``[-0.5, 0.5)``.
    """
    values = np.asarray(delta, dtype=np.float64)
    return values - np.floor(values + 0.5)


def fractional_to_cartesian(positions: ArrayLike, lattice: ArrayLike) -> FloatArray:
    """Convert fractional positions to Cartesian coordinates.

    Parameters
    ----------
    positions : array_like
        Fractional coordinates with final dimension three.
    lattice : array_like
        Row-wise direct lattice matrix.

    Returns
    -------
    ndarray
        Cartesian coordinates in angstrom.
    """
    return np.asarray(positions, dtype=np.float64) @ np.asarray(
        lattice,
        dtype=np.float64,
    )


def cartesian_to_fractional(positions: ArrayLike, lattice: ArrayLike) -> FloatArray:
    """Convert Cartesian positions to fractional coordinates.

    Parameters
    ----------
    positions : array_like
        Cartesian coordinates with final dimension three.
    lattice : array_like
        Row-wise direct lattice matrix.

    Returns
    -------
    ndarray
        Fractional coordinates.

    Raises
    ------
    ValueError
        If the lattice is singular.
    """
    matrix = np.asarray(lattice, dtype=np.float64)
    cell_volume(matrix)
    return np.asarray(positions, dtype=np.float64) @ np.linalg.inv(matrix)
