# -*- coding: utf-8 -*-

"""Spherical coordinates and regular angular sampling grids."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np
from numpy.typing import ArrayLike, NDArray


class Hemisphere(str, Enum):
    """Supported polar domains for regular spherical grids."""

    UPPER = "upper"
    LOWER = "lower"
    FULL = "full"


@dataclass(frozen=True, slots=True)
class SphericalGrid:
    """Regular angular grid on a sphere or hemisphere.

    Parameters
    ----------
    theta, phi : ndarray
        One-dimensional polar and azimuthal coordinates in radians. The
        azimuthal coordinate excludes the endpoint at ``2π``.
    theta_grid, phi_grid : ndarray
        Two-dimensional angular coordinates with shape ``(ntheta, nphi)``.
    directions : ndarray
        Unit Cartesian directions with shape ``(ntheta, nphi, 3)``.
    hemisphere : Hemisphere
        Polar domain represented by the grid.
    """

    theta: NDArray[np.float64]
    phi: NDArray[np.float64]
    theta_grid: NDArray[np.float64]
    phi_grid: NDArray[np.float64]
    directions: NDArray[np.float64]
    hemisphere: Hemisphere

    @property
    def shape(self) -> tuple[int, int]:
        """Return the angular grid shape ``(ntheta, nphi)``."""
        return int(self.theta_grid.shape[0]), int(self.theta_grid.shape[1])

    @property
    def size(self) -> int:
        """Return the number of sampled grid positions."""
        return int(self.theta_grid.size)

    @property
    def flat_directions(self) -> NDArray[np.float64]:
        """Return a read-only ``(n_points, 3)`` view of the directions."""
        view = self.directions.reshape(-1, 3)
        view.setflags(write=False)
        return view

    def reshape_scalar_field(self, values: ArrayLike) -> NDArray[np.float64]:
        """Reshape a scalar field from flat traversal order to the grid.

        Parameters
        ----------
        values : array_like
            Scalar values with shape ``(n_points,)``.

        Returns
        -------
        ndarray
            Read-only array with shape ``(ntheta, nphi)``.

        Raises
        ------
        ValueError
            If the number of values does not match the grid size.
        """
        array = np.asarray(values, dtype=float)
        if array.shape != (self.size,):
            raise ValueError(
                f"Scalar field must have shape ({self.size},), got {array.shape}."
            )
        result = np.array(array.reshape(self.shape), dtype=float, copy=True)
        result.setflags(write=False)
        return result


def create_spherical_grid(
    ntheta: int,
    nphi: int,
    *,
    hemisphere: Hemisphere | str = Hemisphere.FULL,
) -> SphericalGrid:
    """Create a regular spherical grid without duplicating the azimuthal seam.

    Parameters
    ----------
    ntheta : int
        Number of polar samples, including both boundaries of the selected
        polar domain.
    nphi : int
        Number of azimuthal samples from ``0`` to ``2π``. The final endpoint
        is excluded because it duplicates the first meridian.
    hemisphere : Hemisphere or str, optional
        Polar domain: ``"upper"``, ``"lower"`` or ``"full"``.

    Returns
    -------
    SphericalGrid
        Read-only angular coordinates and Cartesian directions.

    Raises
    ------
    ValueError
        If the requested grid dimensions or hemisphere are invalid.
    """
    if isinstance(ntheta, bool) or int(ntheta) != ntheta or int(ntheta) < 2:
        raise ValueError("ntheta must be an integer greater than or equal to 2.")
    if isinstance(nphi, bool) or int(nphi) != nphi or int(nphi) < 3:
        raise ValueError("nphi must be an integer greater than or equal to 3.")

    domain = Hemisphere(hemisphere)
    theta_min, theta_max = {
        Hemisphere.UPPER: (0.0, 0.5 * np.pi),
        Hemisphere.LOWER: (0.5 * np.pi, np.pi),
        Hemisphere.FULL: (0.0, np.pi),
    }[domain]
    theta = np.linspace(
        theta_min,
        theta_max,
        int(ntheta),
        endpoint=True,
        dtype=float,
    )
    phi = np.linspace(
        0.0,
        2.0 * np.pi,
        int(nphi),
        endpoint=False,
        dtype=float,
    )
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")
    directions = spherical_direction(theta_grid, phi_grid)

    arrays = theta, phi, theta_grid, phi_grid, directions
    for array in arrays:
        array.setflags(write=False)

    return SphericalGrid(
        theta=theta,
        phi=phi,
        theta_grid=theta_grid,
        phi_grid=phi_grid,
        directions=directions,
        hemisphere=domain,
    )


def spherical_direction(
    theta: ArrayLike,
    phi: ArrayLike,
) -> NDArray[np.float64]:
    """Convert polar and azimuthal angles to Cartesian unit directions.

    Parameters
    ----------
    theta, phi : array_like
        Broadcast-compatible polar and azimuthal angles in radians.

    Returns
    -------
    ndarray
        Cartesian directions with final axis length three.

    Raises
    ------
    ValueError
        If the angular arrays cannot be broadcast or contain non-finite values.
    """
    polar = np.asarray(theta, dtype=float)
    azimuth = np.asarray(phi, dtype=float)
    if not np.all(np.isfinite(polar)) or not np.all(np.isfinite(azimuth)):
        raise ValueError("Spherical angles must contain finite values.")
    try:
        polar, azimuth = np.broadcast_arrays(polar, azimuth)
    except ValueError as exc:
        raise ValueError("theta and phi must be broadcast-compatible.") from exc

    sine = np.sin(polar)
    cosine = np.cos(polar)
    threshold = 8.0 * np.finfo(float).eps
    sine = np.where(np.abs(sine) <= threshold, 0.0, sine)
    cosine = np.where(np.abs(cosine) <= threshold, 0.0, cosine)
    return np.stack(
        (
            sine * np.cos(azimuth),
            sine * np.sin(azimuth),
            cosine,
        ),
        axis=-1,
    ).astype(float, copy=False)


def cartesian_to_spherical(
    directions: ArrayLike,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Convert Cartesian directions to polar and azimuthal angles.

    Parameters
    ----------
    directions : array_like
        Cartesian vectors with final shape component three.

    Returns
    -------
    theta, phi : tuple of ndarray
        Polar angles in ``[0, π]`` and azimuths in ``[0, 2π)``.

    Raises
    ------
    ValueError
        If the vectors have an invalid shape, contain non-finite values or
        include a zero vector.
    """
    vectors = np.asarray(directions, dtype=float)
    if vectors.ndim < 1 or vectors.shape[-1] != 3:
        raise ValueError("Cartesian directions must have final dimension three.")
    if not np.all(np.isfinite(vectors)):
        raise ValueError("Cartesian directions must contain finite values.")
    norms = np.linalg.norm(vectors, axis=-1)
    if np.any(norms == 0.0):
        raise ValueError("Cartesian directions must be non-zero.")
    unit = vectors / norms[..., np.newaxis]
    theta = np.arccos(np.clip(unit[..., 2], -1.0, 1.0))
    phi = np.mod(np.arctan2(unit[..., 1], unit[..., 0]), 2.0 * np.pi)
    return np.asarray(theta, dtype=float), np.asarray(phi, dtype=float)


def close_periodic_seam(
    values: ArrayLike,
    *,
    axis: int = -1,
) -> NDArray[np.float64]:
    """Append the first azimuthal slice to close a periodic plotting seam.

    Parameters
    ----------
    values : array_like
        Array whose selected axis follows a non-duplicated periodic azimuth.
    axis : int, optional
        Periodic axis.

    Returns
    -------
    ndarray
        Copy with one additional slice equal to the first slice.

    Raises
    ------
    ValueError
        If the input is scalar or the selected axis is empty.
    """
    array = np.asarray(values)
    if array.ndim == 0:
        raise ValueError("A periodic seam cannot be added to a scalar value.")
    normalized_axis = int(axis)
    if normalized_axis < 0:
        normalized_axis += array.ndim
    if normalized_axis < 0 or normalized_axis >= array.ndim:
        raise ValueError("axis is outside the array dimensions.")
    if array.shape[normalized_axis] == 0:
        raise ValueError("The periodic axis must contain at least one value.")
    first = np.take(array, [0], axis=normalized_axis)
    result = np.concatenate((array, first), axis=normalized_axis)
    result.setflags(write=False)
    return result
