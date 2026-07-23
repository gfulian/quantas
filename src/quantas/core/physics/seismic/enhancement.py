# -*- coding: utf-8 -*-

"""Analytical curvature and enhancement of acoustic ray fields."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .group import DirectionalGroupResult, GroupModeResult
from .modes import MODE_INDEX, WaveMode


@dataclass(frozen=True, slots=True)
class EnhancementModeResult:
    """Enhancement solution for one acoustic mode.

    Parameters
    ----------
    group : GroupModeResult
        Phase and group solution associated with the selected mode.
    eigenvalue_hessian : ndarray
        Cartesian Hessian of the Christoffel eigenvalue with shape ``(3, 3)``
        in km^2 s^-2.
    ray_direction_gradient : ndarray
        Cartesian gradient of the unit ray direction with shape ``(3, 3)``.
    area_factor : float
        Geometrical area factor obtained from the cofactor of the ray-direction
        gradient.
    caustic_threshold : float
        Numerical threshold used to identify possible caustics.
    enhancement : float
        Positive enhancement factor :math:`A`. Undefined modes contain
        ``nan`` and exact singularities contain ``inf``.
    log10_enhancement : float
        Base-10 logarithm of :math:`A`.
    valid : bool
        Whether the analytical curvature quantities are defined.
    resolved : bool
        Whether the acoustic mode is uniquely resolved.
    finite : bool
        Whether both enhancement representations are finite.
    caustic_candidate : bool
        Whether the area factor is below the selected numerical threshold.
    """

    group: GroupModeResult
    eigenvalue_hessian: NDArray[np.float64]
    ray_direction_gradient: NDArray[np.float64]
    area_factor: float
    caustic_threshold: float
    enhancement: float
    log10_enhancement: float
    valid: bool
    resolved: bool
    finite: bool
    caustic_candidate: bool


@dataclass(frozen=True, slots=True)
class DirectionalEnhancementResult:
    """Curvature and enhancement results for one wave-normal direction.

    Parameters
    ----------
    group : DirectionalGroupResult
        Pointwise phase and group solution.
    eigenvalue_hessians : ndarray
        Eigenvalue Hessians with shape ``(3, 3, 3)`` in mode order.
    ray_direction_gradients : ndarray
        Gradients of the unit ray directions with shape ``(3, 3, 3)``.
    area_factors : ndarray
        Geometrical area factors with shape ``(3,)``.
    caustic_thresholds : ndarray
        Mode-specific numerical thresholds for possible caustics.
    enhancement : ndarray
        Raw enhancement factors :math:`A` with shape ``(3,)``.
    log10_enhancement : ndarray
        Base-10 logarithms of the enhancement factors.
    valid_mask : ndarray
        Boolean mask identifying defined analytical curvature quantities.
    resolved_mask : ndarray
        Boolean mask identifying uniquely resolved acoustic modes.
    finite_mask : ndarray
        Boolean mask identifying finite enhancement values.
    caustic_candidate_mask : ndarray
        Boolean mask identifying area factors close to zero.
    """

    group: DirectionalGroupResult
    eigenvalue_hessians: NDArray[np.float64]
    ray_direction_gradients: NDArray[np.float64]
    area_factors: NDArray[np.float64]
    caustic_thresholds: NDArray[np.float64]
    enhancement: NDArray[np.float64]
    log10_enhancement: NDArray[np.float64]
    valid_mask: NDArray[np.bool_]
    resolved_mask: NDArray[np.bool_]
    finite_mask: NDArray[np.bool_]
    caustic_candidate_mask: NDArray[np.bool_]

    @property
    def direction(self) -> NDArray[np.float64]:
        """Return the normalized wave-normal direction."""
        return self.group.direction

    def for_mode(self, mode: WaveMode | str) -> EnhancementModeResult:
        """Return the enhancement solution for one acoustic mode.

        Parameters
        ----------
        mode : WaveMode or str
            Mode enum or one of ``"v_s2"``, ``"v_s1"`` and ``"v_p"``.

        Returns
        -------
        EnhancementModeResult
            Selected mode result.

        Raises
        ------
        ValueError
            If ``mode`` is not a supported acoustic mode.
        """
        resolved = WaveMode(mode)
        index = MODE_INDEX[resolved]
        return EnhancementModeResult(
            group=self.group.for_mode(resolved),
            eigenvalue_hessian=self.eigenvalue_hessians[index],
            ray_direction_gradient=self.ray_direction_gradients[index],
            area_factor=float(self.area_factors[index]),
            caustic_threshold=float(self.caustic_thresholds[index]),
            enhancement=float(self.enhancement[index]),
            log10_enhancement=float(self.log10_enhancement[index]),
            valid=bool(self.valid_mask[index]),
            resolved=bool(self.resolved_mask[index]),
            finite=bool(self.finite_mask[index]),
            caustic_candidate=bool(self.caustic_candidate_mask[index]),
        )

    @property
    def is_fully_valid(self) -> bool:
        """Return whether all three curvature solutions are defined."""
        return bool(np.all(self.valid_mask))

    @property
    def is_fully_finite(self) -> bool:
        """Return whether all three enhancement factors are finite."""
        return bool(np.all(self.finite_mask))

    @property
    def has_caustic_candidate(self) -> bool:
        """Return whether at least one mode approaches a caustic."""
        return bool(np.any(self.caustic_candidate_mask))


def calculate_eigenvalue_hessian(
    christoffel: NDArray[np.float64],
    christoffel_gradient: NDArray[np.float64],
    christoffel_hessian: NDArray[np.float64],
    eigenvalue: float,
    polarization: NDArray[np.float64],
    *,
    pseudoinverse_rcond: float,
) -> NDArray[np.float64]:
    """Calculate the analytical Hessian of one Christoffel eigenvalue.

    Parameters
    ----------
    christoffel : ndarray
        Symmetric Christoffel matrix with shape ``(3, 3)``.
    christoffel_gradient : ndarray
        Matrix gradient with shape ``(3, 3, 3)``.
    christoffel_hessian : ndarray
        Matrix Hessian with shape ``(3, 3, 3, 3)``.
    eigenvalue : float
        Selected non-degenerate Christoffel eigenvalue.
    polarization : ndarray
        Normalized eigenvector with shape ``(3,)``.
    pseudoinverse_rcond : float
        Relative cutoff used by :func:`numpy.linalg.pinv`.

    Returns
    -------
    ndarray
        Symmetric eigenvalue Hessian with shape ``(3, 3)`` in km^2 s^-2.
    """
    direct = np.einsum(
        "i,abij,j->ab",
        polarization,
        christoffel_hessian,
        polarization,
        optimize=True,
    )
    shifted = eigenvalue * np.eye(3, dtype=float) - christoffel
    inverse = np.linalg.pinv(shifted, rcond=pseudoinverse_rcond)
    derivative_vectors = np.einsum(
        "aij,j->ai",
        christoffel_gradient,
        polarization,
        optimize=True,
    )
    correction = 2.0 * derivative_vectors @ inverse @ derivative_vectors.T
    result = np.asarray(direct + correction, dtype=float)
    return np.asarray(0.5 * (result + result.T), dtype=float)


def calculate_ray_direction_gradient(
    eigenvalue_hessian: NDArray[np.float64],
    ray_direction: NDArray[np.float64],
    phase_speed: float,
    group_speed: float,
) -> NDArray[np.float64]:
    """Calculate the Cartesian gradient of one unit ray direction.

    Parameters
    ----------
    eigenvalue_hessian : ndarray
        Eigenvalue Hessian with shape ``(3, 3)``.
    ray_direction : ndarray
        Unit group-velocity direction with shape ``(3,)``.
    phase_speed, group_speed : float
        Phase and group speeds in km s^-1.

    Returns
    -------
    ndarray
        Ray-direction gradient with shape ``(3, 3)``.
    """
    projection = np.eye(3, dtype=float) - np.outer(ray_direction, ray_direction)
    return np.asarray(
        projection @ eigenvalue_hessian / (2.0 * phase_speed * group_speed),
        dtype=float,
    )


def calculate_area_factor(
    ray_direction_gradient: NDArray[np.float64],
    wave_normal: NDArray[np.float64],
) -> tuple[float, float]:
    """Calculate the enhancement denominator and its numerical scale.

    Parameters
    ----------
    ray_direction_gradient : ndarray
        Gradient of the unit ray direction with shape ``(3, 3)``.
    wave_normal : ndarray
        Unit wave-normal direction with shape ``(3,)``.

    Returns
    -------
    tuple of float
        Area factor and Frobenius norm of the associated cofactor matrix.
    """
    cofactor = cofactor_matrix(ray_direction_gradient)
    mapped = cofactor @ wave_normal
    return float(np.linalg.norm(mapped)), float(np.linalg.norm(cofactor))


def cofactor_matrix(matrix: NDArray[np.float64]) -> NDArray[np.float64]:
    """Calculate the cofactor matrix of a real 3 x 3 matrix.

    Parameters
    ----------
    matrix : ndarray
        Matrix with shape ``(3, 3)``.

    Returns
    -------
    ndarray
        Cofactor matrix with shape ``(3, 3)``.

    Raises
    ------
    ValueError
        If the input does not have shape ``(3, 3)``.
    """
    value = np.asarray(matrix, dtype=float)
    if value.shape != (3, 3):
        raise ValueError("A cofactor matrix requires an array with shape (3, 3).")

    result = np.empty((3, 3), dtype=float)
    result[0, 0] = value[1, 1] * value[2, 2] - value[1, 2] * value[2, 1]
    result[0, 1] = value[1, 2] * value[2, 0] - value[1, 0] * value[2, 2]
    result[0, 2] = value[1, 0] * value[2, 1] - value[1, 1] * value[2, 0]
    result[1, 0] = value[0, 2] * value[2, 1] - value[0, 1] * value[2, 2]
    result[1, 1] = value[0, 0] * value[2, 2] - value[0, 2] * value[2, 0]
    result[1, 2] = value[0, 1] * value[2, 0] - value[0, 0] * value[2, 1]
    result[2, 0] = value[0, 1] * value[1, 2] - value[0, 2] * value[1, 1]
    result[2, 1] = value[0, 2] * value[1, 0] - value[0, 0] * value[1, 2]
    result[2, 2] = value[0, 0] * value[1, 1] - value[0, 1] * value[1, 0]
    return result
