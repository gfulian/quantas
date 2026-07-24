"""Frozen numerical reference for the Quantas 0.9 SEISMIC implementation.

This module reproduces the numerical formulas used by the historical
``quantas.seismic.utils.seismic_obj.Seismic`` class.  It is deliberately kept
inside the test tree: it is a regression oracle, not production code and not a
backward-compatibility layer.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np
from numpy.typing import ArrayLike, NDArray

_VOIGT_MAP = np.asarray(
    [
        [0, 5, 4],
        [5, 1, 3],
        [4, 3, 2],
    ],
    dtype=int,
)


def reference_stiffness_tensor(stiffness: ArrayLike) -> NDArray[np.float64]:
    """Convert a 6x6 Voigt matrix using the original stiffness mapping.

    Parameters
    ----------
    stiffness : array_like
        Elastic stiffness matrix in GPa.

    Returns
    -------
    numpy.ndarray
        Cartesian fourth-rank stiffness tensor with shape ``(3, 3, 3, 3)``.

    Raises
    ------
    ValueError
        If the input is not a finite symmetric 6x6 matrix.
    """
    matrix = np.asarray(stiffness, dtype=float)
    if matrix.shape != (6, 6):
        raise ValueError("The stiffness matrix must have shape (6, 6).")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("The stiffness matrix must contain finite values.")
    if not np.allclose(matrix, matrix.T, rtol=0.0, atol=1.0e-12):
        raise ValueError("The stiffness matrix must be symmetric.")
    indexes = _VOIGT_MAP
    return matrix[indexes[:, :, None, None], indexes[None, None, :, :]]


class ReferenceWaveMode(str, Enum):
    """Physical wave-mode names in the original ascending eigenvalue order."""

    V_S2 = "v_s2"
    V_S1 = "v_s1"
    V_P = "v_p"


MODE_ORDER = (
    ReferenceWaveMode.V_S2,
    ReferenceWaveMode.V_S1,
    ReferenceWaveMode.V_P,
)

REFERENCE_MODE_NAMES = {
    ReferenceWaveMode.V_S2: "slow_secondary",
    ReferenceWaveMode.V_S1: "fast_secondary",
    ReferenceWaveMode.V_P: "primary",
}


def spherical_direction(theta: float, phi: float) -> NDArray[np.float64]:
    """Return the original Cartesian direction for spherical angles.

    Parameters
    ----------
    theta : float
        Polar angle in radians.
    phi : float
        Azimuth angle in radians.

    Returns
    -------
    numpy.ndarray
        Unit direction with shape ``(3,)``.
    """
    return np.asarray(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ],
        dtype=float,
    )


def original_cartesian_angles(direction: ArrayLike) -> tuple[float, float]:
    """Return angles using the original, reversed ``atan2(x, y)`` azimuth.

    Parameters
    ----------
    direction : array_like
        Cartesian direction.

    Returns
    -------
    tuple of float
        Historical ``(theta, phi)`` angles in radians.
    """
    q = np.asarray(direction, dtype=float)
    q = q / np.linalg.norm(q)
    x, y, z = q
    phi = np.arctan2(x, y)
    theta = np.arctan2(np.sqrt(x * x + y * y), z)
    return float(theta), float(phi)


def physical_cartesian_angles(direction: ArrayLike) -> tuple[float, float]:
    """Return conventional polar and azimuth angles for a Cartesian direction.

    Parameters
    ----------
    direction : array_like
        Cartesian direction.

    Returns
    -------
    tuple of float
        Conventional ``(theta, phi)`` angles in radians.
    """
    q = np.asarray(direction, dtype=float)
    q = q / np.linalg.norm(q)
    x, y, z = q
    phi = np.arctan2(y, x)
    theta = np.arctan2(np.sqrt(x * x + y * y), z)
    return float(theta), float(phi)


def reference_cofactor(matrix: ArrayLike) -> NDArray[np.float64]:
    """Return the 3x3 cofactor matrix following the frozen branch logic.

    The explicit expression remains valid for singular and nearly singular
    matrices and avoids backend-dependent ``inv(A).T * det(A)`` evaluation.

    Parameters
    ----------
    matrix : array_like
        Square matrix with shape ``(3, 3)``.

    Returns
    -------
    numpy.ndarray
        Cofactor matrix.
    """
    value = np.asarray(matrix, dtype=float)

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


@dataclass(frozen=True, slots=True)
class DirectionalReferenceResult:
    """Complete direction-dependent result produced by the original Quantas formulas.

    Parameters
    ----------
    direction : numpy.ndarray
        Normalized wave-normal direction.
    christoffel : numpy.ndarray
        Reduced Christoffel matrix in km2 s-2.
    christoffel_gradient : numpy.ndarray
        Gradient of the reduced Christoffel matrix.
    christoffel_hessian : numpy.ndarray
        Hessian of the reduced Christoffel matrix.
    eigenvalues : numpy.ndarray
        Ascending Christoffel eigenvalues in km2 s-2.
    polarizations : numpy.ndarray
        Row-oriented eigenvectors in ``V_S2, V_S1, V_P`` order.
    phase_speeds : numpy.ndarray
        Phase speeds in km s-1.
    eigenvalue_gradients : numpy.ndarray
        Gradients of the three eigenvalues.
    eigenvalue_hessians : numpy.ndarray
        Hessians of the three eigenvalues.
    group_velocities : numpy.ndarray
        Group-velocity vectors in km s-1.
    group_speeds : numpy.ndarray
        Group-speed magnitudes in km s-1.
    group_directions : numpy.ndarray
        Unit group directions.
    power_flow_angles : numpy.ndarray
        Power-flow angles in radians.
    enhancement : numpy.ndarray
        Raw enhancement factors ``A``.
    """

    direction: NDArray[np.float64]
    christoffel: NDArray[np.float64]
    christoffel_gradient: NDArray[np.float64]
    christoffel_hessian: NDArray[np.float64]
    eigenvalues: NDArray[np.float64]
    polarizations: NDArray[np.float64]
    phase_speeds: NDArray[np.float64]
    eigenvalue_gradients: NDArray[np.float64]
    eigenvalue_hessians: NDArray[np.float64]
    group_velocities: NDArray[np.float64]
    group_speeds: NDArray[np.float64]
    group_directions: NDArray[np.float64]
    power_flow_angles: NDArray[np.float64]
    enhancement: NDArray[np.float64]


class SeismicFormulaReference:
    """Reproduce Quantas 0.9 SEISMIC numerical operations for testing.

    Parameters
    ----------
    stiffness : array_like
        Elastic stiffness matrix in GPa.
    density : float
        Density in kg m-3.

    Raises
    ------
    ValueError
        If the density is non-finite or not positive.
    """

    def __init__(self, stiffness: ArrayLike, density: float) -> None:
        if not np.isfinite(density) or density <= 0.0:
            raise ValueError("Density must be finite and positive.")
        self.stiffness = np.asarray(stiffness, dtype=float)
        stiffness_tensor = reference_stiffness_tensor(self.stiffness)
        self.density = float(density)
        self.reduced_stiffness = np.asarray(
            stiffness_tensor * (1000.0 / self.density),
            dtype=float,
        )
        self.christoffel_hessian = self._christoffel_hessian(self.reduced_stiffness)

    @staticmethod
    def _christoffel_hessian(
        reduced_stiffness: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        result = np.empty((3, 3, 3, 3), dtype=float)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for ell in range(3):
                        result[i, j, k, ell] = (
                            reduced_stiffness[k, i, j, ell]
                            + reduced_stiffness[k, j, i, ell]
                        )
        return result

    def christoffel_matrix(self, direction: ArrayLike) -> NDArray[np.float64]:
        """Return the frozen reference Christoffel matrix.

        Parameters
        ----------
        direction : array_like
            Wave-normal direction. The input is normalized internally.

        Returns
        -------
        numpy.ndarray
            Christoffel matrix in km2 s-2.
        """
        q = np.asarray(direction, dtype=float)
        q = q / np.linalg.norm(q)
        return np.dot(q, np.dot(q, self.reduced_stiffness))

    def christoffel_gradient(self, direction: ArrayLike) -> NDArray[np.float64]:
        """Return the analytical gradient of the Christoffel matrix.

        Parameters
        ----------
        direction : array_like
            Wave-normal direction. The input is normalized internally.

        Returns
        -------
        numpy.ndarray
            Array with shape ``(3, 3, 3)``; the first axis is the derivative
            coordinate.
        """
        q = np.asarray(direction, dtype=float)
        q = q / np.linalg.norm(q)
        symmetrized = self.reduced_stiffness + np.transpose(
            self.reduced_stiffness, (0, 2, 1, 3)
        )
        return np.transpose(np.dot(q, symmetrized), (1, 0, 2))

    def solve(self, direction: ArrayLike) -> DirectionalReferenceResult:
        """Calculate all reference seismic quantities for one direction.

        Parameters
        ----------
        direction : array_like
            Wave-normal direction.

        Returns
        -------
        DirectionalReferenceResult
            Complete frozen numerical result.
        """
        q = np.asarray(direction, dtype=float)
        q = q / np.linalg.norm(q)
        christoffel = np.dot(q, np.dot(q, self.reduced_stiffness))
        gradient = self.christoffel_gradient(q)

        eigenvalues, eigenvectors_columns = np.linalg.eigh(christoffel)
        sorting = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[sorting]
        polarizations = eigenvectors_columns.T[sorting]
        phase_speeds = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues))

        eigenvalue_gradients = np.empty((3, 3), dtype=float)
        group_velocities = np.empty((3, 3), dtype=float)
        for mode in range(3):
            for coordinate in range(3):
                eigenvalue_gradients[mode, coordinate] = np.dot(
                    polarizations[mode],
                    np.dot(gradient[coordinate], polarizations[mode]),
                )
                group_velocities[mode, coordinate] = eigenvalue_gradients[
                    mode, coordinate
                ] / (2.0 * phase_speeds[mode])

        group_speeds = np.linalg.norm(group_velocities, axis=1)
        group_directions = group_velocities / group_speeds[:, None]
        cos_power_flow = np.dot(group_directions, q)
        power_flow_angles = np.arccos(np.around(cos_power_flow, 10))

        eigenvalue_hessians = np.zeros((3, 3, 3), dtype=float)
        for mode in range(3):
            polarization = polarizations[mode]
            eigenvalue_hessians[mode] += np.dot(
                np.dot(self.christoffel_hessian, polarization),
                polarization,
            )
            inverse = np.linalg.pinv(
                eigenvalues[mode] * np.identity(3) - christoffel,
                rcond=1.0e-10,
            )
            derivative_vectors = np.dot(gradient, polarization)
            eigenvalue_hessians[mode] += 2.0 * np.dot(
                np.dot(derivative_vectors, inverse),
                derivative_vectors.T,
            )

        direction_gradients = np.empty((3, 3, 3), dtype=float)
        enhancement = np.empty(3, dtype=float)
        for mode in range(3):
            hessian = eigenvalue_hessians[mode]
            group_velocity = group_velocities[mode]
            group_speed = group_speeds[mode]
            direction_gradients[mode] = hessian / group_speed
            direction_gradients[mode] -= np.outer(
                group_velocity,
                np.dot(hessian, group_velocity),
            ) / (group_speed**3)
            direction_gradients[mode] /= 2.0 * phase_speeds[mode]
            enhancement[mode] = 1.0 / np.linalg.norm(
                np.dot(reference_cofactor(direction_gradients[mode]), q)
            )

        return DirectionalReferenceResult(
            direction=q,
            christoffel=christoffel,
            christoffel_gradient=gradient,
            christoffel_hessian=self.christoffel_hessian,
            eigenvalues=eigenvalues,
            polarizations=polarizations,
            phase_speeds=phase_speeds,
            eigenvalue_gradients=eigenvalue_gradients,
            eigenvalue_hessians=eigenvalue_hessians,
            group_velocities=group_velocities,
            group_speeds=group_speeds,
            group_directions=group_directions,
            power_flow_angles=power_flow_angles,
            enhancement=enhancement,
        )
