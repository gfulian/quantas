"""Experimental batched SEISMIC kernels used only for performance studies."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from tests.reference.seismic_reference import reference_stiffness_tensor


@dataclass(frozen=True, slots=True)
class BatchedSeismicResult:
    """Numerical result from the experimental batched kernel.

    Parameters
    ----------
    directions : numpy.ndarray
        Normalized directions with shape ``(n, 3)``.
    christoffel : numpy.ndarray
        Christoffel matrices with shape ``(n, 3, 3)``.
    eigenvalues : numpy.ndarray
        Ascending eigenvalues with shape ``(n, 3)``.
    polarizations : numpy.ndarray
        Row-oriented polarizations with shape ``(n, 3, 3)``.
    phase_speeds : numpy.ndarray
        Phase speeds with shape ``(n, 3)``.
    eigenvalue_gradients : numpy.ndarray
        Eigenvalue gradients with shape ``(n, 3, 3)``.
    group_velocities : numpy.ndarray
        Group-velocity vectors with shape ``(n, 3, 3)``.
    group_speeds : numpy.ndarray
        Group-speed magnitudes with shape ``(n, 3)``.
    group_directions : numpy.ndarray
        Unit group directions with shape ``(n, 3, 3)``.
    power_flow_angles : numpy.ndarray
        Power-flow angles in radians with shape ``(n, 3)``.
    eigenvalue_hessians : numpy.ndarray or None
        Eigenvalue Hessians when enhancement was requested.
    enhancement : numpy.ndarray or None
        Raw enhancement factors when requested.
    """

    directions: NDArray[np.float64]
    christoffel: NDArray[np.float64]
    eigenvalues: NDArray[np.float64]
    polarizations: NDArray[np.float64]
    phase_speeds: NDArray[np.float64]
    eigenvalue_gradients: NDArray[np.float64]
    group_velocities: NDArray[np.float64]
    group_speeds: NDArray[np.float64]
    group_directions: NDArray[np.float64]
    power_flow_angles: NDArray[np.float64]
    eigenvalue_hessians: NDArray[np.float64] | None
    enhancement: NDArray[np.float64] | None


def _batch_cofactor(matrices: NDArray[np.float64]) -> NDArray[np.float64]:
    """Return explicit cofactors for a stack of 3x3 matrices.

    Parameters
    ----------
    matrices : numpy.ndarray
        Matrix stack with final shape ``(..., 3, 3)``.

    Returns
    -------
    numpy.ndarray
        Cofactor matrices with the same shape.
    """
    a = matrices
    result = np.empty_like(a)
    result[..., 0, 0] = a[..., 1, 1] * a[..., 2, 2] - a[..., 1, 2] * a[..., 2, 1]
    result[..., 0, 1] = a[..., 1, 2] * a[..., 2, 0] - a[..., 1, 0] * a[..., 2, 2]
    result[..., 0, 2] = a[..., 1, 0] * a[..., 2, 1] - a[..., 1, 1] * a[..., 2, 0]
    result[..., 1, 0] = a[..., 0, 2] * a[..., 2, 1] - a[..., 0, 1] * a[..., 2, 2]
    result[..., 1, 1] = a[..., 0, 0] * a[..., 2, 2] - a[..., 0, 2] * a[..., 2, 0]
    result[..., 1, 2] = a[..., 0, 1] * a[..., 2, 0] - a[..., 0, 0] * a[..., 2, 1]
    result[..., 2, 0] = a[..., 0, 1] * a[..., 1, 2] - a[..., 0, 2] * a[..., 1, 1]
    result[..., 2, 1] = a[..., 0, 2] * a[..., 1, 0] - a[..., 0, 0] * a[..., 1, 2]
    result[..., 2, 2] = a[..., 0, 0] * a[..., 1, 1] - a[..., 0, 1] * a[..., 1, 0]
    return result


def _christoffel_hessian(
    reduced_stiffness: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Return the direction-independent Christoffel Hessian.

    Parameters
    ----------
    reduced_stiffness : numpy.ndarray
        Reduced stiffness tensor in km2 s-2.

    Returns
    -------
    numpy.ndarray
        Hessian with shape ``(3, 3, 3, 3)``.
    """
    return np.transpose(reduced_stiffness, (1, 2, 0, 3)) + np.transpose(
        reduced_stiffness, (2, 1, 0, 3)
    )


def solve_batched(
    stiffness: ArrayLike,
    density: float,
    directions: ArrayLike,
    *,
    calculate_enhancement: bool = True,
) -> BatchedSeismicResult:
    """Solve independent Christoffel problems in one NumPy batch.

    Parameters
    ----------
    stiffness : array_like
        Stiffness matrix in GPa.
    density : float
        Density in kg m-3.
    directions : array_like
        Cartesian directions with shape ``(n, 3)``.
    calculate_enhancement : bool, optional
        Whether to evaluate eigenvalue Hessians and raw enhancement factors.

    Returns
    -------
    BatchedSeismicResult
        Batched numerical result.

    Raises
    ------
    ValueError
        If directions do not have shape ``(n, 3)`` or contain zero vectors.
    """
    q = np.asarray(directions, dtype=float)
    if q.ndim != 2 or q.shape[1] != 3:
        raise ValueError("Directions must have shape (n, 3).")
    norms = np.linalg.norm(q, axis=1)
    if np.any(norms == 0.0):
        raise ValueError("Directions must be non-zero.")
    q = q / norms[:, None]

    if not np.isfinite(density) or density <= 0.0:
        raise ValueError("Density must be finite and positive.")
    reduced = np.asarray(
        reference_stiffness_tensor(stiffness) * (1000.0 / density), dtype=float
    )
    christoffel = np.einsum("na,iajb,nb->nij", q, reduced, q, optimize=True)
    eigenvalues, eigenvectors_columns = np.linalg.eigh(christoffel)
    polarizations = np.swapaxes(eigenvectors_columns, 1, 2)
    phase_speeds = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues))

    symmetrized = reduced + np.transpose(reduced, (0, 2, 1, 3))
    christoffel_gradient = np.einsum("nk,iakl->nail", q, symmetrized, optimize=True)
    eigenvalue_gradients = np.einsum(
        "nmi,naij,nmj->nma",
        polarizations,
        christoffel_gradient,
        polarizations,
        optimize=True,
    )
    group_velocities = eigenvalue_gradients / (2.0 * phase_speeds[..., None])
    group_speeds = np.linalg.norm(group_velocities, axis=2)
    group_directions = group_velocities / group_speeds[..., None]
    cos_power_flow = np.einsum("nmi,ni->nm", group_directions, q, optimize=True)
    power_flow_angles = np.arccos(np.around(cos_power_flow, 10))

    eigenvalue_hessians = None
    enhancement = None
    if calculate_enhancement:
        hessian_matrix = _christoffel_hessian(reduced)
        direct = np.einsum(
            "nmi,abij,nmj->nmab",
            polarizations,
            hessian_matrix,
            polarizations,
            optimize=True,
        )
        derivative_vectors = np.einsum(
            "naij,nmj->nmai",
            christoffel_gradient,
            polarizations,
            optimize=True,
        )
        shifted = (
            eigenvalues[..., None, None] * np.eye(3)[None, None, :, :]
            - christoffel[:, None, :, :]
        )
        inverse = np.linalg.pinv(shifted, rcond=1.0e-10)
        correction = 2.0 * np.einsum(
            "nmai,nmij,nmbj->nmab",
            derivative_vectors,
            inverse,
            derivative_vectors,
            optimize=True,
        )
        eigenvalue_hessians = direct + correction

        hessian_times_group = np.einsum(
            "nmij,nmj->nmi", eigenvalue_hessians, group_velocities, optimize=True
        )
        direction_gradients = eigenvalue_hessians / group_speeds[..., None, None]
        direction_gradients -= (
            group_velocities[..., :, None]
            * hessian_times_group[..., None, :]
            / group_speeds[..., None, None] ** 3
        )
        direction_gradients /= 2.0 * phase_speeds[..., None, None]
        cofactors = _batch_cofactor(direction_gradients)
        mapped = np.einsum("nmij,nj->nmi", cofactors, q, optimize=True)
        enhancement = 1.0 / np.linalg.norm(mapped, axis=2)

    return BatchedSeismicResult(
        directions=q,
        christoffel=christoffel,
        eigenvalues=eigenvalues,
        polarizations=polarizations,
        phase_speeds=phase_speeds,
        eigenvalue_gradients=eigenvalue_gradients,
        group_velocities=group_velocities,
        group_speeds=group_speeds,
        group_directions=group_directions,
        power_flow_angles=power_flow_angles,
        eigenvalue_hessians=eigenvalue_hessians,
        enhancement=enhancement,
    )
