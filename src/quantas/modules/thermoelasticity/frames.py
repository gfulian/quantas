# -*- coding: utf-8 -*-

"""Co-rotate sampled elastic tensors into one reference Cartesian frame.

A volume series can contain mathematically equivalent structures whose printed
Cartesian axes differ by a rigid rotation.  Component-wise fitting is meaningful
only after every stiffness matrix is expressed in one common frame.  Quantas
therefore decomposes the lattice deformation gradient relative to the selected
reference state as ``F = R U`` and removes the proper rotation ``R`` from both
lattice and stiffness components while retaining the symmetric stretch ``U``.

The finite-pressure stiffness matrices handled here are the hydrostatic
stress--strain coefficients discussed by Barron and Klein (1965) and Wallace
(1972).  Co-rotation changes only their Cartesian representation; it does not
change the physical stressed state.

References
----------
Canonical citation keys: ``barron_klein_1965`` and ``wallace_1972``.

"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import polar

from quantas.core.physics.elasticity import rotate_voigt_stiffness


FloatArray: TypeAlias = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class ElasticFrameNormalization:
    """Result of co-rotating one elastic point into a reference frame.

    Parameters
    ----------
    lattice : ndarray
        Rotation-free lattice matrix with vectors stored by rows.
    stiffness : ndarray
        Co-rotated ``6 x 6`` Wallace stiffness matrix in GPa.
    rotation_to_reference : ndarray
        Proper active rotation applied to the sampled stiffness tensor.
    removed_rotation_degrees : float
        Magnitude of the rigid rotation removed from the sampled lattice.
    principal_logarithmic_strain : ndarray
        Principal logarithmic stretches retained after co-rotation.
    metadata : dict
        Determinant, volume, and numerical validation diagnostics.
    """

    lattice: FloatArray
    stiffness: FloatArray
    rotation_to_reference: FloatArray
    removed_rotation_degrees: float
    principal_logarithmic_strain: FloatArray
    metadata: dict[str, Any]


def normalize_elastic_frame(
    lattice: ArrayLike,
    stiffness: ArrayLike,
    reference_lattice: ArrayLike,
    *,
    determinant_tolerance: float = 1.0e-10,
) -> ElasticFrameNormalization:
    """Express one sampled lattice and stiffness tensor in a reference frame.

    Parameters
    ----------
    lattice : array_like
        Sampled direct-lattice matrix with vectors stored by rows.
    stiffness : array_like
        Sampled Wallace stiffness matrix in Voigt notation and GPa.
    reference_lattice : array_like
        Reference direct-lattice matrix with vectors stored by rows.
    determinant_tolerance : float, optional
        Absolute tolerance used when validating proper rotations and volume
        preservation.

    Returns
    -------
    ElasticFrameNormalization
        Co-rotated lattice, stiffness, and complete diagnostics.

    Raises
    ------
    ValueError
        If a lattice is singular, the polar factor is improper, or co-rotation
        fails to preserve the sampled volume.

    Notes
    -----
    With row-wise direct-lattice vectors, the deformation gradient is

    ``F = A(V).T @ inv(A_ref.T)``.

    Its right polar decomposition ``F = R U`` separates a proper rigid rotation
    from the symmetric positive-definite stretch.  The returned lattice is
    ``(U @ A_ref.T).T`` and the stiffness is actively rotated by ``R.T``.
    """
    sampled = _validated_lattice(lattice, "lattice")
    reference = _validated_lattice(reference_lattice, "reference_lattice")
    matrix = np.asarray(stiffness, dtype=np.float64)
    if matrix.shape != (6, 6) or np.any(~np.isfinite(matrix)):
        raise ValueError("stiffness must be finite with shape (6, 6)")

    deformation = sampled.T @ np.linalg.inv(reference.T)
    rotation, stretch = polar(deformation, side="right")
    rotation = np.asarray(rotation, dtype=np.float64)
    stretch = np.asarray(stretch, dtype=np.float64)
    determinant = float(np.linalg.det(rotation))
    if not np.isclose(determinant, 1.0, rtol=0.0, atol=determinant_tolerance):
        raise ValueError(
            "elastic lattice frames are related by an improper or singular "
            f"transformation (det R = {determinant:.8g})"
        )
    eigenvalues = np.linalg.eigvalsh(0.5 * (stretch + stretch.T))
    if np.any(eigenvalues <= 0.0) or np.any(~np.isfinite(eigenvalues)):
        raise ValueError("elastic lattice stretch must be positive definite")

    rotation_to_reference = rotation.T
    normalized_lattice = (stretch @ reference.T).T
    normalized_stiffness = rotate_voigt_stiffness(matrix, rotation_to_reference)
    sampled_volume = abs(float(np.linalg.det(sampled)))
    normalized_volume = abs(float(np.linalg.det(normalized_lattice)))
    volume_error = abs(normalized_volume - sampled_volume)
    allowed = max(determinant_tolerance, determinant_tolerance * sampled_volume)
    if volume_error > allowed:
        raise ValueError(
            "elastic frame normalization did not preserve the sampled volume"
        )

    angle = _rotation_angle_degrees(rotation)
    principal_logarithmic_strain = np.log(eigenvalues)
    return ElasticFrameNormalization(
        lattice=np.asarray(normalized_lattice, dtype=np.float64),
        stiffness=np.asarray(normalized_stiffness, dtype=np.float64),
        rotation_to_reference=np.asarray(rotation_to_reference, dtype=np.float64),
        removed_rotation_degrees=angle,
        principal_logarithmic_strain=np.asarray(
            principal_logarithmic_strain, dtype=np.float64
        ),
        metadata={
            "method": "right_polar_decomposition_corotation",
            "rotation_determinant": determinant,
            "sampled_volume_A3": sampled_volume,
            "normalized_volume_A3": normalized_volume,
            "volume_preservation_error_A3": volume_error,
            "maximum_absolute_principal_logarithmic_strain": float(
                np.max(np.abs(principal_logarithmic_strain))
            ),
        },
    )


def maximum_ordered_fractional_displacement(
    reference_fractional: ArrayLike,
    sampled_fractional: ArrayLike,
    reference_lattice: ArrayLike,
) -> float:
    """Return the largest ordered-atom displacement modulo lattice periods.

    Parameters
    ----------
    reference_fractional, sampled_fractional : array_like
        Fractional coordinates with identical ``(natoms, 3)`` shapes and atom
        ordering.
    reference_lattice : array_like
        Reference row-wise lattice used to express the wrapped displacement in
        angstrom.

    Returns
    -------
    float
        Maximum Cartesian displacement in angstrom.

    Raises
    ------
    ValueError
        If coordinate arrays are incompatible.
    """
    reference = np.asarray(reference_fractional, dtype=np.float64)
    sampled = np.asarray(sampled_fractional, dtype=np.float64)
    lattice = _validated_lattice(reference_lattice, "reference_lattice")
    if (
        reference.shape != sampled.shape
        or reference.ndim != 2
        or reference.shape[1] != 3
    ):
        raise ValueError("fractional coordinate arrays must share shape (natoms, 3)")
    if np.any(~np.isfinite(reference)) or np.any(~np.isfinite(sampled)):
        raise ValueError("fractional coordinate arrays must be finite")
    delta = sampled - reference
    delta -= np.rint(delta)
    cartesian = delta @ lattice
    return float(np.max(np.linalg.norm(cartesian, axis=1), initial=0.0))


def _validated_lattice(value: ArrayLike, name: str) -> FloatArray:
    """Return a finite nonsingular row-wise lattice matrix."""
    lattice = np.asarray(value, dtype=np.float64)
    if lattice.shape != (3, 3) or np.any(~np.isfinite(lattice)):
        raise ValueError(f"{name} must be finite with shape (3, 3)")
    if abs(float(np.linalg.det(lattice))) <= np.finfo(np.float64).eps:
        raise ValueError(f"{name} must be nonsingular")
    return lattice


def _rotation_angle_degrees(rotation: FloatArray) -> float:
    """Return the angle of a proper rotation in degrees."""
    cosine = np.clip((float(np.trace(rotation)) - 1.0) / 2.0, -1.0, 1.0)
    return float(np.degrees(np.arccos(cosine)))


__all__ = [
    "ElasticFrameNormalization",
    "maximum_ordered_fractional_displacement",
    "normalize_elastic_frame",
]
