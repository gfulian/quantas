"""Frozen numerical reference for the original Quantas SOEC formulas.

This module reproduces the pointwise tensor conversion and directional
properties of the original ``quantas.soec.utils.soec_obj.SOEC`` class.  It is
kept inside the test tree deliberately: it is a regression oracle, not public
library code and not a backward-compatibility layer.
"""

from __future__ import annotations

import itertools as it

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


def reference_direction_vector(
    theta: float,
    phi: float,
    chi: float | None = None,
) -> NDArray[np.float64]:
    """Return the original Quantas directional vector.

    Parameters
    ----------
    theta, phi : float
        Polar and azimuthal angles in radians.
    chi : float or None, optional
        Transverse rotation angle in radians.

    Returns
    -------
    numpy.ndarray
        Unit Cartesian vector with shape ``(3,)``.
    """
    if chi is None:
        vector = (
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        )
    else:
        vector = (
            np.cos(theta) * np.cos(phi) * np.cos(chi) - np.sin(phi) * np.sin(chi),
            np.cos(theta) * np.sin(phi) * np.cos(chi) + np.cos(phi) * np.sin(chi),
            -np.sin(theta) * np.cos(chi),
        )
    return np.asarray(vector, dtype=float)


def reference_compliance_tensor(
    stiffness: ArrayLike,
) -> NDArray[np.float64]:
    """Build the original Cartesian compliance tensor.

    Parameters
    ----------
    stiffness : array_like
        Symmetric elastic stiffness matrix in Voigt notation, in GPa.

    Returns
    -------
    numpy.ndarray
        Compliance tensor ``S_ijkl`` with shape ``(3, 3, 3, 3)``.

    Raises
    ------
    ValueError
        If the matrix is not finite, symmetric, invertible, or ``6 x 6``.
    """
    matrix = np.asarray(stiffness, dtype=float)
    if matrix.shape != (6, 6):
        raise ValueError("The stiffness matrix must have shape (6, 6).")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("The stiffness matrix must contain finite values.")
    if not np.allclose(matrix, matrix.T, rtol=0.0, atol=1.0e-12):
        raise ValueError("The stiffness matrix must be symmetric.")
    try:
        compliance = np.linalg.inv(matrix)
    except np.linalg.LinAlgError as exc:
        raise ValueError("The stiffness matrix must be invertible.") from exc

    tensor = np.empty((3, 3, 3, 3), dtype=float)
    for i, j, k, ell in it.product((0, 1, 2), repeat=4):
        p = int(_VOIGT_MAP[i, j])
        q = int(_VOIGT_MAP[k, ell])
        coefficient = 1.0 / ((1 + p // 3) * (1 + q // 3))
        tensor[i, j, k, ell] = coefficient * compliance[p, q]
    return tensor


class ElasticityFormulaReference:
    """Reproduce the original Quantas SOEC pointwise formulas.

    Parameters
    ----------
    stiffness : array_like
        Symmetric elastic stiffness matrix in Voigt notation, in GPa.

    Raises
    ------
    ValueError
        If the matrix cannot be converted to a compliance tensor.
    """

    def __init__(self, stiffness: ArrayLike) -> None:
        self.stiffness = np.array(stiffness, dtype=float, copy=True)
        self.compliance = np.linalg.inv(self.stiffness)
        self.compliance_tensor = reference_compliance_tensor(self.stiffness)

    def young_modulus(self, theta: float, phi: float) -> float:
        """Return Young's modulus along one direction.

        Parameters
        ----------
        theta, phi : float
            Polar and azimuthal angles in radians.

        Returns
        -------
        float
            Young's modulus in GPa.
        """
        axis = reference_direction_vector(theta, phi)
        denominator = 0.0
        for i, j, k, ell in it.product((0, 1, 2), repeat=4):
            denominator += (
                self.compliance_tensor[i, j, k, ell]
                * axis[i]
                * axis[j]
                * axis[k]
                * axis[ell]
            )
        return float(1.0 / denominator)

    def linear_compressibility(self, theta: float, phi: float) -> float:
        """Return linear compressibility along one direction.

        Parameters
        ----------
        theta, phi : float
            Polar and azimuthal angles in radians.

        Returns
        -------
        float
            Linear compressibility in TPa^-1.
        """
        axis = reference_direction_vector(theta, phi)
        result = 0.0
        for i, j, k in it.product((0, 1, 2), repeat=3):
            result += axis[i] * axis[j] * self.compliance_tensor[i, j, k, k]
        return float(1000.0 * result)

    def shear_modulus(self, theta: float, phi: float, chi: float) -> float:
        """Return shear modulus for a direction and transverse axis.

        Parameters
        ----------
        theta, phi, chi : float
            Longitudinal and transverse angular coordinates in radians.

        Returns
        -------
        float
            Shear modulus in GPa.
        """
        axis = reference_direction_vector(theta, phi)
        transverse = reference_direction_vector(theta, phi, chi)
        denominator = 0.0
        for i, j, k, ell in it.product((0, 1, 2), repeat=4):
            denominator += (
                axis[i]
                * transverse[j]
                * axis[k]
                * transverse[ell]
                * self.compliance_tensor[i, j, k, ell]
            )
        return float(1.0 / (4.0 * denominator))

    def poisson_ratio(self, theta: float, phi: float, chi: float) -> float:
        """Return Poisson's ratio for a direction and transverse axis.

        Parameters
        ----------
        theta, phi, chi : float
            Longitudinal and transverse angular coordinates in radians.

        Returns
        -------
        float
            Directional Poisson ratio.
        """
        axis = reference_direction_vector(theta, phi)
        transverse = reference_direction_vector(theta, phi, chi)
        numerator = 0.0
        denominator = 0.0
        for i, j, k, ell in it.product((0, 1, 2), repeat=4):
            component = self.compliance_tensor[i, j, k, ell]
            numerator += axis[i] * axis[j] * transverse[k] * transverse[ell] * component
            denominator += axis[i] * axis[j] * axis[k] * axis[ell] * component
        return float(-numerator / denominator)
