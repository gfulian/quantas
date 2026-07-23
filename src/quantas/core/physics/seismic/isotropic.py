# -*- coding: utf-8 -*-

"""Isotropic seismic reference velocities derived from elastic averages."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from quantas.core.physics.elasticity.averages import ElasticAverages


@dataclass(frozen=True, slots=True)
class IsotropicSeismicVelocities:
    """Isotropic shear and compressional wave velocities.

    Parameters
    ----------
    shear : float
        Isotropic shear-wave velocity in km s^-1.
    compressional : float
        Isotropic compressional-wave velocity in km s^-1.
    """

    shear: float
    compressional: float

    def as_array(self) -> np.ndarray:
        """Return velocities in ``Vs, Vs, Vp`` order.

        Returns
        -------
        ndarray
            Velocity vector with shape ``(3,)`` in km s^-1.
        """
        return np.asarray([self.shear, self.shear, self.compressional], dtype=float)


def isotropic_seismic_velocities(
    averages: ElasticAverages,
    density: float,
) -> IsotropicSeismicVelocities:
    """Calculate isotropic seismic velocities from Hill averages.

    Parameters
    ----------
    averages : ElasticAverages
        Voigt-Reuss-Hill elastic averages.
    density : float
        Material density in kg m^-3.

    Returns
    -------
    IsotropicSeismicVelocities
        Isotropic shear and compressional velocities in km s^-1.

    Raises
    ------
    ValueError
        If ``density`` is non-finite or not strictly positive.
    """
    density_value = float(density)
    if not np.isfinite(density_value) or density_value <= 0.0:
        raise ValueError(
            "A finite positive density is required to calculate seismic velocities."
        )

    bulk_modulus = averages.hill.bulk_modulus
    shear_modulus = averages.hill.shear_modulus
    compressional = np.sqrt(
        1000.0 * (bulk_modulus + 4.0 * shear_modulus / 3.0) / density_value
    )
    shear = np.sqrt(1000.0 * shear_modulus / density_value)
    return IsotropicSeismicVelocities(
        shear=float(shear),
        compressional=float(compressional),
    )
