# -*- coding: utf-8 -*-

"""Voigt, Reuss, and Hill averages for elastic tensors."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .tensor import ElasticTensor


@dataclass(frozen=True, slots=True)
class IsotropicElasticProperties:
    """Passive container for one isotropic elastic estimate.

    Parameters
    ----------
    bulk_modulus : float
        Bulk modulus in GPa.
    young_modulus : float
        Young's modulus in GPa.
    shear_modulus : float
        Shear modulus in GPa.
    poisson_ratio : float
        Dimensionless Poisson ratio.
    """

    bulk_modulus: float
    young_modulus: float
    shear_modulus: float
    poisson_ratio: float

    def as_array(self) -> np.ndarray:
        """Return values in ``K, E, G, nu`` order."""
        return np.asarray(
            [
                self.bulk_modulus,
                self.young_modulus,
                self.shear_modulus,
                self.poisson_ratio,
            ],
            dtype=float,
        )


@dataclass(frozen=True, slots=True)
class ElasticAverages:
    """Voigt, Reuss, and Hill isotropic elastic estimates."""

    voigt: IsotropicElasticProperties
    reuss: IsotropicElasticProperties
    hill: IsotropicElasticProperties

    def as_array(self) -> np.ndarray:
        """Return a ``(3, 4)`` array ordered as Voigt, Reuss, and Hill."""
        return np.vstack(
            [
                self.voigt.as_array(),
                self.reuss.as_array(),
                self.hill.as_array(),
            ]
        )


def compute_elastic_averages(tensor: ElasticTensor) -> ElasticAverages:
    """Calculate Voigt, Reuss, and Hill elastic averages.

    Parameters
    ----------
    tensor : ElasticTensor
        Elastic tensor to average.

    Returns
    -------
    ElasticAverages
        Voigt, Reuss, and Hill estimates.
    """
    stiffness = tensor.stiffness
    compliance = tensor.compliance

    kv = float(np.einsum("iijj", tensor.stiffness_tensor) / 9.0)
    kr = float(1.0 / np.einsum("iijj", tensor.compliance_tensor))
    kh = 0.5 * (kv + kr)

    gv = float(
        (
            stiffness[0, 0]
            + stiffness[1, 1]
            + stiffness[2, 2]
            - (stiffness[0, 1] + stiffness[0, 2] + stiffness[1, 2])
            + 3.0 * (stiffness[3, 3] + stiffness[4, 4] + stiffness[5, 5])
        )
        / 15.0
    )

    gr = float(
        15.0
        / (
            4.0 * (compliance[0, 0] + compliance[1, 1] + compliance[2, 2])
            - 4.0 * (compliance[0, 1] + compliance[0, 2] + compliance[1, 2])
            + 3.0 * (compliance[3, 3] + compliance[4, 4] + compliance[5, 5])
        )
    )
    gh = 0.5 * (gv + gr)

    return ElasticAverages(
        voigt=_isotropic_properties(kv, gv),
        reuss=_isotropic_properties(kr, gr),
        hill=_isotropic_properties(kh, gh),
    )


def _isotropic_properties(
    bulk_modulus: float,
    shear_modulus: float,
) -> IsotropicElasticProperties:
    """Build isotropic properties from bulk and shear moduli."""
    young_modulus = 1.0 / (1.0 / (3.0 * shear_modulus) + 1.0 / (9.0 * bulk_modulus))
    poisson_ratio = (
        1.0 - 3.0 * shear_modulus / (3.0 * bulk_modulus + shear_modulus)
    ) / 2.0
    return IsotropicElasticProperties(
        bulk_modulus=float(bulk_modulus),
        young_modulus=float(young_modulus),
        shear_modulus=float(shear_modulus),
        poisson_ratio=float(poisson_ratio),
    )
