# -*- coding: utf-8 -*-

"""Cold Eulerian finite-strain relations for quasi-static elasticity.

The functions in this module implement the volume dependence of hydrostatic
stress--strain coefficients used by the quasi-static thermoelastic workflow.
They are numerical building blocks only and do not read files, fit data, emit
Quantas events, or depend on any frontend.
"""

from __future__ import annotations

from typing import Literal, TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .conventions import VOIGT_PAIRS
from .validation import validate_stiffness_matrix


FloatArray: TypeAlias = NDArray[np.float64]
FiniteStrainOrder = Literal[2, 3]


def eulerian_finite_strain(volume: ArrayLike, reference_volume: float) -> FloatArray:
    r"""Return the isotropic Eulerian finite strain.

    The sign convention is positive on compression,

    .. math::

        f = \frac{1}{2}\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right].

    Parameters
    ----------
    volume : array_like
        Positive volume or volume array.
    reference_volume : float
        Positive reference volume :math:`V_0` in the same units.

    Returns
    -------
    ndarray
        Eulerian finite strain with the shape of ``volume`` and dtype
        ``float64``.

    Raises
    ------
    ValueError
        If a volume is non-finite or non-positive.
    """
    values = np.asarray(volume, dtype=np.float64)
    reference = float(reference_volume)
    if not np.isfinite(reference) or reference <= 0.0:
        raise ValueError("reference_volume must be finite and positive")
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("volume values must be finite and positive")
    return np.asarray(
        0.5 * (np.power(reference / values, 2.0 / 3.0) - 1.0),
        dtype=np.float64,
    )


def wallace_hydrostatic_delta_voigt() -> FloatArray:
    r"""Return Wallace's hydrostatic pre-stress tensor in Voigt form.

    The Cartesian definition is

    .. math::

        \delta^{ij}_{kl} = -\delta_{ij}\delta_{kl}
        -\delta_{il}\delta_{jk}-\delta_{jl}\delta_{ik}.

    Returns
    -------
    ndarray
        Symmetric ``(6, 6)`` matrix using Quantas' stiffness Voigt order
        ``(11, 22, 33, 23, 13, 12)``.
    """
    result: FloatArray = np.zeros((6, 6), dtype=np.float64)
    for first, (i, j) in enumerate(VOIGT_PAIRS):
        for second, (k, ell) in enumerate(VOIGT_PAIRS):
            result[first, second] = -(
                int(i == j) * int(k == ell)
                + int(i == ell) * int(j == k)
                + int(j == ell) * int(i == k)
            )
    return result


def cold_finite_strain_component(
    volume: ArrayLike,
    *,
    reference_volume: float,
    bulk_modulus: float,
    bulk_modulus_derivative: float,
    reference_component: float,
    component_pressure_derivative: float,
    wallace_delta: float,
    order: FiniteStrainOrder = 3,
) -> FloatArray:
    r"""Evaluate one cold finite-strain stiffness component.

    Parameters
    ----------
    volume : array_like
        Positive evaluation volume or array.
    reference_volume, bulk_modulus, bulk_modulus_derivative : float
        Fixed reference EOS parameters ``V0``, ``K0`` and ``Kp``.
    reference_component : float
        Component value ``C0`` at the reference state, in GPa.
    component_pressure_derivative : float
        ``dC/dP`` at the reference state.
    wallace_delta : float
        Matching Voigt entry of Wallace's hydrostatic delta tensor.
    order : {2, 3}, optional
        Finite-strain truncation order.

    Returns
    -------
    ndarray
        Component values with the shape of ``volume``.

    Notes
    -----
    ``wallace_delta`` is an analytical coefficient of the Eulerian
    finite-strain expansion for Wallace stress--strain coefficients. It does
    not apply an additional pressure correction to sampled input tensors.
    Quantas requires CRYSTAL calculations performed with ``PRESSURE`` and fits
    those already corrected coefficients directly.

    References
    ----------
    Canonical citation keys: ``stixrude_lithgow_bertelloni_2005``,
    ``barron_klein_1965``, and ``wallace_1972``.
    """
    if order not in (2, 3):
        raise ValueError("order must be 2 or 3")
    values = np.asarray(volume, dtype=np.float64)
    f = eulerian_finite_strain(values, reference_volume)
    k0 = float(bulk_modulus)
    kp = float(bulk_modulus_derivative)
    c0 = float(reference_component)
    cp = float(component_pressure_derivative)
    delta = float(wallace_delta)
    if not np.all(np.isfinite([k0, kp, c0, cp, delta])) or k0 <= 0.0:
        raise ValueError(
            "finite-strain component parameters must be finite and K0 positive"
        )
    prefactor = np.power(1.0 + 2.0 * f, 2.5)
    polynomial = c0 + (3.0 * k0 * cp - 5.0 * c0) * f
    if order == 3:
        quadratic = 6.0 * k0 * cp - 14.0 * c0 - 1.5 * k0 * delta * (3.0 * kp - 16.0)
        polynomial = polynomial + quadratic * np.square(f)
    return np.asarray(prefactor * polynomial, dtype=np.float64)


def cold_finite_strain_component_jacobian(
    volume: ArrayLike,
    *,
    reference_volume: float,
    bulk_modulus: float,
    bulk_modulus_derivative: float,
    reference_component: float,
    component_pressure_derivative: float,
    wallace_delta: float,
    order: FiniteStrainOrder = 3,
) -> FloatArray:
    r"""Return analytical derivatives of one component prediction.

    The last axis follows ``C0, Cprime, V0, K0, Kprime, V``.

    Parameters are identical to :func:`cold_finite_strain_component`.

    Returns
    -------
    ndarray
        Jacobian with shape ``volume.shape + (6,)``.
    """
    if order not in (2, 3):
        raise ValueError("order must be 2 or 3")
    values = np.asarray(volume, dtype=np.float64)
    f = eulerian_finite_strain(values, reference_volume)
    v0 = float(reference_volume)
    k0 = float(bulk_modulus)
    kp = float(bulk_modulus_derivative)
    c0 = float(reference_component)
    cp = float(component_pressure_derivative)
    delta = float(wallace_delta)
    prefactor = np.power(1.0 + 2.0 * f, 2.5)
    a0 = 1.0 - 5.0 * f
    ap = 3.0 * k0 * f
    constant = np.zeros_like(f, dtype=np.float64)
    da0_df = np.full_like(f, -5.0, dtype=np.float64)
    dap_df = np.full_like(f, 3.0 * k0, dtype=np.float64)
    if order == 3:
        a0 = a0 - 14.0 * np.square(f)
        ap = ap + 6.0 * k0 * np.square(f)
        constant = -1.5 * k0 * delta * (3.0 * kp - 16.0) * np.square(f)
        da0_df = da0_df - 28.0 * f
        dap_df = dap_df + 12.0 * k0 * f
    polynomial = a0 * c0 + ap * cp + constant
    dpoly_df = da0_df * c0 + dap_df * cp
    if order == 3:
        dpoly_df = dpoly_df - 3.0 * k0 * delta * (3.0 * kp - 16.0) * f
    dvalue_df = prefactor * (dpoly_df + 5.0 * polynomial / (1.0 + 2.0 * f))

    dc0 = prefactor * a0
    dcp = prefactor * ap
    dv0 = dvalue_df * (1.0 + 2.0 * f) / (3.0 * v0)
    dk0 = prefactor * (3.0 * cp * f)
    dkp = np.zeros_like(f, dtype=np.float64)
    if order == 3:
        dk0 = prefactor * (
            3.0 * cp * f + (6.0 * cp - 1.5 * delta * (3.0 * kp - 16.0)) * np.square(f)
        )
        dkp = prefactor * (-4.5 * k0 * delta * np.square(f))
    dvolume = -dvalue_df * (1.0 + 2.0 * f) / (3.0 * values)
    return np.stack((dc0, dcp, dv0, dk0, dkp, dvolume), axis=-1).astype(np.float64)


def cold_finite_strain_stiffness(
    volume: ArrayLike,
    *,
    reference_volume: float,
    bulk_modulus: float,
    bulk_modulus_derivative: float,
    reference_stiffness: ArrayLike,
    stiffness_pressure_derivative: ArrayLike,
    order: FiniteStrainOrder = 3,
) -> FloatArray:
    r"""Evaluate the cold finite-strain stiffness matrix at one or more volumes.

    The third-order expression follows the thermodynamically consistent
    Eulerian finite-strain form of Stixrude and Lithgow-Bertelloni.  The input
    ``stiffness_pressure_derivative`` contains :math:`\partial C_{IJ}/\partial
    P` at the reference state and is dimensionless when stiffness and pressure
    are expressed in the same units.

    Parameters
    ----------
    volume : array_like
        Positive evaluation volume or array of volumes.
    reference_volume : float
        Reference volume :math:`V_0`.
    bulk_modulus : float
        Reference isothermal bulk modulus :math:`K_0`, in GPa.
    bulk_modulus_derivative : float
        Reference pressure derivative :math:`K'_0`.
    reference_stiffness : array_like
        Reference stiffness matrix :math:`C^0_{IJ}` with shape ``(6, 6)``, in
        GPa.
    stiffness_pressure_derivative : array_like
        Reference pressure derivatives :math:`C'^0_{IJ}` with shape ``(6, 6)``.
    order : {2, 3}, optional
        Finite-strain truncation order.  Order 3 includes the quadratic strain
        term and is the recommended production formulation.

    Notes
    -----
    The inputs are Wallace stress--strain coefficients at the sampled
    hydrostatic states.  The delta tensor entering the quadratic coefficient is
    part of the finite-strain constitutive relation and is not a second
    pre-stress correction of those inputs.

    Returns
    -------
    ndarray
        Stiffness matrix with shape ``volume.shape + (6, 6)``.  A scalar volume
        therefore returns shape ``(6, 6)``.

    Raises
    ------
    ValueError
        If shapes, orders, or scalar parameters are invalid.
    """
    if order not in (2, 3):
        raise ValueError("order must be 2 or 3")
    k0 = float(bulk_modulus)
    kp0 = float(bulk_modulus_derivative)
    if not np.isfinite(k0) or k0 <= 0.0:
        raise ValueError("bulk_modulus must be finite and positive")
    if not np.isfinite(kp0):
        raise ValueError("bulk_modulus_derivative must be finite")

    c0 = validate_stiffness_matrix(reference_stiffness, copy=True)
    cp = np.asarray(stiffness_pressure_derivative, dtype=np.float64)
    if cp.shape != (6, 6):
        raise ValueError("stiffness_pressure_derivative must have shape (6, 6)")
    if np.any(~np.isfinite(cp)):
        raise ValueError("stiffness_pressure_derivative must be finite")
    if not np.allclose(cp, cp.T, rtol=0.0, atol=1.0e-12):
        raise ValueError("stiffness_pressure_derivative must be symmetric")

    strain = eulerian_finite_strain(volume, reference_volume)
    scalar = strain.ndim == 0
    f = np.atleast_1d(strain)[..., np.newaxis, np.newaxis]
    prefactor = np.power(1.0 + 2.0 * f, 2.5)
    linear = 3.0 * k0 * cp - 5.0 * c0
    polynomial = c0 + linear * f
    if order == 3:
        wallace = wallace_hydrostatic_delta_voigt()
        quadratic = 6.0 * k0 * cp - 14.0 * c0 - 1.5 * k0 * wallace * (3.0 * kp0 - 16.0)
        polynomial = polynomial + quadratic * np.square(f)
    result = np.asarray(prefactor * polynomial, dtype=np.float64)
    return result[0] if scalar else result


__all__ = [
    "FiniteStrainOrder",
    "cold_finite_strain_component",
    "cold_finite_strain_component_jacobian",
    "cold_finite_strain_stiffness",
    "eulerian_finite_strain",
    "wallace_hydrostatic_delta_voigt",
]
