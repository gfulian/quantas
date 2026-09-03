# -*- coding: utf-8 -*-

r"""Acoustic thermodynamics for the sine-wave model of Kieffer.

The model represents each of the three acoustic branches with the dispersion

.. math::

    \nu(k) = \nu_{\max}\sin\left(\frac{\pi k}{2K_{\max}}\right).

This module evaluates the resulting density-of-states integrals from ordinary
cutoff frequencies in hertz.  It contains no elastic, Christoffel, HA, QHA, or
frontend logic.  Energy-like quantities are returned in ``kJ mol^-1`` per mole
of thermodynamic normalization cells; entropy and heat capacity are returned
in ``J mol^-1 K^-1``.
"""

from __future__ import annotations

from collections.abc import Callable
import math

import numpy as np
from scipy import constants as cs
from scipy.integrate import quad

H = cs.Planck
KB = cs.Boltzmann
NA = cs.Avogadro
R = cs.gas_constant

_KJMOL = 1.0e-3
_THETA_MAX = 0.5 * math.pi
_DOS_FACTOR = 3.0 * (2.0 / math.pi) ** 3
_MEAN_SINE = 24.0 * (math.pi - 2.0) / math.pi**3


def _as_temperature_array(temperature: np.ndarray) -> np.ndarray:
    """Return a finite, non-negative one-dimensional temperature array."""
    array = np.asarray(temperature, dtype=np.float64)
    if array.ndim != 1:
        raise ValueError("temperature must be a one-dimensional array")
    if not np.all(np.isfinite(array)) or np.any(array < 0.0):
        raise ValueError("temperature values must be finite and non-negative")
    return array


def _as_cutoff_array(cutoff_frequency: np.ndarray) -> np.ndarray:
    """Return cutoff frequencies with shape ``(3, nvol)``.

    A one-dimensional three-element input represents one volume.  Exactly
    three positive finite branches are required because the Kieffer acoustic
    model represents two quasi-shear and one quasi-longitudinal branch.
    """
    array = np.asarray(cutoff_frequency, dtype=np.float64)
    if array.ndim == 1:
        array = array[:, None]
    if array.ndim != 2 or array.shape[0] != 3:
        raise ValueError("cutoff_frequency must have shape (3,) or (3, nvol)")
    if not np.all(np.isfinite(array)) or np.any(array <= 0.0):
        raise ValueError("cutoff frequencies must be finite and positive")
    return array


def validate_kieffer_inputs(
    temperature: np.ndarray,
    cutoff_frequency: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Validate arrays used by the Kieffer thermodynamic functions.

    Parameters
    ----------
    temperature : array_like
        Temperatures in kelvin with shape ``(nt,)``.
    cutoff_frequency : array_like
        Three ordinary cutoff frequencies in hertz, with shape ``(3,)`` for
        one volume or ``(3, nvol)`` for a volume series.

    Returns
    -------
    tuple of ndarray
        Normalized temperature and cutoff-frequency arrays.

    Raises
    ------
    ValueError
        If a shape, temperature, or cutoff frequency is invalid.
    """
    return _as_temperature_array(temperature), _as_cutoff_array(cutoff_frequency)


def _occupation_energy(x: float) -> float:
    """Return ``x / (exp(x) - 1)`` with stable limiting behaviour."""
    if x == 0.0:
        return 1.0
    if x > 700.0:
        return 0.0
    return x / math.expm1(x)


def _log_bose_factor(x: float) -> float:
    """Return ``log(1 - exp(-x))`` accurately for positive ``x``."""
    if x == 0.0:
        return -math.inf
    return math.log(-math.expm1(-x))


def _heat_capacity_factor(x: float) -> float:
    """Return ``x**2 exp(x) / (exp(x) - 1)**2`` stably."""
    if x == 0.0:
        return 1.0
    if x > 700.0:
        return 0.0
    half = 0.5 * x
    if half > 350.0:
        return 0.0
    return (half / math.sinh(half)) ** 2


def _integrate_branch(
    characteristic_temperature: float,
    temperature: float,
    kernel: Callable[[float], float],
) -> float:
    r"""Integrate one branch after the substitution ``k -> theta``.

    The transformation

    .. math::

        \theta=\arcsin(\nu/\nu_{\max})

    removes the integrable square-root singularity present in the published
    frequency-domain expression.
    """
    xi = characteristic_temperature / temperature

    def integrand(theta: float) -> float:
        if theta == 0.0:
            return 0.0
        return theta * theta * kernel(xi * math.sin(theta))

    value, _error = quad(
        integrand,
        0.0,
        _THETA_MAX,
        epsabs=1.0e-12,
        epsrel=1.0e-10,
        limit=200,
    )
    return float(value)


def _evaluate_finite_temperature(
    temperature: np.ndarray,
    cutoff_frequency: np.ndarray,
    kernel: Callable[[float], float],
    scale: Callable[[float], float],
) -> np.ndarray:
    """Evaluate one Kieffer integral for all temperatures and volumes."""
    temperatures, cutoffs = validate_kieffer_inputs(temperature, cutoff_frequency)
    characteristic = H * cutoffs / KB
    result = np.zeros((temperatures.size, cutoffs.shape[1]), dtype=np.float64)
    for it, current_temperature in enumerate(temperatures):
        if current_temperature == 0.0:
            continue
        for iv in range(cutoffs.shape[1]):
            integral = sum(
                _integrate_branch(value, current_temperature, kernel)
                for value in characteristic[:, iv]
            )
            result[it, iv] = scale(current_temperature) * _DOS_FACTOR * integral
    return result


def zero_point_energy(
    temperature: np.ndarray,
    cutoff_frequency: np.ndarray,
) -> np.ndarray:
    r"""Calculate the Kieffer acoustic zero-point energy.

    For one sine-wave branch,

    .. math::

        \langle\nu\rangle = \nu_{\max}
        \frac{24(\pi-2)}{\pi^3}.

    Parameters
    ----------
    temperature : array_like
        Temperatures in kelvin.  Only the number of values is used.
    cutoff_frequency : array_like
        Ordinary cutoff frequencies in hertz with shape ``(3,)`` or
        ``(3, nvol)``.

    Returns
    -------
    ndarray
        Zero-point energy in ``kJ mol^-1`` with shape ``(nt, nvol)``.

    Raises
    ------
    ValueError
        If the input arrays are invalid.
    """
    temperatures, cutoffs = validate_kieffer_inputs(temperature, cutoff_frequency)
    per_volume = _KJMOL * NA * 0.5 * H * _MEAN_SINE * np.sum(cutoffs, axis=0)
    return np.broadcast_to(
        per_volume[None, :], (temperatures.size, cutoffs.shape[1])
    ).copy()


def thermal_energy(
    temperature: np.ndarray,
    cutoff_frequency: np.ndarray,
) -> np.ndarray:
    r"""Calculate the Kieffer acoustic thermal internal energy.

    Parameters
    ----------
    temperature : array_like
        Temperatures in kelvin with shape ``(nt,)``.
    cutoff_frequency : array_like
        Ordinary cutoff frequencies in hertz with shape ``(3,)`` or
        ``(3, nvol)``.

    Returns
    -------
    ndarray
        ``U(T)-U(0)`` in ``kJ mol^-1`` with shape ``(nt, nvol)``.

    Raises
    ------
    ValueError
        If the input arrays are invalid.
    """
    return _evaluate_finite_temperature(
        temperature,
        cutoff_frequency,
        _occupation_energy,
        lambda value: _KJMOL * R * value,
    )


def entropy(
    temperature: np.ndarray,
    cutoff_frequency: np.ndarray,
) -> np.ndarray:
    r"""Calculate the Kieffer acoustic entropy.

    The implemented kernel is

    .. math::

        \frac{x}{e^x-1}-\ln(1-e^{-x}),

    as printed in Kieffer's equation (27).  In particular, the occupation
    denominator is not squared as it was in the historical Quantas routine.

    Parameters
    ----------
    temperature : array_like
        Temperatures in kelvin with shape ``(nt,)``.
    cutoff_frequency : array_like
        Ordinary cutoff frequencies in hertz with shape ``(3,)`` or
        ``(3, nvol)``.

    Returns
    -------
    ndarray
        Entropy in ``J mol^-1 K^-1`` with shape ``(nt, nvol)``.

    Raises
    ------
    ValueError
        If the input arrays are invalid.
    """
    return _evaluate_finite_temperature(
        temperature,
        cutoff_frequency,
        lambda x: _occupation_energy(x) - _log_bose_factor(x),
        lambda _value: R,
    )


def isochoric_heat_capacity(
    temperature: np.ndarray,
    cutoff_frequency: np.ndarray,
) -> np.ndarray:
    r"""Calculate Kieffer acoustic constant-volume heat capacity.

    Parameters
    ----------
    temperature : array_like
        Temperatures in kelvin with shape ``(nt,)``.
    cutoff_frequency : array_like
        Ordinary cutoff frequencies in hertz with shape ``(3,)`` or
        ``(3, nvol)``.

    Returns
    -------
    ndarray
        ``C_V`` in ``J mol^-1 K^-1`` with shape ``(nt, nvol)``.

    Raises
    ------
    ValueError
        If the input arrays are invalid.
    """
    return _evaluate_finite_temperature(
        temperature,
        cutoff_frequency,
        _heat_capacity_factor,
        lambda _value: R,
    )


def thermal_free_energy(
    temperature: np.ndarray,
    cutoff_frequency: np.ndarray,
) -> np.ndarray:
    r"""Calculate acoustic Helmholtz free energy above its 0 K value.

    This is the quantity evaluated by the historical ``Kieffer.helmholtz``
    method.

    Parameters
    ----------
    temperature : array_like
        Temperatures in kelvin with shape ``(nt,)``.
    cutoff_frequency : array_like
        Ordinary cutoff frequencies in hertz with shape ``(3,)`` or
        ``(3, nvol)``.

    Returns
    -------
    ndarray
        ``F(T)-F(0)`` in ``kJ mol^-1`` with shape ``(nt, nvol)``.

    Raises
    ------
    ValueError
        If the input arrays are invalid.
    """
    return _evaluate_finite_temperature(
        temperature,
        cutoff_frequency,
        _log_bose_factor,
        lambda value: _KJMOL * R * value,
    )


def vibrational_free_energy(
    temperature: np.ndarray,
    cutoff_frequency: np.ndarray,
) -> np.ndarray:
    r"""Calculate total acoustic Helmholtz energy, including zero point.

    Parameters
    ----------
    temperature : array_like
        Temperatures in kelvin with shape ``(nt,)``.
    cutoff_frequency : array_like
        Ordinary cutoff frequencies in hertz with shape ``(3,)`` or
        ``(3, nvol)``.

    Returns
    -------
    ndarray
        ``U_zp + F(T)-F(0)`` in ``kJ mol^-1`` with shape ``(nt, nvol)``.

    Raises
    ------
    ValueError
        If the input arrays are invalid.
    """
    return zero_point_energy(temperature, cutoff_frequency) + thermal_free_energy(
        temperature, cutoff_frequency
    )


__all__ = [
    "entropy",
    "isochoric_heat_capacity",
    "thermal_energy",
    "thermal_free_energy",
    "validate_kieffer_inputs",
    "vibrational_free_energy",
    "zero_point_energy",
]
