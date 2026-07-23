# -*- coding: utf-8 -*-

r"""Statistical thermodynamics of independent harmonic phonons.

The functions in this module operate on phonon frequencies expressed in hertz
and use NumPy vectorization throughout. Frequencies less than or equal to zero
are excluded from the oscillator sums. Energy-like quantities are returned in
``kJ mol^-1``; entropy and heat capacity are returned in
``J mol^-1 K^-1``.
"""

from __future__ import annotations

import numpy as np
from scipy import constants as cs

H = cs.Planck
KB = cs.Boltzmann
NA = cs.Avogadro
_KJMOL = 1.0e-3


def _as_temperature_array(temperature: np.ndarray) -> np.ndarray:
    """Return temperatures as a one-dimensional float array.

    Parameters
    ----------
    temperature : array_like
        Temperature values.

    Returns
    -------
    ndarray
        One-dimensional temperature array.

    Raises
    ------
    ValueError
        If the input is not one-dimensional or contains negative values.
    """
    array = np.asarray(temperature, dtype=np.float64)
    if array.ndim != 1:
        raise ValueError("temperature must be a one-dimensional array")
    if np.any(array < 0.0):
        raise ValueError("temperature values must be non-negative")
    return array


def _as_band_array(band: np.ndarray) -> np.ndarray:
    """Return frequencies as a three-dimensional float array.

    Parameters
    ----------
    band : array_like
        Frequencies with shape ``(nq, nmodes, nvol)``.

    Returns
    -------
    ndarray
        Three-dimensional frequency array.

    Raises
    ------
    ValueError
        If the input does not have three dimensions.
    """
    array = np.asarray(band, dtype=np.float64)
    if array.ndim != 3:
        raise ValueError("band must have shape (nq, nmodes, nvol)")
    return array


def _as_weights_array(weights: np.ndarray, nq: int) -> np.ndarray:
    """Return q-point weights as a one-dimensional float array.

    Parameters
    ----------
    weights : array_like
        Q-point weights.
    nq : int
        Expected number of q-points.

    Returns
    -------
    ndarray
        One-dimensional q-point weight array.

    Raises
    ------
    ValueError
        If the weights have an incompatible shape or a non-positive sum.
    """
    array = np.asarray(weights, dtype=np.float64)
    if array.ndim != 1 or array.shape[0] != nq:
        raise ValueError("weights must have shape (nq,)")
    if float(array.sum()) <= 0.0:
        raise ValueError("the sum of q-point weights must be positive")
    return array


def validate_phonon_inputs(
    temperature: np.ndarray,
    band: np.ndarray,
    weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Validate arrays used by harmonic thermodynamic functions.

    Parameters
    ----------
    temperature : array_like
        Temperature values.
    band : array_like
        Frequencies with shape ``(nq, nmodes, nvol)``.
    weights : array_like
        Q-point weights with shape ``(nq,)``.

    Returns
    -------
    tuple of ndarray
        Normalized temperature, frequency, and weight arrays.

    Raises
    ------
    ValueError
        If any array has an incompatible shape or invalid values.
    """
    temperature_array = _as_temperature_array(temperature)
    band_array = _as_band_array(band)
    weights_array = _as_weights_array(weights, band_array.shape[0])
    return temperature_array, band_array, weights_array


def _weighted_sum(contribution: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Sum oscillator contributions over q-points and modes.

    Parameters
    ----------
    contribution : ndarray
        Array with shape ``(nt, nq, nmodes, nvol)``.
    weights : ndarray
        Q-point weights with shape ``(nq,)``.

    Returns
    -------
    ndarray
        Weighted sum with shape ``(nt, nvol)``.
    """
    return np.sum(contribution * weights[None, :, None, None], axis=(1, 2))


def zero_point_energy(
    temperature: np.ndarray,
    band: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    r"""Calculate the zero-point vibrational energy.

    For independent harmonic oscillators,

    .. math::

        U_{\mathrm{zp}} = \sum_{q\nu} w_q\,\frac{h\nu_{q\nu}}{2}.

    Parameters
    ----------
    temperature : array_like
        Temperature values. Only the number of values is used.
    band : array_like
        Frequencies in hertz with shape ``(nq, nmodes, nvol)``.
    weights : array_like
        Q-point weights with shape ``(nq,)``.

    Returns
    -------
    ndarray
        Zero-point energy in ``kJ mol^-1`` with shape ``(nt, nvol)``.

    Raises
    ------
    ValueError
        If the input arrays are invalid.
    """
    temperature, band, weights = validate_phonon_inputs(temperature, band, weights)
    omega = np.where(band > 0.0, band, 0.0)
    contribution = _KJMOL * NA * H * 0.5 * omega
    contribution = np.broadcast_to(
        contribution[None, :, :, :],
        (temperature.shape[0],) + contribution.shape,
    )
    return _weighted_sum(contribution, weights)


def thermal_energy(
    temperature: np.ndarray,
    band: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    r"""Calculate the thermal vibrational internal energy.

    .. math::

        U_{\mathrm{th}} = \sum_{q\nu} w_q
        \frac{h\nu_{q\nu}}
        {\exp\left(h\nu_{q\nu}/k_{\mathrm B}T\right)-1}.

    Parameters
    ----------
    temperature : array_like
        Temperature values in kelvin.
    band : array_like
        Frequencies in hertz with shape ``(nq, nmodes, nvol)``.
    weights : array_like
        Q-point weights with shape ``(nq,)``.

    Returns
    -------
    ndarray
        Thermal energy in ``kJ mol^-1`` with shape ``(nt, nvol)``.

    Raises
    ------
    ValueError
        If the input arrays are invalid.
    """
    temperature, band, weights = validate_phonon_inputs(temperature, band, weights)
    t = temperature[:, None, None, None]
    omega = band[None, :, :, :]
    valid = (omega > 0.0) & (t != 0.0)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        x = H * omega / (KB * t)
        value = _KJMOL * NA * (H * omega / np.expm1(x))
    contribution = np.where(valid, value, 0.0)
    return _weighted_sum(contribution, weights)


def internal_energy(U0: np.ndarray, Uzp: np.ndarray, Uth: np.ndarray) -> np.ndarray:
    r"""Calculate total internal energy.

    .. math::

        U(T,V) = U_0(V) + U_{\mathrm{zp}}(V) + U_{\mathrm{th}}(T,V).

    Parameters
    ----------
    U0 : array_like
        Static energies with shape ``(nvol,)``.
    Uzp : array_like
        Zero-point energies with shape ``(nvol,)``.
    Uth : array_like
        Thermal energies with shape ``(nt, nvol)``.

    Returns
    -------
    ndarray
        Total internal energy with shape ``(nt, nvol)``.

    Raises
    ------
    ValueError
        If the array shapes are incompatible.
    """
    U0 = np.asarray(U0, dtype=np.float64)
    Uzp = np.asarray(Uzp, dtype=np.float64)
    Uth = np.asarray(Uth, dtype=np.float64)
    if U0.ndim != 1 or Uzp.ndim != 1 or Uth.ndim != 2:
        raise ValueError("U0 and Uzp must be 1D, while Uth must be 2D")
    if U0.shape[0] != Uzp.shape[0] or U0.shape[0] != Uth.shape[1]:
        raise ValueError("energy arrays have incompatible volume dimensions")
    return U0[None, :] + Uzp[None, :] + Uth


def entropy(
    temperature: np.ndarray,
    band: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    r"""Calculate harmonic vibrational entropy.

    With :math:`x=h\nu/(k_{\mathrm B}T)`, the oscillator entropy is

    .. math::

        S = k_{\mathrm B}\left[\frac{x}{e^x-1}
        - \ln\left(1-e^{-x}\right)\right].

    Parameters
    ----------
    temperature : array_like
        Temperature values in kelvin.
    band : array_like
        Frequencies in hertz with shape ``(nq, nmodes, nvol)``.
    weights : array_like
        Q-point weights with shape ``(nq,)``.

    Returns
    -------
    ndarray
        Entropy in ``J mol^-1 K^-1`` with shape ``(nt, nvol)``.

    Raises
    ------
    ValueError
        If the input arrays are invalid.
    """
    temperature, band, weights = validate_phonon_inputs(temperature, band, weights)
    t = temperature[:, None, None, None]
    omega = band[None, :, :, :]
    valid = (omega > 0.0) & (t != 0.0)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        x = H * omega / (KB * t)
        occupation = 1.0 / np.expm1(x)
        logarithm = np.log(1.0 - np.exp(-x))
        value = NA * KB * (occupation * x - logarithm)
    contribution = np.where(valid, value, 0.0)
    return _weighted_sum(contribution, weights)


def vibrational_free_energy(
    temperature: np.ndarray,
    band: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    r"""Calculate harmonic vibrational Helmholtz free energy.

    .. math::

        F_{\mathrm{vib}} = \sum_{q\nu} w_q\left[
        \frac{h\nu_{q\nu}}{2}
        + k_{\mathrm B}T\ln\left(1-e^{-h\nu_{q\nu}/k_{\mathrm B}T}\right)
        \right].

    At zero temperature the logarithmic term vanishes and the expression
    reduces to the zero-point energy.

    Parameters
    ----------
    temperature : array_like
        Temperature values in kelvin.
    band : array_like
        Frequencies in hertz with shape ``(nq, nmodes, nvol)``.
    weights : array_like
        Q-point weights with shape ``(nq,)``.

    Returns
    -------
    ndarray
        Vibrational Helmholtz free energy in ``kJ mol^-1`` with shape
        ``(nt, nvol)``.

    Raises
    ------
    ValueError
        If the input arrays are invalid.
    """
    temperature, band, weights = validate_phonon_inputs(temperature, band, weights)
    t = temperature[:, None, None, None]
    omega = band[None, :, :, :]
    valid = omega > 0.0
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        x = H * omega / (KB * t)
        finite_t = _KJMOL * NA * (0.5 * H * omega + KB * t * np.log(1.0 - np.exp(-x)))
        zero_t = _KJMOL * NA * (0.5 * H * omega)
        value = np.where(t == 0.0, zero_t, finite_t)
    contribution = np.where(valid, value, 0.0)
    return _weighted_sum(contribution, weights)


def free_energy(U0: np.ndarray, Fvib: np.ndarray) -> np.ndarray:
    r"""Calculate total Helmholtz free energy.

    .. math::

        F(T,V) = U_0(V) + F_{\mathrm{vib}}(T,V).

    Parameters
    ----------
    U0 : array_like
        Static energies with shape ``(nvol,)``.
    Fvib : array_like
        Vibrational free energies with shape ``(nt, nvol)``.

    Returns
    -------
    ndarray
        Total Helmholtz free energy with shape ``(nt, nvol)``.

    Raises
    ------
    ValueError
        If the array shapes are incompatible.
    """
    U0 = np.asarray(U0, dtype=np.float64)
    Fvib = np.asarray(Fvib, dtype=np.float64)
    if U0.ndim != 1 or Fvib.ndim != 2:
        raise ValueError("U0 must be 1D and Fvib must be 2D")
    if U0.shape[0] != Fvib.shape[1]:
        raise ValueError("energy arrays have incompatible volume dimensions")
    return U0[None, :] + Fvib


def isochoric_heat_capacity(
    temperature: np.ndarray,
    band: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    r"""Calculate harmonic isochoric heat capacity.

    .. math::

        C_V = \sum_{q\nu} w_q k_{\mathrm B}
        \frac{x^2e^x}{(e^x-1)^2},
        \qquad x=\frac{h\nu_{q\nu}}{k_{\mathrm B}T}.

    Parameters
    ----------
    temperature : array_like
        Temperature values in kelvin.
    band : array_like
        Frequencies in hertz with shape ``(nq, nmodes, nvol)``.
    weights : array_like
        Q-point weights with shape ``(nq,)``.

    Returns
    -------
    ndarray
        Isochoric heat capacity in ``J mol^-1 K^-1`` with shape
        ``(nt, nvol)``.

    Raises
    ------
    ValueError
        If the input arrays are invalid.
    """
    temperature, band, weights = validate_phonon_inputs(temperature, band, weights)
    t = temperature[:, None, None, None]
    omega = band[None, :, :, :]
    valid = (omega > 0.0) & (t != 0.0)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        x = H * omega / (KB * t)
        exponential = np.exp(x)
        n = (x / np.expm1(x)) ** 2
        value = NA * KB * exponential * n
    return _weighted_sum(np.where(valid & (n != 0.0), value, 0.0), weights)
