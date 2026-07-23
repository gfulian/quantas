from __future__ import annotations

import math

import numpy as np
from scipy import constants as cs

from quantas.core.physics.thermodynamics import (
    entropy,
    free_energy,
    internal_energy,
    isochoric_heat_capacity,
    thermal_energy,
    vibrational_free_energy,
    zero_point_energy,
)

H = cs.Planck
KB = cs.Boltzmann
NA = cs.Avogadro


def _reference_zero_point(temperature, band, weights):
    result = np.zeros((len(temperature), band.shape[2]), dtype=np.float64)
    for i, _ in enumerate(temperature):
        for j in range(band.shape[0]):
            for k in range(band.shape[1]):
                for n in range(band.shape[2]):
                    omega = band[j, k, n]
                    if omega > 0.0:
                        result[i, n] += 1.0e-3 * NA * H * (0.5 * omega) * weights[j]
    return result


def _reference_thermal(temperature, band, weights):
    result = np.zeros((len(temperature), band.shape[2]), dtype=np.float64)
    for i, temp in enumerate(temperature):
        for j in range(band.shape[0]):
            for k in range(band.shape[1]):
                for n in range(band.shape[2]):
                    omega = band[j, k, n]
                    if omega > 0.0 and temp != 0.0:
                        x = H * omega / (KB * temp)
                        result[i, n] += (
                            1.0e-3 * NA * (H * omega / math.expm1(x)) * weights[j]
                        )
    return result


def _reference_entropy(temperature, band, weights):
    result = np.zeros((len(temperature), band.shape[2]), dtype=np.float64)
    for i, temp in enumerate(temperature):
        for j in range(band.shape[0]):
            for k in range(band.shape[1]):
                for n in range(band.shape[2]):
                    omega = band[j, k, n]
                    if omega > 0.0 and temp != 0.0:
                        x = H * omega / (KB * temp)
                        occ = 1.0 / math.expm1(x)
                        log_term = math.log(1.0 - math.exp(-x))
                        result[i, n] += NA * KB * (occ * x - log_term) * weights[j]
    return result


def _reference_fvib(temperature, band, weights):
    result = np.zeros((len(temperature), band.shape[2]), dtype=np.float64)
    for i, temp in enumerate(temperature):
        for j in range(band.shape[0]):
            for k in range(band.shape[1]):
                for n in range(band.shape[2]):
                    omega = band[j, k, n]
                    if omega <= 0.0:
                        continue
                    hw = H * omega
                    if temp == 0.0:
                        value = 1.0e-3 * NA * (0.5 * hw)
                    else:
                        x = H * omega / (KB * temp)
                        kt = KB * temp
                        value = (
                            1.0e-3 * NA * (0.5 * hw + kt * math.log(1.0 - math.exp(-x)))
                        )
                    result[i, n] += value * weights[j]
    return result


def _reference_cv(temperature, band, weights):
    result = np.zeros((len(temperature), band.shape[2]), dtype=np.float64)
    for i, temp in enumerate(temperature):
        for j in range(band.shape[0]):
            for k in range(band.shape[1]):
                for n in range(band.shape[2]):
                    omega = band[j, k, n]
                    if omega > 0.0 and temp != 0.0:
                        x = H * omega / (KB * temp)
                        exp_x = math.exp(x)
                        nvalue = (x / math.expm1(x)) ** 2
                        if nvalue != 0.0:
                            result[i, n] += NA * KB * exp_x * nvalue * weights[j]
    return result


def test_phonon_thermodynamics_match_explicit_oscillator_sums():
    temperature = np.array([0.0, 100.0, 300.0], dtype=np.float64)
    band = np.array(
        [
            [[0.0, -1.0], [1.0e12, 1.2e12], [2.0e12, 2.4e12]],
            [[1.5e12, 1.8e12], [3.0e12, 3.6e12], [4.5e12, 5.4e12]],
        ],
        dtype=np.float64,
    )
    weights = np.array([0.25, 0.75], dtype=np.float64)

    np.testing.assert_allclose(
        zero_point_energy(temperature, band, weights),
        _reference_zero_point(temperature, band, weights),
        rtol=1.0e-13,
        atol=1.0e-13,
    )
    np.testing.assert_allclose(
        thermal_energy(temperature, band, weights),
        _reference_thermal(temperature, band, weights),
        rtol=1.0e-13,
        atol=1.0e-13,
    )
    np.testing.assert_allclose(
        entropy(temperature, band, weights),
        _reference_entropy(temperature, band, weights),
        rtol=1.0e-13,
        atol=1.0e-13,
    )
    np.testing.assert_allclose(
        vibrational_free_energy(temperature, band, weights),
        _reference_fvib(temperature, band, weights),
        rtol=1.0e-13,
        atol=1.0e-13,
    )
    np.testing.assert_allclose(
        isochoric_heat_capacity(temperature, band, weights),
        _reference_cv(temperature, band, weights),
        rtol=1.0e-13,
        atol=1.0e-13,
    )


def test_total_energy_helpers_match_expected_broadcasting():
    U0 = np.array([-10.0, -9.5], dtype=np.float64)
    Uzp = np.array([0.1, 0.2], dtype=np.float64)
    Uth = np.array([[0.0, 0.0], [0.3, 0.4]], dtype=np.float64)
    Fvib = np.array([[0.1, 0.2], [-0.2, -0.1]], dtype=np.float64)

    np.testing.assert_allclose(
        internal_energy(U0, Uzp, Uth),
        np.array([[-9.9, -9.3], [-9.6, -8.9]], dtype=np.float64),
    )
    np.testing.assert_allclose(
        free_energy(U0, Fvib),
        np.array([[-9.9, -9.3], [-10.2, -9.6]], dtype=np.float64),
    )


def test_shape_validation_raises_clear_errors():
    temperature = np.array([300.0], dtype=np.float64)
    band = np.ones((1, 3, 1), dtype=np.float64)

    with np.testing.assert_raises(ValueError):
        zero_point_energy(temperature[:, None], band, np.ones(1))

    with np.testing.assert_raises(ValueError):
        zero_point_energy(temperature, band[0], np.ones(1))

    with np.testing.assert_raises(ValueError):
        zero_point_energy(temperature, band, np.ones(2))
