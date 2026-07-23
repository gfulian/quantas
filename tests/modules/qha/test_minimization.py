from __future__ import annotations

import numpy as np

from quantas.core.physics.units import pressure_to_energy
from quantas.core.physics.eos import FittedEnergyEOS
from quantas.core.math.fitting import FitStatus
from quantas.modules.qha.core.minimization import (
    MinimumStatus,
    VolumeRangeStatus,
    classify_volume,
    fine_grid,
    minimize_eos,
    minimize_polynomial,
    numerical_bulk_modulus,
    polynomial_fit,
    polynomial_stationary_points,
    pressure_shifted_free_energy,
)


def test_fine_grid_returns_centered_relative_factors():
    grid, index = fine_grid(npoints=5, separation=0.05)

    np.testing.assert_allclose(grid, np.array([-0.001, -0.0005, 0.0, 0.0005, 0.001]))
    assert index == 2


def test_fine_grid_rejects_even_number_of_points():
    with np.testing.assert_raises(ValueError):
        fine_grid(npoints=4)


def test_pressure_shifted_free_energy_adds_pv_term():
    volume = np.array([10.0, 11.0, 12.0])
    free_energy = np.array([1.0, 2.0, 3.0])

    shifted = pressure_shifted_free_energy(free_energy, volume, 0.2)

    np.testing.assert_allclose(shifted, np.array([3.0, 4.2, 5.4]))


def test_classify_volume_marks_inside_boundary_and_outside_points():
    volume = np.linspace(10.0, 20.0, 11)

    assert classify_volume(15.0, volume) is VolumeRangeStatus.INSIDE
    assert classify_volume(10.2, volume) is VolumeRangeStatus.NEAR_BOUNDARY
    assert classify_volume(21.0, volume) is VolumeRangeStatus.OUTSIDE


def test_polynomial_fit_returns_coefficients_and_diagnostics():
    volume = np.linspace(8.0, 12.0, 9)
    free_energy = 3.0 + 0.5 * (volume - 10.0) ** 2

    result = polynomial_fit(volume, free_energy, degree=2)

    assert result.success
    assert result.status is FitStatus.SUCCESS
    assert result.metadata["model"] == "polynomial"
    assert result.parameters is not None
    np.testing.assert_allclose(
        np.polynomial.polynomial.polyval(volume, result.parameters),
        free_energy,
        atol=1.0e-12,
    )
    assert result.rmse is not None and result.rmse < 1.0e-12


def test_polynomial_stationary_points_returns_minima_only():
    # f(V) = (V - 10)^2, coefficients in ascending order.
    parameters = np.array([100.0, -20.0, 1.0])

    minima = polynomial_stationary_points(parameters)

    np.testing.assert_allclose(minima, np.array([10.0]))


def test_minimize_polynomial_finds_pressure_shifted_minimum():
    volume = np.linspace(8.0, 12.0, 9)
    free_energy = (volume - 10.0) ** 2

    result = minimize_polynomial(
        volume, free_energy, pressure_energy_density=0.5, degree=2
    )

    assert result.success
    assert result.status is MinimumStatus.SUCCESS
    assert result.method == "polynomial"
    assert result.range_status is VolumeRangeStatus.INSIDE
    assert result.volume is not None
    np.testing.assert_allclose(result.volume, 9.75, atol=1.0e-10)
    assert result.fit is not None


def test_minimize_polynomial_warns_for_extrapolated_minimum():
    volume = np.linspace(8.0, 12.0, 9)
    free_energy = (volume - 10.0) ** 2

    result = minimize_polynomial(
        volume, free_energy, pressure_energy_density=8.0, degree=2
    )

    assert result.success
    assert result.range_status is VolumeRangeStatus.OUTSIDE
    assert result.warnings


def test_minimize_eos_finds_reference_volume_and_bulk_modulus():
    from quantas.core.physics.eos import EnergyEOS

    eos = EnergyEOS()
    volume = np.linspace(68.0, 76.0, 11)
    pars = np.array([-100.0, 0.55, 4.2, 72.0])
    energy = eos.birchmurnaghan(volume, *pars)

    result = minimize_eos(volume, energy, eos="birchmurnaghan")

    assert result.success, result.message
    assert result.method == "eos"
    assert result.range_status is VolumeRangeStatus.INSIDE
    assert result.volume is not None
    assert result.bulk_modulus is not None
    assert result.bulk_modulus_derivative is not None
    np.testing.assert_allclose(result.volume, pars[3], atol=1.0e-8)
    np.testing.assert_allclose(result.bulk_modulus, pars[1], atol=1.0e-8)
    np.testing.assert_allclose(result.bulk_modulus_derivative, pars[2], atol=1.0e-8)
    assert result.fit is not None and result.fit.success


def test_minimize_eos_returns_failed_result_when_fit_cannot_start():
    volume = np.array([1.0, 2.0, 3.0, 4.0])
    energy = np.array([3.0, 2.0, 2.0, 3.0])

    result = minimize_eos(volume, energy, eos="BM")

    assert not result.success
    assert result.status is MinimumStatus.FAILED
    assert result.fit is not None


def test_minimize_eos_evaluates_nonzero_pressure_from_unshifted_fit():
    from quantas.core.physics.eos import EnergyEOS

    eos = EnergyEOS()
    volume = np.linspace(68.0, 76.0, 11)
    k0_gpa = 160.0
    k0_native = float(pressure_to_energy(k0_gpa, "eV", "A", "GPa"))
    pars = np.array([-100.0, k0_native, 4.2, 72.0])
    energy = eos.birchmurnaghan(volume, *pars)
    pressure = 20.0
    pressure_energy_density = float(pressure_to_energy(pressure, "eV", "A", "GPa"))

    result = minimize_eos(
        volume,
        energy,
        pressure_energy_density=pressure_energy_density,
        eos="BM",
        energy_unit="eV",
        volume_unit="A",
        pressure_unit="GPa",
    )
    expected = FittedEnergyEOS(
        "BM",
        [-100.0, k0_gpa, 4.2, 72.0],
        sampled_volumes=volume,
    ).state_at_pressure(pressure)

    assert result.success
    np.testing.assert_allclose(result.volume, expected.volume, rtol=1.0e-9)
    np.testing.assert_allclose(
        result.bulk_modulus,
        expected.bulk_modulus,
        rtol=1.0e-9,
    )
    np.testing.assert_allclose(
        result.bulk_modulus_derivative,
        expected.bulk_modulus_derivative,
        rtol=1.0e-9,
    )


def test_numerical_bulk_modulus_returns_expected_quadratic_value():
    volume = np.linspace(8.0, 12.0, 9)
    free_energy = 3.0 + 0.5 * (volume - 10.0) ** 2
    index = int(np.argmin(np.abs(volume - 10.0)))

    kt, kp, fit = numerical_bulk_modulus(volume, free_energy, index=index, degree=2)

    assert fit.success
    np.testing.assert_allclose(kt, 10.0, atol=1.0e-10)
    assert np.isfinite(kp)
