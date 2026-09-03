"""Scientific characterization of Kieffer acoustic thermodynamics."""

from __future__ import annotations

import numpy as np
from scipy import constants as cs

from quantas.core.physics.thermodynamics import (
    kieffer_entropy,
    kieffer_isochoric_heat_capacity,
    kieffer_thermal_energy,
    kieffer_thermal_free_energy,
    kieffer_vibrational_free_energy,
    kieffer_zero_point_energy,
    validate_kieffer_inputs,
)


CUTOFFS = np.array([3.0e12, 4.0e12, 5.0e12], dtype=np.float64)


def test_historical_helmholtz_and_cv_reference() -> None:
    """The corrected core preserves legacy quantities that match the paper."""
    temperature = np.array([298.15], dtype=np.float64)

    # Frozen outputs from the historical Kieffer class.  Its Helmholtz method
    # returns J mol^-1 despite not documenting the unit.
    legacy_helmholtz_j_mol = -6521.4320263502705
    legacy_cv_j_mol_k = 24.237038581996995

    np.testing.assert_allclose(
        kieffer_thermal_free_energy(temperature, CUTOFFS)[0, 0],
        1.0e-3 * legacy_helmholtz_j_mol,
        rtol=2.0e-12,
    )
    np.testing.assert_allclose(
        kieffer_isochoric_heat_capacity(temperature, CUTOFFS)[0, 0],
        legacy_cv_j_mol_k,
        rtol=2.0e-12,
    )


def test_historical_entropy_defect_is_corrected() -> None:
    """Freeze and expose the squared-denominator error in legacy entropy."""
    temperature = np.array([298.15], dtype=np.float64)
    legacy_entropy_j_mol_k = 50.62987574670751
    corrected = kieffer_entropy(temperature, CUTOFFS)[0, 0]

    assert not np.isclose(corrected, legacy_entropy_j_mol_k, rtol=1.0e-3)
    np.testing.assert_allclose(corrected, 40.437182509075214, rtol=2.0e-11)


def test_thermodynamic_identities_hold() -> None:
    """Entropy and heat capacity agree with derivatives of free energy."""
    temperature = 298.15
    step = 1.0e-2
    grid = np.array([temperature - step, temperature, temperature + step])
    free = kieffer_thermal_free_energy(grid, CUTOFFS)[:, 0]
    entropy_from_free = -1.0e3 * (free[2] - free[0]) / (2.0 * step)

    entropy_grid = kieffer_entropy(grid, CUTOFFS)[:, 0]
    cv_from_entropy = temperature * (entropy_grid[2] - entropy_grid[0]) / (2.0 * step)

    np.testing.assert_allclose(
        kieffer_entropy(np.array([temperature]), CUTOFFS)[0, 0],
        entropy_from_free,
        rtol=2.0e-8,
    )
    np.testing.assert_allclose(
        kieffer_isochoric_heat_capacity(np.array([temperature]), CUTOFFS)[0, 0],
        cv_from_entropy,
        rtol=2.0e-8,
    )


def test_energy_identity_holds() -> None:
    """Thermal internal energy satisfies ``U = F + T S``."""
    temperature = np.array([50.0, 298.15, 1000.0], dtype=np.float64)
    thermal = kieffer_thermal_energy(temperature, CUTOFFS)
    free = kieffer_thermal_free_energy(temperature, CUTOFFS)
    entropy = kieffer_entropy(temperature, CUTOFFS)
    np.testing.assert_allclose(
        thermal,
        free + 1.0e-3 * temperature[:, None] * entropy,
        rtol=2.0e-12,
        atol=1.0e-12,
    )


def test_zero_temperature_and_zero_point_limits() -> None:
    """Thermal functions vanish at 0 K while total Helmholtz retains ZPE."""
    temperature = np.array([0.0], dtype=np.float64)
    mean_sine = 24.0 * (np.pi - 2.0) / np.pi**3
    expected_zpe = 1.0e-3 * cs.Avogadro * 0.5 * cs.Planck * mean_sine * CUTOFFS.sum()

    np.testing.assert_allclose(kieffer_thermal_energy(temperature, CUTOFFS), 0.0)
    np.testing.assert_allclose(kieffer_entropy(temperature, CUTOFFS), 0.0)
    np.testing.assert_allclose(
        kieffer_isochoric_heat_capacity(temperature, CUTOFFS), 0.0
    )
    np.testing.assert_allclose(kieffer_thermal_free_energy(temperature, CUTOFFS), 0.0)
    np.testing.assert_allclose(
        kieffer_zero_point_energy(temperature, CUTOFFS), expected_zpe
    )
    np.testing.assert_allclose(
        kieffer_vibrational_free_energy(temperature, CUTOFFS), expected_zpe
    )


def test_high_temperature_heat_capacity_tends_to_three_r() -> None:
    """The three acoustic branches approach the Dulong-Petit branch limit."""
    result = kieffer_isochoric_heat_capacity(np.array([1.0e7]), CUTOFFS)[0, 0]
    np.testing.assert_allclose(result, 3.0 * cs.gas_constant, rtol=1.0e-8)


def test_volume_series_preserves_shape_and_float64() -> None:
    """Three cutoff rows support an arbitrary volume axis."""
    cutoffs = np.column_stack((CUTOFFS, 1.1 * CUTOFFS))
    result = kieffer_entropy(np.array([100.0, 300.0]), cutoffs)
    assert result.shape == (2, 2)
    assert result.dtype == np.dtype("float64")


def test_invalid_inputs_raise_clear_errors() -> None:
    """Reject invalid dimensions, temperatures, and cutoff frequencies."""
    for temperature, cutoffs in (
        (np.array([[300.0]]), CUTOFFS),
        (np.array([-1.0]), CUTOFFS),
        (np.array([np.nan]), CUTOFFS),
        (np.array([300.0]), np.array([1.0e12, 2.0e12])),
        (np.array([300.0]), np.array([1.0e12, 2.0e12, 0.0])),
        (np.array([300.0]), np.array([1.0e12, 2.0e12, np.inf])),
    ):
        with np.testing.assert_raises(ValueError):
            validate_kieffer_inputs(temperature, cutoffs)
