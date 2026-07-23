"""Tests for QHA thermal-expansion calculation methods."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.units import (
    N,
    convert_energy_per_temperature,
    convert_pressure,
    convert_volume,
    pressure_to_energy,
)
from quantas.modules.qha.models import QHAOptions, QHAResult
from quantas.modules.qha.thermodynamics import (
    THERMAL_EXPANSION_SOURCE_CODES,
    calculate_mixed_thermal_expansion,
    calculate_mode_thermal_expansion,
    refresh_thermal_expansion_dependents,
)


def test_mixed_derivative_recovers_analytic_alpha() -> None:
    """Mixed derivatives must evaluate pressure at fixed volume."""
    temperature = np.array([0.0, 100.0, 200.0, 300.0])
    bulk_modulus = 120.0
    expected_alpha = 2.5e-5
    result = QHAResult(
        temperature=temperature,
        pressure=np.array([0.0]),
        equilibrium_volume=np.array([[10.0], [10.1], [10.2], [10.3]]),
        isothermal_bulk_modulus=np.full((4, 1), bulk_modulus),
        valid_mask=np.ones((4, 1), dtype=bool),
    )
    sampled_volume = np.linspace(9.0, 12.0, 6)
    entropy_slope = float(
        pressure_to_energy(
            bulk_modulus * expected_alpha,
            "Ha",
            "A",
            "GPa",
        )
    )
    sampled_entropy = np.vstack(
        [entropy_slope * sampled_volume + 0.01 * t for t in temperature]
    )
    options = QHAOptions(
        temperature_min=0.0,
        temperature_max=300.0,
        temperature_step=100.0,
        thermal_expansion_method="mixed_derivative",
    )

    calculate_mixed_thermal_expansion(result, sampled_volume, sampled_entropy, options)
    refresh_thermal_expansion_dependents(result, options)

    assert result.thermal_expansion_mixed is not None
    assert result.thermal_expansion is not None
    np.testing.assert_allclose(
        result.thermal_expansion_mixed,
        expected_alpha,
        rtol=1.0e-9,
        atol=2.0e-14,
    )
    np.testing.assert_allclose(result.thermal_expansion, expected_alpha)
    assert result.thermal_expansion_source is not None
    assert np.all(
        result.thermal_expansion_source
        == THERMAL_EXPANSION_SOURCE_CODES["mixed_derivative"]
    )


def test_mixed_derivative_uses_numerical_fallback_locally() -> None:
    """Unavailable mixed derivatives must not invalidate numerical alpha."""
    temperature = np.array([0.0, 100.0, 200.0])
    volume = np.array([[10.0], [10.1], [10.2]])
    result = QHAResult(
        temperature=temperature,
        pressure=np.array([0.0]),
        equilibrium_volume=volume,
        isothermal_bulk_modulus=np.full((3, 1), 100.0),
        valid_mask=np.ones((3, 1), dtype=bool),
    )
    options = QHAOptions(
        temperature_min=0.0,
        temperature_max=200.0,
        temperature_step=100.0,
        thermal_expansion_method="mixed_derivative",
    )

    result.thermal_expansion_mixed = np.full_like(volume, np.nan)
    refresh_thermal_expansion_dependents(result, options)

    assert result.thermal_expansion is not None
    assert np.all(np.isfinite(result.thermal_expansion))
    assert result.thermal_expansion_source is not None
    assert np.all(
        result.thermal_expansion_source
        == THERMAL_EXPANSION_SOURCE_CODES["numerical_fallback"]
    )


def test_mode_gruneisen_thermal_expansion_inverts_definition() -> None:
    """Mode-based alpha must invert gamma = alpha K V / Cv."""
    expected_alpha = 3.0e-5
    gamma = 1.25
    volume = np.array([[108.0]])
    bulk_modulus = np.array([[95.0]])
    volume_m3 = np.asarray(convert_volume(volume, "A", "m"), dtype=float)
    bulk_modulus_pa = np.asarray(
        convert_pressure(bulk_modulus, "GPa", "Pa"), dtype=float
    )
    cv_j_mol_k = expected_alpha * bulk_modulus_pa * volume_m3 * N / gamma
    cv_native = np.asarray(
        convert_energy_per_temperature(
            cv_j_mol_k,
            "J mol^-1 K^-1",
            "Ha cell^-1 K^-1",
        ),
        dtype=float,
    )
    result = QHAResult(
        temperature=np.array([300.0]),
        pressure=np.array([0.0]),
        equilibrium_volume=volume,
        isothermal_bulk_modulus=bulk_modulus,
        isochoric_heat_capacity=cv_native,
        mode_weighted_gruneisen=np.array([[gamma]]),
        valid_mask=np.ones((1, 1), dtype=bool),
    )
    options = QHAOptions(
        temperature_min=300.0,
        temperature_max=300.0,
        thermal_expansion_method="mode_gruneisen",
    )

    calculate_mode_thermal_expansion(result, options)
    result.thermal_expansion_numerical = np.array([[1.0e-5]])
    refresh_thermal_expansion_dependents(result, options)

    assert result.thermal_expansion_mode is not None
    np.testing.assert_allclose(
        result.thermal_expansion_mode,
        expected_alpha,
        rtol=1.0e-12,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(result.thermal_expansion, expected_alpha)
    assert result.thermal_expansion_source is not None
    assert (
        int(result.thermal_expansion_source[0, 0])
        == (THERMAL_EXPANSION_SOURCE_CODES["mode_gruneisen"])
    )


def test_numerical_method_preserves_existing_gradient_definition() -> None:
    """The fallback numerical method must retain the current implementation."""
    temperature = np.array([0.0, 10.0, 20.0])
    volume = np.array([[10.0], [10.1], [10.4]])
    expected = np.gradient(volume[:, 0], temperature) / volume[:, 0]
    result = QHAResult(
        temperature=temperature,
        pressure=np.array([0.0]),
        equilibrium_volume=volume,
        valid_mask=np.ones((3, 1), dtype=bool),
    )
    options = QHAOptions(
        temperature_min=0.0,
        temperature_max=20.0,
        temperature_step=10.0,
        thermal_expansion_method="numerical",
    )

    refresh_thermal_expansion_dependents(result, options)

    assert result.thermal_expansion_numerical is not None
    np.testing.assert_allclose(result.thermal_expansion_numerical[:, 0], expected)
    np.testing.assert_allclose(result.thermal_expansion[:, 0], expected)


def test_mode_gruneisen_method_requires_frequency_scheme() -> None:
    """The mode-based method is not defined for the TD scheme."""
    options = QHAOptions(
        scheme="td",
        thermal_expansion_method="mode_gruneisen",
    )
    with pytest.raises(ValueError, match="requires scheme='freq'"):
        options.validate()
