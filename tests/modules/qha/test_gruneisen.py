"""Tests for QHA Grüneisen-parameter utilities."""

from __future__ import annotations

import numpy as np

from quantas.modules.qha.core.gruneisen import (
    GruneisenStatus,
    gruneisen_from_power_law,
    mode_gruneisen,
    thermal_gruneisen,
)


def test_mode_gruneisen_recovers_constant_power_law() -> None:
    """Mode Grüneisen values match exact power-law frequencies."""
    volumes = np.array([70.0, 75.0, 80.0, 85.0, 90.0])
    gamma_values = np.array([1.2, 1.8, 2.4])
    frequencies = np.vstack(
        [
            gruneisen_from_power_law(volumes, gamma, 500.0 + i * 100.0, 80.0)
            for i, gamma in enumerate(gamma_values)
        ]
    )

    result = mode_gruneisen(frequencies, volumes, degree=1)

    assert result.success
    assert result.status == GruneisenStatus.SUCCESS
    assert result.gamma is not None
    np.testing.assert_allclose(
        result.gamma,
        np.repeat(gamma_values[:, None], volumes.size, axis=1),
        rtol=1.0e-12,
        atol=1.0e-12,
    )


def test_mode_gruneisen_handles_multidimensional_mode_arrays() -> None:
    """The last axis is interpreted as the volume coordinate."""
    volumes = np.array([8.0, 9.0, 10.0, 11.0])
    frequencies = np.empty((2, 3, volumes.size))
    expected = np.empty((2, 3, volumes.size))
    for i in range(2):
        for j in range(3):
            gamma = 0.5 + 0.1 * i + 0.2 * j
            frequencies[i, j] = gruneisen_from_power_law(
                volumes, gamma, 100.0 + 10.0 * j, 10.0
            )
            expected[i, j] = gamma

    result = mode_gruneisen(frequencies, volumes, degree=1)

    assert result.success
    assert result.gamma is not None
    assert result.gamma.shape == frequencies.shape
    np.testing.assert_allclose(result.gamma, expected, rtol=1.0e-12, atol=1.0e-12)


def test_mode_gruneisen_rejects_nonpositive_frequencies_by_default() -> None:
    """Logarithmic derivatives require positive frequencies."""
    volumes = np.array([10.0, 11.0, 12.0])
    frequencies = np.array([[100.0, 95.0, 90.0], [0.0, 20.0, 25.0]])

    result = mode_gruneisen(frequencies, volumes, degree=1)

    assert not result.success
    assert result.status == GruneisenStatus.INVALID_INPUT
    assert "positive" in result.warnings[0]


def test_mode_gruneisen_can_skip_nonpositive_modes() -> None:
    """Invalid modes can be reported while preserving usable modes."""
    volumes = np.array([10.0, 11.0, 12.0])
    frequencies = np.array([[100.0, 95.0, 90.0], [0.0, 20.0, 25.0]])

    result = mode_gruneisen(frequencies, volumes, degree=1, allow_nonpositive=True)

    assert not result.success
    assert result.status == GruneisenStatus.PARTIAL
    assert result.gamma is not None
    assert np.all(np.isfinite(result.gamma[0]))
    assert np.all(np.isnan(result.gamma[1]))
    assert result.failed_modes == [(1,)]


def test_thermal_gruneisen_uses_heat_capacity_weights() -> None:
    """Thermal Grüneisen parameters are weighted along the mode axis."""
    gamma = np.array([[1.0, 2.0], [3.0, 4.0]])
    cv = np.array([[2.0, 1.0], [1.0, 3.0]])

    result = thermal_gruneisen(gamma, cv, axis=0)

    assert result.success
    assert result.gamma is not None
    expected = np.array([(1.0 * 2.0 + 3.0 * 1.0) / 3.0, (2.0 * 1.0 + 4.0 * 3.0) / 4.0])
    np.testing.assert_allclose(result.gamma, expected)


def test_thermal_gruneisen_rejects_zero_weight_sum() -> None:
    """The weighted average is invalid if all weights vanish."""
    gamma = np.array([[1.0, 2.0], [3.0, 4.0]])
    cv = np.zeros_like(gamma)

    result = thermal_gruneisen(gamma, cv, axis=0)

    assert not result.success
    assert result.status == GruneisenStatus.INVALID_INPUT
