"""Tests for anisotropic isothermal-to-adiabatic conversion."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from quantas.core.physics.elasticity import adiabatic_stiffness_field


def _isotropic_stiffness(bulk: float, shear: float) -> np.ndarray:
    c11 = bulk + 4.0 * shear / 3.0
    c12 = bulk - 2.0 * shear / 3.0
    matrix = np.zeros((6, 6), dtype=np.float64)
    matrix[:3, :3] = c12
    np.fill_diagonal(matrix[:3, :3], c11)
    np.fill_diagonal(matrix[3:, 3:], shear)
    return matrix


def test_isotropic_limit_changes_only_bulk_response() -> None:
    """The tensor identity reproduces the scalar K_S-K_T relation."""
    bulk = 120.0
    shear = 70.0
    temperature = 1000.0
    volume = 1.0e-28
    cv = 4.0e-22
    alpha_v = 3.0e-5
    alpha = np.eye(3, dtype=np.float64) * (alpha_v / 3.0)
    isothermal = _isotropic_stiffness(bulk, shear)

    result = adiabatic_stiffness_field(
        isothermal,
        temperature,
        volume,
        cv,
        alpha,
    )

    expected_delta_k_pa = temperature * volume / cv * (bulk * 1.0e9 * alpha_v) ** 2
    expected_delta_k = expected_delta_k_pa / 1.0e9
    adiabatic = result.stiffness
    calculated_bulk = (adiabatic[0, 0] + 2.0 * adiabatic[0, 1]) / 3.0
    np.testing.assert_allclose(calculated_bulk, bulk + expected_delta_k, rtol=1e-13)
    np.testing.assert_allclose(adiabatic[3:, 3:], isothermal[3:, 3:])
    assert result.valid_mask
    assert not result.invalid_mask


def test_zero_kelvin_limit_does_not_require_cv_or_expansion() -> None:
    """At exactly zero kelvin the adiabatic and isothermal tensors coincide."""
    stiffness = _isotropic_stiffness(100.0, 60.0)
    result = adiabatic_stiffness_field(
        stiffness,
        0.0,
        1.0e-28,
        None,
        None,
    )
    np.testing.assert_array_equal(result.stiffness, stiffness)
    np.testing.assert_array_equal(result.correction, np.zeros((6, 6)))
    assert result.valid_mask


def test_zero_kelvin_uncertainty_equals_isothermal_uncertainty() -> None:
    """The exact zero-temperature limit must not double-count stiffness error."""
    stiffness = _isotropic_stiffness(100.0, 60.0)
    sigma = np.full((6, 6), 0.1, dtype=np.float64)
    result = adiabatic_stiffness_field(
        stiffness,
        0.0,
        1.0e-28,
        None,
        None,
        sigma_stiffness_isothermal=sigma,
    )
    assert result.sigma_stiffness is not None
    np.testing.assert_array_equal(result.sigma_stiffness, sigma)


def test_nonzero_temperature_missing_inputs_is_explicitly_invalid() -> None:
    """Missing thermodynamic fields never cause a silent C_S=C_T fallback."""
    stiffness = _isotropic_stiffness(100.0, 60.0)
    result = adiabatic_stiffness_field(
        stiffness,
        300.0,
        1.0e-28,
        None,
        None,
    )
    assert result.invalid_mask
    assert np.all(np.isnan(result.stiffness))


def test_first_order_uncertainty_is_finite_and_nonnegative() -> None:
    """Available independent input errors produce a finite sigma field."""
    stiffness = _isotropic_stiffness(120.0, 70.0)
    alpha = np.eye(3, dtype=np.float64) * 1.0e-5
    result = adiabatic_stiffness_field(
        stiffness,
        800.0,
        1.0e-28,
        4.0e-22,
        alpha,
        sigma_stiffness_isothermal=np.full((6, 6), 0.1),
        sigma_volume_m3=1.0e-31,
        sigma_heat_capacity_j_per_k=1.0e-24,
        sigma_thermal_expansion_tensor=np.full((3, 3), 1.0e-7),
    )
    assert result.sigma_stiffness is not None
    assert np.all(np.isfinite(result.sigma_stiffness))
    assert np.all(result.sigma_stiffness >= 0.0)


@pytest.mark.parametrize("parameter", ("volume", "heat_capacity", "alpha11"))
def test_scalar_uncertainty_jacobians_match_central_differences(
    parameter: str,
) -> None:
    """Analytical volume, heat-capacity, and expansion Jacobians are correct."""
    stiffness = _isotropic_stiffness(120.0, 70.0)
    temperature = 800.0
    volume = 1.0e-28
    heat_capacity = 4.0e-22
    alpha = np.eye(3, dtype=np.float64) * 1.0e-5
    sigma_volume = None
    sigma_heat_capacity = None
    sigma_alpha = None
    if parameter == "volume":
        sigma = 1.0e-31
        step = 1.0e-32
        sigma_volume = sigma
    elif parameter == "heat_capacity":
        sigma = 1.0e-24
        step = 1.0e-25
        sigma_heat_capacity = sigma
    else:
        sigma = 1.0e-7
        step = 1.0e-8
        sigma_alpha = np.zeros((3, 3), dtype=np.float64)
        sigma_alpha[0, 0] = sigma

    propagated = adiabatic_stiffness_field(
        stiffness,
        temperature,
        volume,
        heat_capacity,
        alpha,
        sigma_volume_m3=sigma_volume,
        sigma_heat_capacity_j_per_k=sigma_heat_capacity,
        sigma_thermal_expansion_tensor=sigma_alpha,
    )
    assert propagated.sigma_stiffness is not None

    volume_plus = volume_minus = volume
    cv_plus = cv_minus = heat_capacity
    alpha_plus = alpha.copy()
    alpha_minus = alpha.copy()
    if parameter == "volume":
        volume_plus += step
        volume_minus -= step
    elif parameter == "heat_capacity":
        cv_plus += step
        cv_minus -= step
    else:
        alpha_plus[0, 0] += step
        alpha_minus[0, 0] -= step
    plus = adiabatic_stiffness_field(
        stiffness,
        temperature,
        volume_plus,
        cv_plus,
        alpha_plus,
    ).stiffness
    minus = adiabatic_stiffness_field(
        stiffness,
        temperature,
        volume_minus,
        cv_minus,
        alpha_minus,
    ).stiffness
    numerical_sigma = np.abs((plus - minus) / (2.0 * step)) * sigma
    np.testing.assert_allclose(
        propagated.sigma_stiffness,
        numerical_sigma,
        rtol=2.0e-6,
        atol=1.0e-12,
    )


def test_mixed_zero_temperature_field_emits_no_divide_warning() -> None:
    """Zero-K rows are masked before thermodynamic division is evaluated."""
    stiffness = np.broadcast_to(_isotropic_stiffness(100.0, 60.0), (2, 6, 6)).copy()
    temperature = np.asarray([0.0, 300.0], dtype=np.float64)
    volume = np.full(2, 1.0e-28, dtype=np.float64)
    cv = np.asarray([0.0, 4.0e-22], dtype=np.float64)
    alpha = np.broadcast_to(np.eye(3, dtype=np.float64) * 1.0e-5, (2, 3, 3)).copy()

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = adiabatic_stiffness_field(
            stiffness,
            temperature,
            volume,
            cv,
            alpha,
            sigma_volume_m3=np.asarray([0.0, 1.0e-31]),
            sigma_heat_capacity_j_per_k=np.asarray([0.0, 1.0e-24]),
        )

    assert np.all(result.valid_mask)
    np.testing.assert_array_equal(result.stiffness[0], stiffness[0])
