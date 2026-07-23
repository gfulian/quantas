"""Tests for cold quasi-static finite-strain elasticity relations."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import (
    cold_finite_strain_stiffness,
    eulerian_finite_strain,
    wallace_hydrostatic_delta_voigt,
)


def test_eulerian_finite_strain_reference_and_compression() -> None:
    """The reference strain is zero and compression is positive."""
    values = eulerian_finite_strain([100.0, 90.0, 110.0], 100.0)
    assert values.dtype == np.float64
    assert values[0] == pytest.approx(0.0)
    assert values[1] > 0.0
    assert values[2] < 0.0


def test_wallace_delta_voigt_known_components() -> None:
    """Wallace's hydrostatic tensor has the expected normal and shear terms."""
    delta = wallace_hydrostatic_delta_voigt()
    assert delta.shape == (6, 6)
    assert np.allclose(delta, delta.T)
    assert delta[0, 0] == pytest.approx(-3.0)
    assert delta[0, 1] == pytest.approx(-1.0)
    assert delta[3, 3] == pytest.approx(-1.0)
    assert delta[0, 3] == pytest.approx(0.0)


def test_cold_finite_strain_returns_reference_stiffness_at_v0() -> None:
    """The cold finite-strain model exactly recovers C0 at zero strain."""
    c0 = np.diag([200.0, 200.0, 200.0, 80.0, 80.0, 80.0])
    c0[0, 1] = c0[1, 0] = 100.0
    c0[0, 2] = c0[2, 0] = 100.0
    c0[1, 2] = c0[2, 1] = 100.0
    cp = np.full((6, 6), 0.5)
    cp = 0.5 * (cp + cp.T)
    result = cold_finite_strain_stiffness(
        100.0,
        reference_volume=100.0,
        bulk_modulus=100.0,
        bulk_modulus_derivative=4.0,
        reference_stiffness=c0,
        stiffness_pressure_derivative=cp,
        order=3,
    )
    assert result.shape == (6, 6)
    assert np.allclose(result, c0)


def test_cold_finite_strain_vectorizes_over_volumes() -> None:
    """Several volumes are evaluated without changing matrix axes."""
    c0 = np.eye(6) * 100.0
    cp = np.eye(6)
    result = cold_finite_strain_stiffness(
        [90.0, 100.0, 110.0],
        reference_volume=100.0,
        bulk_modulus=120.0,
        bulk_modulus_derivative=4.5,
        reference_stiffness=c0,
        stiffness_pressure_derivative=cp,
    )
    assert result.shape == (3, 6, 6)
    assert np.allclose(result[1], c0)


def test_component_jacobian_matches_central_finite_differences() -> None:
    """The analytical QSA Jacobian supports covariance propagation."""
    from quantas.core.physics.elasticity import (
        cold_finite_strain_component,
        cold_finite_strain_component_jacobian,
    )

    volumes = np.asarray([94.0, 101.0, 108.0])
    parameters = np.asarray([100.0, 120.0, 4.5, 210.0, 5.0, 0.0])
    analytic = cold_finite_strain_component_jacobian(
        volumes,
        reference_volume=parameters[0],
        bulk_modulus=parameters[1],
        bulk_modulus_derivative=parameters[2],
        reference_component=parameters[3],
        component_pressure_derivative=parameters[4],
        wallace_delta=-3.0,
        order=3,
    )
    numerical = np.empty_like(analytic)
    for index in range(6):
        scale = max(abs(parameters[index]), 1.0)
        step = 1.0e-6 * scale
        plus = parameters.copy()
        minus = parameters.copy()
        if index < 5:
            plus[index] += step
            minus[index] -= step
            value_plus = cold_finite_strain_component(
                volumes,
                reference_volume=plus[0],
                bulk_modulus=plus[1],
                bulk_modulus_derivative=plus[2],
                reference_component=plus[3],
                component_pressure_derivative=plus[4],
                wallace_delta=-3.0,
                order=3,
            )
            value_minus = cold_finite_strain_component(
                volumes,
                reference_volume=minus[0],
                bulk_modulus=minus[1],
                bulk_modulus_derivative=minus[2],
                reference_component=minus[3],
                component_pressure_derivative=minus[4],
                wallace_delta=-3.0,
                order=3,
            )
        else:
            value_plus = cold_finite_strain_component(
                volumes + step,
                reference_volume=parameters[0],
                bulk_modulus=parameters[1],
                bulk_modulus_derivative=parameters[2],
                reference_component=parameters[3],
                component_pressure_derivative=parameters[4],
                wallace_delta=-3.0,
                order=3,
            )
            value_minus = cold_finite_strain_component(
                volumes - step,
                reference_volume=parameters[0],
                bulk_modulus=parameters[1],
                bulk_modulus_derivative=parameters[2],
                reference_component=parameters[3],
                component_pressure_derivative=parameters[4],
                wallace_delta=-3.0,
                order=3,
            )
        numerical[:, index] = (value_plus - value_minus) / (2.0 * step)
    numerical = numerical[:, [3, 4, 0, 1, 2, 5]]
    np.testing.assert_allclose(analytic, numerical, rtol=3.0e-6, atol=3.0e-7)
