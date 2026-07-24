from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.eos import EnergyEOS
from quantas.core.math.fitting import FitQuality, FitStatus


def _reference_murnaghan(volume, E0, K0, KP, V0):
    return (
        E0
        + K0 * volume / KP * (((V0 / volume) ** KP) / (KP - 1) + 1)
        - V0 * K0 / (KP - 1)
    )


def _reference_birchmurnaghan(volume, E0, K0, KP, V0):
    eta = (V0 / volume) ** (1.0 / 3.0)
    return E0 + (9.0 * K0 * V0 / 16.0) * (
        KP * (eta**2 - 1) ** 3 + ((eta**2 - 1) ** 2) * (6.0 - 4.0 * eta**2)
    )


def _reference_poirier_tarantola(volume, E0, K0, KP, V0):
    eta = (volume / V0) ** (1.0 / 3.0)
    rho = -3.0 * np.log(eta)
    return E0 + (K0 * V0 * rho**2) / 6.0 * (3.0 + rho * (KP - 2.0))


def _reference_vinet(volume, E0, K0, KP, V0):
    eta = (volume / V0) ** (1.0 / 3.0)
    return E0 + 2.0 * K0 * V0 / (KP - 1.0) ** 2 * (
        2.0
        - (5.0 + 3.0 * KP * (eta - 1.0) - 3.0 * eta)
        * np.exp(-3.0 * (KP - 1.0) * (eta - 1.0) / 2.0)
    )


def test_energy_eos_formulas_match_reference_values():
    eos = EnergyEOS()
    volume = np.array([68.0, 70.0, 72.0, 74.0, 76.0], dtype=np.float64)
    pars = (-100.0, 0.55, 4.2, 72.0)

    np.testing.assert_allclose(
        eos.murnaghan(volume, *pars),
        _reference_murnaghan(volume, *pars),
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        eos.birchmurnaghan(volume, *pars),
        _reference_birchmurnaghan(volume, *pars),
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        eos.poirier_tarantola(volume, *pars),
        _reference_poirier_tarantola(volume, *pars),
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        eos.vinet(volume, *pars), _reference_vinet(volume, *pars), rtol=0.0, atol=0.0
    )


def test_energy_eos_fits_return_structured_diagnostics():
    eos = EnergyEOS()
    volume = np.linspace(68.0, 76.0, 11)
    expected = np.array([-100.0, 0.55, 4.2, 72.0], dtype=np.float64)

    for tag in ["murnaghan", "birchmurnaghan", "poirier-tarantola", "vinet"]:
        energy = eos.function(tag)(volume, *expected)
        result = eos.fit(tag, volume, energy)
        assert result.success, result.message
        assert result.status is FitStatus.SUCCESS
        assert result.metadata["model"] == eos.canonical_name(tag)
        assert result.parameters is not None
        np.testing.assert_allclose(
            result.parameters, expected, rtol=1.0e-9, atol=1.0e-9
        )
        assert result.errors is not None
        assert result.errors.shape == result.parameters.shape
        assert result.residuals is not None
        assert result.rmse is not None


def test_energy_eos_fit_parameters_preserves_tuple_access():
    eos = EnergyEOS()
    volume = np.linspace(68.0, 76.0, 11)
    expected = np.array([-100.0, 0.55, 4.2, 72.0], dtype=np.float64)
    energy = eos.birchmurnaghan(volume, *expected)

    fitted, errors = eos.fit_parameters("BM", volume, energy)

    np.testing.assert_allclose(fitted, expected, rtol=1.0e-9, atol=1.0e-9)
    assert errors.shape == fitted.shape


def test_energy_fit_returns_failure_for_invalid_data():
    eos = EnergyEOS()
    result = eos.fit("BM", [1.0, 2.0], [1.0])

    assert not result.success
    assert result.status is FitStatus.INVALID_INPUT


def test_energy_fit_accepts_exactly_determined_four_point_bm3():
    eos = EnergyEOS()
    volume = np.array([68.0, 70.5, 73.0, 76.0], dtype=np.float64)
    expected = np.array([-100.0, 0.55, 4.2, 72.0], dtype=np.float64)
    energy = eos.evaluate("BM3", volume, expected)

    result = eos.fit("BM3", volume, energy)

    assert result.success, result.message
    assert result.status is FitStatus.SUCCESS
    assert result.quality is FitQuality.POOR
    assert result.dof == 0
    assert result.covariance is None
    assert result.errors is None
    assert result.metadata["covariance_status"] == (
        "unavailable_zero_degrees_of_freedom"
    )
    assert any("exactly determined" in warning for warning in result.warnings)
    np.testing.assert_allclose(result.parameters, expected, rtol=1.0e-7, atol=1.0e-8)


def test_energy_fit_rejects_four_points_for_five_parameter_bm4():
    eos = EnergyEOS()
    volume = np.array([68.0, 70.5, 73.0, 76.0], dtype=np.float64)
    parameters = np.array([-100.0, 0.55, 4.2, -0.04, 72.0], dtype=np.float64)
    energy = eos.evaluate("BM4", volume, parameters)

    result = eos.fit("BM4", volume, energy)

    assert not result.success
    assert result.status is FitStatus.INVALID_INPUT
    assert "smaller than the number of parameters" in result.message


def test_energy_eos_pressure_is_zero_at_reference_volume():
    eos = EnergyEOS()
    pars = np.array([-100.0, 0.55, 4.2, 72.0], dtype=np.float64)
    volume = np.array([pars[3]], dtype=np.float64)

    for tag in ["M", "BM", "PT", "V"]:
        np.testing.assert_allclose(
            eos.pressure(tag, pars, volume), np.array([0.0]), atol=1.0e-14
        )


def test_energy_eos_rejects_unknown_tag():
    eos = EnergyEOS()
    with np.testing.assert_raises(ValueError):
        eos.function("unknown")


@pytest.mark.parametrize(
    ("tag", "parameters"),
    [
        ("M", [-100.0, 0.55, 4.2, 72.0]),
        ("BM2", [-100.0, 0.55, 72.0]),
        ("BM3", [-100.0, 0.55, 4.2, 72.0]),
        ("BM4", [-100.0, 0.55, 4.2, -0.04, 72.0]),
        ("PT2", [-100.0, 0.55, 72.0]),
        ("PT3", [-100.0, 0.55, 4.2, 72.0]),
        ("PT4", [-100.0, 0.55, 4.2, -0.04, 72.0]),
        ("V2", [-100.0, 0.55, 72.0]),
        ("V3", [-100.0, 0.55, 4.2, 72.0]),
    ],
)
def test_every_integrated_order_differentiates_to_matching_pressure(tag, parameters):
    eos = EnergyEOS()
    volume = np.linspace(67.0, 77.0, 9)
    step = 1.0e-5

    energy_plus = eos.evaluate(tag, volume + step, parameters)
    energy_minus = eos.evaluate(tag, volume - step, parameters)
    numerical = -(energy_plus - energy_minus) / (2.0 * step)
    pressure = eos.pressure(tag, parameters, volume)

    np.testing.assert_allclose(pressure, numerical, rtol=3.0e-7, atol=2.0e-9)


@pytest.mark.parametrize(
    ("tag", "expected_names"),
    [
        ("BM2", ["E0", "K0", "V0"]),
        ("BM3", ["E0", "K0", "KP", "V0"]),
        ("BM4", ["E0", "K0", "KP", "KPP", "V0"]),
    ],
)
def test_birch_murnaghan_fit_uses_order_dependent_free_parameters(tag, expected_names):
    eos = EnergyEOS()
    volume = np.linspace(68.0, 76.0, 15)
    full = {
        "E0": -100.0,
        "K0": 0.55,
        "KP": 4.2,
        "KPP": -0.04,
        "V0": 72.0,
    }
    energy = eos.evaluate(tag, volume, full)

    result = eos.fit(tag, volume, energy, p0=full)

    assert result.success, result.message
    assert result.metadata["eos_tag"] == tag
    assert result.metadata["parameter_order"] == expected_names
    assert len(result.parameters) == len(expected_names)
    assert result.metadata["resolved_parameter_order"] == [
        "E0",
        "K0",
        "KP",
        "KPP",
        "V0",
    ]
