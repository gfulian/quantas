from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.eos import PressureEOS
from quantas.core.physics.eos import EnergyEOS

EOS_NAMES = ["murnaghan", "birchmurnaghan", "poirier-tarantola", "vinet"]
PARAMETERS = np.array([-100.0, 0.55, 4.2, 72.0], dtype=np.float64)


@pytest.mark.parametrize("name", EOS_NAMES)
def test_pressure_eos_recovers_reference_parameters(name: str):
    eos = PressureEOS()
    volume = PARAMETERS[3]

    pressure = float(np.asarray(eos.pressure(name, PARAMETERS, volume)))
    bulk_modulus = float(np.asarray(eos.bulk_modulus(name, PARAMETERS, volume)))
    derivative = float(
        np.asarray(eos.bulk_modulus_derivative(name, PARAMETERS, volume))
    )

    assert pressure == pytest.approx(0.0, abs=1.0e-12)
    assert bulk_modulus == pytest.approx(PARAMETERS[1], rel=1.0e-13)
    assert derivative == pytest.approx(PARAMETERS[2], rel=1.0e-13)


@pytest.mark.parametrize("name", EOS_NAMES)
def test_pressure_is_negative_energy_derivative(name: str):
    energy_eos = EnergyEOS()
    pressure_eos = PressureEOS()
    volume = np.linspace(66.0, 76.0, 9)
    step = 1.0e-5

    function = energy_eos.function(name)
    numerical_pressure = -(
        function(volume + step, *PARAMETERS) - function(volume - step, *PARAMETERS)
    ) / (2.0 * step)
    pressure = pressure_eos.pressure(name, PARAMETERS, volume)

    np.testing.assert_allclose(pressure, numerical_pressure, rtol=2.0e-7, atol=2.0e-10)


@pytest.mark.parametrize("name", EOS_NAMES)
def test_bulk_modulus_matches_pressure_volume_derivative(name: str):
    eos = PressureEOS()
    volume = np.linspace(66.0, 75.0, 7)
    step = 1.0e-5

    p_plus = eos.pressure(name, PARAMETERS, volume + step)
    p_minus = eos.pressure(name, PARAMETERS, volume - step)
    numerical = -volume * (p_plus - p_minus) / (2.0 * step)
    analytical = eos.bulk_modulus(name, PARAMETERS, volume)

    np.testing.assert_allclose(analytical, numerical, rtol=2.0e-8, atol=2.0e-10)


@pytest.mark.parametrize("name", EOS_NAMES)
def test_bulk_modulus_derivative_matches_numerical_result(name: str):
    eos = PressureEOS()
    volume = np.linspace(66.0, 75.0, 7)
    step = 1.0e-5

    k_plus = eos.bulk_modulus(name, PARAMETERS, volume + step)
    k_minus = eos.bulk_modulus(name, PARAMETERS, volume - step)
    p_plus = eos.pressure(name, PARAMETERS, volume + step)
    p_minus = eos.pressure(name, PARAMETERS, volume - step)
    numerical = (k_plus - k_minus) / (p_plus - p_minus)
    analytical = eos.bulk_modulus_derivative(name, PARAMETERS, volume)

    np.testing.assert_allclose(analytical, numerical, rtol=2.0e-8, atol=2.0e-9)


def test_pressure_eos_accepts_pressure_parameter_order():
    eos = PressureEOS()
    pressure_parameters = PARAMETERS[1:]
    volume = np.array([68.0, 72.0, 75.0])

    expected = eos.pressure("BM", PARAMETERS, volume)
    calculated = eos.pressure("BM", pressure_parameters, volume)

    np.testing.assert_allclose(calculated, expected, rtol=0.0, atol=0.0)


def test_energy_eos_delegates_pressure_to_shared_pressure_model():
    energy_eos = EnergyEOS()
    pressure_eos = PressureEOS()
    volume = np.array([68.0, 72.0, 75.0])

    for name in EOS_NAMES:
        expected = pressure_eos.pressure(name, PARAMETERS, volume)
        calculated = energy_eos.pressure(name, PARAMETERS, volume)
        np.testing.assert_allclose(calculated, expected, rtol=0.0, atol=0.0)


@pytest.mark.parametrize(
    ("tag", "parameters"),
    [
        ("M", [160.0, 4.6, 20.0]),
        ("BM2", [160.0, 20.0]),
        ("BM3", [160.0, 4.6, 20.0]),
        ("BM4", [160.0, 4.6, -0.025, 20.0]),
        ("PT2", [160.0, 20.0]),
        ("PT3", [160.0, 4.6, 20.0]),
        ("PT4", [160.0, 4.6, -0.025, 20.0]),
        ("V2", [160.0, 20.0]),
        ("V3", [160.0, 4.6, 20.0]),
        ("T2", [160.0, 20.0]),
        ("T3", [160.0, 4.6, 20.0]),
        ("T4", [160.0, 4.6, -0.025, 20.0]),
    ],
)
def test_every_pressure_order_recovers_reference_k_kp_and_kpp(tag, parameters):
    eos = PressureEOS()
    v0 = 20.0
    pressure = float(eos.pressure(tag, parameters, v0))
    k0 = float(eos.bulk_modulus(tag, parameters, v0))
    kp = float(eos.bulk_modulus_derivative(tag, parameters, v0))
    kpp = float(eos.bulk_modulus_second_derivative(tag, parameters, v0))

    from quantas.core.physics.eos.parameters import resolve_pressure_parameters

    resolved = resolve_pressure_parameters(tag, parameters)
    assert pressure == pytest.approx(0.0, abs=1.0e-12)
    assert k0 == pytest.approx(resolved.K0, rel=1.0e-12)
    assert kp == pytest.approx(resolved.KP, rel=1.0e-12)
    assert kpp == pytest.approx(resolved.KPP, rel=1.0e-10, abs=1.0e-12)


@pytest.mark.parametrize(
    ("tag", "parameters"),
    [
        ("BM2", [160.0, 20.0]),
        ("BM3", [160.0, 4.6, 20.0]),
        ("BM4", [160.0, 4.6, -0.025, 20.0]),
        ("PT2", [160.0, 20.0]),
        ("PT3", [160.0, 4.6, 20.0]),
        ("PT4", [160.0, 4.6, -0.025, 20.0]),
        ("V2", [160.0, 20.0]),
        ("V3", [160.0, 4.6, 20.0]),
        ("T2", [160.0, 20.0]),
        ("T3", [160.0, 4.6, 20.0]),
        ("T4", [160.0, 4.6, -0.025, 20.0]),
    ],
)
def test_second_bulk_modulus_derivative_matches_numerical_dkprime_dp(tag, parameters):
    eos = PressureEOS()
    volume = np.linspace(18.0, 22.0, 7)
    step = 1.0e-5
    kp_plus = eos.bulk_modulus_derivative(tag, parameters, volume + step)
    kp_minus = eos.bulk_modulus_derivative(tag, parameters, volume - step)
    p_plus = eos.pressure(tag, parameters, volume + step)
    p_minus = eos.pressure(tag, parameters, volume - step)
    numerical = (kp_plus - kp_minus) / (p_plus - p_minus)
    analytical = eos.bulk_modulus_second_derivative(tag, parameters, volume)

    np.testing.assert_allclose(analytical, numerical, rtol=2.0e-6, atol=2.0e-9)


@pytest.mark.parametrize(
    ("tag", "energy_parameters"),
    [
        ("BM2", [-100.0, 160.0, 20.0]),
        ("BM3", [-100.0, 160.0, 4.6, 20.0]),
        ("BM4", [-100.0, 160.0, 4.6, -0.025, 20.0]),
        ("PT2", [-100.0, 160.0, 20.0]),
        ("V2", [-100.0, 160.0, 20.0]),
    ],
)
def test_pressure_eos_accepts_matching_free_energy_parameter_vector(
    tag, energy_parameters
):
    energy_eos = EnergyEOS()
    pressure_eos = PressureEOS()
    volume = np.array([18.5, 20.0, 21.5])

    expected = energy_eos.pressure(tag, energy_parameters, volume)
    calculated = pressure_eos.pressure(tag, energy_parameters, volume)

    np.testing.assert_allclose(calculated, expected, rtol=0.0, atol=0.0)
