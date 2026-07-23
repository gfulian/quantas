"""Scientific tests for compositional P-V-T equations of state."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.eos import (
    PVTCouplingFamily,
    PVTEOS,
    PVTModel,
    PressureEOS,
    available_eos_models,
    available_pvt_couplings,
    parse_pvt_coupling,
)


PRESSURE_PARAMETERS = {"K0": 160.0, "KP": 4.2, "KPP": -0.02, "V0": 100.0}
THERMAL_PARAMETERS = {
    "V0": 100.0,
    "temperature_ref": 300.0,
    "alpha0": 3.0e-5,
    "alpha1": 1.0e-8,
}


def _model(coupling: str) -> PVTModel:
    return PVTModel("BM3", coupling, "berman:quadratic")


def test_registry_and_model_validation() -> None:
    assert available_pvt_couplings() == tuple(PVTCouplingFamily)
    assert parse_pvt_coupling("dKdt") is PVTCouplingFamily.LINEAR_BULK_MODULUS
    assert parse_pvt_coupling("AG") is PVTCouplingFamily.ANDERSON_GRUNEISEN
    assert parse_pvt_coupling("Pth") is PVTCouplingFamily.THERMAL_PRESSURE
    with pytest.raises(ValueError, match="requires a temperature_model"):
        PVTModel("BM3", "linear")
    with pytest.raises(ValueError, match="must be omitted"):
        PVTModel("BM3", "thermal-pressure", "berman")


@pytest.mark.parametrize("coupling", ["linear", "anderson-gruneisen"])
def test_reference_isotherm_and_zero_pressure_volume(coupling: str) -> None:
    eos = PVTEOS()
    model = _model(coupling)
    coupling_parameters = {"dK0_dT": -0.02} if coupling == "linear" else {"delta": 4.0}
    volume = np.asarray([92.0, 100.0, 108.0])
    expected = PressureEOS().pressure("BM3", PRESSURE_PARAMETERS, volume)
    actual = eos.pressure(
        model,
        PRESSURE_PARAMETERS,
        THERMAL_PARAMETERS,
        coupling_parameters,
        volume,
        300.0,
    )
    np.testing.assert_allclose(actual, expected, rtol=2.0e-13, atol=2.0e-13)
    temperatures = np.asarray([300.0, 600.0, 900.0])
    v0 = eos.reference_volume(
        model,
        PRESSURE_PARAMETERS,
        THERMAL_PARAMETERS,
        coupling_parameters,
        temperatures,
    )
    zero = eos.pressure(
        model,
        PRESSURE_PARAMETERS,
        THERMAL_PARAMETERS,
        coupling_parameters,
        v0,
        temperatures,
    )
    np.testing.assert_allclose(zero, 0.0, atol=2.0e-12)


def test_linear_bulk_modulus_and_temperature_derivative() -> None:
    eos = PVTEOS()
    temperatures = np.asarray([300.0, 500.0, 900.0])
    k0 = eos.zero_pressure_bulk_modulus(
        _model("linear"),
        PRESSURE_PARAMETERS,
        THERMAL_PARAMETERS,
        {"dK0_dT": -0.02},
        temperatures,
    )
    np.testing.assert_allclose(k0, 160.0 - 0.02 * (temperatures - 300.0))
    derivative = eos.zero_pressure_bulk_modulus_temperature_derivative(
        _model("linear"),
        PRESSURE_PARAMETERS,
        THERMAL_PARAMETERS,
        {"dK0_dT": -0.02},
        temperatures,
    )
    np.testing.assert_allclose(derivative, -0.02)


def test_anderson_gruneisen_identity_and_derivative() -> None:
    eos = PVTEOS()
    model = _model("ag")
    temperatures = np.asarray([350.0, 600.0, 900.0])
    delta = 4.3
    v0 = eos.reference_volume(
        model, PRESSURE_PARAMETERS, THERMAL_PARAMETERS, {"delta": delta}, temperatures
    )
    k0 = eos.zero_pressure_bulk_modulus(
        model, PRESSURE_PARAMETERS, THERMAL_PARAMETERS, {"delta": delta}, temperatures
    )
    np.testing.assert_allclose(k0, 160.0 * (100.0 / v0) ** delta)
    alpha = eos.expansion_coefficient_zero_pressure(
        model, PRESSURE_PARAMETERS, THERMAL_PARAMETERS, {"delta": delta}, temperatures
    )
    derivative = eos.zero_pressure_bulk_modulus_temperature_derivative(
        model, PRESSURE_PARAMETERS, THERMAL_PARAMETERS, {"delta": delta}, temperatures
    )
    np.testing.assert_allclose(derivative, -delta * alpha * k0)


def test_thermal_pressure_reference_and_alpha_identity() -> None:
    eos = PVTEOS()
    model = PVTModel("BM3", "thermal-pressure")
    coupling = {"temperature_ref": 300.0, "alpha_ref": 3.0e-5, "theta_e": 500.0}
    pth = eos.thermal_pressure(PRESSURE_PARAMETERS, coupling, [0.0, 300.0, 800.0])
    assert pth[0] < 0.0
    assert pth[1] == pytest.approx(0.0, abs=1.0e-15)
    assert pth[2] > 0.0
    derivative = eos.thermal_pressure_temperature_derivative(
        PRESSURE_PARAMETERS, coupling, [300.0]
    )[0]
    assert derivative == pytest.approx(3.0e-5 * 160.0, rel=2.0e-14)
    alpha = eos.expansion_coefficient_zero_pressure(
        model, PRESSURE_PARAMETERS, None, coupling, [300.0]
    )[0]
    assert alpha == pytest.approx(3.0e-5, rel=2.0e-12)


def test_pressure_volume_round_trip_for_all_couplings() -> None:
    eos = PVTEOS()
    pressure = np.asarray([0.0, 2.0, 8.0, 15.0])
    temperature = np.asarray([300.0, 500.0, 700.0, 900.0])
    cases = [
        (_model("linear"), THERMAL_PARAMETERS, {"dK0_dT": -0.02}),
        (_model("ag"), THERMAL_PARAMETERS, {"delta": 4.0}),
        (
            PVTModel("BM3", "thermal-pressure"),
            None,
            {"temperature_ref": 300.0, "alpha_ref": 3.0e-5, "theta_e": 500.0},
        ),
    ]
    for model, thermal, coupling in cases:
        volume = eos.volume(
            model, PRESSURE_PARAMETERS, thermal, coupling, pressure, temperature
        )
        calculated = eos.pressure(
            model, PRESSURE_PARAMETERS, thermal, coupling, volume, temperature
        )
        np.testing.assert_allclose(calculated, pressure, rtol=1.0e-10, atol=2.0e-10)


@pytest.mark.parametrize(
    "pressure_model", available_eos_models(), ids=lambda item: item.tag
)
def test_linear_coupling_supports_every_pressure_eos(pressure_model) -> None:
    eos = PVTEOS()
    parameters = {"K0": 160.0, "KP": 4.2, "KPP": -0.02, "V0": 100.0}
    model = PVTModel(pressure_model, "linear", "berman:linear")
    pressure = eos.pressure(
        model,
        parameters,
        {"V0": 100.0, "temperature_ref": 300.0, "alpha0": 3.0e-5},
        {"dK0_dT": -0.02},
        [95.0, 100.0, 105.0],
        [300.0, 500.0, 700.0],
    )
    assert pressure.shape == (3,)
    assert np.all(np.isfinite(pressure))


def test_reference_volume_mismatch_is_rejected() -> None:
    with pytest.raises(ValueError, match="same reference volume"):
        PVTEOS().pressure(
            _model("linear"),
            PRESSURE_PARAMETERS,
            {**THERMAL_PARAMETERS, "V0": 101.0},
            {"dK0_dT": -0.02},
            [100.0],
            [300.0],
        )


def test_linear_coupling_invalid_bulk_modulus_reports_temperature_and_units() -> None:
    eos = PVTEOS()
    with pytest.raises(ValueError) as caught:
        eos.zero_pressure_bulk_modulus(
            _model("linear"),
            PRESSURE_PARAMETERS,
            THERMAL_PARAMETERS,
            {"dK0_dT": -0.8},
            [300.0, 1065.0],
        )

    message = str(caught.value)
    assert "T=1065 K" in message
    assert "K0(T)=" in message
    assert "GPa" in message
    assert "dK0_dT=-0.8 GPa K^-1" in message
    assert "Tref=300 K" in message
