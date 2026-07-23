from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.eos import EOSState, EnergyEOS, FittedEnergyEOS

EOS_NAMES = ["murnaghan", "birchmurnaghan", "poirier-tarantola", "vinet"]
PARAMETERS = np.array([-100.0, 0.55, 4.2, 72.0], dtype=np.float64)
SAMPLED_VOLUMES = np.linspace(66.0, 78.0, 13)


@pytest.mark.parametrize("name", EOS_NAMES)
def test_fitted_eos_returns_reference_state_at_zero_pressure(name: str):
    model = FittedEnergyEOS(name, PARAMETERS, sampled_volumes=SAMPLED_VOLUMES)

    state = model.state_at_pressure(0.0)

    assert isinstance(state, EOSState)
    assert state.volume == pytest.approx(PARAMETERS[3], rel=1.0e-12)
    assert state.bulk_modulus == pytest.approx(PARAMETERS[1], rel=1.0e-12)
    assert state.bulk_modulus_derivative == pytest.approx(PARAMETERS[2], rel=1.0e-12)
    assert not state.extrapolated
    assert state.metadata["eos"] == model.eos
    assert state.metadata["eos_tag"] == model.eos_tag
    assert state.metadata["eos_order"] == model.eos_order


@pytest.mark.parametrize("name", EOS_NAMES)
def test_volume_pressure_round_trip(name: str):
    model = FittedEnergyEOS(name, PARAMETERS, sampled_volumes=SAMPLED_VOLUMES)
    volumes = np.array([67.0, 69.0, 72.0, 75.0, 77.0])

    pressures = model.pressure(volumes)
    recovered = np.array(
        [model.volume_at_pressure(float(value)) for value in pressures]
    )

    np.testing.assert_allclose(recovered, volumes, rtol=2.0e-11, atol=2.0e-10)


def test_murnaghan_state_matches_closed_form_pressure_relations():
    model = FittedEnergyEOS("murnaghan", PARAMETERS, sampled_volumes=SAMPLED_VOLUMES)
    pressures = np.array([0.00, 0.02, 0.05, 0.10])

    states = model.states_at_pressures(pressures)
    volumes = np.array([state.volume for state in states])
    moduli = np.array([state.bulk_modulus for state in states])
    derivatives = np.array([state.bulk_modulus_derivative for state in states])

    k0 = PARAMETERS[1]
    kp = PARAMETERS[2]
    v0 = PARAMETERS[3]
    expected_volumes = v0 * (1.0 + kp * pressures / k0) ** (-1.0 / kp)

    np.testing.assert_allclose(volumes, expected_volumes, rtol=1.0e-11, atol=1.0e-11)
    np.testing.assert_allclose(moduli, k0 + kp * pressures, rtol=1.0e-11, atol=1.0e-11)
    np.testing.assert_allclose(
        derivatives, np.full_like(pressures, kp), rtol=0.0, atol=0.0
    )


def test_state_marks_volume_outside_sampled_fit_interval():
    model = FittedEnergyEOS(
        "BM", PARAMETERS, sampled_volumes=np.linspace(69.0, 75.0, 7)
    )

    state = model.state_at_pressure(0.15)

    assert state.volume < 69.0
    assert state.extrapolated


def test_model_keeps_copies_of_parameters_and_covariance():
    covariance = np.eye(4) * 0.01
    model = FittedEnergyEOS("V", PARAMETERS, covariance=covariance)

    parameters = model.parameters
    stored_covariance = model.covariance
    parameters[1] = 999.0
    assert stored_covariance is not None
    stored_covariance[0, 0] = 999.0

    assert model.reference_bulk_modulus == pytest.approx(PARAMETERS[1])
    assert model.covariance is not None
    assert model.covariance[0, 0] == pytest.approx(0.01)


def test_model_rejects_invalid_parameters_and_covariance():
    with pytest.raises(ValueError, match="E0, K0, KP, V0"):
        FittedEnergyEOS("BM", [0.5, 4.0, 70.0])
    with pytest.raises(ValueError, match="K0 must be positive"):
        FittedEnergyEOS("BM", [-100.0, -0.5, 4.0, 70.0])
    with pytest.raises(ValueError, match="shape"):
        FittedEnergyEOS("BM", PARAMETERS, covariance=np.eye(3))


def test_model_rejects_unreachable_tensile_pressure():
    model = FittedEnergyEOS("murnaghan", PARAMETERS)
    limiting_pressure = -PARAMETERS[1] / PARAMETERS[2]

    with pytest.raises(ValueError, match="could not be bracketed"):
        model.volume_at_pressure(limiting_pressure - 0.01)


@pytest.mark.parametrize(
    ("tag", "parameters", "covariance_shape"),
    [
        ("BM2", [-100.0, 0.55, 72.0], (3, 3)),
        ("BM3", [-100.0, 0.55, 4.2, 72.0], (4, 4)),
        ("BM4", [-100.0, 0.55, 4.2, -0.04, 72.0], (5, 5)),
    ],
)
def test_fit_accepts_order_dependent_shapes(
    tag, parameters, covariance_shape
):
    model = FittedEnergyEOS(
        tag, parameters, covariance=np.eye(len(parameters)) * 1.0e-8
    )

    assert model.eos_tag == tag
    assert model.covariance.shape == covariance_shape
    state = model.state_at_pressure(0.01, uncertainty_method="covariance")
    assert state.uncertainty is not None
    assert np.all(np.isfinite(state.uncertainty.covariance))


def test_bm2_resolves_kp_and_kpp_while_bm4_keeps_fitted_kpp():
    bm2 = FittedEnergyEOS("BM2", [-100.0, 0.55, 72.0])
    bm4 = FittedEnergyEOS("BM4", [-100.0, 0.55, 4.2, -0.04, 72.0])

    assert bm2.reference_bulk_modulus_derivative == pytest.approx(4.0)
    assert bm2.reference_bulk_modulus_second_derivative == pytest.approx(
        -35.0 / (9.0 * 0.55)
    )
    assert bm4.reference_bulk_modulus_second_derivative == pytest.approx(-0.04)


def test_fitted_bm2_exposes_complete_parameter_uncertainties() -> None:
    covariance = np.diag([1.0e-8, 4.0e-6, 1.0e-4])
    model = FittedEnergyEOS(
        "BM2",
        [-100.0, 0.55, 72.0],
        covariance=covariance,
    )

    resolved_covariance = model.resolved_covariance
    resolved_errors = model.resolved_errors

    assert resolved_covariance is not None
    assert resolved_errors is not None
    assert resolved_covariance.shape == (5, 5)
    assert resolved_errors.shape == (5,)
    assert resolved_errors[2] == pytest.approx(0.0)
    assert resolved_errors[3] > 0.0


@pytest.mark.parametrize(
    ("tag", "parameters", "expected_kp"),
    [
        ("M", [-100.0, 0.55, 4.2, 72.0], 4.2),
        ("BM2", [-100.0, 0.55, 72.0], 4.0),
        ("BM3", [-100.0, 0.55, 4.2, 72.0], 4.2),
        ("BM4", [-100.0, 0.55, 4.2, -0.04, 72.0], 4.2),
        ("PT2", [-100.0, 0.55, 72.0], 2.0),
        ("PT3", [-100.0, 0.55, 4.2, 72.0], 4.2),
        ("PT4", [-100.0, 0.55, 4.2, -0.04, 72.0], 4.2),
        ("V2", [-100.0, 0.55, 72.0], 1.0),
        ("V3", [-100.0, 0.55, 4.2, 72.0], 4.2),
    ],
)
def test_integrated_and_pressure_eos_share_the_same_reference_state(
    tag: str,
    parameters: list[float],
    expected_kp: float,
) -> None:
    model = FittedEnergyEOS(tag, parameters)

    state = model.state_at_pressure(0.0)

    assert state.volume == pytest.approx(72.0, rel=1.0e-12)
    assert state.bulk_modulus == pytest.approx(0.55, rel=1.0e-12)
    assert state.bulk_modulus_derivative == pytest.approx(expected_kp, rel=1.0e-12)


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
def test_pressure_solution_matches_direct_enthalpy_minimization(tag, parameters):
    from scipy.optimize import minimize_scalar

    pressure = 0.04
    energy_eos = EnergyEOS()
    fitted = FittedEnergyEOS(tag, parameters)

    direct = minimize_scalar(
        lambda volume: float(
            energy_eos.evaluate(tag, volume, parameters) + pressure * volume
        ),
        bounds=(55.0, 90.0),
        method="bounded",
        options={"xatol": 1.0e-12},
    )

    assert direct.success
    solved = fitted.volume_at_pressure(pressure)
    assert solved == pytest.approx(direct.x, rel=0.0, abs=3.0e-6)
