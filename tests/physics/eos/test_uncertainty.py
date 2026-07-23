from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.eos import EOSStateUncertainty, FittedEnergyEOS

PARAMETERS = np.array([-100.0, 0.55, 4.2, 72.0], dtype=np.float64)
SAMPLED_VOLUMES = np.linspace(66.0, 78.0, 13)


def _murnaghan_jacobian(pressure: float) -> np.ndarray:
    _, k0, kp, v0 = PARAMETERS
    factor = 1.0 + kp * pressure / k0
    volume = v0 * factor ** (-1.0 / kp)

    dvolume_dk0 = volume * pressure / (factor * k0**2)
    dvolume_dkp = volume * (np.log(factor) / kp**2 - pressure / (kp * k0 * factor))
    dvolume_dv0 = volume / v0

    return np.array(
        [
            [0.0, dvolume_dk0, dvolume_dkp, dvolume_dv0],
            [0.0, 1.0, pressure, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ],
        dtype=np.float64,
    )


def test_state_jacobian_matches_murnaghan_closed_form():
    model = FittedEnergyEOS("murnaghan", PARAMETERS, sampled_volumes=SAMPLED_VOLUMES)
    pressure = 0.05

    calculated = model.state_jacobian_at_pressure(pressure)
    expected = _murnaghan_jacobian(pressure)

    np.testing.assert_allclose(calculated, expected, rtol=2.0e-6, atol=2.0e-9)


def test_covariance_propagation_matches_linear_transformation():
    covariance = np.array(
        [
            [4.0e-6, 0.0, 0.0, 0.0],
            [0.0, 2.5e-5, -4.0e-5, 1.0e-4],
            [0.0, -4.0e-5, 1.6e-3, -2.0e-3],
            [0.0, 1.0e-4, -2.0e-3, 4.0e-2],
        ],
        dtype=np.float64,
    )
    model = FittedEnergyEOS(
        "murnaghan",
        PARAMETERS,
        sampled_volumes=SAMPLED_VOLUMES,
        covariance=covariance,
    )
    pressure = 0.05

    state = model.state_at_pressure(pressure, uncertainty_method="covariance")

    assert isinstance(state.uncertainty, EOSStateUncertainty)
    assert state.uncertainty.method == "covariance"
    expected_jacobian = _murnaghan_jacobian(pressure)
    expected_covariance = expected_jacobian @ covariance @ expected_jacobian.T
    np.testing.assert_allclose(
        state.uncertainty.covariance,
        expected_covariance,
        rtol=3.0e-6,
        atol=2.0e-10,
    )
    assert state.sigma_volume == pytest.approx(np.sqrt(expected_covariance[0, 0]))
    assert state.sigma_bulk_modulus == pytest.approx(np.sqrt(expected_covariance[1, 1]))
    assert state.sigma_bulk_modulus_derivative == pytest.approx(
        np.sqrt(expected_covariance[2, 2])
    )


def test_energy_offset_uncertainty_does_not_affect_pressure_state():
    covariance = np.zeros((4, 4), dtype=np.float64)
    covariance[0, 0] = 100.0
    model = FittedEnergyEOS("BM", PARAMETERS, covariance=covariance)

    state = model.state_at_pressure(0.05, uncertainty_method="covariance")

    assert state.uncertainty is not None
    np.testing.assert_allclose(state.uncertainty.covariance, np.zeros((3, 3)))
    np.testing.assert_allclose(state.uncertainty.standard_deviations, np.zeros(3))


def test_uncertainty_requires_parameter_covariance():
    model = FittedEnergyEOS("V", PARAMETERS)

    with pytest.raises(ValueError, match="covariance"):
        model.state_at_pressure(0.05, uncertainty_method="covariance")
    with pytest.raises(ValueError, match="covariance"):
        model.state_at_pressure(0.05, uncertainty_method="montecarlo")


def test_covariance_validation_rejects_asymmetric_and_indefinite_matrices():
    asymmetric = np.eye(4)
    asymmetric[0, 1] = 0.1
    with pytest.raises(ValueError, match="symmetric"):
        FittedEnergyEOS("BM", PARAMETERS, covariance=asymmetric)

    indefinite = np.eye(4)
    indefinite[0, 1] = indefinite[1, 0] = 2.0
    with pytest.raises(ValueError, match="positive semidefinite"):
        FittedEnergyEOS("BM", PARAMETERS, covariance=indefinite)


def test_monte_carlo_matches_delta_method():
    covariance = np.diag([1.0e-8, 4.0e-6, 4.0e-4, 1.0e-2])
    model = FittedEnergyEOS(
        "murnaghan",
        PARAMETERS,
        sampled_volumes=SAMPLED_VOLUMES,
        covariance=covariance,
    )

    covariance_state = model.state_at_pressure(0.05, uncertainty_method="covariance")
    monte_carlo_1 = model.state_at_pressure(
        0.05,
        uncertainty_method="montecarlo",
        monte_carlo_samples=3000,
        random_state=1234,
    )
    monte_carlo_2 = model.state_at_pressure(
        0.05,
        uncertainty_method="montecarlo",
        monte_carlo_samples=3000,
        random_state=1234,
    )

    assert monte_carlo_1.uncertainty is not None
    assert monte_carlo_2.uncertainty is not None
    assert monte_carlo_1.uncertainty.method == "montecarlo"
    np.testing.assert_allclose(
        monte_carlo_1.uncertainty.covariance,
        monte_carlo_2.uncertainty.covariance,
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        monte_carlo_1.uncertainty.standard_deviations,
        covariance_state.uncertainty.standard_deviations,
        rtol=0.12,
        atol=2.0e-5,
    )
    assert monte_carlo_1.uncertainty.n_successful == 3000
    assert monte_carlo_1.uncertainty.confidence_intervals is not None


@pytest.mark.parametrize(
    "name",
    ["murnaghan", "birchmurnaghan", "poirier-tarantola", "vinet"],
)
def test_covariance_uncertainty_is_finite_for_all_supported_eos(name: str):
    covariance = np.diag([1.0e-8, 4.0e-6, 4.0e-4, 1.0e-2])
    model = FittedEnergyEOS(
        name,
        PARAMETERS,
        sampled_volumes=SAMPLED_VOLUMES,
        covariance=covariance,
    )

    state = model.state_at_pressure(0.05, uncertainty_method="covariance")

    assert state.uncertainty is not None
    assert np.all(np.isfinite(state.uncertainty.covariance))
    assert np.all(state.uncertainty.standard_deviations >= 0.0)


def test_state_dictionary_serializes_uncertainty_arrays():
    covariance = np.diag([0.0, 1.0e-5, 1.0e-3, 1.0e-2])
    model = FittedEnergyEOS("PT", PARAMETERS, covariance=covariance)

    state = model.state_at_pressure(0.02, uncertainty_method="covariance")
    serialized = state.as_dict()

    assert serialized["uncertainty"]["method"] == "covariance"
    assert len(serialized["uncertainty"]["covariance"]) == 3
    assert len(serialized["uncertainty"]["standard_deviations"]) == 3
