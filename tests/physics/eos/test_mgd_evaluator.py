"""Numerical validation of Mie--Grüneisen--Debye thermal pressure."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad

from quantas.core.physics.eos import (
    MGDNormalization,
    MGDParameters,
    MGDThermalPressure,
    PVTEOS,
    PVTModel,
    PressureEOS,
    debye_function_3,
)
from quantas.core.physics.units import N


def _quadrature_debye_3(argument: float) -> float:
    if argument == 0.0:
        return 1.0
    integral, _ = quad(
        lambda value: value**3 / np.expm1(value),
        0.0,
        argument,
        epsabs=1.0e-13,
        epsrel=1.0e-13,
        limit=300,
    )
    return 3.0 * integral / argument**3


@pytest.mark.parametrize(
    "argument",
    (0.0, 1.0e-8, 1.0e-4, 0.1, 1.0, 2.0, 3.0, 10.0, 50.0, 100.0),
)
def test_debye_function_matches_independent_quadrature(argument: float) -> None:
    calculated = float(debye_function_3(argument))
    expected = _quadrature_debye_3(argument)
    assert calculated == pytest.approx(expected, rel=2.0e-12, abs=2.0e-14)


def test_debye_function_validates_its_domain() -> None:
    with pytest.raises(ValueError, match="cannot be negative"):
        debye_function_3(-1.0)
    with pytest.raises(ValueError, match="must be finite"):
        debye_function_3(float("nan"))


def test_full_mgd_reference_identity_and_zero_temperature_limit() -> None:
    evaluator = MGDThermalPressure()
    parameters = MGDParameters(295.0, 459.0, 1.547, 0.94)
    normalization = MGDNormalization.cell(formula="NaF", formula_units_per_cell=4)
    volumes = np.array([140.0, 150.0, 160.0])

    reference = evaluator.pressure(
        "mgd", parameters, normalization, 150.0, volumes, 295.0
    )
    zero_temperature = evaluator.pressure(
        "mgd", parameters, normalization, 150.0, volumes, 0.0
    )

    assert np.array_equal(reference, np.zeros(3, dtype=np.float64))
    assert np.all(np.isfinite(zero_temperature))
    assert np.all(zero_temperature < 0.0)


def test_cell_and_molar_normalizations_are_equivalent() -> None:
    evaluator = MGDThermalPressure()
    parameters = MGDParameters(295.0, 459.0, 1.547, 0.94)
    cell = MGDNormalization.cell(formula="NaF", formula_units_per_cell=4)
    molar = MGDNormalization.molar_formula_unit(formula="NaF")
    reference_cell = 150.0
    volume_cell = np.array([142.0, 150.0, 158.0])
    reference_molar = N * reference_cell * 1.0e-24 / 4.0
    volume_molar = N * volume_cell * 1.0e-24 / 4.0

    cell_pressure = evaluator.pressure(
        "mgd",
        parameters,
        cell,
        reference_cell,
        volume_cell,
        900.0,
    )
    molar_pressure = evaluator.pressure(
        "mgd",
        parameters,
        molar,
        reference_molar,
        volume_molar,
        900.0,
    )

    assert molar_pressure == pytest.approx(cell_pressure, rel=2.0e-14, abs=2.0e-14)


def test_full_mgd_q_to_zero_limit_is_continuous() -> None:
    evaluator = MGDThermalPressure()
    normalization = MGDNormalization.cell(atoms_per_cell=8)
    exact = MGDParameters(295.0, 459.0, 1.547, 0.0)
    near = MGDParameters(295.0, 459.0, 1.547, 1.0e-12)
    volumes = np.array([120.0, 145.0, 150.0, 175.0])

    assert evaluator.gamma("mgd", exact, 150.0, volumes) == pytest.approx(
        evaluator.gamma("mgd", near, 150.0, volumes), rel=3.0e-13
    )
    assert evaluator.debye_temperature("mgd", exact, 150.0, volumes) == pytest.approx(
        evaluator.debye_temperature("mgd", near, 150.0, volumes),
        rel=3.0e-13,
    )
    assert evaluator.pressure(
        "mgd", exact, normalization, 150.0, volumes, 900.0
    ) == pytest.approx(
        evaluator.pressure("mgd", near, normalization, 150.0, volumes, 900.0),
        rel=3.0e-13,
        abs=1.0e-14,
    )


def test_q_compromise_has_parallel_isochors() -> None:
    evaluator = MGDThermalPressure()
    parameters = MGDParameters(295.0, 459.0, 1.547)
    normalization = MGDNormalization.cell(atoms_per_cell=8)
    volumes = np.array([120.0, 150.0, 180.0])

    gamma = evaluator.gamma("mgd:q-compromise", parameters, 150.0, volumes)
    theta = evaluator.debye_temperature("mgd:q-compromise", parameters, 150.0, volumes)
    pressure = evaluator.pressure(
        "mgd:q-compromise",
        parameters,
        normalization,
        150.0,
        volumes,
        900.0,
    )
    derivative = evaluator.volume_derivative(
        "mgd:q-compromise",
        parameters,
        normalization,
        150.0,
        volumes,
        900.0,
    )

    assert gamma == pytest.approx(parameters.gamma0 * volumes / 150.0)
    assert np.array_equal(theta, np.full(3, parameters.theta_d0))
    assert pressure == pytest.approx(np.full(3, pressure[0]), rel=0.0, abs=0.0)
    assert np.array_equal(derivative, np.zeros(3, dtype=np.float64))


def test_mgd_analytical_state_derivatives_match_finite_differences() -> None:
    evaluator = MGDThermalPressure()
    parameters = MGDParameters(295.0, 459.0, 1.547, 0.94)
    normalization = MGDNormalization.cell(atoms_per_cell=8)
    volume = 145.0
    temperature = 800.0
    volume_step = 1.0e-4
    temperature_step = 1.0e-3

    numerical_volume = (
        float(
            evaluator.pressure(
                "mgd",
                parameters,
                normalization,
                150.0,
                volume + volume_step,
                temperature,
            )
        )
        - float(
            evaluator.pressure(
                "mgd",
                parameters,
                normalization,
                150.0,
                volume - volume_step,
                temperature,
            )
        )
    ) / (2.0 * volume_step)
    numerical_temperature = (
        float(
            evaluator.pressure(
                "mgd",
                parameters,
                normalization,
                150.0,
                volume,
                temperature + temperature_step,
            )
        )
        - float(
            evaluator.pressure(
                "mgd",
                parameters,
                normalization,
                150.0,
                volume,
                temperature - temperature_step,
            )
        )
    ) / (2.0 * temperature_step)

    assert float(
        evaluator.volume_derivative(
            "mgd", parameters, normalization, 150.0, volume, temperature
        )
    ) == pytest.approx(numerical_volume, rel=2.0e-8, abs=1.0e-11)
    assert float(
        evaluator.temperature_derivative(
            "mgd", parameters, normalization, 150.0, volume, temperature
        )
    ) == pytest.approx(numerical_temperature, rel=2.0e-8, abs=1.0e-11)


@pytest.mark.parametrize("pressure_model", ("BM3", "Vinet"))
def test_pvt_mgd_composition_and_volume_inversion(pressure_model: str) -> None:
    normalization = MGDNormalization.cell(formula="NaF", formula_units_per_cell=4)
    model = PVTModel(
        pressure_model,
        "thermal-pressure",
        thermal_pressure_model="mgd",
        mgd_normalization=normalization,
    )
    pressure_parameters = {"K0": 48.0, "KP": 4.5, "KPP": -0.1, "V0": 150.0}
    coupling_parameters = {
        "temperature_ref": 295.0,
        "theta_d0": 459.0,
        "gamma0": 1.547,
        "q": 0.94,
    }
    evaluator = PVTEOS()
    volume = np.array([132.0, 145.0, 150.0, 160.0])
    temperature = np.array([1000.0, 800.0, 295.0, 500.0])

    pressure = evaluator.pressure(
        model,
        pressure_parameters,
        None,
        coupling_parameters,
        volume,
        temperature,
    )
    recovered = evaluator.volume(
        model,
        pressure_parameters,
        None,
        coupling_parameters,
        pressure,
        temperature,
    )

    assert recovered == pytest.approx(volume, rel=2.0e-11, abs=2.0e-10)
    expected_tag = "BM3" if pressure_model == "BM3" else "V3"
    assert model.tag == (f"{expected_tag}+mie-gruneisen-debye:full+thermal-pressure")


def test_reference_isotherm_reduces_to_static_pressure_eos() -> None:
    normalization = MGDNormalization.cell(atoms_per_cell=8)
    model = PVTModel(
        "BM3",
        "thermal-pressure",
        thermal_pressure_model="mgd",
        mgd_normalization=normalization,
    )
    pressure_parameters = {"K0": 48.0, "KP": 4.5, "KPP": -0.1, "V0": 150.0}
    coupling_parameters = {
        "temperature_ref": 295.0,
        "theta_d0": 459.0,
        "gamma0": 1.547,
        "q": 0.94,
    }
    volumes = np.array([135.0, 150.0, 165.0])

    pvt = PVTEOS().pressure(
        model,
        pressure_parameters,
        None,
        coupling_parameters,
        volumes,
        295.0,
    )
    static = PressureEOS().pressure("BM3", pressure_parameters, volumes)

    assert pvt == pytest.approx(static, rel=0.0, abs=0.0)


def test_pvt_model_requires_mgd_normalization() -> None:
    with pytest.raises(ValueError, match="requires mgd_normalization"):
        PVTModel(
            "BM3",
            "thermal-pressure",
            thermal_pressure_model="mgd",
        )
    with pytest.raises(ValueError, match="only valid"):
        PVTModel(
            "BM3",
            "thermal-pressure",
            mgd_normalization=MGDNormalization.cell(atoms_per_cell=8),
        )


def test_mgd_fitting_adapter_exposes_full_parameter_contract() -> None:
    from quantas.modules.eos import PVTEOSFitModel

    model = PVTModel(
        "BM3",
        "thermal-pressure",
        thermal_pressure_model="mgd",
        mgd_normalization=MGDNormalization.cell(atoms_per_cell=8),
    )
    adapter = PVTEOSFitModel(model)
    assert adapter.parameter_names == (
        "K0",
        "KP",
        "KPP",
        "V0",
        "temperature_ref",
        "theta_d0",
        "gamma0",
        "q",
    )


def test_composed_bulk_modulus_matches_total_pressure_derivative() -> None:
    normalization = MGDNormalization.cell(atoms_per_cell=8)
    model = PVTModel(
        "BM3",
        "thermal-pressure",
        thermal_pressure_model="mgd",
        mgd_normalization=normalization,
    )
    pressure_parameters = {"K0": 48.0, "KP": 4.5, "KPP": -0.1, "V0": 150.0}
    coupling_parameters = {
        "temperature_ref": 295.0,
        "theta_d0": 459.0,
        "gamma0": 1.547,
        "q": 0.94,
    }
    evaluator = PVTEOS()
    volume = 145.0
    temperature = 800.0
    step = 1.0e-4
    numerical = (
        -volume
        * (
            float(
                evaluator.pressure(
                    model,
                    pressure_parameters,
                    None,
                    coupling_parameters,
                    volume + step,
                    temperature,
                )
            )
            - float(
                evaluator.pressure(
                    model,
                    pressure_parameters,
                    None,
                    coupling_parameters,
                    volume - step,
                    temperature,
                )
            )
        )
        / (2.0 * step)
    )
    analytical = float(
        evaluator.bulk_modulus(
            model,
            pressure_parameters,
            None,
            coupling_parameters,
            volume,
            temperature,
        )
    )

    assert analytical == pytest.approx(numerical, rel=2.0e-8, abs=2.0e-8)
