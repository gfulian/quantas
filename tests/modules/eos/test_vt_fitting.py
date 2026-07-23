"""Tests for frontend-neutral volumetric and axial V-T fitting."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.math.fitting import (
    EffectiveVarianceOptions,
    FitMethod,
    ODROptions,
    OLSOptions,
    ParameterState,
    WLSOptions,
)
from quantas.core.physics.eos import (
    TemperatureEOS,
    TemperatureEOSFamily,
    TemperatureEOSVariant,
    available_temperature_eos_models,
    parse_temperature_eos_model,
)
from quantas.modules.eos import (
    EOSDataset,
    EOSFitDomain,
    EOSFitOptions,
    EOSFitRequest,
    EOSFitter,
    ParameterConstraint,
    TemperatureEOSFitModel,
    build_temperature_parameter_map,
    read_eos_input,
)

_DATA = Path(__file__).with_name("data")


def _parameters(model) -> dict[str, float]:
    common = {"V0": 100.0, "temperature_ref": 298.15}
    if model.family is TemperatureEOSFamily.BERMAN:
        return {**common, "alpha0": 3.0e-5, "alpha1": 1.0e-8}
    if model.family is TemperatureEOSFamily.FEI:
        return {
            **common,
            "alpha0": 2.5e-5,
            "alpha1": 1.0e-8,
            "alpha2": 0.08,
        }
    if model.family is TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL:
        return {**common, "alpha0": 3.0e-5, "alpha1": 5.0e-5}
    if model.family is TemperatureEOSFamily.SALJE:
        return {
            "V0": 99.5,
            "temperature_ref": 0.0,
            "p1": 1.0e-4,
            "theta_sat": 260.0,
        }
    return {
        **common,
        "alpha_ref": 3.0e-5,
        "theta_e": 500.0,
        "kp": 4.0,
    }


def _synthetic(model, *, axial: bool = False) -> tuple[EOSDataset, dict[str, float]]:
    model = parse_temperature_eos_model(model)
    temperature = np.linspace(120.0, 900.0, 48, dtype=np.float64)
    parameters = _parameters(model)
    value = TemperatureEOS().value(model, parameters, temperature)
    value = value * (1.0 + 2.0e-7 * np.sin(np.arange(value.size)))
    if axial:
        length = np.cbrt(value)
        dataset = EOSDataset(
            jobname=f"synthetic-axial-{model.tag}",
            columns={
                "temperature": temperature,
                "sigma_temperature": np.full(temperature.shape, 0.2),
                "a": length,
                "sigma_a": np.full(temperature.shape, 2.0e-5),
            },
            units={"temperature": "K", "a": "angstrom"},
        )
    else:
        dataset = EOSDataset(
            jobname=f"synthetic-{model.tag}",
            columns={
                "temperature": temperature,
                "sigma_temperature": np.full(temperature.shape, 0.2),
                "volume": value,
                "sigma_volume": np.full(temperature.shape, 2.0e-3),
            },
            units={"temperature": "K", "volume": "angstrom^3"},
        )
    return dataset, parameters


@pytest.mark.parametrize(
    "model", available_temperature_eos_models(), ids=lambda item: item.tag
)
def test_ols_fits_every_temperature_model(model) -> None:
    dataset, expected = _synthetic(model)
    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model=model,
            domain=EOSFitDomain.VOLUME_TEMPERATURE,
        ),
    )

    assert result.fit.success
    assert result.fit.method is FitMethod.OLS
    assert result.parameter_values["V0"] == pytest.approx(expected["V0"], rel=2.0e-4)
    assert result.predictions["volume"].shape == dataset.column("volume").shape
    assert result.predictions["expansion_coefficient"].shape == (dataset.npoints,)
    assert result.metadata["relationship"] == "volume(temperature)"


def test_temperature_parameter_map_tracks_fixed_and_implied_parameters() -> None:
    dataset, _ = _synthetic(
        next(
            model
            for model in available_temperature_eos_models()
            if model.family is TemperatureEOSFamily.BERMAN
            and model.variant is TemperatureEOSVariant.LINEAR
        )
    )
    mapping = build_temperature_parameter_map(
        "berman:linear",
        dataset.column("temperature"),
        dataset.column("volume"),
    )

    assert mapping.names == ("V0", "temperature_ref", "alpha0", "alpha1")
    assert mapping.states == (
        ParameterState.FREE,
        ParameterState.FIXED,
        ParameterState.FREE,
        ParameterState.IMPLIED,
    )
    assert mapping.expand(mapping.initial_free_values()).as_mapping()["alpha1"] == 0.0


def test_khp_defaults_keep_poorly_identified_parameters_fixed() -> None:
    dataset, _ = _synthetic("khp")
    mapping = build_temperature_parameter_map(
        "khp", dataset.column("temperature"), dataset.column("volume")
    )

    states = dict(zip(mapping.names, mapping.states, strict=True))
    assert states["temperature_ref"] is ParameterState.FIXED
    assert states["theta_e"] is ParameterState.FIXED
    assert states["kp"] is ParameterState.FIXED
    assert mapping.free_names == ("V0", "alpha_ref")


def test_reference_temperature_can_be_changed_but_not_refined() -> None:
    dataset, _ = _synthetic("berman")
    request = EOSFitRequest(
        model="berman",
        domain="vt",
        constraints=(ParameterConstraint.fixed("TREF", 400.0),),
    )
    result = EOSFitter().fit(dataset, request)

    assert result.fit.success
    assert result.parameter_values["temperature_ref"] == 400.0
    with pytest.raises(ValueError, match="must be fixed"):
        build_temperature_parameter_map(
            "berman",
            dataset.column("temperature"),
            dataset.column("volume"),
            (ParameterConstraint.free("temperature_ref", 300.0),),
        )


def test_effective_variance_uses_temperature_derivative() -> None:
    dataset, _ = _synthetic("berman:quadratic")
    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="berman:quadratic",
            domain="vt",
            options=EOSFitOptions(
                solver_options=EffectiveVarianceOptions(max_iterations=20)
            ),
        ),
    )

    assert result.fit.success
    diagnostics = result.fit.diagnostics
    assert diagnostics is not None
    derivative = np.asarray(diagnostics.metadata["derivative_x"])
    expected = TemperatureEOSFitModel("berman:quadratic", "volume").derivative_x(
        dataset.column("temperature"), result.fit.parameters
    )
    np.testing.assert_allclose(derivative, expected, rtol=2.0e-10)
    sigma = np.sqrt(
        dataset.column("sigma_volume") ** 2
        + (expected * dataset.column("sigma_temperature")) ** 2
    )
    np.testing.assert_allclose(
        diagnostics.metadata["effective_sigma"], sigma, rtol=2.0e-9
    )


@pytest.mark.parametrize(
    "solver",
    [
        OLSOptions(),
        WLSOptions(),
        EffectiveVarianceOptions(max_iterations=25),
        ODROptions(max_iterations=100),
    ],
    ids=lambda item: item.method.value,
)
def test_axial_vt_uses_cubed_length_with_all_solvers(solver) -> None:
    dataset, expected = _synthetic("berman:quadratic", axial=True)
    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="berman:quadratic",
            domain="vt",
            target="a",
            options=EOSFitOptions(solver_options=solver),
        ),
    )

    assert result.fit.success
    assert result.parameter_values["L0"] == pytest.approx(
        np.cbrt(expected["V0"]), rel=1.0e-4
    )
    assert result.derived["length_ref"] == pytest.approx(result.parameter_values["L0"])
    np.testing.assert_allclose(
        result.predictions["linear_expansion_coefficient"],
        result.predictions["auxiliary_expansion_coefficient"] / 3.0,
    )
    assert result.metadata["coefficient_space"] == "auxiliary_cubed_length"


def test_axial_vt_matches_direct_fit_of_cubed_coordinate() -> None:
    axial_dataset, _ = _synthetic("fei:linear", axial=True)
    q = axial_dataset.column("a") ** 3
    sigma_q = 3.0 * axial_dataset.column("a") ** 2 * axial_dataset.column("sigma_a")
    volume_dataset = EOSDataset(
        columns={
            "temperature": axial_dataset.column("temperature"),
            "sigma_temperature": axial_dataset.column("sigma_temperature"),
            "volume": q,
            "sigma_volume": sigma_q,
        }
    )
    options = EOSFitOptions(solver_options=EffectiveVarianceOptions(max_iterations=25))
    axial = EOSFitter().fit(
        axial_dataset,
        EOSFitRequest(model="fei:linear", domain="vt", target="a", options=options),
    )
    volume = EOSFitter().fit(
        volume_dataset,
        EOSFitRequest(model="fei:linear", domain="vt", options=options),
    )

    assert axial.fit.parameters is not None
    assert volume.fit.parameters is not None
    converted = axial.fit.parameters.copy()
    converted[0] = converted[0] ** 3
    np.testing.assert_allclose(converted, volume.fit.parameters, rtol=2e-10)
    assert axial.fit.covariance is not None
    assert volume.fit.covariance is not None
    jacobian = np.eye(axial.fit.parameters.size)
    jacobian[0, 0] = 3.0 * axial.fit.parameters[0] ** 2
    np.testing.assert_allclose(
        jacobian @ axial.fit.covariance @ jacobian.T,
        volume.fit.covariance,
        rtol=1e-6,
    )


def test_real_triclinic_file_runs_berman_quadratic_ols() -> None:
    dataset = read_eos_input(_DATA / "T_triclinic.dat")
    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(model="berman:quadratic", domain="vt"),
    )

    assert result.fit.success
    assert result.parameter_values["V0"] == pytest.approx(669.2920234)
    assert result.parameter_values["alpha0"] == pytest.approx(1.5505181e-5)
    assert result.parameter_values["alpha1"] == pytest.approx(2.2329229e-8)
    assert result.fit.rmse == pytest.approx(0.07233945)


@pytest.mark.parametrize(
    ("solver", "expected_k"),
    [
        (WLSOptions(), 1.0005401182),
        (EffectiveVarianceOptions(max_iterations=30), 1.0005769060),
        (ODROptions(max_iterations=100), 1.0005664413),
    ],
    ids=lambda item: item.method.value if hasattr(item, "method") else str(item),
)
def test_real_rutile_volume_supports_weighted_vt_methods(solver, expected_k) -> None:
    dataset = read_eos_input(_DATA / "rutile.dat")
    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="berman:quadratic",
            domain="vt",
            options=EOSFitOptions(solver_options=solver),
        ),
    )

    assert result.fit.success
    assert result.parameter_values["V0"] == pytest.approx(expected_k, rel=2e-8)
    assert result.metadata["is_isobaric"] is True
    assert result.metadata["reference_pressure"] == pytest.approx(1.0e-4)
    assert result.metadata["input_target_scale"] == "V/V0"


def test_real_rutile_axial_fit_reports_original_length_predictions() -> None:
    dataset = read_eos_input(_DATA / "rutile.dat")
    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="berman:quadratic",
            domain="vt",
            target="a",
            options=EOSFitOptions(
                solver_options=EffectiveVarianceOptions(max_iterations=30)
            ),
        ),
    )

    assert result.fit.success
    assert result.predictions["a"].shape == (dataset.npoints,)
    assert result.parameter_values["L0"] == pytest.approx(np.cbrt(1.0004836931))
    assert result.derived["length_ref"] == pytest.approx(result.parameter_values["L0"])
    assert result.metadata["input_target_scale"] == "L/L0"


def test_vt_rejects_constant_temperature_before_solver() -> None:
    dataset = EOSDataset(
        columns={
            "temperature": np.full(5, 298.15),
            "volume": np.linspace(99.0, 101.0, 5),
        }
    )
    with pytest.raises(ValueError, match="temperature column is constant"):
        EOSFitter().fit(
            dataset,
            EOSFitRequest(model="berman", domain="vt"),
        )
