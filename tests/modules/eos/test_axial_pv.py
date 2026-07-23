"""Tests for rigorous cubed-length axial EOS fitting."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.math.fitting import (
    EffectiveVarianceOptions,
    FitMethod,
    OLSOptions,
    ODROptions,
    ParameterState,
    WLSOptions,
)
from quantas.core.physics.eos import (
    PressureEOS,
    available_eos_models,
    implied_kp,
    implied_kpp,
)
from quantas.modules.eos import (
    AxialEOSFitModel,
    EOSDataset,
    EOSFitOptions,
    EOSFitRequest,
    EOSFitter,
    ParameterConstraint,
    axial_to_volume_parameters,
    build_axial_parameter_map,
    estimate_axial_parameters,
    read_eos_input,
)

_DATA = Path(__file__).with_name("data")


def _synthetic_axial_dataset(model_tag: str) -> tuple[EOSDataset, dict[str, float]]:
    model = next(item for item in available_eos_models() if item.tag == model_tag)
    auxiliary_kp = implied_kp(model) if model.order == 2 else 4.5
    if auxiliary_kp is None:
        auxiliary_kp = 4.5
    auxiliary_kpp = (
        implied_kpp(model, 150.0, auxiliary_kp) if model.order != 4 else -0.02
    )
    linear = {
        "M0": 450.0,
        "MP": 3.0 * auxiliary_kp,
        "MPP": 3.0 * auxiliary_kpp,
        "L0": 5.0,
    }
    axis = np.linspace(np.cbrt(0.82 * linear["L0"] ** 3), linear["L0"], 24)
    pressure = PressureEOS().pressure(
        model,
        axial_to_volume_parameters(linear),
        axis**3,
    )
    pressure = pressure + 5.0e-4 * np.sin(np.arange(axis.size))
    dataset = EOSDataset(
        jobname=f"synthetic-axial-{model_tag}",
        columns={
            "pressure": pressure,
            "sigma_pressure": np.full(axis.shape, 0.01),
            "a": axis,
            "sigma_a": np.full(axis.shape, 0.002),
        },
        units={"pressure": "GPa", "a": "angstrom"},
    )
    return dataset, linear


def test_axial_estimate_uses_cubed_length_and_returns_linear_parameters() -> None:
    dataset, expected = _synthetic_axial_dataset("BM3")

    estimate = estimate_axial_parameters(
        dataset.column("a"), dataset.column("pressure")
    )

    assert estimate["M0"] == pytest.approx(expected["M0"], rel=0.15)
    assert estimate["L0"] == pytest.approx(expected["L0"], rel=0.01)
    assert estimate["MP"] == 12.0


@pytest.mark.parametrize("model", available_eos_models(), ids=lambda item: item.tag)
def test_ols_fits_every_axial_eos_family_and_order(model) -> None:
    dataset, expected = _synthetic_axial_dataset(model.tag)

    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(model=model, target="a"),
    )

    assert result.fit.success
    assert result.fit.method is FitMethod.OLS
    assert result.fit.parameter_names == ("M0", "MP", "MPP", "L0")
    assert result.parameter_values["M0"] == pytest.approx(expected["M0"], rel=1e-3)
    assert result.parameter_values["L0"] == pytest.approx(expected["L0"], rel=1e-4)
    assert result.metadata["coordinate_transform"] == "q=x^3"
    assert result.metadata["relationship"] == "pressure(a)"


def test_second_order_axial_parameters_use_linear_implied_values() -> None:
    dataset, _ = _synthetic_axial_dataset("BM2")
    mapping = build_axial_parameter_map(
        "BM2", dataset.column("a"), dataset.column("pressure")
    )

    assert mapping.names == ("M0", "MP", "MPP", "L0")
    assert mapping.free_names == ("M0", "L0")
    assert mapping.states == (
        ParameterState.FREE,
        ParameterState.IMPLIED,
        ParameterState.IMPLIED,
        ParameterState.FREE,
    )
    resolved = mapping.expand(mapping.initial_free_values()).as_mapping()
    assert resolved["MP"] == 12.0
    auxiliary = axial_to_volume_parameters(resolved)
    assert auxiliary["KP"] == 4.0
    assert resolved["MPP"] == pytest.approx(3.0 * auxiliary["KPP"])


@pytest.mark.parametrize(
    "solver_options",
    (OLSOptions(), WLSOptions(), EffectiveVarianceOptions(), ODROptions()),
    ids=lambda item: item.method.value,
)
def test_axial_and_auxiliary_volume_fits_are_equivalent(solver_options) -> None:
    axial_dataset, _ = _synthetic_axial_dataset("BM3")
    axis = axial_dataset.column("a")
    sigma_axis = axial_dataset.column("sigma_a")
    volume_dataset = EOSDataset(
        columns={
            "pressure": axial_dataset.column("pressure"),
            "sigma_pressure": axial_dataset.column("sigma_pressure"),
            "volume": axis**3,
            "sigma_volume": 3.0 * axis**2 * sigma_axis,
        }
    )
    options = EOSFitOptions(solver_options=solver_options)

    axial = EOSFitter().fit(
        axial_dataset,
        EOSFitRequest(model="BM3", target="a", options=options),
    )
    volume = EOSFitter().fit(
        volume_dataset,
        EOSFitRequest(model="BM3", options=options),
    )

    assert axial.fit.success and volume.fit.success
    converted = axial_to_volume_parameters(axial.parameter_values)
    for name in ("K0", "KP", "KPP", "V0"):
        assert converted[name] == pytest.approx(
            volume.parameter_values[name], rel=2.0e-8, abs=2.0e-8
        )


def test_fixed_linear_parameter_is_not_optimized() -> None:
    dataset, _ = _synthetic_axial_dataset("BM3")
    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="BM3",
            target="a",
            constraints=(ParameterConstraint.fixed("MP", 13.5),),
        ),
    )

    assert result.fit.success
    assert result.parameter_values["MP"] == 13.5
    index = result.fit.parameter_names.index("MP")
    assert result.fit.parameter_states[index] is ParameterState.FIXED
    assert result.fit.errors is not None
    assert result.fit.errors[index] == 0.0


def test_effective_variance_reports_length_space_derivatives() -> None:
    dataset, _ = _synthetic_axial_dataset("BM3")
    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="BM3",
            target="a",
            options=EOSFitOptions(solver_options=EffectiveVarianceOptions()),
        ),
    )

    assert result.fit.success
    assert result.fit.diagnostics is not None
    metadata = result.fit.diagnostics.metadata
    derivative = np.asarray(metadata["derivative_length"])
    modulus = np.asarray(metadata["linear_modulus"])
    axis = dataset.column("a")
    np.testing.assert_allclose(derivative, -modulus / axis, rtol=2.0e-12)
    projected = derivative * dataset.column("sigma_a")
    np.testing.assert_allclose(
        np.abs(projected),
        np.abs(np.asarray(metadata["projected_sigma_x"])),
        rtol=2.0e-12,
    )


def test_odr_retains_cubed_coordinate_and_reports_length_corrections() -> None:
    dataset, _ = _synthetic_axial_dataset("BM3")
    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="BM3",
            target="a",
            options=EOSFitOptions(solver_options=ODROptions()),
        ),
    )

    assert result.fit.success
    assert result.fit.diagnostics is not None
    diagnostics = result.fit.diagnostics
    assert diagnostics.x_corrections is not None
    metadata = diagnostics.metadata
    assert metadata["solver_coordinate"] == "a^3"
    adjusted = np.asarray(metadata["adjusted_length"])
    correction = np.asarray(metadata["length_corrections"])
    np.testing.assert_allclose(adjusted - dataset.column("a"), correction)
    assert result.fit.metadata["coordinate_transform"] == "q=x^3"


@pytest.mark.parametrize("target", ("a", "b", "c"))
@pytest.mark.parametrize(
    "solver_options",
    (OLSOptions(), WLSOptions(), EffectiveVarianceOptions(), ODROptions()),
    ids=lambda item: item.method.value,
)
def test_real_topaz_axes_fit_with_every_solver(target, solver_options) -> None:
    dataset = read_eos_input(_DATA / "PV_topaz.dat")

    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="BM3",
            target=target,
            options=EOSFitOptions(solver_options=solver_options),
        ),
    )

    assert result.fit.success
    assert result.parameter_values["M0"] > 0.0
    assert result.parameter_values["L0"] == pytest.approx(
        dataset.column(target)[0], rel=0.02
    )
    assert result.predictions["pressure"].shape == dataset.column("pressure").shape


def test_axial_model_evaluates_linear_parameters_through_shared_pressure_core() -> None:
    model = AxialEOSFitModel("BM3", "a")
    parameters = np.asarray([450.0, 13.5, -0.06, 5.0])
    axis = np.asarray([5.0, 4.9, 4.8])

    calculated = model.evaluate(axis**3, parameters)
    expected = PressureEOS().pressure(
        "BM3", axial_to_volume_parameters(parameters), axis**3
    )

    np.testing.assert_allclose(calculated, expected, rtol=0.0, atol=0.0)
