"""Tests for the first volumetric P-V Quantas EOS API vertical."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.math.fitting import (
    EffectiveVarianceOptions,
    FitMethod,
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
    EOSDataset,
    EOSFitOptions,
    EOSFitRequest,
    EOSFitter,
    ParameterConstraint,
    build_pressure_parameter_map,
    estimate_pressure_parameters,
    read_eos_input,
)


def _synthetic_dataset(model_tag: str) -> tuple[EOSDataset, dict[str, float]]:
    model = next(item for item in available_eos_models() if item.tag == model_tag)
    kp = implied_kp(model) if model.order == 2 else 4.5
    if kp is None:
        kp = 4.5
    kpp = implied_kpp(model, 150.0, kp) if model.order != 4 else -0.02
    parameters = {"K0": 150.0, "KP": kp, "KPP": kpp, "V0": 100.0}
    volume = np.linspace(82.0, 100.0, 24, dtype=np.float64)
    pressure = PressureEOS().pressure(model, parameters, volume)
    pressure = pressure + 5.0e-4 * np.sin(np.arange(volume.size))
    dataset = EOSDataset(
        jobname=f"synthetic-{model_tag}",
        columns={
            "pressure": pressure,
            "sigma_pressure": np.full(volume.shape, 0.01),
            "volume": volume,
            "sigma_volume": np.full(volume.shape, 0.02),
        },
    )
    return dataset, parameters


def test_pressure_estimate_recovers_physical_starting_values() -> None:
    dataset, expected = _synthetic_dataset("BM3")

    estimate = estimate_pressure_parameters(
        dataset.column("volume"), dataset.column("pressure")
    )

    assert estimate["K0"] == pytest.approx(expected["K0"], rel=0.15)
    assert estimate["V0"] == pytest.approx(expected["V0"], rel=0.01)
    assert estimate["KP"] == 4.0


@pytest.mark.parametrize("model", available_eos_models(), ids=lambda item: item.tag)
def test_ols_fits_every_supported_pressure_eos_family_and_order(model) -> None:
    dataset, expected = _synthetic_dataset(model.tag)

    result = EOSFitter().fit(dataset, EOSFitRequest(model=model))

    assert result.fit.success is True
    assert result.fit.method is FitMethod.OLS
    assert result.parameter_values["K0"] == pytest.approx(expected["K0"], rel=1.0e-3)
    assert result.parameter_values["V0"] == pytest.approx(expected["V0"], rel=1.0e-4)
    assert result.fit.rmse is not None and result.fit.rmse < 1.0e-3
    assert result.predictions["pressure"].shape == dataset.column("pressure").shape
    assert result.metadata["residual_definition"] == (
        "observed_pressure-calculated_pressure"
    )


def test_model_order_is_reflected_in_free_and_implied_parameters() -> None:
    dataset, _ = _synthetic_dataset("BM2")
    mapping = build_pressure_parameter_map(
        "BM2", dataset.column("volume"), dataset.column("pressure")
    )

    assert mapping.names == ("K0", "KP", "KPP", "V0")
    assert mapping.free_names == ("K0", "V0")
    assert mapping.states == (
        ParameterState.FREE,
        ParameterState.IMPLIED,
        ParameterState.IMPLIED,
        ParameterState.FREE,
    )
    resolved = mapping.expand(mapping.initial_free_values()).as_mapping()
    assert resolved["KP"] == 4.0
    assert resolved["KPP"] == pytest.approx(
        implied_kpp(
            next(m for m in available_eos_models() if m.tag == "BM2"),
            resolved["K0"],
            4.0,
        )
    )


def test_fixed_parameter_is_not_optimized_and_has_zero_uncertainty() -> None:
    dataset, _ = _synthetic_dataset("BM3")
    request = EOSFitRequest(
        model="BM3",
        constraints=(ParameterConstraint.fixed("KP", 4.5),),
    )

    result = EOSFitter().fit(dataset, request)

    assert result.fit.success is True
    assert result.parameter_values["KP"] == 4.5
    kp_index = result.fit.parameter_names.index("KP")
    assert result.fit.parameter_states[kp_index] is ParameterState.FIXED
    assert result.fit.errors is not None
    assert result.fit.errors[kp_index] == 0.0


def test_wls_uses_pressure_uncertainties_and_absolute_sigma() -> None:
    dataset, _ = _synthetic_dataset("V3")
    request = EOSFitRequest(
        model="V3",
        options=EOSFitOptions(solver_options=WLSOptions(absolute_sigma=True)),
    )

    result = EOSFitter().fit(dataset, request)

    assert result.fit.success is True
    assert result.fit.diagnostics is not None
    assert result.fit.diagnostics.weighted is True
    assert result.fit.diagnostics.standardized_residuals is not None
    assert result.fit.diagnostics.metadata["absolute_sigma"] is True


def test_ols_covariance_is_scaled_from_residual_variance() -> None:
    dataset, _ = _synthetic_dataset("BM3")

    result = EOSFitter().fit(dataset, EOSFitRequest(model="BM3"))

    assert result.fit.success is True
    assert result.fit.diagnostics is not None
    assert result.fit.diagnostics.metadata["absolute_sigma"] is False
    assert result.fit.errors is not None
    assert result.fit.errors[result.fit.parameter_names.index("K0")] < 1.0


def test_mask_is_applied_before_estimation_and_fitting() -> None:
    dataset, _ = _synthetic_dataset("T3")
    mask = np.ones(dataset.npoints, dtype=np.bool_)
    mask[-3:] = False

    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(model="T3", mask=mask),
    )

    assert result.fit.success is True
    assert result.fit.n_points == int(np.count_nonzero(mask))
    assert result.metadata["selected_mask"] == mask.tolist()
    assert result.predictions["pressure"].size == dataset.npoints


def test_current_vertical_accepts_axial_and_odr_targets() -> None:
    dataset, _ = _synthetic_dataset("BM3")
    axis = np.cbrt(dataset.column("volume"))
    dataset.columns["a"] = axis
    dataset.columns["sigma_a"] = np.full(axis.shape, 0.001)

    axial = EOSFitter().fit(dataset, EOSFitRequest(model="BM3", target="a"))
    odr = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="BM3",
            options=EOSFitOptions(solver_options=ODROptions()),
        ),
    )

    assert axial.fit.success
    assert axial.fit.parameter_names == ("M0", "MP", "MPP", "L0")
    assert odr.fit.success
    assert odr.fit.method is FitMethod.ODR


def test_real_quartz_file_runs_ols_and_pressure_wls() -> None:
    dataset = read_eos_input(Path(__file__).with_name("data") / "PV_quartz.dat")

    ols = EOSFitter().fit(dataset, EOSFitRequest(model="BM3"))
    wls = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="BM3",
            options=EOSFitOptions(solver_options=WLSOptions()),
        ),
    )

    assert ols.fit.success is True
    assert ols.parameter_values["K0"] == pytest.approx(37.28542, rel=2.0e-5)
    assert ols.parameter_values["KP"] == pytest.approx(5.93353, rel=2.0e-5)
    assert ols.parameter_values["V0"] == pytest.approx(112.96753, rel=2.0e-6)
    assert ols.fit.rmse is not None and ols.fit.rmse < 0.012

    assert wls.fit.success is True
    assert wls.parameter_values["K0"] == pytest.approx(37.14549, rel=2.0e-5)
    assert wls.parameter_values["KP"] == pytest.approx(5.97666, rel=2.0e-5)
    assert wls.parameter_values["V0"] == pytest.approx(112.98100, rel=2.0e-7)
    assert wls.fit.diagnostics is not None
    assert wls.fit.diagnostics.weighted is True


def test_real_topaz_volume_runs_without_discarding_axial_columns() -> None:
    dataset = read_eos_input(Path(__file__).with_name("data") / "PV_topaz.dat")

    result = EOSFitter().fit(dataset, EOSFitRequest(model="BM3"))

    assert result.fit.success is True
    assert set(("a", "b", "c")).issubset(dataset.columns)
    assert 120.0 < result.parameter_values["K0"] < 250.0
    assert 340.0 < result.parameter_values["V0"] < 350.0


@pytest.mark.parametrize("model", available_eos_models(), ids=lambda item: item.tag)
def test_effective_variance_fits_every_supported_pressure_model(model) -> None:
    dataset, expected = _synthetic_dataset(model.tag)

    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model=model,
            options=EOSFitOptions(
                solver_options=EffectiveVarianceOptions(max_iterations=25),
            ),
        ),
    )

    assert result.fit.success is True
    assert result.parameter_values["K0"] == pytest.approx(expected["K0"], rel=1.0e-3)
    assert result.parameter_values["V0"] == pytest.approx(expected["V0"], rel=1.0e-4)
    assert result.fit.diagnostics is not None
    assert result.fit.diagnostics.n_iterations is not None


def test_effective_variance_uses_pressure_and_volume_uncertainties():
    dataset, expected = _synthetic_dataset("BM3")
    request = EOSFitRequest(
        model="BM3",
        options=EOSFitOptions(
            solver_options=EffectiveVarianceOptions(max_iterations=20),
        ),
    )

    result = EOSFitter().fit(dataset, request)

    assert result.fit.success is True
    assert result.fit.method is FitMethod.EFFECTIVE_VARIANCE
    assert result.parameter_values["K0"] == pytest.approx(expected["K0"], rel=1.0e-3)
    assert result.parameter_values["V0"] == pytest.approx(expected["V0"], rel=1.0e-4)
    assert result.fit.diagnostics is not None
    metadata = result.fit.diagnostics.metadata
    volume = dataset.column("volume")
    physical = result.parameter_values
    bulk = PressureEOS().bulk_modulus("BM3", physical, volume)
    expected_sigma = np.sqrt(
        dataset.column("sigma_pressure") ** 2
        + (bulk / volume * dataset.column("sigma_volume")) ** 2
    )
    np.testing.assert_allclose(metadata["effective_sigma"], expected_sigma, rtol=2.0e-9)
    np.testing.assert_allclose(metadata["derivative_x"], -bulk / volume, rtol=2.0e-12)


def test_real_quartz_file_runs_effective_variance():
    dataset = read_eos_input(Path(__file__).with_name("data") / "PV_quartz.dat")

    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="BM3",
            options=EOSFitOptions(solver_options=EffectiveVarianceOptions()),
        ),
    )

    assert result.fit.success is True
    assert result.parameter_values["K0"] == pytest.approx(37.12593, rel=3.0e-5)
    assert result.parameter_values["KP"] == pytest.approx(5.98826, rel=3.0e-5)
    assert result.parameter_values["V0"] == pytest.approx(112.98088, rel=2.0e-6)
    diagnostics = result.fit.diagnostics
    assert diagnostics is not None
    assert diagnostics.n_iterations is not None
    assert 1 <= diagnostics.n_iterations <= 25
    assert diagnostics.stop_reason == "effective parameters and weights converged"
    assert diagnostics.metadata["reweighting_iterations"] == diagnostics.n_iterations
    assert diagnostics.reduced_chi_square is not None
