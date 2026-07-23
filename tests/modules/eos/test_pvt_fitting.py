"""Global and controlled fitting tests for P-V-T EOS models."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.math.fitting import (
    EffectiveVarianceOptions,
    ODROptions,
    OLSOptions,
    ParameterState,
    WLSOptions,
)
from quantas.core.physics.eos import PVTEOS, PVTModel
from quantas.modules.eos import (
    EOSDataset,
    EOSFitDomain,
    EOSFitOptions,
    EOSFitRequest,
    EOSFitter,
    ParameterConstraint,
    build_pvt_parameter_map,
    pvt_parameter_names,
)


PRESSURE = {"K0": 160.0, "KP": 4.2, "KPP": -0.02, "V0": 100.0}
THERMAL = {
    "V0": 100.0,
    "temperature_ref": 300.0,
    "alpha0": 3.0e-5,
    "alpha1": 1.0e-8,
}
LINEAR = {"dK0_dT": -0.02}


def _dataset(
    model: PVTModel,
    thermal: dict[str, float] | None,
    coupling: dict[str, float],
    *,
    with_sigma: bool = True,
) -> EOSDataset:
    temperature = np.repeat(np.asarray([300.0, 450.0, 600.0, 750.0, 900.0]), 8)
    pressure = np.tile(np.linspace(0.0, 14.0, 8), 5)
    volume = PVTEOS().volume(model, PRESSURE, thermal, coupling, pressure, temperature)
    columns = {
        "pressure": pressure,
        "volume": volume,
        "temperature": temperature,
    }
    if with_sigma:
        columns.update(
            {
                "sigma_pressure": np.full(pressure.shape, 0.02),
                "sigma_volume": np.full(pressure.shape, 0.002),
                "sigma_temperature": np.full(pressure.shape, 0.2),
            }
        )
    return EOSDataset(
        jobname="synthetic PVT",
        columns=columns,
        units={"pressure": "GPa", "volume": "angstrom^3", "temperature": "K"},
    )


def _initial_constraints() -> tuple[ParameterConstraint, ...]:
    return (
        ParameterConstraint.free("K0", 158.0),
        ParameterConstraint.free("KP", 4.0),
        ParameterConstraint.free("V0", 99.8),
        ParameterConstraint.fixed("temperature_ref", 300.0),
        ParameterConstraint.free("alpha0", 2.8e-5),
        ParameterConstraint.free("alpha1", 0.8e-8),
        ParameterConstraint.free("dK0_dT", -0.018),
    )


@pytest.mark.parametrize(
    "solver",
    [
        OLSOptions(max_iterations=5000),
        WLSOptions(max_iterations=5000),
        EffectiveVarianceOptions(max_iterations=30, inner_max_iterations=5000),
        ODROptions(max_iterations=100),
    ],
    ids=lambda item: type(item).__name__,
)
def test_global_linear_pvt_fit_recovers_synthetic_parameters(solver) -> None:
    model = PVTModel("BM3", "linear", "berman:quadratic")
    result = EOSFitter().fit(
        _dataset(model, THERMAL, LINEAR),
        EOSFitRequest(
            model=model,
            domain=EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE,
            constraints=_initial_constraints(),
            options=EOSFitOptions(solver_options=solver),
        ),
    )
    assert result.fit.success, result.fit.message
    assert result.parameter_values["K0"] == pytest.approx(160.0, rel=2.0e-7)
    assert result.parameter_values["KP"] == pytest.approx(4.2, rel=2.0e-7)
    assert result.parameter_values["V0"] == pytest.approx(100.0, rel=2.0e-8)
    assert result.parameter_values["alpha0"] == pytest.approx(3.0e-5, rel=2.0e-6)
    assert result.parameter_values["alpha1"] == pytest.approx(1.0e-8, rel=3.0e-5)
    assert result.parameter_values["dK0_dT"] == pytest.approx(-0.02, rel=2.0e-6)
    assert result.fit.parameter_states[2] is ParameterState.IMPLIED
    assert result.predictions["pressure"].shape == (40,)
    assert result.predictions["bulk_modulus"].shape == (40,)
    assert result.predictions["reference_volume"].shape == (40,)
    assert result.predictions["zero_pressure_bulk_modulus"].shape == (40,)
    assert result.predictions["zero_pressure_expansion_coefficient"].shape == (40,)
    assert result.predictions["zero_pressure_dK0_dT"].shape == (40,)
    assert result.metadata["fit_mode"] == "global"
    assert result.metadata["controlled_fit"] is False


def test_anderson_gruneisen_global_fit() -> None:
    model = PVTModel("BM3", "ag", "berman:quadratic")
    coupling = {"delta": 4.4}
    constraints = tuple(
        constraint
        for constraint in _initial_constraints()
        if constraint.name != "dK0_dT"
    ) + (ParameterConstraint.free("delta", 4.0),)
    result = EOSFitter().fit(
        _dataset(model, THERMAL, coupling),
        EOSFitRequest(
            model=model,
            domain="pvt",
            constraints=constraints,
            options=EOSFitOptions(solver_options=EffectiveVarianceOptions()),
        ),
    )
    assert result.fit.success
    assert result.parameter_values["delta"] == pytest.approx(4.4, rel=3.0e-5)
    assert result.derived["dK0_dT_at_reference"] < 0.0


def test_thermal_pressure_global_fit() -> None:
    model = PVTModel("BM3", "thermal-pressure")
    coupling = {"temperature_ref": 300.0, "alpha_ref": 3.0e-5, "theta_e": 500.0}
    dataset = _dataset(model, None, coupling)
    constraints = (
        ParameterConstraint.free("K0", 158.0),
        ParameterConstraint.free("KP", 4.0),
        ParameterConstraint.free("V0", 99.8),
        ParameterConstraint.fixed("temperature_ref", 300.0),
        ParameterConstraint.free("alpha_ref", 2.8e-5),
        ParameterConstraint.fixed("theta_e", 500.0),
    )
    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model=model,
            domain="pvt",
            constraints=constraints,
            options=EOSFitOptions(solver_options=WLSOptions(max_iterations=5000)),
        ),
    )
    assert result.fit.success
    assert result.parameter_values["K0"] == pytest.approx(160.0, rel=1.0e-7)
    assert result.parameter_values["alpha_ref"] == pytest.approx(3.0e-5, rel=1.0e-7)
    assert result.derived["zero_pressure_expansion_at_reference"] == pytest.approx(
        3.0e-5, rel=1.0e-7
    )
    assert result.predictions["thermal_pressure"].shape == (40,)


def test_controlled_fit_can_fix_pv_and_vt_parameters() -> None:
    model = PVTModel("BM3", "linear", "berman:quadratic")
    fixed = (
        ParameterConstraint.fixed("K0", 160.0),
        ParameterConstraint.fixed("KP", 4.2),
        ParameterConstraint.fixed("V0", 100.0),
        ParameterConstraint.fixed("temperature_ref", 300.0),
        ParameterConstraint.fixed("alpha0", 3.0e-5),
        ParameterConstraint.fixed("alpha1", 1.0e-8),
        ParameterConstraint.free("dK0_dT", -0.015),
    )
    result = EOSFitter().fit(
        _dataset(model, THERMAL, LINEAR),
        EOSFitRequest(
            model=model,
            domain="pvt",
            constraints=fixed,
            options=EOSFitOptions(solver_options=OLSOptions()),
        ),
    )
    assert result.fit.success
    assert result.parameter_values["dK0_dT"] == pytest.approx(-0.02, rel=2.0e-8)
    assert result.metadata["controlled_fit"] is True
    assert result.metadata["fit_mode"] == "controlled"
    assert set(result.metadata["controlled_parameters"]) == {
        "K0",
        "KP",
        "V0",
        "alpha0",
        "alpha1",
    }
    assert result.fit.n_parameters == 1


def test_pvt_effective_variance_uses_volume_and_temperature_terms() -> None:
    model = PVTModel("BM3", "linear", "berman:quadratic")
    result = EOSFitter().fit(
        _dataset(model, THERMAL, LINEAR),
        EOSFitRequest(
            model=model,
            domain="pvt",
            constraints=_initial_constraints(),
            options=EOSFitOptions(solver_options=EffectiveVarianceOptions()),
        ),
    )
    assert result.fit.success
    diagnostics = result.fit.diagnostics
    assert diagnostics is not None
    derivative = np.asarray(diagnostics.metadata["derivative_x"])
    projected = np.asarray(diagnostics.metadata["projected_sigma_x"])
    assert derivative.shape == (2, 40)
    assert projected.shape == (2, 40)
    variance = np.asarray(diagnostics.metadata["variance_x_component"])
    assert variance.shape == (2, 40)
    assert np.all(np.sum(variance, axis=0) > 0.0)


def test_parameter_map_shares_reference_volume_and_kp_with_khp() -> None:
    model = PVTModel("BM3", "linear", "khp")
    dataset = _dataset(PVTModel("BM3", "linear", "berman:quadratic"), THERMAL, LINEAR)
    mapping = build_pvt_parameter_map(
        model,
        dataset.column("volume"),
        dataset.column("temperature"),
        dataset.column("pressure"),
    )
    assert mapping.names.count("V0") == 1
    assert mapping.names.count("KP") == 1
    assert "kp" not in mapping.names


def test_pvt_request_and_data_validation() -> None:
    with pytest.raises(ValueError, match="PVTModel"):
        EOSFitRequest(model="BM3", domain="pvt")
    model = PVTModel("BM3", "linear", "berman:quadratic")
    dataset = _dataset(model, THERMAL, LINEAR)
    constant = EOSDataset(
        columns={
            **dataset.columns,
            "temperature": np.full(dataset.npoints, 300.0),
        }
    )
    with pytest.raises(ValueError, match="temperature column is constant"):
        EOSFitter().fit(constant, EOSFitRequest(model=model, domain="pvt"))


def test_pvt_odr_requires_positive_sigma_for_both_coordinates() -> None:
    model = PVTModel("BM3", "linear", "berman:quadratic")
    dataset = _dataset(model, THERMAL, LINEAR, with_sigma=False)
    result = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model=model,
            domain="pvt",
            constraints=_initial_constraints(),
            options=EOSFitOptions(solver_options=ODROptions()),
        ),
    )
    assert not result.fit.success
    assert "sigma_x" in result.fit.message


def test_pvt_parameter_names_are_stable() -> None:
    assert pvt_parameter_names(PVTModel("BM3", "linear", "berman:linear")) == (
        "K0",
        "KP",
        "KPP",
        "V0",
        "temperature_ref",
        "alpha0",
        "alpha1",
        "dK0_dT",
    )


def _mgd_dataset(model: PVTModel) -> EOSDataset:
    """Return a broad synthetic NaF-like P-V-T table for MGD fitting."""
    temperature = np.repeat(
        np.asarray([100.0, 200.0, 295.0, 450.0, 700.0, 1000.0, 1200.0]),
        9,
    )
    pressure = np.tile(np.linspace(0.0, 24.0, 9), 7)
    pressure_parameters = {"K0": 48.0, "KP": 4.5, "KPP": -0.1, "V0": 150.0}
    coupling_parameters = {
        "temperature_ref": 295.0,
        "theta_d0": 459.0,
        "gamma0": 1.547,
        "q": 0.94,
    }
    if model.thermal_pressure_spec.mgd_variant.value == "q-compromise":
        coupling_parameters.pop("q")
    volume = PVTEOS().volume(
        model,
        pressure_parameters,
        None,
        coupling_parameters,
        pressure,
        temperature,
    )
    return EOSDataset(
        jobname="synthetic MGD PVT",
        columns={
            "pressure": pressure,
            "volume": volume,
            "temperature": temperature,
            "sigma_pressure": np.full(pressure.shape, 0.02),
            "sigma_volume": np.full(pressure.shape, 0.002),
            "sigma_temperature": np.full(pressure.shape, 0.2),
        },
        units={"pressure": "GPa", "volume": "angstrom^3", "temperature": "K"},
    )


def test_full_mgd_global_fit_recovers_synthetic_parameters() -> None:
    from quantas.core.physics.eos import MGDNormalization

    model = PVTModel(
        "BM3",
        "thermal-pressure",
        thermal_pressure_model="mgd",
        mgd_normalization=MGDNormalization.cell(atoms_per_cell=8),
    )
    constraints = (
        ParameterConstraint.free("K0", 47.0),
        ParameterConstraint.free("KP", 4.3),
        ParameterConstraint.free("V0", 149.5),
        ParameterConstraint.fixed("temperature_ref", 295.0),
        ParameterConstraint.free("theta_d0", 440.0, lower_bound=100.0),
        ParameterConstraint.free("gamma0", 1.4),
        ParameterConstraint.free("q", 0.8),
    )
    result = EOSFitter().fit(
        _mgd_dataset(model),
        EOSFitRequest(
            model=model,
            domain="pvt",
            constraints=constraints,
            options=EOSFitOptions(solver_options=OLSOptions(max_iterations=10000)),
        ),
    )
    assert result.fit.success, result.fit.message
    assert result.parameter_values["K0"] == pytest.approx(48.0, rel=2.0e-6)
    assert result.parameter_values["KP"] == pytest.approx(4.5, rel=2.0e-6)
    assert result.parameter_values["V0"] == pytest.approx(150.0, rel=2.0e-7)
    assert result.parameter_values["theta_d0"] == pytest.approx(459.0, rel=2.0e-5)
    assert result.parameter_values["gamma0"] == pytest.approx(1.547, rel=2.0e-6)
    assert result.parameter_values["q"] == pytest.approx(0.94, rel=2.0e-5)
    assert result.predictions["thermal_pressure"].shape == (63,)
    assert result.metadata["pvt_model"]["thermal_pressure_model"]["variant"] == "full"


def test_q_compromise_mgd_fit_has_no_q_parameter() -> None:
    from quantas.core.physics.eos import MGDNormalization

    model = PVTModel(
        "BM3",
        "thermal-pressure",
        thermal_pressure_model="mgd:q-compromise",
        mgd_normalization=MGDNormalization.cell(atoms_per_cell=8),
    )
    constraints = (
        ParameterConstraint.free("K0", 47.0),
        ParameterConstraint.free("KP", 4.3),
        ParameterConstraint.free("V0", 149.5),
        ParameterConstraint.fixed("temperature_ref", 295.0),
        ParameterConstraint.free("theta_d0", 440.0),
        ParameterConstraint.free("gamma0", 1.4),
    )
    result = EOSFitter().fit(
        _mgd_dataset(model),
        EOSFitRequest(
            model=model,
            domain="pvt",
            constraints=constraints,
            options=EOSFitOptions(solver_options=OLSOptions(max_iterations=10000)),
        ),
    )
    assert result.fit.success, result.fit.message
    assert "q" not in result.fit.parameter_names
    assert result.parameter_values["theta_d0"] == pytest.approx(459.0, rel=2.0e-5)
    assert result.parameter_values["gamma0"] == pytest.approx(1.547, rel=2.0e-6)


def test_mgd_effective_variance_projects_volume_and_temperature_errors() -> None:
    """Use analytical MGD coordinate derivatives in effective variance."""
    from quantas.core.physics.eos import MGDNormalization

    model = PVTModel(
        "BM3",
        "thermal-pressure",
        thermal_pressure_model="mgd",
        mgd_normalization=MGDNormalization.cell(atoms_per_cell=8),
    )
    result = EOSFitter().fit(
        _mgd_dataset(model),
        EOSFitRequest(
            model=model,
            domain="pvt",
            constraints=(
                ParameterConstraint.free("K0", 47.0),
                ParameterConstraint.free("KP", 4.3),
                ParameterConstraint.free("V0", 149.5),
                ParameterConstraint.fixed("temperature_ref", 295.0),
                ParameterConstraint.fixed("theta_d0", 459.0),
                ParameterConstraint.free("gamma0", 1.4),
                ParameterConstraint.fixed("q", 0.94),
            ),
            options=EOSFitOptions(
                solver_options=EffectiveVarianceOptions(
                    max_iterations=30,
                    inner_max_iterations=5000,
                )
            ),
        ),
    )
    assert result.fit.success, result.fit.message
    diagnostics = result.fit.diagnostics
    assert diagnostics is not None
    derivative = np.asarray(diagnostics.metadata["derivative_x"])
    variance = np.asarray(diagnostics.metadata["variance_x_component"])
    assert derivative.shape == (2, 63)
    assert variance.shape == (2, 63)
    assert np.all(np.isfinite(derivative))
    assert np.all(np.sum(variance, axis=0) > 0.0)
    assert result.parameter_values["gamma0"] == pytest.approx(1.547, rel=3e-6)
