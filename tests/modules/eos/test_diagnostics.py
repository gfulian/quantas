from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.math.fitting import EffectiveVarianceOptions, OLSOptions
from quantas.modules.eos import (
    EOSArchive,
    EOSCalculator,
    EOSDiagnostics,
    EOSFitOptions,
    EOSFitRequest,
    EOSFitter,
    read_eos_input,
)

DATA = Path(__file__).with_name("data")


def _quartz_record(tmp_path: Path):
    dataset = read_eos_input(DATA / "PV_quartz.dat")
    request = EOSFitRequest(
        model="BM3",
        domain="pv",
        target="volume",
        options=EOSFitOptions(
            solver_options=EffectiveVarianceOptions(max_iterations=50)
        ),
    )
    result = EOSFitter().fit(dataset, request)
    archive_path = tmp_path / "quartz.hdf5"
    with EOSArchive.create(archive_path, dataset=dataset) as archive:
        record = archive.store_fit(1, request, result)
    return archive_path, dataset, record


def test_pressure_diagnostics_include_residuals_and_finite_strain(tmp_path):
    archive_path, dataset, record = _quartz_record(tmp_path)

    diagnostic = EOSDiagnostics.from_archive(archive_path).build()

    assert diagnostic.record_id == record.record_id
    assert diagnostic.nrows == dataset.npoints
    assert "finite_strain" in diagnostic.columns
    assert "normalized_pressure" in diagnostic.columns
    assert diagnostic.metadata["normalized_pressure"]["strain_family"] == "eulerian"
    np.testing.assert_allclose(
        diagnostic.columns["residual"],
        diagnostic.columns["observed_pressure"]
        - diagnostic.columns["calculated_pressure"],
    )
    assert (
        np.count_nonzero(np.isfinite(diagnostic.columns["standardized_residual"]))
        == dataset.npoints
    )


def test_pv_calculator_evaluates_pressure_grid_and_covariance(tmp_path):
    archive_path, _, record = _quartz_record(tmp_path)
    calculator = EOSCalculator.from_archive(archive_path)

    calculation = calculator.calculate(pressure=np.array([0.0, 2.0, 5.0]))

    assert calculation.record_id == record.record_id
    assert calculation.nrows == 3
    np.testing.assert_allclose(calculation.columns["pressure"], [0.0, 2.0, 5.0])
    assert np.all(calculation.columns["volume"] > 0.0)
    assert np.all(calculation.columns["bulk_modulus"] > 0.0)
    assert "volume" in calculation.uncertainties
    assert "pressure" not in calculation.uncertainties
    assert np.all(calculation.uncertainties["volume"] >= 0.0)


def test_pv_calculator_round_trip_volume_and_pressure(tmp_path):
    archive_path, _, _ = _quartz_record(tmp_path)
    calculator = EOSCalculator.from_archive(archive_path)

    from_pressure = calculator.calculate(
        pressure=np.array([1.0, 4.0]), propagate_uncertainty=False
    )
    from_volume = calculator.calculate(
        volume=from_pressure.columns["volume"], propagate_uncertainty=False
    )

    np.testing.assert_allclose(from_volume.columns["pressure"], [1.0, 4.0], rtol=1e-9)

    propagated = calculator.calculate(volume=from_pressure.columns["volume"])
    assert "volume" not in propagated.uncertainties
    assert "pressure" in propagated.uncertainties


def test_calculator_requires_slot_for_multiple_results(tmp_path):
    dataset = read_eos_input(DATA / "PV_topaz.dat")
    archive_path = tmp_path / "topaz.hdf5"
    with EOSArchive.create(archive_path, dataset=dataset) as archive:
        for target in ("a", "b"):
            request = EOSFitRequest(
                model="BM3",
                domain="pv",
                target=target,
                options=EOSFitOptions(solver_options=OLSOptions()),
            )
            result = EOSFitter().fit(dataset, request)
            archive.store_fit(1, request, result)

    with pytest.raises(ValueError, match="exactly one accepted"):
        EOSCalculator.from_archive(archive_path)
    calculator = EOSCalculator.from_archive(archive_path, slot="pv/a")
    result = calculator.calculate(pressure=[0.0, 1.0], propagate_uncertainty=False)
    assert "a" in result.columns
    assert "linear_modulus" in result.columns


def test_unsupported_normalized_pressure_still_returns_residual_diagnostics(tmp_path):
    dataset = read_eos_input(DATA / "PV_quartz.dat")
    request = EOSFitRequest(
        model="M",
        domain="pv",
        target="volume",
        options=EOSFitOptions(solver_options=OLSOptions()),
    )
    fit = EOSFitter().fit(dataset, request)
    path = tmp_path / "murnaghan.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        archive.store_fit(1, request, fit)

    diagnostic = EOSDiagnostics.from_archive(path).build()

    assert "residual" in diagnostic.columns
    assert "normalized_pressure" not in diagnostic.columns
    assert diagnostic.metadata["normalized_pressure"]["available"] is False
    assert any("not defined" in warning for warning in diagnostic.warnings)


def test_vt_calculator_returns_value_alpha_and_derivative(tmp_path):
    from quantas.core.physics.eos import TemperatureEOS
    from quantas.modules.eos import EOSDataset, ParameterConstraint

    temperature = np.linspace(200.0, 800.0, 30)
    parameters = {
        "V0": 100.0,
        "temperature_ref": 300.0,
        "alpha0": 3.0e-5,
        "alpha1": 1.0e-8,
    }
    volume = TemperatureEOS().value("berman:quadratic", parameters, temperature)
    dataset = EOSDataset(
        columns={"temperature": temperature, "volume": volume},
        units={"temperature": "K", "volume": "angstrom^3"},
    )
    request = EOSFitRequest(
        model="berman:quadratic",
        domain="vt",
        constraints=(ParameterConstraint.fixed("temperature_ref", 300.0),),
    )
    fit = EOSFitter().fit(dataset, request)
    path = tmp_path / "vt.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        archive.store_fit(1, request, fit)

    result = EOSCalculator.from_archive(path).calculate(
        temperature=[300.0, 600.0], propagate_uncertainty=False
    )

    assert result.columns["volume"][0] == pytest.approx(100.0, rel=1e-8)
    assert np.all(result.columns["expansion_coefficient"] > 0.0)
    propagated = EOSCalculator.from_archive(path).calculate(temperature=[300.0, 600.0])
    assert "temperature" not in propagated.uncertainties
    assert "volume" in propagated.uncertainties
    np.testing.assert_allclose(
        result.columns["temperature_derivative"],
        result.columns["expansion_coefficient"] * result.columns["volume"],
    )


def test_pvt_calculator_round_trip_and_derived_properties(tmp_path):
    from quantas.core.physics.eos import PVTEOS, PVTModel
    from quantas.modules.eos import EOSDataset, ParameterConstraint

    model = PVTModel("BM3", "linear", "berman:quadratic")
    pressure_parameters = {"K0": 160.0, "KP": 4.2, "KPP": -0.02, "V0": 100.0}
    thermal_parameters = {
        "V0": 100.0,
        "temperature_ref": 300.0,
        "alpha0": 3.0e-5,
        "alpha1": 1.0e-8,
    }
    coupling_parameters = {"dK0_dT": -0.02}
    temperature = np.repeat([300.0, 600.0, 900.0], 6)
    pressure = np.tile(np.linspace(0.0, 10.0, 6), 3)
    volume = PVTEOS().volume(
        model,
        pressure_parameters,
        thermal_parameters,
        coupling_parameters,
        pressure,
        temperature,
    )
    dataset = EOSDataset(
        columns={"pressure": pressure, "temperature": temperature, "volume": volume},
        units={"pressure": "GPa", "temperature": "K", "volume": "angstrom^3"},
    )
    constraints = (
        ParameterConstraint.free("K0", 158.0),
        ParameterConstraint.free("KP", 4.0),
        ParameterConstraint.free("V0", 99.8),
        ParameterConstraint.fixed("temperature_ref", 300.0),
        ParameterConstraint.free("alpha0", 2.8e-5),
        ParameterConstraint.free("alpha1", 0.8e-8),
        ParameterConstraint.free("dK0_dT", -0.018),
    )
    request = EOSFitRequest(
        model=model,
        domain="pvt",
        constraints=constraints,
        options=EOSFitOptions(solver_options=OLSOptions(max_iterations=5000)),
    )
    fit = EOSFitter().fit(dataset, request)
    path = tmp_path / "pvt.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        archive.store_fit(1, request, fit)

    calculator = EOSCalculator.from_archive(path)
    calculated = calculator.calculate(
        pressure=[0.0, 5.0],
        temperature=[300.0, 600.0],
        propagate_uncertainty=False,
    )
    inverse = calculator.calculate(
        volume=calculated.columns["volume"],
        temperature=[300.0, 600.0],
        propagate_uncertainty=False,
    )

    np.testing.assert_allclose(inverse.columns["pressure"], [0.0, 5.0], atol=1e-8)
    assert np.all(calculated.columns["bulk_modulus"] > 0.0)
    assert np.all(np.isfinite(calculated.columns["bulk_modulus_derivative"]))
    assert np.all(np.isfinite(calculated.columns["bulk_modulus_second_derivative"]))
    assert np.all(calculated.columns["expansion_coefficient"] > 0.0)
    assert "thermal_pressure" not in calculated.columns
    propagated = calculator.calculate(
        pressure=[0.0, 5.0],
        temperature=[300.0, 600.0],
    )
    assert "pressure" not in propagated.uncertainties
    assert "temperature" not in propagated.uncertainties
    assert "volume" in propagated.uncertainties
    np.testing.assert_allclose(
        calculated.columns["zero_pressure_dK0_dT"], -0.02, rtol=2e-5
    )


def test_diagnostics_preserve_excluded_rows_and_groups(tmp_path):
    from quantas.modules.eos import EOSDataset

    volume = np.linspace(100.0, 92.0, 10)
    pressure = np.linspace(0.0, 8.0, 10)
    default_mask = np.ones(10, dtype=np.bool_)
    default_mask[[2, 7]] = False
    groups = np.repeat([1, 2], 5)
    dataset = EOSDataset(
        columns={"pressure": pressure, "volume": volume},
        units={"pressure": "GPa", "volume": "angstrom^3"},
        default_mask=default_mask,
        groups=groups,
    )
    request = EOSFitRequest(
        model="BM3",
        domain="pv",
        options=EOSFitOptions(solver_options=OLSOptions()),
    )
    fit = EOSFitter().fit(dataset, request)
    request = fit.request
    path = tmp_path / "selected.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        archive.store_fit(1, request, fit)

    diagnostic = EOSDiagnostics.from_archive(path).build()

    np.testing.assert_array_equal(
        diagnostic.columns["included"].astype(np.bool_),
        default_mask,
    )
    np.testing.assert_array_equal(diagnostic.columns["group"], groups)
    assert diagnostic.nrows == dataset.npoints
    assert np.all(np.isnan(diagnostic.columns["standardized_residual"][~default_mask]))
    assert np.all(np.isfinite(diagnostic.columns["residual"]))
    assert diagnostic.metadata["selected_observations"] == int(
        np.count_nonzero(default_mask)
    )


def test_axial_diagnostics_use_cubed_length_convention(tmp_path):
    dataset = read_eos_input(DATA / "PV_topaz.dat")
    request = EOSFitRequest(
        model="BM3",
        domain="pv",
        target="a",
        options=EOSFitOptions(solver_options=OLSOptions()),
    )
    fit = EOSFitter().fit(dataset, request)
    path = tmp_path / "axial.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        archive.store_fit(1, request, fit)

    diagnostic = EOSDiagnostics.from_archive(path).build()
    calculation = EOSCalculator.from_archive(path).calculate(
        pressure=[0.0],
        propagate_uncertainty=False,
    )

    assert diagnostic.metadata["normalized_pressure"]["linear_eos"] is True
    assert any("M/3" in warning for warning in diagnostic.warnings)
    assert calculation.columns["linear_modulus"][0] == pytest.approx(
        fit.parameter_values["M0"],
        rel=1.0e-9,
    )


def test_mgd_calculator_reports_thermal_pressure_only_for_thermal_coupling(tmp_path):
    from quantas.core.physics.eos import MGDNormalization, PVTEOS, PVTModel
    from quantas.modules.eos import EOSDataset, ParameterConstraint

    model = PVTModel(
        "BM3",
        "thermal-pressure",
        thermal_pressure_model="mgd",
        mgd_normalization=MGDNormalization.cell(atoms_per_cell=8),
    )
    pressure_parameters = {"K0": 48.0, "KP": 4.5, "KPP": -0.02, "V0": 150.0}
    coupling_parameters = {
        "temperature_ref": 295.0,
        "theta_d0": 459.0,
        "gamma0": 1.547,
        "q": 0.94,
    }
    temperature = np.repeat([295.0, 600.0, 900.0], 5)
    pressure = np.tile(np.linspace(0.0, 8.0, 5), 3)
    volume = PVTEOS().volume(
        model,
        pressure_parameters,
        None,
        coupling_parameters,
        pressure,
        temperature,
    )
    dataset = EOSDataset(
        columns={"pressure": pressure, "temperature": temperature, "volume": volume},
        units={"pressure": "GPa", "temperature": "K", "volume": "angstrom^3"},
    )
    request = EOSFitRequest(
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
        options=EOSFitOptions(solver_options=OLSOptions(max_iterations=5000)),
    )
    fit = EOSFitter().fit(dataset, request)
    path = tmp_path / "mgd.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        archive.store_fit(1, request, fit)

    calculation = EOSCalculator.from_archive(path).calculate(
        pressure=[0.0, 5.0],
        temperature=[295.0, 900.0],
        propagate_uncertainty=False,
    )

    assert "thermal_pressure" in calculation.columns
    assert calculation.columns["thermal_pressure"][0] == pytest.approx(0.0, abs=1e-12)
    assert calculation.columns["thermal_pressure"][1] > 0.0
