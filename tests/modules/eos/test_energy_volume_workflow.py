"""Public Energy-EOS workflow characterization tests."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner
import numpy as np
import pytest

from quantas.api import eos as public_eos
from quantas.cli.main import main
from quantas.core.math.fitting import OLSOptions, WLSOptions
from quantas.core.physics.eos import EnergyEOS, PressureEOS
from quantas.core.physics.units import (
    convert_energy,
    energy_to_pressure,
    pressure_to_energy,
)
from quantas.modules.eos import (
    EOSArchive,
    EOSBatchJob,
    EOSBatchPlan,
    EOSBatchWorkflow,
    EOSCalculator,
    EOSDataset,
    EOSDiagnostics,
    EOSFitOptions,
    EOSFitRecord,
    EOSFitRequest,
    EOSFitter,
    EOSPlotter,
    EOSReportDetail,
    EOSReportOptions,
    build_eos_batch_report,
    infer_result_slots,
    parse_eos_spec,
    resolve_eos_spec,
)


def _energy_dataset(
    model: str = "BM3",
    *,
    include_pressure: bool = False,
    include_sigma: bool = False,
) -> EOSDataset:
    """Return a DFT-scale synthetic E-V dataset in public workflow units."""
    volume = np.linspace(17.2, 20.8, 13, dtype=np.float64)
    k0_internal = float(
        pressure_to_energy(180.0, "Ha", "angstrom", "GPa")
    )
    parameters = {
        "E0": -275.0,
        "K0": k0_internal,
        "KP": 4.2,
        "V0": 19.0,
    }
    energy = EnergyEOS().evaluate(model, volume, parameters)
    columns: dict[str, np.ndarray] = {
        "volume": volume,
        "energy": np.asarray(energy, dtype=np.float64),
    }
    units = {"volume": "angstrom^3", "energy": "Ha"}
    if include_pressure:
        pressure_internal = PressureEOS().pressure(model, parameters, volume)
        columns["pressure"] = np.asarray(
            energy_to_pressure(
                pressure_internal,
                "Ha",
                "angstrom",
                "GPa",
            ),
            dtype=np.float64,
        )
        units["pressure"] = "GPa"
    if include_sigma:
        columns["sigma_energy"] = np.full(volume.size, 2.0e-8, dtype=np.float64)
        units["sigma_energy"] = "Ha"
    return EOSDataset(
        jobname=f"synthetic {model} energy EOS",
        columns=columns,
        units=units,
    )


def _fit_energy_dataset(dataset: EOSDataset, model: str = "BM3"):
    request = EOSFitRequest(
        model=model,
        domain="ev",
        target="energy",
        options=EOSFitOptions(solver_options=OLSOptions(max_iterations=5000)),
    )
    result = EOSFitter().fit(dataset, request)
    assert result.fit.success, result.fit.message
    return request, result


def test_public_api_executes_energy_fit() -> None:
    dataset = _energy_dataset()
    request = public_eos.FitRequest(
        model="BM3",
        domain="ev",
        target="energy",
        options=public_eos.FitOptions(solver_options=public_eos.OLSOptions()),
    )

    result = public_eos.fit(dataset, request)

    assert result.fit.success, result.fit.message
    assert result.parameter_values["V0"] == pytest.approx(19.0, abs=2.0e-5)
    assert result.parameter_values["K0"] == pytest.approx(180.0, rel=2.0e-5)


@pytest.mark.parametrize("model", ["BM3", "T3", "SJ"])
def test_public_energy_fit_recovers_dft_scale_parameters(model: str) -> None:
    dataset = _energy_dataset(model)
    _, result = _fit_energy_dataset(dataset, model)

    values = result.parameter_values
    assert values["E0"] == pytest.approx(-275.0, abs=2.0e-8)
    assert values["V0"] == pytest.approx(19.0, abs=2.0e-5)
    assert values["K0"] == pytest.approx(180.0, rel=2.0e-5)
    assert values["KP"] == pytest.approx(4.2, rel=2.0e-5)
    assert result.metadata["pressure_relation"] == "P(V)=-dE/dV"
    assert set(result.predictions) == {
        "energy",
        "pressure",
        "bulk_modulus",
        "bulk_modulus_derivative",
        "bulk_modulus_second_derivative",
    }


@pytest.mark.parametrize("energy_unit", ["Ha", "eV", "Ry"])
def test_energy_fit_retains_dataset_energy_unit(energy_unit: str) -> None:
    """Equivalent E-V datasets fit directly in their declared energy unit."""
    reference = _energy_dataset("BM3")
    energy = np.asarray(
        convert_energy(reference.column("energy"), "Ha", energy_unit),
        dtype=np.float64,
    )
    dataset = EOSDataset(
        jobname=f"BM3 in {energy_unit}",
        columns={"volume": reference.column("volume"), "energy": energy},
        units={"volume": "angstrom^3", "energy": energy_unit},
    )

    _, result = _fit_energy_dataset(dataset, "BM3")

    assert result.parameter_values["E0"] == pytest.approx(
        float(convert_energy(-275.0, "Ha", energy_unit)), rel=2.0e-10
    )
    assert result.parameter_values["V0"] == pytest.approx(19.0, abs=2.0e-5)
    assert result.parameter_values["K0"] == pytest.approx(180.0, rel=2.0e-5)
    assert result.parameter_values["KP"] == pytest.approx(4.2, rel=2.0e-5)
    np.testing.assert_allclose(
        result.predictions["energy"],
        energy,
        rtol=2.0e-8,
        atol=2.0e-8 * max(1.0, float(np.max(np.abs(energy)))),
    )


def test_energy_fit_uncertainty_scales_with_energy_unit() -> None:
    """E0 uncertainty and energy covariance retain the dataset unit scale."""
    base = _energy_dataset("BM3", include_sigma=True)
    perturbation = np.linspace(-1.0, 1.0, base.npoints) * 1.0e-6

    results = {}
    for energy_unit in ("Ha", "eV"):
        energy = np.asarray(
            convert_energy(base.column("energy") + perturbation, "Ha", energy_unit),
            dtype=np.float64,
        )
        sigma = np.asarray(
            convert_energy(base.column("sigma_energy"), "Ha", energy_unit),
            dtype=np.float64,
        )
        dataset = EOSDataset(
            jobname=f"BM3 uncertainty in {energy_unit}",
            columns={
                "volume": base.column("volume"),
                "energy": energy,
                "sigma_energy": sigma,
            },
            units={
                "volume": "angstrom^3",
                "energy": energy_unit,
                "sigma_energy": energy_unit,
            },
        )
        request = EOSFitRequest(
            model="BM3",
            domain="ev",
            target="energy",
            options=EOSFitOptions(
                solver_options=WLSOptions(max_iterations=5000)
            ),
        )
        result = EOSFitter().fit(dataset, request)
        assert result.fit.success, result.fit.message
        results[energy_unit] = result

    factor = float(convert_energy(1.0, "Ha", "eV"))
    ha = results["Ha"].fit
    ev = results["eV"].fit
    # Covariance comes from a nonlinear numerical Jacobian.  The physical unit
    # scaling is exact, while a few-per-mille numerical variation is observed
    # across supported Python/NumPy/SciPy combinations.  The variance inherits
    # roughly twice the relative variation of the standard error.
    assert ev.errors[0] == pytest.approx(ha.errors[0] * factor, rel=2.0e-3)
    assert ev.covariance[0, 0] == pytest.approx(
        ha.covariance[0, 0] * factor**2, rel=4.0e-3
    )
    assert ev.errors[1] == pytest.approx(ha.errors[1], rel=2.0e-2)
    assert ev.errors[2] == pytest.approx(ha.errors[2], rel=2.0e-2)
    assert ev.errors[4] == pytest.approx(ha.errors[4], rel=2.0e-2)


def test_energy_wls_uses_explicit_sigma_energy() -> None:
    dataset = _energy_dataset(include_sigma=True)
    request = EOSFitRequest(
        model="BM3",
        domain="ev",
        target="energy",
        options=EOSFitOptions(solver_options=WLSOptions(max_iterations=5000)),
    )

    result = EOSFitter().fit(dataset, request)

    assert result.fit.success, result.fit.message
    assert result.parameter_values["K0"] == pytest.approx(180.0, rel=2.0e-5)


def test_energy_diagnostics_compare_optional_source_pressure() -> None:
    dataset = _energy_dataset(include_pressure=True)
    request, result = _fit_energy_dataset(dataset)
    record = EOSFitRecord(
        record_id=1,
        dataset_id=1,
        request=request,
        result=result,
    )
    diagnostic = EOSDiagnostics(record, dataset).build()

    np.testing.assert_allclose(
        diagnostic.columns["residual"],
        diagnostic.columns["observed_energy"]
        - diagnostic.columns["calculated_energy"],
    )
    np.testing.assert_allclose(
        diagnostic.columns["pressure_difference"],
        diagnostic.columns["source_pressure"] - diagnostic.columns["eos_pressure"],
        atol=2.0e-5,
    )
    assert diagnostic.units["residual"] == "Ha"
    assert diagnostic.units["eos_pressure"] == "GPa"


def test_energy_archive_preserves_non_hartree_unit(tmp_path: Path) -> None:
    """HDF5, diagnostics, and post-fit calculation retain eV semantics."""
    reference = _energy_dataset("BM3")
    energy = np.asarray(
        convert_energy(reference.column("energy"), "Ha", "eV"),
        dtype=np.float64,
    )
    dataset = EOSDataset(
        jobname="BM3 eV archive",
        columns={"volume": reference.column("volume"), "energy": energy},
        units={"volume": "angstrom^3", "energy": "eV"},
    )
    request, result = _fit_energy_dataset(dataset)
    path = tmp_path / "energy-ev.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        archive.store_fit(1, request, result)

    calculated = EOSCalculator.from_archive(path).calculate(
        volume=[19.0],
        propagate_uncertainty=False,
    )
    diagnostic = EOSDiagnostics.from_archive(path).build()

    assert calculated.units["energy"] == "eV"
    assert calculated.columns["energy"][0] == pytest.approx(
        float(convert_energy(-275.0, "Ha", "eV")), rel=2.0e-10
    )
    assert diagnostic.units["observed_energy"] == "eV"
    assert diagnostic.units["calculated_energy"] == "eV"
    assert diagnostic.units["residual"] == "eV"


def test_energy_archive_calculator_and_plot_round_trip(tmp_path: Path) -> None:
    dataset = _energy_dataset(include_pressure=True)
    request, result = _fit_energy_dataset(dataset)
    path = tmp_path / "energy.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        stored = archive.store_fit(1, request, result)

    assert {slot.key for slot in infer_result_slots(dataset)} == {
        "ev/energy",
        "pv/volume",
    }
    calculator = EOSCalculator.from_archive(path)
    at_volume = calculator.calculate(
        volume=[19.0],
        propagate_uncertainty=False,
    )
    assert at_volume.record_id == stored.record_id
    assert at_volume.columns["pressure"][0] == pytest.approx(0.0, abs=2.0e-5)
    assert at_volume.columns["energy"][0] == pytest.approx(-275.0, abs=2.0e-8)
    assert at_volume.columns["bulk_modulus"][0] == pytest.approx(180.0, rel=2.0e-5)

    inverse = calculator.calculate(pressure=[0.0], propagate_uncertainty=False)
    assert inverse.columns["volume"][0] == pytest.approx(19.0, abs=2.0e-5)

    diagnostic = EOSDiagnostics.from_archive(path).build()
    assert "observed_energy" in diagnostic.columns
    assert "eos_pressure" in diagnostic.columns

    plotter = EOSPlotter.from_archive(path)
    assert plotter.available_plot_types()[:3] == ("fit", "pressure", "residuals")
    collection = plotter.build(("fit", "pressure", "residuals"))
    assert {plot.key for plot in collection.plots} == {
        "fit",
        "pressure",
        "residuals_vs_volume",
    }

    listed = CliRunner().invoke(
        main,
        ["eos", "plot", str(path), "--slot", "ev/energy", "--list-plots"],
    )
    assert listed.exit_code == 0, listed.output
    assert "pressure" in listed.output


def test_energy_cli_runs_without_explicit_fit_target(tmp_path: Path) -> None:
    dataset = _energy_dataset()
    source = tmp_path / "energy.dat"
    with source.open("w", encoding="utf-8") as handle:
        handle.write("JOB synthetic energy EOS\n")
        handle.write("UNITS V=angstrom^3 E=Ha\n")
        handle.write("FORMAT V E\n")
        handle.write("DATA\n")
        for volume, energy in zip(
            dataset.column("volume"), dataset.column("energy"), strict=True
        ):
            handle.write(f"{volume:.12f} {energy:.16f}\n")

    archive = tmp_path / "energy.hdf5"
    report = tmp_path / "energy.log"
    result = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(source),
            "--domain",
            "ev",
            "--ev-eos",
            "BM3",
            "--output",
            str(archive),
            "--report",
            str(report),
            "--force",
        ],
    )

    assert result.exit_code == 0, result.output
    assert archive.exists()
    assert "Energy-volume" in result.output
    assert "Energy" in result.output
    report_text = report.read_text(encoding="utf-8")
    assert "Ab initio energy" in report_text
    assert "EOS energy" in report_text
    assert "EOS pressure" in report_text


def test_energy_spec_defaults_resolve_all_to_energy() -> None:
    dataset = _energy_dataset()
    document = parse_eos_spec(
        """# QUANTAS EOS SPEC 1

[defaults.ev]
model = SJ
solver = ols

[job static-energy]
domain = ev
targets = all
"""
    )

    resolved = resolve_eos_spec(document, dataset)

    assert len(resolved.plan.jobs) == 1
    request = resolved.plan.jobs[0].request
    assert request.domain.value == "ev"
    assert request.target == "energy"
    assert request.model.tag == "SJ"


def test_energy_extended_report_data(
    tmp_path: Path,
) -> None:
    dataset = _energy_dataset(include_sigma=True)
    request = EOSFitRequest(
        model="BM3",
        domain="ev",
        target="energy",
        options=EOSFitOptions(solver_options=WLSOptions(max_iterations=5000)),
    )
    plan = EOSBatchPlan(jobs=(EOSBatchJob(request, job_id="static-energy"),))
    archive = tmp_path / "energy-report.hdf5"

    result = EOSBatchWorkflow().run(dataset, plan, archive)
    tables = build_eos_batch_report(
        result,
        EOSReportOptions(
            detail=EOSReportDetail.EXTENDED,
            show_uncertainties=True,
            max_data_rows=3,
        ),
    )
    by_title = {table.title: table for table in tables}

    assert by_title["EOS input data"].columns == ["Volume", "Energy"]
    assert by_title["EOS input standard uncertainties"].columns == ["sigma(Energy)"]
    observed = by_title["Observed and calculated EOS data"]
    assert observed.columns == [
        "Volume",
        "Ab initio energy",
        "EOS energy",
        "Residual",
        "EOS pressure",
    ]
    assert observed.metadata["column_formats"] == [
        "eos_structural",
        "energy",
        "energy",
        "eos_residual",
        "eos_pressure",
    ]
    assert observed.metadata["column_units"] == [
        "angstrom^3",
        "Ha",
        "Ha",
        "Ha",
        "GPa",
    ]
    assert len(observed.rows) == dataset.npoints
    assert observed.rows[0][4] == pytest.approx(
        result.jobs[0].result.predictions["pressure"][0]
    )
