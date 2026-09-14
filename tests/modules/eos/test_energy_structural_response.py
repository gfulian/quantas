"""Theoretical axial response derived from public EnergyEOS fits."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.math.fitting import OLSOptions
from quantas.core.physics.eos import EnergyEOS
from quantas.core.physics.units import pressure_to_energy
from quantas.modules.eos import (
    EOSArchive,
    EOSCalculator,
    EOSDataset,
    EOSFitOptions,
    EOSFitRequest,
    EOSFitter,
)


def _energy_parameters() -> dict[str, float]:
    return {
        "E0": -275.0,
        "K0": float(pressure_to_energy(180.0, "Ha", "angstrom", "GPa")),
        "KP": 4.2,
        "V0": 19.0,
    }


def _structural_dataset(
    *,
    model: str = "BM3",
    tetragonal_exponent: float | None = None,
    perturb_energy: bool = False,
) -> EOSDataset:
    volume = np.linspace(17.2, 20.8, 13, dtype=np.float64)
    energy = np.asarray(
        EnergyEOS().evaluate(model, volume, _energy_parameters()),
        dtype=np.float64,
    )
    if perturb_energy:
        energy = energy + np.linspace(-1.0, 1.0, volume.size) * 2.0e-7
    if tetragonal_exponent is None:
        a = np.cbrt(volume)
        b = a.copy()
        c = a.copy()
        system = "cubic"
        space_group = 221
        symbol = "Pm-3m"
    else:
        exponent = float(tetragonal_exponent)
        a0 = 3.0
        a = a0 * (volume / 19.0) ** exponent
        b = a.copy()
        c = volume / a**2
        system = "tetragonal"
        space_group = 123
        symbol = "P4/mmm"
    ninety = np.full(volume.size, 90.0, dtype=np.float64)
    return EOSDataset(
        jobname="synthetic structural EnergyEOS",
        columns={
            "volume": volume,
            "energy": energy,
            "a": a,
            "b": b,
            "c": c,
            "alpha": ninety,
            "beta": ninety,
            "gamma": ninety,
        },
        units={
            "volume": "angstrom^3",
            "energy": "Ha",
            "a": "angstrom",
            "b": "angstrom",
            "c": "angstrom",
            "alpha": "degree",
            "beta": "degree",
            "gamma": "degree",
        },
        metadata={
            "crystal_reference": "crystallographic",
            "crystal_system": system,
            "space_group_number": space_group,
            "space_group_symbol": symbol,
        },
    )


def _fit(dataset: EOSDataset, *, model: str = "BM3", axial_model: str | None = None):
    request = EOSFitRequest(
        model=model,
        axial_model=axial_model,
        domain="ev",
        target="energy",
        options=EOSFitOptions(solver_options=OLSOptions(max_iterations=5000)),
    )
    result = EOSFitter().fit(dataset, request)
    assert result.fit.success, result.fit.message
    return request, result


def test_cubic_primary_response_is_exact() -> None:
    dataset = _structural_dataset()
    _, result = _fit(dataset)

    assert result.derived["a0"] == pytest.approx(19.0 ** (1.0 / 3.0), rel=2.0e-5)
    assert result.derived["eta_a"] == pytest.approx(1.0 / 3.0, abs=1.0e-14)
    assert result.derived["M_a"] == pytest.approx(540.0, rel=3.0e-5)
    assert result.metadata["structural_response"]["independent_axes"] == ["a"]
    assert result.metadata["structural_response"]["uncertainty_sources"][
        "cross_covariance_energy_structural"
    ] is False
    assert "structural_a" in result.predictions
    assert "eta_a" in result.predictions
    assert "M_a" in result.predictions


def test_tetragonal_response_has_independent_axes() -> None:
    exponent = 0.20
    dataset = _structural_dataset(tetragonal_exponent=exponent)
    _, result = _fit(dataset)

    assert result.derived["eta_a"] == pytest.approx(exponent, abs=2.0e-6)
    assert result.derived["eta_c"] == pytest.approx(1.0 - 2.0 * exponent, abs=2.0e-6)
    assert result.derived["M_a"] == pytest.approx(180.0 / exponent, rel=3.0e-5)
    assert result.derived["M_c"] == pytest.approx(
        180.0 / (1.0 - 2.0 * exponent), rel=3.0e-5
    )
    assert result.metadata["structural_response"]["independent_axes"] == ["a", "c"]


def test_sjeos_supports_primary_structural_response() -> None:
    dataset = _structural_dataset(model="SJ")
    _, result = _fit(dataset, model="SJ")

    assert result.derived["eta_a"] == pytest.approx(1.0 / 3.0, abs=1.0e-14)
    assert result.derived["M_a"] == pytest.approx(3.0 * result.parameter_values["K0"])


def test_secondary_axial_fit_uses_derived_pressure_wls() -> None:
    dataset = _structural_dataset(model="SJ", perturb_energy=True)
    request, result = _fit(dataset, model="SJ", axial_model="BM3")

    assert request.axial_model is not None
    secondary = result.metadata["secondary_axial_fits"]
    assert secondary["model"]["tag"] == "BM3"
    assert secondary["method"] == "diagonal_wls_marginal_pressure_uncertainties"
    assert np.asarray(secondary["pressure_covariance"]).shape == (13, 13)
    axial = secondary["fits"]["a"]
    assert axial["parameter_values"]["M0"] == pytest.approx(
        result.derived["M_a"], rel=5.0e-3
    )
    assert result.derived["axial_a_M0"] == pytest.approx(
        axial["parameter_values"]["M0"]
    )


def test_structural_energy_calculator_round_trip(tmp_path: Path) -> None:
    dataset = _structural_dataset(perturb_energy=True)
    request, result = _fit(dataset)
    archive_path = tmp_path / "structural-energy.hdf5"
    with EOSArchive.create(archive_path, dataset=dataset) as archive:
        archive.store_fit(1, request, result)

    calculator = EOSCalculator.from_archive(archive_path)
    calculation = calculator.calculate(volume=[19.0])

    assert calculation.columns["structural_a"][0] == pytest.approx(
        19.0 ** (1.0 / 3.0), rel=2.0e-5
    )
    assert calculation.columns["eta_a"][0] == pytest.approx(1.0 / 3.0)
    assert calculation.columns["M_a"][0] == pytest.approx(
        3.0 * calculation.columns["bulk_modulus"][0]
    )
    assert calculation.units["M_a"] == "GPa"
    assert "M_a" in calculation.uncertainties
    assert calculation.metadata["structural_response"]["basis"] == "sampled"


def test_request_round_trip_retains_axial_model(tmp_path: Path) -> None:
    dataset = _structural_dataset(perturb_energy=True)
    request, result = _fit(dataset, model="BM3", axial_model="BM2")
    archive_path = tmp_path / "secondary-axial.hdf5"
    with EOSArchive.create(archive_path, dataset=dataset) as archive:
        stored = archive.store_fit(1, request, result)
        restored = archive.record(stored.record_id)

    assert restored.request.axial_model is not None
    assert restored.request.axial_model.tag == "BM2"
    assert restored.result.metadata["secondary_axial_fits"]["model"]["tag"] == "BM2"


def test_energy_spec_accepts_secondary_axial_model() -> None:
    from quantas.modules.eos import parse_eos_spec, resolve_eos_spec

    dataset = _structural_dataset()
    document = parse_eos_spec(
        """# QUANTAS EOS SPEC 1

[defaults.ev]
model = SJ
axial_model = BM3

[job structural]
domain = ev
targets = energy
"""
    )

    resolved = resolve_eos_spec(document, dataset)
    request = resolved.plan.jobs[0].request
    assert request.model.tag == "SJ"
    assert request.axial_model is not None
    assert request.axial_model.tag == "BM3"


def test_axial_model_validation_checks_domain_and_capability() -> None:
    with pytest.raises(ValueError, match="only for E-V"):
        EOSFitRequest(model="BM3", domain="pv", target="volume", axial_model="BM3")
    with pytest.raises(ValueError, match="not exposed for axial pressure fitting"):
        EOSFitRequest(model="BM3", domain="ev", target="energy", axial_model="SJ")


def test_mgo_crystallographic_normalization_preserves_bulk() -> None:
    root = Path(__file__).resolve().parents[3]
    from quantas.api import eos as eos_api

    primitive = eos_api.read_input(root / "examples/eos/EV_mgo_pbe.dat")
    crystallographic = eos_api.read_input(
        root / "examples/eos/EV_mgo_pbe_crystallographic.dat"
    )
    request = eos_api.FitRequest(model="BM3", domain="ev", target="energy")
    primitive_result = eos_api.fit(primitive, request)
    crystallographic_result = eos_api.fit(crystallographic, request)

    assert primitive_result.fit.success
    assert crystallographic_result.fit.success
    assert crystallographic_result.parameter_values["V0"] == pytest.approx(
        4.0 * primitive_result.parameter_values["V0"], rel=1.0e-12
    )
    assert crystallographic_result.parameter_values["E0"] == pytest.approx(
        4.0 * primitive_result.parameter_values["E0"], rel=1.0e-12
    )
    assert crystallographic_result.parameter_values["K0"] == pytest.approx(
        primitive_result.parameter_values["K0"], rel=1.0e-12
    )
    assert crystallographic_result.parameter_values["KP"] == pytest.approx(
        primitive_result.parameter_values["KP"], rel=1.0e-12
    )
    volumes = crystallographic.column("volume")
    energies = crystallographic.column("energy")
    a_values = crystallographic.column("a")
    reference_index = int(np.argmin(energies))
    expected_a0 = float(a_values[reference_index]) * (
        crystallographic_result.parameter_values["V0"]
        / float(volumes[reference_index])
    ) ** (1.0 / 3.0)
    assert crystallographic_result.derived["a0"] == pytest.approx(
        expected_a0, rel=1.0e-12
    )
    assert crystallographic_result.derived["a0"] == pytest.approx(
        4.22221, abs=5.0e-6
    )
    assert crystallographic_result.derived["eta_a"] == pytest.approx(1.0 / 3.0)
    assert crystallographic_result.derived["M_a"] == pytest.approx(
        3.0 * crystallographic_result.parameter_values["K0"], rel=1.0e-12
    )
    assert crystallographic_result.metadata["structural_response"][
        "independent_axes"
    ] == ["a"]


def test_structural_response_report_tables_are_explicit() -> None:
    from quantas.modules.eos.report import (
        eos_energy_structural_response_table,
        eos_secondary_axial_table,
    )

    dataset = _structural_dataset(model="SJ", perturb_energy=True)
    _, result = _fit(dataset, model="SJ", axial_model="BM3")

    primary = eos_energy_structural_response_table(result)
    secondary = eos_secondary_axial_table(result)
    assert primary is not None
    assert primary.title == "E-V structural response"
    assert primary.title.isascii()
    assert all(str(row[0]).isascii() for row in primary.rows)
    assert any(row[0] == "a0" for row in primary.rows)
    assert any(row[0] == "eta_a" for row in primary.rows)
    assert any(row[0] == "M_a" for row in primary.rows)
    assert any(row[3] == "angstrom" for row in primary.rows)
    assert any(row[3] == "-" for row in primary.rows)
    assert secondary is not None
    assert secondary.title == "Secondary axial EOS fits"
    assert any(row[0] == "a" and row[1] == "M0" for row in secondary.rows)


def test_cli_secondary_axial_fit_persists_requested_model(tmp_path: Path) -> None:
    from click.testing import CliRunner

    from quantas.cli.main import main

    dataset = _structural_dataset(model="SJ", perturb_energy=True)
    source = tmp_path / "structural-ev.dat"
    lines = [
        "JOB structural EnergyEOS CLI",
        "SYSTEM cubic",
        "SPACE_GROUP_NUMBER 221",
        "SPACE_GROUP_SYMBOL Pm-3m",
        "CRYSTAL_REFERENCE crystallographic",
        "CELL_MULTIPLICITY 1",
        "UNITS V=angstrom^3 A=angstrom B=angstrom C=angstrom "
        "ALPHA=degree BETA=degree GAMMA=degree E=Ha",
        "FORMAT V A B C ALPHA BETA GAMMA E",
        "DATA",
    ]
    for index in range(dataset.npoints):
        lines.append(
            " ".join(
                format(float(dataset.column(name)[index]), ".15g")
                for name in (
                    "volume",
                    "a",
                    "b",
                    "c",
                    "alpha",
                    "beta",
                    "gamma",
                    "energy",
                )
            )
        )
    source.write_text("\n".join(lines) + "\n", encoding="utf-8")
    archive_path = tmp_path / "structural-ev.hdf5"
    report_path = tmp_path / "structural-ev.log"

    completed = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(source),
            "--domain",
            "ev",
            "--ev-eos",
            "SJ",
            "--axial-eos",
            "BM3",
            "--output",
            str(archive_path),
            "--report",
            str(report_path),
            "--quiet",
        ],
    )

    assert completed.exit_code == 0, completed.output
    with EOSArchive(archive_path) as archive:
        record = archive.accepted("ev/energy")
        assert record is not None
        assert record.request.axial_model is not None
        assert record.request.axial_model.tag == "BM3"
        assert record.result.metadata["secondary_axial_fits"]["model"]["tag"] == "BM3"
    report = report_path.read_text(encoding="utf-8")
    assert "E-V structural response" in report
    assert "Secondary axial EOS fits" in report


def test_cli_rejects_secondary_axial_option_outside_ev(tmp_path: Path) -> None:
    from click.testing import CliRunner

    from quantas.cli.main import main

    root = Path(__file__).resolve().parent
    completed = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(root / "data/PV_quartz.dat"),
            "--axial-eos",
            "BM3",
            "--output",
            str(tmp_path / "bad.hdf5"),
            "--report",
            str(tmp_path / "bad.log"),
        ],
    )

    assert completed.exit_code != 0
    assert "--axial-eos is valid only with --domain ev" in completed.output
