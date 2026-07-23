"""Click adapter for ``quantas eos run``."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner
import numpy as np
import pytest

from quantas.cli.main import main
from quantas.core.physics.eos import MGDNormalization, PVTEOS, PVTModel
from quantas.modules.eos import EOSArchive, EOSSlotStatus

DATA = Path(__file__).with_name("data")


def test_eos_run_help_groups_options_by_responsibility() -> None:
    result = CliRunner().invoke(main, ["eos", "run", "--help"])

    assert result.exit_code == 0
    for heading in (
        "Scientific model:",
        "Units:",
        "Fit selection:",
        "P-V model:",
        "V-T model:",
        "P-V-T coupling:",
        "MGD normalization:",
        "Parameter constraints:",
        "Numerical solver:",
        "Output and reporting:",
    ):
        assert heading in result.output
    assert "--specfile" not in result.output
    assert "-v, --verbosity" in result.output
    assert "--debug" not in result.output


def test_cli_debug_shows_solver_trace_and_persists_failure_metadata(
    tmp_path: Path,
) -> None:
    archive_path = tmp_path / "failed-debug.hdf5"
    result = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(DATA / "PV_quartz.dat"),
            "--solver",
            "ols",
            "--max-iterations",
            "1",
            "--verbosity",
            "debug",
            "--output",
            str(archive_path),
            "--report",
            str(tmp_path / "failed-debug.log"),
        ],
    )

    assert result.exit_code != 0
    assert "EOS solver debug summary" in result.output
    assert "Termination category" in result.output
    assert "iteration_limit" in result.output
    assert "EOS solver numerical ranges" in result.output
    assert "EOS solver parameter path" in result.output
    assert "Solver model-evaluation trace" in result.output
    assert "Rerun the EOS batch with -v debug" not in result.output
    with EOSArchive(archive_path) as archive:
        record = archive.record(1)
        assert record.result.fit.diagnostics is not None
        metadata = record.result.fit.diagnostics.metadata
        assert metadata["termination_category"] == "iteration_limit"
        assert metadata["detailed_trace_requested"] is True
        assert metadata["evaluation_trace"]


def test_failed_run_without_debug_suggests_detailed_diagnostics(
    tmp_path: Path,
) -> None:
    result = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(DATA / "PV_quartz.dat"),
            "--solver",
            "ols",
            "--max-iterations",
            "1",
            "--output",
            str(tmp_path / "failed-normal.hdf5"),
            "--report",
            str(tmp_path / "failed-normal.log"),
        ],
    )

    assert result.exit_code != 0
    assert "rerun with -v debug for detailed solver diagnostics" in result.output
    assert "EOS solver debug summary" not in result.output


def test_quiet_run_writes_archive_and_plain_report(tmp_path: Path) -> None:
    archive_path = tmp_path / "quartz.hdf5"
    report_path = tmp_path / "quartz.txt"
    result = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(DATA / "PV_quartz.dat"),
            "--solver",
            "effective-variance",
            "--output",
            str(archive_path),
            "--report",
            str(report_path),
            "--show-uncertainties",
            "--max-data-rows",
            "2",
            "--quiet",
        ],
    )

    assert result.exit_code == 0, result.output
    assert result.output == ""
    assert archive_path.is_file()
    text = report_path.read_text(encoding="utf-8")
    assert "EOS input data" in text
    assert "EOS parameters" in text
    assert "EOS input standard uncertainties" in text
    assert "\x1b[" not in text
    with EOSArchive(archive_path) as archive:
        assert archive.slot_state("pv/volume").status is EOSSlotStatus.ACCEPTED
        assert (
            archive.batch_manifest()["jobs"][0]["request"]["options"]["method"]
            == "effective_variance"
        )


def test_fit_all_executes_every_available_topaz_target(tmp_path: Path) -> None:
    archive_path = tmp_path / "topaz.hdf5"
    result = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(DATA / "PV_topaz.dat"),
            "--fit",
            "all",
            "--output",
            str(archive_path),
            "--report",
            str(tmp_path / "topaz.log"),
            "--quiet",
        ],
    )

    assert result.exit_code == 0, result.output
    with EOSArchive(archive_path) as archive:
        assert archive.record_ids == (1, 2, 3, 4)
        assert {
            state.slot.key
            for state in archive.slots()
            if state.status is EOSSlotStatus.ACCEPTED
        } == {"pv/volume", "pv/a", "pv/b", "pv/c"}


def test_existing_archive_requires_force(tmp_path: Path) -> None:
    archive_path = tmp_path / "existing.hdf5"
    archive_path.write_bytes(b"existing")
    result = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(DATA / "PV_quartz.dat"),
            "--output",
            str(archive_path),
            "--report",
            str(tmp_path / "existing.log"),
            "--quiet",
        ],
    )

    assert result.exit_code != 0
    assert "use --force" in result.output
    assert archive_path.read_bytes() == b"existing"


def test_cli_fixed_parameter_is_stored_in_request(tmp_path: Path) -> None:
    archive_path = tmp_path / "fixed.hdf5"
    result = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(DATA / "PV_quartz.dat"),
            "--fix",
            "KP=4",
            "--output",
            str(archive_path),
            "--report",
            str(tmp_path / "fixed.log"),
            "--quiet",
        ],
    )

    assert result.exit_code == 0, result.output
    with EOSArchive(archive_path) as archive:
        record = archive.record(1)
        constraints = record.request.as_dict()["constraints"]
        assert len(constraints) == 1
        assert constraints[0]["name"] == "KP"
        assert constraints[0]["state"] == "fixed"
        assert constraints[0]["value"] == 4.0


def test_cli_applies_dataset_selection_and_spec_groups(tmp_path: Path) -> None:
    data_path = tmp_path / "grouped.dat"
    data_path.write_text(
        "\n".join(
            (
                "TITLE Grouped P-V data",
                "SYSTEM tetragonal",
                "FORMAT P V GROUP USE",
                "0.0 100.0 1 1",
                "1.0  99.2 1 1 *",
                "2.0  98.5 2 0",
                "3.0  97.8 2 1",
                "4.0  97.2 2 1",
                "5.0  96.7 2 1",
            )
        )
        + "\n",
        encoding="utf-8",
    )
    spec_path = tmp_path / "grouped.spec"
    spec_path.write_text(
        "\n".join(
            (
                "# QUANTAS EOS SPEC 1",
                "[job group-2]",
                "domain = pv",
                "targets = volume",
                "model = BM2",
                "groups = 2",
                "solver = ols",
            )
        )
        + "\n",
        encoding="utf-8",
    )
    archive_path = tmp_path / "grouped.hdf5"

    result = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(data_path),
            "--spec",
            str(spec_path),
            "--output",
            str(archive_path),
            "--report",
            str(tmp_path / "grouped.log"),
            "--quiet",
        ],
    )

    assert result.exit_code == 0, result.output
    with EOSArchive(archive_path) as archive:
        dataset = archive.dataset(1)
        record = archive.record(1)
        np.testing.assert_array_equal(
            dataset.default_mask, [True, False, False, True, True, True]
        )
        np.testing.assert_array_equal(
            record.request.mask, [False, False, False, True, True, True]
        )
        assert record.request.metadata["selection"]["selected"] == 3
        assert dataset.crystal_system.value == "tetragonal"


def _write_mgd_pvt_data(path: Path) -> Path:
    """Write one exact synthetic P-V-T table in the public text format."""
    model = PVTModel(
        "BM3",
        "thermal-pressure",
        thermal_pressure_model="mgd",
        mgd_normalization=MGDNormalization.cell(atoms_per_cell=8),
    )
    temperature = np.repeat(np.asarray([100.0, 295.0, 600.0, 1000.0]), 7)
    pressure = np.tile(np.linspace(0.0, 18.0, 7), 4)
    volume = PVTEOS().volume(
        model,
        {"K0": 48.0, "KP": 4.5, "KPP": -0.1, "V0": 150.0},
        None,
        {
            "temperature_ref": 295.0,
            "theta_d0": 459.0,
            "gamma0": 1.547,
            "q": 0.94,
        },
        pressure,
        temperature,
    )
    lines = ["TITLE Synthetic MGD PVT", "FORMAT 1 T P V"]
    lines.extend(
        f"{t:.8f} {p_value:.12f} {v:.12f}"
        for t, p_value, v in zip(temperature, pressure, volume, strict=True)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_cli_executes_volume_only_mgd_fit(tmp_path: Path) -> None:
    data_path = _write_mgd_pvt_data(tmp_path / "mgd.dat")
    archive_path = tmp_path / "mgd.hdf5"
    result = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(data_path),
            "--domain",
            "pvt",
            "--pvt-coupling",
            "thermal-pressure",
            "--thermal-pressure-model",
            "mgd",
            "--formula",
            "NaF",
            "--formula-units-per-cell",
            "4",
            "--initial",
            "K0=47",
            "--initial",
            "KP=4.3",
            "--initial",
            "V0=149.5",
            "--fix",
            "temperature_ref=295",
            "--initial",
            "theta_d0=440",
            "--bound",
            "theta_d0=100:1000",
            "--initial",
            "gamma0=1.4",
            "--initial",
            "q=0.8",
            "--max-iterations",
            "10000",
            "--output",
            str(archive_path),
            "--report",
            str(tmp_path / "mgd.log"),
            "--quiet",
        ],
    )

    assert result.exit_code == 0, result.output
    with EOSArchive(archive_path) as archive:
        record = archive.record(1)
        assert record.request.model.tag.endswith(
            "mie-gruneisen-debye:full+thermal-pressure"
        )
        assert record.result.parameter_values["theta_d0"] == pytest.approx(
            459.0, rel=3e-5
        )
        assert record.result.parameter_values["gamma0"] == pytest.approx(
            1.547, rel=3e-6
        )
        assert record.result.parameter_values["q"] == pytest.approx(0.94, rel=3e-5)


def test_cli_rejects_mgd_options_outside_pvt_domain(tmp_path: Path) -> None:
    result = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(DATA / "PV_quartz.dat"),
            "--thermal-pressure-model",
            "mgd",
            "--atoms-per-cell",
            "8",
            "--dry-run",
            "--report",
            str(tmp_path / "invalid-mgd.log"),
        ],
    )
    assert result.exit_code != 0
    assert "require --domain pvt" in result.output


def test_eos_diagnose_writes_csv(tmp_path: Path) -> None:
    archive_path = tmp_path / "diagnose.hdf5"
    run = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(DATA / "PV_quartz.dat"),
            "--pv-eos",
            "BM",
            "--pv-order",
            "3",
            "--output",
            str(archive_path),
            "--report",
            str(tmp_path / "diagnose.log"),
            "--quiet",
        ],
    )
    assert run.exit_code == 0, run.output
    csv_path = tmp_path / "diagnostics.csv"

    result = CliRunner().invoke(
        main,
        ["eos", "diagnose", str(archive_path), "--output", str(csv_path)],
    )

    assert result.exit_code == 0, result.output
    assert "EOS diagnostics summary" in result.output
    assert "Strain family" in result.output
    assert csv_path.is_file()
    header = csv_path.read_text(encoding="utf-8").splitlines()[0]
    assert "normalized_pressure [GPa]" in header


def test_eos_calculate_pressure_range_writes_csv(tmp_path: Path) -> None:
    archive_path = tmp_path / "calculate.hdf5"
    run = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(DATA / "PV_quartz.dat"),
            "--output",
            str(archive_path),
            "--report",
            str(tmp_path / "calculate.log"),
            "--quiet",
        ],
    )
    assert run.exit_code == 0, run.output
    csv_path = tmp_path / "states.csv"

    result = CliRunner().invoke(
        main,
        [
            "eos",
            "calculate",
            str(archive_path),
            "--pressure-range",
            "0:2:1",
            "--output",
            str(csv_path),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "EOS calculator summary" in result.output
    assert "EOS calculated properties" in result.output
    assert csv_path.is_file()
    lines = csv_path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 4
    assert "sigma_volume" in lines[0]


def test_eos_calculate_help_documents_grid_and_uncertainty_options() -> None:
    result = CliRunner().invoke(main, ["eos", "calculate", "--help"])

    assert result.exit_code == 0
    for heading in (
        "Record selection:",
        "State coordinates:",
        "Uncertainty propagation:",
        "Output and reporting:",
    ):
        assert heading in result.output
    assert "--pressure-range" in result.output
    assert "--pairwise" in result.output
