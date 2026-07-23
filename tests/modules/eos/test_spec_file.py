"""Strict EOS specification files and CLI dry-run integration."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner
import pytest

from quantas.cli.main import main
from quantas.core.math.fitting import FitMethod
from quantas.core.physics.eos import PVTModel
from quantas.modules.eos import (
    EOSArchive,
    EOSReportDetail,
    EOSSpecError,
    EOSSlotStatus,
    parse_eos_spec,
    read_eos_input,
    read_eos_spec,
    resolve_eos_spec,
)

DATA = Path(__file__).with_name("data")


def _write_spec(path: Path, body: str) -> Path:
    path.write_text("# QUANTAS EOS SPEC 1\n\n" + body.strip() + "\n", encoding="utf-8")
    return path


def test_spec_signature_and_unknown_key_are_strict() -> None:
    with pytest.raises(EOSSpecError, match="first non-empty line"):
        parse_eos_spec("[job volume]\ndomain = pv\ntargets = volume\n")

    text = """# QUANTAS EOS SPEC 1

[job volume]
domain = pv
targets = volume
solvert = odr
"""
    with pytest.raises(EOSSpecError) as caught:
        parse_eos_spec(text, source="bad.spec")
    message = str(caught.value)
    assert "bad.spec:6" in message
    assert "unknown key 'solvert'" in message


def test_spec_resolves_defaults_multi_target_and_presentation(tmp_path: Path) -> None:
    spec_path = _write_spec(
        tmp_path / "topaz.anything",
        """
[metadata]
title = Topaz mixed analysis

[input]
pressure_unit = GPa
length_unit = angstrom

[batch]
failure_policy = continue

[defaults]
solver = effective-variance
covariance_scaling = inflate-only
accept = yes

[defaults.pv]
model = BM3

[presentation]
detail = extended
show_uncertainties = yes
max_data_rows = 7

[job volume]
domain = pv
targets = volume
model = BM2
fix.KP = 4.0

[job axes]
domain = pv
targets = a, b, c
solver = ols
model = PT3
""",
    )
    document = read_eos_spec(spec_path)
    dataset = read_eos_input(DATA / "PV_topaz.dat", **document.input_options.as_dict())
    resolved = resolve_eos_spec(document, dataset)

    assert document.metadata["title"] == "Topaz mixed analysis"
    assert resolved.plan.failure_policy.value == "continue"
    assert tuple(job.job_id for job in resolved.plan.jobs) == (
        "volume",
        "axes-a",
        "axes-b",
        "axes-c",
    )
    assert tuple(job.request.model.tag for job in resolved.plan.jobs) == (
        "BM2",
        "PT3",
        "PT3",
        "PT3",
    )
    assert resolved.plan.jobs[0].request.options.method is FitMethod.EFFECTIVE_VARIANCE
    assert resolved.plan.jobs[1].request.options.method is FitMethod.OLS
    constraint = resolved.plan.jobs[0].request.constraints[0]
    assert constraint.name == "KP"
    assert constraint.value == 4.0
    assert resolved.report_options.detail is EOSReportDetail.EXTENDED
    assert resolved.report_options.show_uncertainties is True
    assert resolved.report_options.max_data_rows == 7


def test_spec_rejects_unknown_solver_and_incompatible_solver_option(
    tmp_path: Path,
) -> None:
    dataset = read_eos_input(DATA / "PV_quartz.dat")
    unknown = read_eos_spec(
        _write_spec(
            tmp_path / "unknown",
            """
[job volume]
domain = pv
targets = volume
solver = odl
""",
        )
    )
    with pytest.raises(EOSSpecError, match="unknown solver 'odl'"):
        resolve_eos_spec(unknown, dataset)

    incompatible = read_eos_spec(
        _write_spec(
            tmp_path / "incompatible",
            """
[job volume]
domain = pv
targets = volume
solver = ols
odr_ndigit = 12
""",
        )
    )
    with pytest.raises(EOSSpecError, match="odr_ndigit is valid only"):
        resolve_eos_spec(incompatible, dataset)


def test_spec_rejects_constraint_ambiguity_and_unknown_parameter(
    tmp_path: Path,
) -> None:
    dataset = read_eos_input(DATA / "PV_quartz.dat")
    conflict = read_eos_spec(
        _write_spec(
            tmp_path / "conflict.spec",
            """
[job volume]
domain = pv
targets = volume
fix.KP = 4
initial.KP = 5
""",
        )
    )
    with pytest.raises(EOSSpecError, match="fixed parameter 'KP'"):
        resolve_eos_spec(conflict, dataset)

    unknown = read_eos_spec(
        _write_spec(
            tmp_path / "unknown-parameter.spec",
            """
[job volume]
domain = pv
targets = volume
initial.banana = 1
""",
        )
    )
    with pytest.raises(EOSSpecError, match="parameter 'banana' is not available"):
        resolve_eos_spec(unknown, dataset)


def test_spec_thermal_pressure_rejects_independent_vt_model(tmp_path: Path) -> None:
    dataset = read_eos_input(DATA / "PV_quartz.dat")
    document = read_eos_spec(
        _write_spec(
            tmp_path / "pvt.spec",
            """
[job pvt]
domain = pvt
targets = volume
pv_model = BM3
vt_model = berman:quadratic
coupling = thermal-pressure
""",
        )
    )
    with pytest.raises(EOSSpecError, match="remove vt_model"):
        resolve_eos_spec(document, dataset)


def test_spec_resolves_mgd_volume_only_fit(tmp_path: Path) -> None:
    """Resolve MGD normalization and refinable parameters without K_S data."""
    dataset = read_eos_input(DATA / "PV_quartz.dat")
    # Resolution validates the scientific contract; the existing PV-only file
    # is sufficient here because no numerical fit is executed in this test.
    document = read_eos_spec(
        _write_spec(
            tmp_path / "mgd.spec",
            """
[job mgd]
domain = pvt
targets = volume
pv_model = BM3
coupling = thermal-pressure
thermal_pressure_model = mgd
volume_basis = cell
formula = NaF
formula_units_per_cell = 4
fix.temperature_ref = 295
initial.theta_d0 = 450
bound.theta_d0 = 100 : 1000
initial.gamma0 = 1.5
initial.q = 1.0
solver = ols
""",
        )
    )
    resolved = resolve_eos_spec(document, dataset)
    request = resolved.plan.jobs[0].request
    assert isinstance(request.model, PVTModel)
    assert request.model.thermal_pressure_spec is not None
    assert request.model.thermal_pressure_spec.tag == "mie-gruneisen-debye:full"
    normalization = request.model.mgd_normalization_spec
    assert normalization is not None
    assert normalization.atoms_per_cell == pytest.approx(8.0)
    assert normalization.formula == "NaF"
    constraints = {item.name: item for item in request.constraints}
    assert constraints["temperature_ref"].value == pytest.approx(295.0)
    assert constraints["theta_d0"].initial_value == pytest.approx(450.0)
    assert constraints["gamma0"].initial_value == pytest.approx(1.5)
    assert constraints["q"].initial_value == pytest.approx(1.0)


def test_spec_mgd_requires_cell_atom_normalization(tmp_path: Path) -> None:
    dataset = read_eos_input(DATA / "PV_quartz.dat")
    document = read_eos_spec(
        _write_spec(
            tmp_path / "bad-mgd.spec",
            """
[job mgd]
domain = pvt
targets = volume
coupling = thermal-pressure
thermal_pressure_model = mgd
""",
        )
    )
    with pytest.raises(EOSSpecError, match="atoms_per_cell|formula"):
        resolve_eos_spec(document, dataset)


def test_cli_spec_dry_run_creates_no_archive_and_prints_resolved_plan(
    tmp_path: Path,
) -> None:
    spec_path = _write_spec(
        tmp_path / "analysis.txt",
        """
[presentation]
detail = short
max_data_rows = 2

[job volume]
domain = pv
targets = volume
model = BM2
solver = effective-variance
""",
    )
    archive_path = tmp_path / "dry.hdf5"
    report_path = tmp_path / "dry.txt"
    result = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(DATA / "PV_quartz.dat"),
            "--spec",
            str(spec_path),
            "--dry-run",
            "--output",
            str(archive_path),
            "--report",
            str(report_path),
            "--quiet",
        ],
    )

    assert result.exit_code == 0, result.output
    assert not archive_path.exists()
    text = report_path.read_text(encoding="utf-8")
    assert "EOS specification:" in text
    assert "EOS batch plan" in text
    assert "BM2" in text
    assert "effective_variance" in text
    assert "no fit was executed" in text


def test_cli_spec_executes_mixed_plan_and_persists_manifest(tmp_path: Path) -> None:
    spec_path = _write_spec(
        tmp_path / "mixed",
        """
[metadata]
title = Mixed topaz plan

[defaults.pv]
model = BM3

[job volume]
domain = pv
targets = volume
model = BM2
solver = effective-variance

[job axes]
domain = pv
targets = a, b
solver = ols
""",
    )
    archive_path = tmp_path / "topaz.hdf5"
    result = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(DATA / "PV_topaz.dat"),
            "--spec",
            str(spec_path),
            "--output",
            str(archive_path),
            "--report",
            str(tmp_path / "topaz.log"),
            "--quiet",
        ],
    )

    assert result.exit_code == 0, result.output
    with EOSArchive(archive_path) as archive:
        assert archive.record_ids == (1, 2, 3)
        assert archive.slot_state("pv/volume").status is EOSSlotStatus.ACCEPTED
        assert archive.slot_state("pv/a").status is EOSSlotStatus.ACCEPTED
        assert archive.slot_state("pv/b").status is EOSSlotStatus.ACCEPTED
        manifest = archive.batch_manifest()
        assert manifest["metadata"]["source"] == "spec"
        assert manifest["metadata"]["title"] == "Mixed topaz plan"
        assert [job["job_id"] for job in manifest["jobs"]] == [
            "volume",
            "axes-a",
            "axes-b",
        ]


def test_cli_spec_rejects_scientific_cli_overrides(tmp_path: Path) -> None:
    spec_path = _write_spec(
        tmp_path / "analysis.spec",
        """
[job volume]
domain = pv
targets = volume
""",
    )
    result = CliRunner().invoke(
        main,
        [
            "eos",
            "run",
            str(DATA / "PV_quartz.dat"),
            "--spec",
            str(spec_path),
            "--solver",
            "ols",
            "--dry-run",
            "--report",
            str(tmp_path / "invalid-override.log"),
        ],
    )

    assert result.exit_code != 0
    assert "--spec is the authority" in result.output
    assert "--solver" in result.output


def test_spec_uses_v0_for_volume_and_l0_for_axes(tmp_path: Path) -> None:
    dataset = read_eos_input(DATA / "rutile.dat")
    document = read_eos_spec(
        _write_spec(
            tmp_path / "vt-names.spec",
            """
[job volume]
domain = vt
targets = volume
model = salje
initial.V0 = 0.9956

[job axis]
domain = vt
targets = a
model = salje
initial.L0 = 0.9987
""",
        )
    )
    resolved = resolve_eos_spec(document, dataset)

    volume_constraint = resolved.plan.jobs[0].request.constraints[0]
    axis_constraint = resolved.plan.jobs[1].request.constraints[0]
    assert volume_constraint.name == "V0"
    assert axis_constraint.name == "L0"

    for name in ("value_ref", "X0"):
        invalid = read_eos_spec(
            _write_spec(
                tmp_path / f"invalid-{name}.spec",
                f"""
[job axis]
domain = vt
targets = a
model = salje
initial.{name} = 1.0
""",
            )
        )
        with pytest.raises(EOSSpecError, match="not available|use V0|use L0"):
            resolve_eos_spec(invalid, dataset)
