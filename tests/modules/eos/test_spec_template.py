"""Generation and CLI delivery of the EOS specification template."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner
import pytest

from quantas.cli.main import main
from quantas.core.math.fitting import FitMethod
from quantas.modules.eos import (
    EOS_SPEC_SIGNATURE,
    EOS_SPEC_TEMPLATE_FILENAME,
    eos_spec_template,
    parse_eos_spec,
    read_eos_input,
    resolve_eos_spec,
    write_eos_spec_template,
)

DATA = Path(__file__).with_name("data")


def test_template_is_complete_commented_and_resolvable() -> None:
    text = eos_spec_template()

    assert text.startswith(EOS_SPEC_SIGNATURE + "\n")
    assert text.endswith("\n")
    for section in (
        "[metadata]",
        "[input]",
        "[batch]",
        "[defaults]",
        "[defaults.pv]",
        "[defaults.vt]",
        "[defaults.pvt]",
        "[presentation]",
        "[job volume-example]",
    ):
        assert section in text
    for option in (
        "pressure_unit",
        "length_unit",
        "temperature_unit",
        "failure_policy",
        "solver",
        "covariance_scaling",
        "max_iterations",
        "ftol",
        "xtol",
        "gtol",
        "inner_max_iterations",
        "odr_difference",
        "odr_ndigit",
        "accept",
        "replace_accepted",
        "domain",
        "targets",
        "model",
        "pv_model",
        "vt_model",
        "coupling",
        "note",
        "fix.KP",
        "initial.K0",
        "bound.K0",
        "detail",
        "show_uncertainties",
        "max_data_rows",
    ):
        assert option in text

    root = Path(__file__).parents[3]
    distributed = root / "examples" / "eos" / "eos.spec"
    download = root / "docs" / "source" / "_downloads" / "eos.spec"
    assert distributed.read_text(encoding="utf-8") == text
    assert download.read_text(encoding="utf-8") == text

    document = parse_eos_spec(text, source="generated-eos.spec")
    dataset = read_eos_input(DATA / "PV_quartz.dat")
    resolved = resolve_eos_spec(document, dataset)

    assert len(resolved.plan.jobs) == 1
    job = resolved.plan.jobs[0]
    assert job.job_id == "volume-example"
    assert job.request.model.tag == "BM3"
    assert job.request.options.method is FitMethod.OLS


def test_write_template_refuses_implicit_overwrite(tmp_path: Path) -> None:
    destination = tmp_path / "analysis.any-extension"

    returned = write_eos_spec_template(destination)
    original = destination.read_text(encoding="utf-8")

    assert returned == destination
    assert original == eos_spec_template()
    with pytest.raises(FileExistsError):
        write_eos_spec_template(destination)

    destination.write_text("obsolete", encoding="utf-8")
    write_eos_spec_template(destination, overwrite=True)
    assert destination.read_text(encoding="utf-8") == original


def test_cli_generates_default_template_file() -> None:
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(main, ["eos", "spec-template"])

        assert result.exit_code == 0, result.output
        destination = Path(EOS_SPEC_TEMPLATE_FILENAME)
        assert destination.is_file()
        assert destination.read_text(encoding="utf-8") == eos_spec_template()
        assert "EOS specification template written to" in result.output
        assert "\x1b[" not in result.output

        report = Path("resolved.txt")
        dry_run = runner.invoke(
            main,
            [
                "eos",
                "run",
                str(DATA / "PV_quartz.dat"),
                "--spec",
                str(destination),
                "--dry-run",
                "--quiet",
                "--report",
                str(report),
            ],
        )
        assert dry_run.exit_code == 0, dry_run.output
        assert not Path("PV_quartz_EOS.hdf5").exists()
        resolved = report.read_text(encoding="utf-8")
        assert "volume-example" in resolved
        assert "BM3" in resolved
        assert "no fit was executed" in resolved


def test_cli_custom_path_requires_force_for_replacement(tmp_path: Path) -> None:
    destination = tmp_path / "eos-template.dat"
    runner = CliRunner()

    first = runner.invoke(main, ["eos", "spec-template", str(destination)])
    second = runner.invoke(main, ["eos", "spec-template", str(destination)])
    destination.write_text("obsolete", encoding="utf-8")
    forced = runner.invoke(
        main,
        ["eos", "spec-template", str(destination), "--force"],
    )

    assert first.exit_code == 0, first.output
    assert second.exit_code != 0
    assert "use --force" in second.output
    assert forced.exit_code == 0, forced.output
    assert destination.read_text(encoding="utf-8") == eos_spec_template()


def test_cli_template_help_documents_output_and_force() -> None:
    result = CliRunner().invoke(main, ["eos", "spec-template", "--help"])

    assert result.exit_code == 0
    assert "OUTPUT" in result.output
    assert "eos.spec" in result.output
    assert "Output and reporting:" in result.output
    assert "--force" in result.output
