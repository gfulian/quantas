"""CLI discovery and validation for equation-of-state model tags."""

from __future__ import annotations

import click
from click.testing import CliRunner
import pytest

from quantas.cli.eos import eos
from quantas.cli.main import main
from quantas.cli.shell_completion import PowerShellComplete
from quantas.cli.eos_model_type import ENERGY_EOS_MODEL, EOSModelParamType
from quantas.cli.qha import qha


def test_eos_type_normalizes_compact_and_long_aliases() -> None:
    """EOS-valued options normalize aliases to stable compact tags."""
    assert ENERGY_EOS_MODEL.convert("bm3", None, None) == "BM3"
    assert ENERGY_EOS_MODEL.convert("birch-murnaghan3", None, None) == "BM3"
    assert ENERGY_EOS_MODEL.convert("NS4", None, None) == "PT4"
    assert ENERGY_EOS_MODEL.convert("natural-strain4", None, None) == "PT4"
    assert ENERGY_EOS_MODEL.convert("SJEOS", None, None) == "SJ"


def test_eos_type_reports_actionable_model_errors() -> None:
    """Invalid tags report family-aware recovery guidance."""
    with pytest.raises(click.BadParameter, match=r"Available models: BM2, BM3, BM4"):
        ENERGY_EOS_MODEL.convert("BM5", None, None)
    with pytest.raises(click.BadParameter, match=r"does not define selectable EOS orders"):
        ENERGY_EOS_MODEL.convert("SJ3", None, None)
    with pytest.raises(click.BadParameter, match=r"quantas eos show-models"):
        ENERGY_EOS_MODEL.convert("unknown", None, None)


def test_pressure_type_rejects_energy_only_sjeos() -> None:
    """Capability validation distinguishes E(V) from direct P-V fitting."""
    pressure_type = EOSModelParamType(require_pressure_fit=True)
    with pytest.raises(click.BadParameter, match=r"not available for direct P-V fitting"):
        pressure_type.convert("SJ", None, None)


def test_eos_type_completion_exposes_historical_aliases() -> None:
    """Shell completion includes concise canonical and historical tags."""
    ns_items = ENERGY_EOS_MODEL.shell_complete(None, None, "NS")  # type: ignore[arg-type]
    sj_items = ENERGY_EOS_MODEL.shell_complete(None, None, "SJ")  # type: ignore[arg-type]
    assert [item.value for item in ns_items] == ["NS", "NS2", "NS3", "NS4"]
    assert [item.value for item in sj_items] == ["SJ", "SJEOS"]
    assert all(item.help for item in ns_items + sj_items)


def test_show_models_lists_all_domains_and_selected_sections() -> None:
    """EOS discovery shows every domain and can select several sections."""
    runner = CliRunner()
    full = runner.invoke(eos, ["show-models"])
    pv = runner.invoke(eos, ["show-models", "--domain", "pv"])
    ev_vt = runner.invoke(
        eos,
        ["show-models", "--domain", "ev", "--domain", "vt"],
    )
    pvt = runner.invoke(eos, ["show-models", "--domain", "pvt"])

    assert full.exit_code == pv.exit_code == ev_vt.exit_code == pvt.exit_code == 0
    for domain in ("P-V", "E-V", "V-T", "P-V-T"):
        assert domain in full.output
    assert "E-V integrated energy models" in full.output
    assert "V-T thermal-expansion models" in full.output
    assert "P-V-T coupling models" in full.output
    assert "P-V-T thermal-pressure components" in full.output
    assert " SJ " not in pv.output
    assert " SJ " in ev_vt.output
    assert "BERMAN:linear" in ev_vt.output
    assert "P-V isothermal models" not in ev_vt.output
    assert "linear-bulk-modulus" in pvt.output
    assert "mie-gruneisen-debye:full" in pvt.output


def test_powershell_completion_backend_is_registered_and_rendered() -> None:
    """PowerShell receives a native Click-backed argument completer."""
    from click.shell_completion import get_completion_class

    assert get_completion_class("powershell") is PowerShellComplete
    result = CliRunner().invoke(main, ["completion", "powershell"])
    assert result.exit_code == 0, result.output
    assert "Register-ArgumentCompleter -Native -CommandName" in result.output
    assert "powershell_complete" in result.output
    assert "CompletionCompleters]::CompleteFilename" in result.output

    completer = PowerShellComplete(main, {}, "quantas", "_QUANTAS_COMPLETE")
    items = completer.get_completions(
        ["qha", "run", "input.yaml", "--eos"],
        "SJ",
    )
    assert [item.value for item in items] == ["SJ", "SJEOS"]

    completion = CliRunner().invoke(
        main,
        [],
        prog_name="quantas",
        env={
            "_QUANTAS_COMPLETE": "powershell_complete",
            "QUANTAS_COMPLETE_ARGS": (
                '["quantas","qha","run","input.yaml","--eos"]'
            ),
            "QUANTAS_COMPLETE_WORD": "SJ",
        },
    )
    assert completion.exit_code == 0, completion.output
    assert '"value":"SJ"' in completion.output
    assert '"value":"SJEOS"' in completion.output


def test_qha_eos_help_uses_compact_model_metavar() -> None:
    """QHA help no longer renders the complete EOS catalogue inline."""
    result = CliRunner().invoke(qha, ["run", "--help"])
    assert result.exit_code == 0
    assert "--eos MODEL" in result.output
    assert "BM2|BM3|BM4" not in result.output
