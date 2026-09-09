"""CLI discovery and validation for equation-of-state model tags."""

from __future__ import annotations

import click
from click.testing import CliRunner
import pytest

from quantas.cli.eos import eos
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


def test_show_models_lists_capabilities_and_filters_domains() -> None:
    """The EOS catalogue can be restricted to compatible scientific domains."""
    runner = CliRunner()
    full = runner.invoke(eos, ["show-models"])
    pv = runner.invoke(eos, ["show-models", "--domain", "pv"])
    ev = runner.invoke(eos, ["show-models", "--domain", "ev"])
    common = runner.invoke(
        eos,
        ["show-models", "--domain", "pv", "--domain", "ev"],
    )

    assert full.exit_code == pv.exit_code == ev.exit_code == common.exit_code == 0
    assert "Available equations of state" in full.output
    assert "P-V fit" in full.output
    assert "Stabilized jellium (SJEOS)" in full.output
    assert " SJ " not in pv.output
    assert " SJ " in ev.output
    assert " SJ " not in common.output
    assert " BM3 " in pv.output
    assert " BM3 " in ev.output
    assert " BM3 " in common.output


def test_qha_eos_help_uses_compact_model_metavar() -> None:
    """QHA help no longer renders the complete EOS catalogue inline."""
    result = CliRunner().invoke(qha, ["run", "--help"])
    assert result.exit_code == 0
    assert "--eos MODEL" in result.output
    assert "BM2|BM3|BM4" not in result.output
