"""Structural checks for the public interoperability documentation."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / "docs/source/workflows/interoperability.rst"
API_PAGE = ROOT / "docs/source/api/interoperability.rst"
EXAMPLE = ROOT / "examples/interoperability/workflow_api.py"
CLI_EXAMPLE = ROOT / "examples/interoperability/workflow_cli.sh"
DOWNLOAD_DIR = ROOT / "docs/source/_downloads/interoperability"


def test_interoperability_workflow_is_complete() -> None:
    """The workflow page must cover contracts, both frontends, and failures."""
    text = WORKFLOW.read_text(encoding="utf-8")
    required = (
        "QHA to Thermoelasticity",
        "Thermoelasticity to one material state",
        "State to Elasticity",
        "State to SEISMIC",
        "Complete CLI example",
        "Complete Python API example",
        "CLI/API equivalence",
        "Provenance and persistence",
        "Failure modes and diagnostics",
        "Performance and reuse",
        "Boundaries of the current implementation",
        "qha_to_thermoelastic_context",
        "thermoelastic_to_seismic",
        "--tensor-condition adiabatic",
        "--extrapolation fail",
    )
    for item in required:
        assert item in text
    assert "Work in progress" not in text


def test_examples_use_only_public_api() -> None:
    """The Python example must import through quantas.api only."""
    text = EXAMPLE.read_text(encoding="utf-8")
    assert "from quantas.api import" in text
    assert "quantas.modules" not in text
    assert "quantas.core" not in text
    assert "quantas.models" not in text


def test_downloaded_examples_match_curated_examples() -> None:
    """Documentation downloads must remain byte-identical to examples."""
    assert (DOWNLOAD_DIR / EXAMPLE.name).read_bytes() == EXAMPLE.read_bytes()
    assert (DOWNLOAD_DIR / CLI_EXAMPLE.name).read_bytes() == CLI_EXAMPLE.read_bytes()


def test_api_reference_links_to_complete_example() -> None:
    """The API page must direct readers to the workflow and executable script."""
    text = API_PAGE.read_text(encoding="utf-8")
    assert ":doc:`../workflows/interoperability`" in text
    assert "workflow_api.py" in text
