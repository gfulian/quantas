"""Structural checks for the EOS implementation/workflow chapter."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PAGE = ROOT / "docs" / "source" / "workflows" / "eos.rst"


def _text() -> str:
    return PAGE.read_text(encoding="utf-8")


def test_eos_workflow_archive_semantics() -> None:
    text = _text()
    for phrase in (
        "Why EOS uses a different command-line workflow",
        "immutable fit attempts",
        "accepted record",
        "candidate record",
        "replace_accepted",
        "Persistent archives and immutable history",
    ):
        assert phrase in text


def test_eos_workflow_domains_and_batch() -> None:
    text = _text()
    for phrase in (
        "P--V implementation choices",
        "V--T implementation choices",
        "P--V--T implementation choices",
        "Ordinary least squares (OLS)",
        "Weighted least squares (WLS)",
        "Iterative effective variance",
        "Orthogonal distance regression (ODR)",
        "QUANTAS EOS SPEC 1",
        "--dry-run",
    ):
        assert phrase in text


def test_eos_workflow_postfit_limits() -> None:
    text = _text()
    for phrase in (
        "Post-fit diagnostics",
        "Property calculation",
        "Plotting implementation",
        "Fit success versus scientific adequacy",
        "Failed fits and last iterates",
        "EOS has no batching control analogous to SEISMIC",
        "Common interpretation errors",
    ):
        assert phrase in text
