"""Regression tests for compact deterministic plain-text reports."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.api import qha, rendering
from quantas.models import InputData, ReportTable, ResultData, ResultMetadata
from quantas.modules.qha.models import QHAResult
from quantas.references.render import render_citation_notice


pytestmark = pytest.mark.infrastructure


def test_multiline_cells_use_visible_width() -> None:
    """A multiline cell must not widen a column by its total character count."""
    table = ReportTable(
        title="Multiline",
        columns=["Property", "Value"],
        rows=[["Status", "first line\nsecond line"], ["Count", 2]],
    )

    text = rendering.render_table(table)

    assert "first line" in text
    assert "second line" in text
    assert max(len(line) for line in text.splitlines()) < 80
    assert all(line == line.rstrip() for line in text.splitlines())


def test_qha_api_report_summarizes_structure() -> None:
    """The public QHA report must not dump the complete structural mapping."""
    structure = {
        "representation": "primitive",
        "orientation": "crystal",
        "reference_index": 6,
        "atomic_numbers": np.array([12, 8]),
        "volume_series": {
            "volume": np.arange(11.0),
            "lattice": np.zeros((11, 3, 3)),
            "fractional_positions": np.zeros((11, 2, 3)),
        },
        "normalization": {
            "basis": "primitive",
            "source_basis": "supercell",
            "repetitions": 8,
        },
        "symmetry": {
            "international_symbol": "Fm-3m",
            "space_group_number": 225,
        },
        "reconstruction": [{"status": "exact"} for _ in range(11)],
    }
    result = ResultData(
        metadata=ResultMetadata(module="qha", method="quasi-harmonic"),
        input_data=InputData(
            source="mgo.yaml",
            data={
                "jobname": "MgO",
                "natoms": 2,
                "has_structure": True,
                "structure": structure,
                "mode_continuity": "assumed",
            },
        ),
        results={"qha": QHAResult(jobname="MgO")},
    )

    text = rendering.render_tables(qha.build_report(result))

    assert "Structure representation" in text
    assert "primitive from supercell (8 repetitions)" in text
    assert "Fm-3m (225)" in text
    assert "fractional_positions" not in text
    assert max(len(line) for line in text.splitlines()) < 160
    assert not any(len(line) > 500 and not line.strip() for line in text.splitlines())


def test_citation_notice_preserves_academic_credit_message() -> None:
    """The standard citation notice must remain readable and typo-free."""
    text = render_citation_notice(["quantas_2022"])

    assert text.startswith("_" * 80)
    assert text.endswith("_" * 80)
    assert "Scientific recognition" in text
    assert "sscientific" not in text
    assert not text.endswith("\n")
