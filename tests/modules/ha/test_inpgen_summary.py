"""Tests for concise HA/QHA input-generation summaries."""

from __future__ import annotations

from quantas.modules.ha.io.inpgen import (
    phonon_input_summary_table,
    phonon_qpoint_sampling_table,
)


def _sample_input() -> dict:
    """Return a compact normalized phonon input mapping for report tests."""
    return {
        "natom": 2,
        "formula_units": 1,
        "qpoints": 2,
        "q_position_source": "crystal-dispersion-table",
        "volume": [18.1, 18.5],
        "phonon": [
            {
                "q-position": [0.0, 0.0, 0.0],
                "weight": 1.0,
                "band": [{"frequency": [0.0, 0.0]}] * 6,
            },
            {
                "q-position": [0.0, 0.0, 0.5],
                "weight": 2.0,
                "band": [{"frequency": [1.0, 2.0]}] * 6,
            },
        ],
        "mode_continuity": "verified",
        "provenance": {
            "energy": {
                "selected_quantity": "total_energy",
                "states": [
                    {"corrections": ["DISP"]},
                    {"corrections": ["DISP"]},
                ],
            }
        },
    }


def test_phonon_input_summary_exposes_energy_and_qpoint_provenance() -> None:
    """The standard summary reports the scientifically selected input data."""
    table = phonon_input_summary_table(
        _sample_input(),
        interface="crystal",
        source_count=2,
        eigenvectors_available=True,
    )

    values = {str(row[0]): row[1] for row in table.rows}
    assert values["Q-position source"] == "crystal-dispersion-table"
    assert values["Energy quantity"] == "total_energy"
    assert values["Energy corrections"] == "DISP"
    assert values["Formula units"] == 1


def test_phonon_qpoint_sampling_table_is_compact() -> None:
    """The q-point preview preserves coordinates and truncates long meshes."""
    data = _sample_input()
    data["phonon"] = data["phonon"] * 8

    table = phonon_qpoint_sampling_table(data, max_rows=3)

    assert table is not None
    assert len(table.rows) == 3
    assert table.rows[0] == [1, 0.0, 0.0, 0.0, 1.0]
    assert table.rows[1] == [2, 0.0, 0.0, 0.5, 2.0]
    assert "13 additional q-point(s)" in table.metadata["notes"][1]


def test_phonon_qpoint_sampling_table_skips_unavailable_positions() -> None:
    """No coordinate table is emitted when the source has no reliable q labels."""
    data = _sample_input()
    for qpoint in data["phonon"]:
        qpoint["q-position"] = None

    assert phonon_qpoint_sampling_table(data) is None
