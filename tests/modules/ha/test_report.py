"""Tests for the HA report."""

from __future__ import annotations

from pathlib import Path

from quantas.modules.ha.api import read_ha_input, run_ha
from quantas.modules.ha.models import HAOptions
from quantas.modules.ha.report import (
    all_tables,
    input_table,
    options_table,
    static_data_table,
    thermodynamic_summary_table,
    thermodynamic_table,
)


DATA = Path(__file__).parent / "data/mgo_b3lyp_qha.yaml"


def test_report_tables_are_frontend_neutral_objects():
    input_data = read_ha_input(DATA)
    options = HAOptions(
        temperature_min=0.0, temperature_max=100.0, temperature_step=100.0
    )
    result = run_ha(input_data, options=options).results["ha"]

    table = input_table(input_data)
    assert table.title == "Input data"
    assert table.columns == ["Property", "Value"]
    assert any(row[0] == "Job name" for row in table.rows)
    assert ["Formula units per cell (Z)", 1] in table.rows
    assert ["Atoms per formula unit", 2.0] in table.rows

    opts = options_table(options)
    assert opts.title == "Selected options"
    assert any(row == ["Energy unit", options.energy_unit] for row in opts.rows)

    summary = thermodynamic_summary_table(result)
    assert summary.title == "Thermodynamic properties"
    assert any(row[0] == "F" and row[2] == "yes" for row in summary.rows)

    tables = all_tables(input_data, options, result)
    assert [item.title for item in tables] == [
        "Input data",
        "Selected options",
        "Static and zero-point energies",
        "Thermodynamic properties",
        "Thermal energy",
        "Internal energy",
        "Entropy",
        "Vibrational Helmholtz free energy",
        "Helmholtz free energy",
        "Isochoric heat capacity",
    ]


def test_thermodynamic_table_accepts_symbol_and_can_truncate():
    input_data = read_ha_input(DATA)
    options = HAOptions(
        temperature_min=0.0, temperature_max=200.0, temperature_step=100.0
    )
    result = run_ha(input_data, options=options).results["ha"]

    table = thermodynamic_table(result, "F", max_rows=2)

    assert table.title == "Helmholtz free energy"
    assert table.columns[0] == "T"
    assert len(table.columns) == result.volume.shape[0] + 1
    assert len(table.rows) == 2
    assert table.metadata["symbol"] == "F"
    assert table.metadata["truncated"] is True


def test_multivolume_zero_point_table_is_volume_resolved():
    input_data = read_ha_input(DATA)
    options = HAOptions(
        temperature_min=0.0, temperature_max=200.0, temperature_step=100.0
    )
    result = run_ha(input_data, options=options).results["ha"]

    table = thermodynamic_table(result, "Uzp")

    assert table.columns == ["V", "Uzp"]
    assert len(table.rows) == input_data.nvol
    assert table.metadata["independent_variable"] == "volume"
    assert table.metadata["temperature_independent"] is True
    assert table.metadata["total_rows"] == input_data.nvol
    assert table.rows[0][0] == result.volume[0]
    assert table.rows[-1][0] == result.volume[-1]


def test_multivolume_report_preserves_temperature_volume_matrix():
    input_data = read_ha_input(DATA)
    options = HAOptions(
        temperature_min=0.0, temperature_max=200.0, temperature_step=100.0
    )
    result = run_ha(input_data, options=options).results["ha"]

    static = static_data_table(result)
    free = thermodynamic_table(result, "F", row_indices=(0, -1))

    assert static.columns == ["V", "U0", "Uzp"]
    assert len(static.rows) == input_data.nvol
    assert free.columns[0] == "T"
    assert len(free.columns) == input_data.nvol + 1
    assert len(free.rows) == 2
    assert all(len(row) == input_data.nvol + 1 for row in free.rows)
