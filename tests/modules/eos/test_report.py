"""Tests for stable EOS names and deterministic report formatting."""

from __future__ import annotations

from pathlib import Path

from quantas.core.math.fitting import OLSOptions
from quantas.modules.eos import EOSFitOptions, EOSFitRequest, EOSFitter, read_eos_input
from quantas.modules.eos.report import (
    eos_data_table,
    eos_parameter_table,
    eos_uncertainty_table,
)
from quantas.renderers.tables import format_numeric, render_table

_DATA = Path(__file__).with_name("data")


def test_eos_numeric_profiles_preserve_small_values() -> None:
    assert format_numeric(0.0001, "eos_pressure") == "0.0001"
    assert format_numeric(303.0, "eos_temperature") == "303.00"
    assert format_numeric(1.0, "eos_structural") == "1.000000"
    assert format_numeric(1.0e-6, "eos_pressure_uncertainty") == "1.000000e-06"
    assert format_numeric(3.0, "eos_temperature_uncertainty") == "3.00"
    assert format_numeric(9.957938e-6, "eos_parameter") == "9.957938e-06"
    assert format_numeric(188.4286764, "eos_parameter") == "188.428676"
    assert format_numeric(1.3628324e-14, "eos_covariance") == "1.362832e-14"
    assert format_numeric(0.79516817, "eos_correlation") == "0.795168"


def test_input_tables_use_quantity_specific_profiles() -> None:
    dataset = read_eos_input(_DATA / "rutile.dat")
    data_text = render_table(eos_data_table(dataset, max_rows=1))
    sigma_text = render_table(eos_uncertainty_table(dataset, max_rows=1))

    assert "0.0001" in data_text
    assert "303.00" in data_text
    assert "1.000000" in data_text
    assert "1.000000e-06" in sigma_text
    assert "3.00" in sigma_text
    assert "0.000500" in sigma_text


def test_vt_reports_use_v0_or_physical_l0_and_clear_delta_label() -> None:
    dataset = read_eos_input(_DATA / "rutile.dat")
    options = EOSFitOptions(solver_options=OLSOptions())

    volume = EOSFitter().fit(
        dataset,
        EOSFitRequest(model="salje", domain="vt", target="volume", options=options),
    )
    axis = EOSFitter().fit(
        dataset,
        EOSFitRequest(model="salje", domain="vt", target="a", options=options),
    )

    volume_text = render_table(eos_parameter_table(volume))
    axis_text = render_table(eos_parameter_table(axis))

    assert "V0" in volume_text
    assert "L0" not in volume_text
    assert "L0" in axis_text
    assert "V0" not in axis_text
    assert "Final - initial" in volume_text
    assert "Shift" not in volume_text
    assert "value_ref" not in volume_text + axis_text
    assert "X0" not in volume_text + axis_text
    assert "e-06" in volume_text + axis_text
