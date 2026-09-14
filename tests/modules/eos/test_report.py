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


def test_text_presentation_uses_portable_ascii_scientific_notation() -> None:
    """Terminal/report labels and units avoid Unicode-dependent notation."""
    from quantas.modules.eos.presentation import (
        domain_label,
        format_text_unit,
        parameter_label,
    )

    assert format_text_unit("angstrom^3") == "angstrom^3"
    assert format_text_unit("Å") == "angstrom"
    assert format_text_unit("GPa^-1") == "GPa^-1"
    assert format_text_unit("(GPa)^-1") == "GPa^-1"
    assert format_text_unit("GPa⁻¹") == "GPa^-1"
    assert format_text_unit("J mol⁻¹") == "J/mol"
    assert parameter_label("KP") == "K'"
    assert parameter_label("KPP") == "K''"
    assert parameter_label("alpha0") == "alpha0"
    assert parameter_label("delta") == "delta"
    assert domain_label("ev") == "Energy-volume"
    values = (
        format_text_unit("Å"),
        format_text_unit("GPa⁻¹"),
        parameter_label("KP"),
        parameter_label("alpha0"),
        domain_label("pvt"),
    )
    assert all(value is None or value.isascii() for value in values)


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


def test_non_ev_configuration_omits_secondary_axial_formulation() -> None:
    """Secondary axial settings are reported only for the E-V domain."""
    from quantas.modules.eos.batch import EOSBatchJob, EOSBatchJobResult
    from quantas.modules.eos.report import (
        eos_job_configuration_table,
        eos_requested_fit_table,
    )

    dataset = read_eos_input(_DATA / "PV_quartz.dat")
    request = EOSFitRequest(
        model="BM3",
        domain="pv",
        target="volume",
        options=EOSFitOptions(solver_options=OLSOptions()),
    )
    result = EOSFitter().fit(dataset, request)
    job = EOSBatchJob(request=request, job_id="pv-volume")
    completed = EOSBatchJobResult(
        job_id="pv-volume",
        request=request,
        result=result,
        record_id=1,
        accepted=True,
    )

    requested = eos_requested_fit_table(job, 1)
    configured = eos_job_configuration_table(completed)
    assert "Secondary axial formulation" not in {row[0] for row in requested.rows}
    assert "Secondary axial formulation" not in {row[0] for row in configured.rows}


def test_eos_run_help_uses_ev_specific_model_option() -> None:
    """The E-V model option is namespaced consistently with other EOS domains."""
    from click.testing import CliRunner

    from quantas.cli.main import main

    result = CliRunner().invoke(main, ["eos", "run", "--help"])

    assert result.exit_code == 0
    assert "--ev-eos MODEL" in result.output
    assert "--eos MODEL" not in result.output
    assert "secondary axial P(l^3) fits" in result.output
    assert "Angel-style" not in result.output
