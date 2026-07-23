"""Stable EOS application facade, capability matrix, and report summaries."""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from quantas import api as quantas_api
from quantas.api import eos
from quantas.modules.eos import MODULE_CONTRACT
from quantas.modules.eos import (
    EOSBatchJob,
    EOSBatchPlan,
    EOSBatchWorkflow,
    EOSReportOptions,
    build_eos_batch_result_tables,
)
from quantas.modules.eos.domains.pv import PressureEOSFitModel as DomainPVModel
from quantas.modules.eos.fitting import PressureEOSFitModel as FacadePVModel


DATA = Path(__file__).with_name("data")


def test_public_facade_declares_ev_as_core_only() -> None:
    assert MODULE_CONTRACT.name == "eos"
    assert MODULE_CONTRACT.archive_schema_version == "1.1"
    assert MODULE_CONTRACT.capability("pv").fitting

    energy = eos.domain_capability("ev")
    assert energy.status is eos.CapabilityStatus.CORE_ONLY
    assert not energy.fitting
    assert "QHA" in energy.note


def test_public_eos_namespace_is_the_single_application_facade() -> None:
    assert quantas_api.eos is eos
    assert eos.fit is not None
    assert eos.open_archive is not None


def test_public_eos_solver_and_pvt_contracts_are_self_contained() -> None:
    """Advanced requests can be constructed without implementation imports."""
    options = eos.FitOptions(
        solver_options=eos.default_solver_options(eos.FitMethod.WLS)
    )
    normalization = eos.MGDNormalization.cell(
        formula="NaF",
        formula_units_per_cell=4.0,
    )
    model = eos.PVTModel(
        pressure_model="BM3",
        coupling="thermal-pressure",
        thermal_pressure_model="mgd",
        mgd_normalization=normalization,
    )

    assert isinstance(options.solver_options, eos.WLSOptions)
    assert model.mgd_normalization is normalization
    assert "SolverOptions | None" in str(inspect.signature(eos.FitOptions))
    assert "FitOptions | None" not in str(inspect.signature(eos.FitOptions))


def test_fitting_facade_reexports_domain_model() -> None:
    assert FacadePVModel is DomainPVModel


def test_stable_fit_facade_executes_pressure_volume_request() -> None:
    request = eos.FitRequest(
        model="BM3",
        domain="pv",
        target="volume",
        options=eos.FitOptions(solver_options=eos.OLSOptions()),
    )
    result = eos.fit(DATA / "PV_quartz.dat", request)

    assert result.fit.success
    values = dict(zip(result.fit.parameter_names, result.fit.parameters, strict=True))
    assert values["V0"] == pytest.approx(112.98, rel=1.0e-3)


def test_batch_summary_is_sorted_and_reports_human_units(tmp_path: Path) -> None:
    jobs = (
        EOSBatchJob(
            eos.FitRequest(
                model="BM3",
                domain="pv",
                target="volume",
                options=eos.FitOptions(solver_options=eos.OLSOptions()),
            ),
            job_id="volume",
        ),
        EOSBatchJob(
            eos.FitRequest(
                model="BM3",
                domain="pv",
                target="a",
                options=eos.FitOptions(solver_options=eos.OLSOptions()),
            ),
            job_id="axis-a",
        ),
    )
    result = EOSBatchWorkflow().run(
        DATA / "PV_topaz.dat",
        EOSBatchPlan(jobs=jobs),
        tmp_path / "topaz.hdf5",
    )

    tables = build_eos_batch_result_tables(result, EOSReportOptions())
    by_title = {table.title: table for table in tables}
    fit_summary = by_title["EOS batch fit summary"]
    parameter_summary = by_title["EOS batch parameter summary"]

    assert fit_summary.columns[:4] == [
        "Domain",
        "Quantity",
        "Formulation",
        "Solver",
    ]
    assert [row[1] for row in fit_summary.rows] == [
        "Cell parameter a",
        "Volume",
    ]
    assert "BM3" in fit_summary.rows[0][2]
    assert fit_summary.rows[0][3] == "Ordinary least squares"
    assert any(row[-1] == "Å" for row in parameter_summary.rows)
    assert any(row[-1] == "GPa" for row in parameter_summary.rows)
