"""Declarative EOS batch workflow and frontend-neutral reports."""

from __future__ import annotations

from pathlib import Path

import pytest

from quantas.core.math.fitting import OLSOptions, WLSOptions
from quantas.modules.eos import (
    EOSArchive,
    EOSBatchFailurePolicy,
    EOSBatchJob,
    EOSBatchPlan,
    EOSBatchWorkflow,
    EOSFitDomain,
    EOSFitOptions,
    EOSFitRequest,
    EOSReportDetail,
    EOSReportOptions,
    EOSSlotStatus,
    build_eos_batch_report,
    read_eos_input,
)

DATA = Path(__file__).with_name("data")


def _request(
    *,
    target: str = "volume",
    domain: EOSFitDomain = EOSFitDomain.PRESSURE_VOLUME,
    solver_options: OLSOptions | WLSOptions | None = None,
) -> EOSFitRequest:
    model = "BM3" if domain is EOSFitDomain.PRESSURE_VOLUME else "berman:quadratic"
    return EOSFitRequest(
        model=model,
        target=target,
        domain=domain,
        options=EOSFitOptions(solver_options=solver_options or OLSOptions()),
    )


def test_empty_plan_creates_self_contained_archive_and_manifest(tmp_path: Path) -> None:
    destination = tmp_path / "empty.hdf5"
    result = EOSBatchWorkflow().run(
        DATA / "PV_topaz.dat",
        EOSBatchPlan(metadata={"purpose": "input-only"}),
        destination,
    )

    assert result.successful
    assert result.jobs == ()
    with EOSArchive(destination) as archive:
        assert archive.record_ids == ()
        assert archive.batch_manifest() == {
            "failure_policy": "stop",
            "jobs": [],
            "metadata": {"purpose": "input-only"},
        }
        states = {state.slot.key: state.status for state in archive.slots()}
        assert states == {
            "pv/a": EOSSlotStatus.NOT_PROCESSED,
            "pv/b": EOSSlotStatus.NOT_PROCESSED,
            "pv/c": EOSSlotStatus.NOT_PROCESSED,
            "pv/volume": EOSSlotStatus.NOT_PROCESSED,
        }


def test_ordered_jobs_are_persisted_and_accepted_explicitly(tmp_path: Path) -> None:
    destination = tmp_path / "topaz.hdf5"
    plan = EOSBatchPlan(
        jobs=(
            EOSBatchJob(_request(target="volume"), job_id="volume"),
            EOSBatchJob(_request(target="b"), job_id="axis-b"),
        )
    )

    result = EOSBatchWorkflow().run(DATA / "PV_topaz.dat", plan, destination)

    assert result.successful
    assert [job.record_id for job in result.jobs] == [1, 2]
    assert all(job.accepted for job in result.jobs)
    with EOSArchive(destination) as archive:
        assert archive.slot_state("pv/volume").status is EOSSlotStatus.ACCEPTED
        assert archive.slot_state("pv/b").status is EOSSlotStatus.ACCEPTED
        assert archive.slot_state("pv/a").status is EOSSlotStatus.NOT_PROCESSED
        assert archive.accepted("pv/volume").record_id == 1
        assert archive.accepted("pv/b").record_id == 2


def test_duplicate_accepted_slot_requires_explicit_replacement() -> None:
    first = EOSBatchJob(_request(), job_id="bm3")
    second = EOSBatchJob(_request(), job_id="repeat")
    with pytest.raises(ValueError, match="accepts more than one result"):
        EOSBatchPlan(jobs=(first, second))

    plan = EOSBatchPlan(
        jobs=(first, EOSBatchJob(_request(), job_id="repeat", replace_accepted=True))
    )
    assert len(plan.jobs) == 2


def test_continue_policy_persists_failed_attempt_and_later_success(
    tmp_path: Path,
) -> None:
    destination = tmp_path / "continue.hdf5"
    failed = EOSBatchJob(
        _request(
            domain=EOSFitDomain.VOLUME_TEMPERATURE,
            solver_options=WLSOptions(),
        ),
        job_id="weighted-without-sigma",
        accept=False,
    )
    success = EOSBatchJob(
        _request(domain=EOSFitDomain.VOLUME_TEMPERATURE),
        job_id="unweighted",
    )
    plan = EOSBatchPlan(
        jobs=(failed, success),
        failure_policy=EOSBatchFailurePolicy.CONTINUE,
    )

    result = EOSBatchWorkflow().run(DATA / "T_triclinic.dat", plan, destination)

    assert len(result.jobs) == 2
    assert not result.jobs[0].successful
    assert result.jobs[1].successful
    with EOSArchive(destination) as archive:
        assert archive.record_ids == (1, 2)
        state = archive.slot_state("vt/volume")
        assert state.accepted_record_id == 2
        assert state.attempted_record_ids == (1, 2)


def test_reports_separate_data_uncertainties_and_extended_diagnostics(
    tmp_path: Path,
) -> None:
    destination = tmp_path / "quartz.hdf5"
    plan = EOSBatchPlan(
        jobs=(EOSBatchJob(_request(solver_options=WLSOptions()), job_id="quartz"),)
    )
    result = EOSBatchWorkflow().run(DATA / "PV_quartz.dat", plan, destination)

    tables = build_eos_batch_report(
        result,
        EOSReportOptions(
            detail=EOSReportDetail.EXTENDED,
            show_uncertainties=True,
            max_data_rows=3,
        ),
    )

    by_title = {table.title: table for table in tables}
    assert by_title["EOS input data"].columns == ["Pressure", "Volume"]
    assert len(by_title["EOS input data"].rows) == 3
    assert by_title["EOS input standard uncertainties"].columns[0] == "σ(Pressure)"
    assert "Observed and calculated EOS data" in by_title
    assert "Parameter covariance" in by_title
    assert "Parameter correlation" in by_title
    assert "EOS batch fit summary" in by_title
    assert "EOS batch parameter summary" in by_title


def test_explicit_unit_overrides_normalize_values_and_preserve_raw_data(
    tmp_path: Path,
) -> None:
    source = tmp_path / "units.anything"
    source.write_text(
        "FORMAT P T V A\n10.0 25.0 1.0 2.0\n",
        encoding="utf-8",
    )

    dataset = read_eos_input(
        source,
        pressure_unit="kbar",
        temperature_unit="C",
        length_unit="bohr",
    )

    assert dataset.column("pressure")[0] == pytest.approx(1.0)
    assert dataset.column("temperature")[0] == pytest.approx(298.15)
    assert dataset.column("a")[0] == pytest.approx(1.0583544218)
    assert dataset.column("volume")[0] == pytest.approx(0.1481847115)
    assert dataset.raw_columns["pressure"][0] == pytest.approx(10.0)
    assert dataset.raw_units["volume"] == "bohr^3"
    assert dataset.units["volume"] == "angstrom^3"
