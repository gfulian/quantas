"""Archive inspection and frontend-neutral EOS session controller."""

from __future__ import annotations

from pathlib import Path
from shutil import copyfile

import h5py
import numpy as np
import pytest

from quantas.core.events import EventLevel, ListObserver
from quantas.core.math.fitting import OLSOptions, WLSOptions
from quantas.modules.eos import (
    EOSArchive,
    EOSDataset,
    EOSFitter,
    EOSFitOptions,
    EOSFitRequest,
    EOSFitResult,
    EOSRecordDisposition,
    EOSSession,
    EOSSlotStatus,
    EOSStateEventType,
    read_eos_input,
)

DATA = Path(__file__).with_name("data")


class _RaisingFitter(EOSFitter):
    """Test double that raises one exception before numerical execution."""

    def __init__(self, exc: Exception) -> None:
        self.exc = exc

    def fit(self, dataset: EOSDataset, request: EOSFitRequest) -> EOSFitResult:
        raise self.exc


def _request(model: str) -> EOSFitRequest:
    return EOSFitRequest(
        model=model,
        options=EOSFitOptions(solver_options=WLSOptions()),
    )


def test_session_stores_attempt_before_decision_and_resumes(tmp_path: Path) -> None:
    path = tmp_path / "interactive.hdf5"
    source = tmp_path / "PV_quartz.dat"
    copyfile(DATA / "PV_quartz.dat", source)
    with EOSSession.create(source, path) as session:
        first = session.fit(_request("BM3"), note="first inspection")
        assert first.record_id == 1
        assert session.inspect_record(first.record_id).disposition is (
            EOSRecordDisposition.UNDECIDED
        )
        assert session.archive.slot_state("pv/volume").status is EOSSlotStatus.ATTEMPTED
        assert session.archive.accepted("pv/volume") is None

    source.unlink()
    with EOSSession.resume(path) as session:
        assert session.inspect_record(1).record.result.parameter_values
        session.accept(1, note="reference fit")
        second = session.fit(_request("BM2"))
        assert second.record_id == 2
        assert session.archive.event_ids == tuple(
            range(1, len(session.archive.event_ids) + 1)
        )
        assert session.inspect_record(1).disposition is EOSRecordDisposition.ACCEPTED
        assert session.archive.slot_state("pv/volume").accepted_record_id == 1


def test_session_reuses_compatible_values_by_name_and_records_lineage(
    tmp_path: Path,
) -> None:
    path = tmp_path / "lineage.hdf5"
    with EOSSession.create(DATA / "PV_quartz.dat", path) as session:
        first = session.fit(_request("BM3"))
        session.accept(first.record_id)
        second = session.fit(
            _request("BM2"),
            initial_guess_record_id=first.record_id,
        )
        copied = {item.name: item.initial_value for item in second.request.constraints}
        assert copied["V0"] == pytest.approx(first.result.parameter_values["V0"])
        assert copied["K0"] == pytest.approx(first.result.parameter_values["K0"])
        assert "KP" not in copied
        assert second.parent_record_id == first.record_id
        source_view = session.inspect_record(first.record_id)
        assert source_view.used_as_initial_guess
        assert source_view.child_record_ids == (second.record_id,)
        assert EOSStateEventType.RECORD_USED_AS_INITIAL_GUESS in {
            event.event_type for event in source_view.events
        }

        third = session.fit(
            _request("BM3"),
            initial_guess_record_id=second.record_id,
        )
        third_guesses = {
            item.name: item.initial_value for item in third.request.constraints
        }
        assert third_guesses["KP"] == pytest.approx(
            second.result.parameter_values["KP"]
        )


def test_acceptance_can_be_revoked_without_erasing_history(tmp_path: Path) -> None:
    path = tmp_path / "decisions.hdf5"
    with EOSSession.create(DATA / "PV_quartz.dat", path) as session:
        first = session.fit(_request("BM3"))
        second = session.fit(_request("BM2"))
        session.accept(first.record_id)
        with pytest.raises(ValueError, match="revoke acceptance"):
            session.reject(first.record_id)

        session.accept(second.record_id)
        assert session.inspect_record(first.record_id).disposition is (
            EOSRecordDisposition.SUPERSEDED
        )
        assert session.inspect_record(second.record_id).disposition is (
            EOSRecordDisposition.ACCEPTED
        )

        state = session.unaccept("pv/volume", note="model comparison reopened")
        assert state.status is EOSSlotStatus.ATTEMPTED
        assert state.accepted_record_id is None
        assert session.inspect_record(second.record_id).disposition is (
            EOSRecordDisposition.WITHDRAWN
        )
        session.reject(second.record_id, note="residual trend")
        assert session.inspect_record(second.record_id).disposition is (
            EOSRecordDisposition.REJECTED
        )

    with h5py.File(path, "r") as h5:
        assert "results/accepted/pv/volume" not in h5
        assert "session/records/000001" in h5
        assert "session/records/000002" in h5


def test_archive_inspection_reports_slots_records_and_size(tmp_path: Path) -> None:
    path = tmp_path / "inspection.hdf5"
    with EOSSession.create(DATA / "PV_quartz.dat", path) as session:
        record = session.fit(_request("BM3"))
        session.mark_candidate(record.record_id, note="compare later")
        inspection = session.inspect()
        slot = next(
            item for item in inspection.slots if item.state.slot.key == "pv/volume"
        )
        assert inspection.record_count == 1
        assert inspection.event_count >= 2
        assert inspection.size.size_bytes > 0
        assert slot.last is not None
        assert slot.last.disposition is EOSRecordDisposition.CANDIDATE
        assert "compare later" in slot.last.notes


def test_session_persists_numerical_validation_failures(tmp_path: Path) -> None:
    dataset = EOSDataset(
        jobname="unweighted synthetic data",
        columns={
            "pressure": np.linspace(0.0, 5.0, 8),
            "volume": np.linspace(100.0, 96.0, 8),
        },
        units={"pressure": "GPa", "volume": "angstrom^3"},
    )
    request = EOSFitRequest(
        model="BM3",
        options=EOSFitOptions(solver_options=WLSOptions()),
    )
    path = tmp_path / "failed.hdf5"
    with EOSSession.create(dataset, path) as session:
        record = session.fit(request)
        assert not record.successful
        assert record.result.fit.status.value == "invalid_input"
        assert session.archive.record_ids == (record.record_id,)
        assert session.inspect_record(record.record_id).disposition is (
            EOSRecordDisposition.UNDECIDED
        )


def test_session_classifies_documented_request_exception_as_invalid_input(
    tmp_path: Path,
) -> None:
    """Documented request errors remain persistent INVALID_INPUT attempts."""
    path = tmp_path / "expected_error.hdf5"
    fitter = _RaisingFitter(ValueError("invalid synthetic request"))
    with EOSSession.create(DATA / "PV_quartz.dat", path, fitter=fitter) as session:
        record = session.fit(_request("BM3"))

    assert record.result.fit.status.value == "invalid_input"
    assert record.result.metadata["workflow_exception"] == "ValueError"
    assert (
        record.result.metadata["workflow_failure_classification"]
        == "request_error"
    )


def test_session_propagates_unexpected_workflow_exception(
    tmp_path: Path,
) -> None:
    """Unexpected workflow exceptions must remain visible to the caller."""
    path = tmp_path / "unexpected_error.hdf5"
    fitter = _RaisingFitter(RuntimeError("synthetic backend failure"))
    with EOSSession.create(DATA / "PV_quartz.dat", path, fitter=fitter) as session:
        with pytest.raises(RuntimeError, match="synthetic backend failure"):
            session.fit(_request("BM3"))
        assert session.archive.record_ids == ()


def test_archive_size_warning_is_advisory_and_emitted_once(tmp_path: Path) -> None:
    observer = ListObserver()
    path = tmp_path / "size_warning.hdf5"
    with EOSSession.create(
        DATA / "PV_quartz.dat",
        path,
        observer=observer,
        archive_size_warning_mib=1.0e-9,
    ) as session:
        session.fit(
            EOSFitRequest(
                model="BM3",
                options=EOSFitOptions(solver_options=OLSOptions()),
            )
        )
        warnings = [
            event
            for event in observer.events
            if event.level is EventLevel.WARNING
            and event.data.get("kind") == "archive_size_warning"
        ]
        assert len(warnings) == 1
        assert session.archive_size().exceeds_warning_threshold


def test_batch_style_archive_creation_still_refuses_implicit_overwrite(
    tmp_path: Path,
) -> None:
    path = tmp_path / "protected.hdf5"
    dataset = read_eos_input(DATA / "PV_quartz.dat")
    with EOSArchive.create(path, dataset=dataset):
        pass
    with pytest.raises(FileExistsError):
        EOSArchive.create(path, dataset=dataset)
