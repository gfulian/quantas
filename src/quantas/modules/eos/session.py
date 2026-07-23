# -*- coding: utf-8 -*-

"""Frontend-neutral persistent EOS fitting sessions.

An :class:`EOSSession` coordinates the shared numerical fitter, a writable
native EOS archive, and Quantas observers. It contains no Click, Rich, terminal,
or plotting dependencies. Interactive and graphical frontends remain
responsible only for collecting choices and rendering returned contracts.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np

from quantas.core.events import Event, EventLevel, EventRecord, NullObserver
from quantas.core.math.fitting import FitResult, FitStatus, ParameterState

from .api import EOSFitter
from .archive import EOSArchive, infer_result_slots
from .history import EOSFitRecord, EOSResultSlot, EOSSlotState, EOSStateEvent
from .inspection import (
    EOSArchiveInspection,
    EOSArchiveSizeInfo,
    EOSRecordInspection,
    EOSSlotInspection,
)
from .io import read_eos_input
from .models import (
    EOSDataset,
    EOSFitRequest,
    EOSFitResult,
    ParameterConstraint,
)


class EOSSession:
    """Coordinate one persistent, resumable EOS analysis session.

    Parameters
    ----------
    archive : EOSArchive
        Writable native EOS archive. The session does not assume ownership
        unless ``owns_archive`` is true.
    dataset_id : int, optional
        Embedded dataset used for new fitting attempts.
    fitter : EOSFitter or None, optional
        Shared scientific fitting service.
    observer : callable or None, optional
        Quantas event observer used by a CLI, GUI, notebook, or test.
    archive_size_warning_mib : float or None, optional
        Advisory archive-size threshold in mebibytes. The default is 100 MiB.
        ``None`` disables warnings. Growth is never blocked.
    owns_archive : bool, optional
        Close the archive when the session closes.

    Raises
    ------
    OSError
        If the archive is read-only.
    KeyError
        If ``dataset_id`` is not embedded in the archive.
    ValueError
        If the size-warning threshold is not positive.

    Notes
    -----
    A session stores every attempted fit before any acceptance, rejection, or
    candidate decision. Reopening the HDF5 file therefore resumes both the
    scientific history and the current accepted state without requiring the
    original text input file.
    """

    def __init__(
        self,
        archive: EOSArchive,
        *,
        dataset_id: int = 1,
        fitter: EOSFitter | None = None,
        observer: Callable[[Event], None] | None = None,
        archive_size_warning_mib: float | None = 100.0,
        owns_archive: bool = False,
    ) -> None:
        if not isinstance(archive, EOSArchive):
            raise TypeError("EOSSession archive must be an EOSArchive")
        if not archive.writable:
            raise OSError("EOSSession requires a writable EOS archive")
        self.archive = archive
        self.dataset_id = int(dataset_id)
        self._dataset = archive.dataset(self.dataset_id)
        self.fitter = EOSFitter() if fitter is None else fitter
        self.observer = NullObserver() if observer is None else observer
        self.archive_size_warning_mib = _normalize_warning_threshold(
            archive_size_warning_mib
        )
        self.owns_archive = bool(owns_archive)
        self._events: list[EventRecord] = []
        self._closed = False
        self._size_warning_emitted = False
        self._check_archive_size()

    @classmethod
    def create(
        cls,
        source: str | Path | EOSDataset,
        archive_path: str | Path,
        *,
        creator: str = "quantas eos interactive",
        overwrite: bool = False,
        pressure_unit: str | None = None,
        length_unit: str | None = None,
        temperature_unit: str | None = None,
        fitter: EOSFitter | None = None,
        observer: Callable[[Event], None] | None = None,
        archive_size_warning_mib: float | None = 100.0,
    ) -> EOSSession:
        """Create a new persistent EOS session.

        Parameters
        ----------
        source : str, Path, or EOSDataset
            Keyword-directed text input or a normalized in-memory dataset.
        archive_path : str or Path
            Destination native EOS HDF5 file.
        creator : str, optional
            Archive creation provenance.
        overwrite : bool, optional
            Replace an existing archive. Interactive frontends should normally
            leave this false and offer explicit resume instead.
        pressure_unit, length_unit, temperature_unit : str or None, optional
            Input-unit overrides used only for text input.
        fitter, observer, archive_size_warning_mib : optional
            Services and advisory size threshold passed to the session.

        Returns
        -------
        EOSSession
            New writable session owning its archive handle.

        Raises
        ------
        FileExistsError
            If the destination exists and ``overwrite`` is false.
        OSError, ValueError
            If the input cannot be read or the archive cannot be created.
        """
        dataset = (
            source
            if isinstance(source, EOSDataset)
            else read_eos_input(
                source,
                pressure_unit=pressure_unit,
                length_unit=length_unit,
                temperature_unit=temperature_unit,
            )
        )
        archive = EOSArchive.create(
            archive_path,
            dataset=dataset,
            creator=creator,
            overwrite=overwrite,
        )
        try:
            session = cls(
                archive,
                dataset_id=1,
                fitter=fitter,
                observer=observer,
                archive_size_warning_mib=archive_size_warning_mib,
                owns_archive=True,
            )
        except Exception:
            archive.close()
            raise
        session._emit(
            f"Created EOS session archive at {archive.path}",
            data={"kind": "session_created", "path": str(archive.path)},
        )
        return session

    @classmethod
    def resume(
        cls,
        archive_path: str | Path,
        *,
        dataset_id: int = 1,
        fitter: EOSFitter | None = None,
        observer: Callable[[Event], None] | None = None,
        archive_size_warning_mib: float | None = 100.0,
    ) -> EOSSession:
        """Resume a writable session from an existing native EOS archive.

        Parameters
        ----------
        archive_path : str or Path
            Existing native EOS HDF5 archive.
        dataset_id : int, optional
            Embedded dataset used for subsequent attempts.
        fitter, observer, archive_size_warning_mib : optional
            Services and advisory size threshold passed to the session.

        Returns
        -------
        EOSSession
            Resumed session owning the opened archive handle.

        Raises
        ------
        OSError
            If the file cannot be opened for update.
        ValueError
            If the file is not a supported EOS archive.
        KeyError
            If ``dataset_id`` is absent.
        """
        archive = EOSArchive(archive_path, mode="r+")
        try:
            session = cls(
                archive,
                dataset_id=dataset_id,
                fitter=fitter,
                observer=observer,
                archive_size_warning_mib=archive_size_warning_mib,
                owns_archive=True,
            )
        except Exception:
            archive.close()
            raise
        session._emit(
            f"Resumed EOS session archive at {archive.path}",
            data={
                "kind": "session_resumed",
                "path": str(archive.path),
                "dataset_id": dataset_id,
            },
        )
        return session

    @property
    def dataset(self) -> EOSDataset:
        """Return the normalized embedded dataset used by the session."""
        return self._dataset

    @property
    def events(self) -> tuple[EventRecord, ...]:
        """Return meaningful operational events emitted by this controller."""
        return tuple(self._events)

    @property
    def available_slots(self) -> tuple[EOSResultSlot, ...]:
        """Return result slots scientifically supported by the session dataset."""
        return infer_result_slots(self.dataset)

    def close(self) -> None:
        """Flush the archive and close an owned handle."""
        if self._closed:
            return
        self.archive.flush()
        if self.owns_archive:
            self.archive.close()
        self._closed = True

    def __enter__(self) -> EOSSession:
        """Return the open session in a context manager."""
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        """Close the session when leaving a context manager."""
        self.close()

    def prepare_request_from_record(
        self,
        request: EOSFitRequest,
        source_record_id: int,
    ) -> EOSFitRequest:
        """Build a request using compatible fitted values as initial guesses.

        Parameters
        ----------
        request : EOSFitRequest
            Target model and explicit user constraints. Explicit constraints
            always take precedence over reused values.
        source_record_id : int
            Successful earlier record in the same result slot.

        Returns
        -------
        EOSFitRequest
            Detached request containing model-aware initial-value overrides.

        Raises
        ------
        ValueError
            If the source failed, addresses another slot, or belongs to another
            embedded dataset.

        Notes
        -----
        Matching is performed by stable parameter name, never by vector
        position. Only parameters free in the target model are copied. Fixed,
        implied, derived, missing, or out-of-bound source values are ignored.
        """
        source = self._session_record(source_record_id)
        target_slot = EOSResultSlot.from_request(request)
        if source.slot != target_slot:
            raise ValueError(
                "initial-guess source must address the same EOS result slot"
            )
        if not source.successful:
            raise ValueError("initial-guess source must be a successful EOS fit record")

        definitions = self.fitter.parameter_definitions(self.dataset, request)
        existing = {constraint.name: constraint for constraint in request.constraints}
        source_values = source.result.parameter_values
        copied: list[ParameterConstraint] = []
        copied_names: list[str] = []
        for definition in definitions:
            if definition.name in existing:
                continue
            if definition.state is not ParameterState.FREE:
                continue
            if definition.name not in source_values:
                continue
            value = float(source_values[definition.name])
            if not np.isfinite(value):
                continue
            if not definition.lower_bound <= value <= definition.upper_bound:
                continue
            copied.append(
                ParameterConstraint.free(
                    definition.name,
                    value,
                    lower_bound=definition.lower_bound,
                    upper_bound=definition.upper_bound,
                    unit=definition.unit,
                    description=definition.description,
                    metadata={
                        **dict(definition.metadata),
                        "initial_guess_record_id": source.record_id,
                    },
                )
            )
            copied_names.append(definition.name)

        metadata = {
            **request.metadata,
            "initial_guess": {
                "source_record_id": source.record_id,
                "copied_parameters": copied_names,
            },
        }
        return replace(
            request,
            constraints=(*request.constraints, *copied),
            metadata=metadata,
        )

    def fit(
        self,
        request: EOSFitRequest,
        *,
        initial_guess_record_id: int | None = None,
        note: str | None = None,
        provenance: dict[str, Any] | None = None,
    ) -> EOSFitRecord:
        """Execute and immediately persist one EOS fitting attempt.

        Parameters
        ----------
        request : EOSFitRequest
            Complete target model and numerical strategy.
        initial_guess_record_id : int or None, optional
            Earlier compatible record whose fitted values initialize the new
            target model by parameter name.
        note : str or None, optional
            Creation note stored with the immutable record.
        provenance : dict or None, optional
            Additional passive caller provenance.

        Returns
        -------
        EOSFitRecord
            Immutable stored record. No acceptance decision is made.

        Raises
        ------
        ValueError
            If the requested slot is unavailable or the initial-guess source is
            incompatible.

        Notes
        -----
        Numerical and request-validation failures produced after a valid slot
        is selected are converted to failed :class:`EOSFitResult` objects and
        persisted. The frontend may then inspect or reject them like any other
        attempt.
        """
        self._require_open()
        slot = EOSResultSlot.from_request(request)
        if slot not in self.available_slots:
            raise ValueError(
                f"dataset {self.dataset_id} does not provide EOS result slot "
                f"{slot.key!r}"
            )
        actual_request = request
        if initial_guess_record_id is not None:
            actual_request = self.prepare_request_from_record(
                request,
                initial_guess_record_id,
            )
        self._emit(
            f"Starting EOS fit for {slot.key} with "
            f"{getattr(actual_request.model, 'tag', actual_request.model)}",
            data={"kind": "fit_start", "request": actual_request},
        )
        result = self._execute_fit(actual_request)
        record = self.archive.append_fit(
            self.dataset_id,
            actual_request,
            result,
            parent_record_id=initial_guess_record_id,
            note=note,
            provenance={"workflow": "interactive-session", **dict(provenance or {})},
        )
        if initial_guess_record_id is not None:
            self.archive.mark_used_as_initial_guess(
                initial_guess_record_id,
                child_record_id=record.record_id,
            )
        level = EventLevel.RESULT if record.successful else EventLevel.WARNING
        self._emit(
            f"Stored EOS fit record #{record.record_id}: {result.fit.status.value}",
            level=level,
            data={"kind": "fit_record", "record": record},
        )
        self._check_archive_size()
        return record

    def accept(self, record_id: int, *, note: str | None = None) -> EOSSlotState:
        """Accept one successful session record as the current result."""
        record = self._session_record(record_id)
        state = self.archive.accept(record.record_id, note=note)
        self._emit(
            f"Accepted EOS fit record #{record.record_id} for {record.slot.key}",
            level=EventLevel.RESULT,
            data={"kind": "record_accepted", "record_id": record.record_id},
        )
        self._check_archive_size()
        return state

    def unaccept(
        self,
        slot: str | EOSResultSlot,
        *,
        note: str | None = None,
    ) -> EOSSlotState:
        """Revoke the current accepted result for one slot."""
        resolved = EOSResultSlot.parse(slot)
        previous = self.archive.slot_state(resolved).accepted_record_id
        state = self.archive.unaccept(resolved, note=note)
        self._emit(
            f"Revoked accepted EOS result for {resolved.key}",
            level=EventLevel.WARNING,
            data={
                "kind": "record_unaccepted",
                "record_id": previous,
                "slot": resolved.key,
            },
        )
        self._check_archive_size()
        return state

    def reject(self, record_id: int, *, note: str | None = None) -> EOSStateEvent:
        """Reject one stored record without deleting it."""
        record = self._session_record(record_id)
        event = self.archive.reject(record.record_id, note=note)
        self._emit(
            f"Rejected EOS fit record #{record.record_id}",
            level=EventLevel.WARNING,
            data={"kind": "record_rejected", "record_id": record.record_id},
        )
        self._check_archive_size()
        return event

    def mark_candidate(
        self,
        record_id: int,
        *,
        note: str | None = None,
    ) -> EOSStateEvent:
        """Bookmark one stored record for later comparison."""
        record = self._session_record(record_id)
        event = self.archive.mark_candidate(record.record_id, note=note)
        self._emit(
            f"Marked EOS fit record #{record.record_id} as a candidate",
            data={"kind": "record_candidate", "record_id": record.record_id},
        )
        self._check_archive_size()
        return event

    def add_note(self, record_id: int, note: str) -> EOSStateEvent:
        """Append a scientific note to one session record."""
        record = self._session_record(record_id)
        event = self.archive.add_note(record.record_id, note)
        self._emit(
            f"Added note to EOS fit record #{record.record_id}",
            data={"kind": "record_note", "record_id": record.record_id},
        )
        self._check_archive_size()
        return event

    def inspect_record(self, record_id: int) -> EOSRecordInspection:
        """Return the derived inspection view for one session record."""
        self._session_record(record_id)
        return self.archive.inspect_record(record_id)

    def inspect_slot(self, slot: str | EOSResultSlot) -> EOSSlotInspection:
        """Return the derived history and current state of one slot."""
        return self.archive.inspect_slot(slot)

    def inspect(self) -> EOSArchiveInspection:
        """Return a complete frontend-neutral archive inspection snapshot."""
        return self.archive.inspect(warning_threshold_mib=self.archive_size_warning_mib)

    def archive_size(self) -> EOSArchiveSizeInfo:
        """Return the current archive size and advisory warning state."""
        return self.archive.size_info(
            warning_threshold_mib=self.archive_size_warning_mib
        )

    def _session_record(self, record_id: int) -> EOSFitRecord:
        """Return one record after checking the selected embedded dataset."""
        self._require_open()
        record = self.archive.record(record_id)
        if record.dataset_id != self.dataset_id:
            raise ValueError(
                f"EOS record {record_id} belongs to dataset {record.dataset_id}, "
                f"not session dataset {self.dataset_id}"
            )
        return record

    def _execute_fit(self, request: EOSFitRequest) -> EOSFitResult:
        """Convert workflow exceptions into persistent failed fit results."""
        try:
            return self.fitter.fit(self.dataset, request)
        except Exception as exc:
            fit = FitResult.failed(
                str(exc),
                status=FitStatus.INVALID_INPUT,
                method=(
                    None
                    if request.options.solver_options is None
                    else request.options.solver_options.method
                ),
            )
            return EOSFitResult(
                request=request,
                fit=fit,
                warnings=[str(exc)],
                metadata={"workflow_exception": type(exc).__name__},
            )

    def _check_archive_size(self) -> EOSArchiveSizeInfo:
        """Emit one warning after the advisory size threshold is reached."""
        info = self.archive_size()
        if info.exceeds_warning_threshold and not self._size_warning_emitted:
            threshold = info.warning_threshold_mib
            self._emit(
                "EOS archive size is "
                f"{info.size_mib:.2f} MiB and has reached the advisory "
                f"threshold of {threshold:.2f} MiB",
                level=EventLevel.WARNING,
                data={"kind": "archive_size_warning", "size": info.as_dict()},
            )
            self._size_warning_emitted = True
        return info

    def _emit(
        self,
        message: str,
        *,
        level: EventLevel = EventLevel.INFO,
        progress: float | None = None,
        data: dict[str, Any] | None = None,
    ) -> None:
        """Emit one operational Quantas event and retain meaningful records."""
        event = Event(message=message, level=level, progress=progress, data=data or {})
        if level is not EventLevel.PROGRESS:
            self._events.append(EventRecord.from_event(event))
        self.observer(event)

    def _require_open(self) -> None:
        """Reject operations after the controller has been closed."""
        if self._closed:
            raise OSError("EOS session is closed")


def _normalize_warning_threshold(value: float | None) -> float | None:
    """Return a positive advisory archive-size threshold or ``None``."""
    if value is None:
        return None
    threshold = float(value)
    if not np.isfinite(threshold) or threshold <= 0.0:
        raise ValueError("EOS archive size warning threshold must be positive")
    return threshold


__all__ = ["EOSSession"]
