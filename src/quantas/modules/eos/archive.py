# -*- coding: utf-8 -*-

"""Active interface for persistent native EOS HDF5 archives."""

from __future__ import annotations

from os import replace
import json
from pathlib import Path
from typing import Any

import h5py

from .hdf5 import (
    initialize_eos_archive,
    read_eos_dataset,
    read_fit_record,
    read_slot_state,
    read_state_event,
    sorted_numeric_children,
    validate_eos_archive,
    write_accepted_result,
    write_eos_dataset,
    write_fit_record,
    write_slot_state,
    write_state_event,
)
from .history import (
    EOSFitRecord,
    EOSResultSlot,
    EOSSlotState,
    EOSSlotStatus,
    EOSStateEvent,
    EOSStateEventType,
)
from .inspection import (
    EOSArchiveInspection,
    EOSArchiveSizeInfo,
    EOSRecordDisposition,
    EOSRecordInspection,
    EOSSlotInspection,
)
from .models import EOSDataset, EOSFitRequest, EOSFitResult


class EOSArchive:
    """Read and update one native EOS HDF5 archive.

    The archive is independent of batch, interactive, CLI, and GUI workflows.
    Every numerical attempt is stored as an immutable record. Acceptance,
    rejection, notes, and use as an initial guess are append-only state events,
    while ``session/current`` and ``results/accepted`` provide compact mutable
    views of the current scientific state.

    Parameters
    ----------
    path : str or Path
        Existing native EOS archive.
    mode : {"r", "r+"}, optional
        Read-only or writable access.

    Raises
    ------
    ValueError
        If the file is not a supported native EOS archive.
    """

    def __init__(self, path: str | Path, mode: str = "r") -> None:
        if mode not in {"r", "r+"}:
            raise ValueError("EOSArchive mode must be 'r' or 'r+'")
        self.path = Path(path)
        self.mode = mode
        self._h5 = h5py.File(self.path, mode)
        try:
            validate_eos_archive(self._h5)
        except Exception:
            self._h5.close()
            raise

    @classmethod
    def create(
        cls,
        path: str | Path,
        *,
        dataset: EOSDataset | None = None,
        creator: str | None = None,
        overwrite: bool = False,
    ) -> EOSArchive:
        """Create a new archive atomically and return it opened for writing.

        Parameters
        ----------
        path : str or Path
            Destination HDF5 file.
        dataset : EOSDataset or None, optional
            Optional first input dataset. Compatible result slots are
            registered immediately with ``not_processed`` status.
        creator : str or None, optional
            Human or workflow identifier.
        overwrite : bool, optional
            Replace an existing destination.

        Returns
        -------
        EOSArchive
            Writable archive.
        """
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and not overwrite:
            raise FileExistsError(destination)
        temporary = destination.with_name(f".{destination.name}.tmp")
        if temporary.exists():
            temporary.unlink()
        try:
            with h5py.File(temporary, "w") as h5:
                initialize_eos_archive(h5, creator=creator)
            replace(temporary, destination)
        except Exception:
            if temporary.exists():
                temporary.unlink()
            raise
        archive = cls(destination, mode="r+")
        if dataset is not None:
            archive.add_dataset(dataset)
        return archive

    @property
    def writable(self) -> bool:
        """Return whether mutation operations are allowed."""
        return self.mode == "r+"

    def close(self) -> None:
        """Flush and close the underlying HDF5 file."""
        if self._h5:
            if self.writable:
                self._h5.flush()
            self._h5.close()

    def flush(self) -> None:
        """Flush pending HDF5 changes to disk.

        Raises
        ------
        OSError
            If the archive handle is already closed.
        """
        if not self._h5:
            raise OSError("EOS archive is closed")
        self._h5.flush()

    def __enter__(self) -> EOSArchive:
        """Return the open archive in a context manager."""
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        """Close the archive when leaving a context manager."""
        self.close()

    @property
    def dataset_ids(self) -> tuple[int, ...]:
        """Return all dataset identifiers in creation order."""
        group = self._h5["input/datasets"]
        return tuple(int(name) for name in sorted_numeric_children(group))

    @property
    def record_ids(self) -> tuple[int, ...]:
        """Return all immutable fit-record identifiers."""
        group = self._h5["session/records"]
        return tuple(int(name) for name in sorted_numeric_children(group))

    @property
    def event_ids(self) -> tuple[int, ...]:
        """Return all append-only state-event identifiers."""
        group = self._h5["session/state_events"]
        return tuple(int(name) for name in sorted_numeric_children(group))

    def add_dataset(self, dataset: EOSDataset) -> int:
        """Append a complete input dataset and register available result slots.

        Parameters
        ----------
        dataset : EOSDataset
            Raw and normalized EOS data.

        Returns
        -------
        int
            Monotonic dataset identifier.
        """
        self._require_writable()
        session = self._h5["session"]
        dataset_id = int(session.attrs["next_dataset_id"])
        root = self._h5["input/datasets"]
        temporary_name = f"_pending_{dataset_id:06d}"
        final_name = f"{dataset_id:06d}"
        pending = root.create_group(temporary_name)
        try:
            write_eos_dataset(pending, dataset_id, dataset)
            root.move(temporary_name, final_name)
            session.attrs["next_dataset_id"] = dataset_id + 1
            for slot in infer_result_slots(dataset):
                self.register_slot(slot)
            self._h5.flush()
        except Exception:
            if temporary_name in root:
                del root[temporary_name]
            raise
        return dataset_id

    def dataset(self, dataset_id: int = 1) -> EOSDataset:
        """Read one archived input dataset."""
        path = f"input/datasets/{int(dataset_id):06d}"
        if path not in self._h5:
            raise KeyError(f"EOS archive has no dataset {dataset_id}")
        return read_eos_dataset(self._h5[path])

    def register_slot(self, slot: str | EOSResultSlot) -> EOSSlotState:
        """Ensure that a result slot exists with an explicit empty state."""
        self._require_writable()
        resolved = EOSResultSlot.parse(slot)
        group = self._slot_group(resolved, create=True)
        if "status" not in group.attrs:
            write_slot_state(group, EOSSlotState(resolved))
            self._h5.flush()
        return read_slot_state(group)

    def slots(self) -> tuple[EOSSlotState, ...]:
        """Return all registered result slots in stable key order."""
        root = self._h5["session/current"]
        states: list[EOSSlotState] = []
        for domain_name in sorted(root):
            domain = root[domain_name]
            for target in sorted(domain):
                states.append(read_slot_state(domain[target]))
        return tuple(states)

    def slot_state(self, slot: str | EOSResultSlot) -> EOSSlotState:
        """Return the current state of one registered result slot.

        Raises
        ------
        KeyError
            If the archive has no such scientifically available slot. This
            distinguishes an unavailable property from an available but
            ``not_processed`` property.
        """
        resolved = EOSResultSlot.parse(slot)
        path = f"session/current/{resolved.key}"
        if path not in self._h5:
            raise KeyError(
                f"EOS archive has no registered result slot {resolved.key!r}"
            )
        return read_slot_state(self._h5[path])

    def append_fit(
        self,
        dataset_id: int,
        request: EOSFitRequest,
        result: EOSFitResult,
        *,
        parent_record_id: int | None = None,
        note: str | None = None,
        provenance: dict[str, Any] | None = None,
        accept: bool = False,
    ) -> EOSFitRecord:
        """Append one immutable fit attempt and optionally accept it.

        Failed and invalid-input fit results are persisted exactly like
        successful attempts, but they cannot be accepted.
        """
        self._require_writable()
        if int(dataset_id) not in self.dataset_ids:
            raise KeyError(f"EOS archive has no dataset {dataset_id}")
        if result.request.as_dict() != request.as_dict():
            raise ValueError("request and result.request must be identical")
        slot = EOSResultSlot.from_request(request)
        self.register_slot(slot)
        if parent_record_id is not None:
            parent = self.record(parent_record_id)
            if parent.slot != slot:
                raise ValueError("parent EOS record must address the same result slot")

        session = self._h5["session"]
        record_id = int(session.attrs["next_record_id"])
        record = EOSFitRecord(
            record_id=record_id,
            dataset_id=int(dataset_id),
            request=request,
            result=result,
            parent_record_id=parent_record_id,
            note=note,
            provenance=dict(provenance or {}),
        )
        root = self._h5["session/records"]
        temporary_name = f"_pending_{record_id:06d}"
        final_name = f"{record_id:06d}"
        pending = root.create_group(temporary_name)
        try:
            write_fit_record(pending, record)
            root.move(temporary_name, final_name)
            session.attrs["next_record_id"] = record_id + 1
            state = self.slot_state(slot)
            attempts = (*state.attempted_record_ids, record_id)
            updated = EOSSlotState(
                slot=slot,
                status=(
                    EOSSlotStatus.ACCEPTED
                    if state.accepted_record_id is not None
                    else EOSSlotStatus.ATTEMPTED
                ),
                accepted_record_id=state.accepted_record_id,
                last_record_id=record_id,
                attempted_record_ids=attempts,
            )
            write_slot_state(self._slot_group(slot, create=True), updated)
            self._append_event(
                EOSStateEventType.RECORD_CREATED,
                record_id=record_id,
                slot=slot,
                note=note,
                metadata={"fit_status": result.fit.status.value},
            )
            self._h5.flush()
        except Exception:
            if temporary_name in root:
                del root[temporary_name]
            raise
        if accept:
            self.accept(record_id)
        return record

    def store_fit(
        self,
        dataset_id: int,
        request: EOSFitRequest,
        result: EOSFitResult,
        *,
        parent_record_id: int | None = None,
        note: str | None = None,
        provenance: dict[str, Any] | None = None,
    ) -> EOSFitRecord:
        """Append and accept one fit, convenient for Python and batch use."""
        return self.append_fit(
            dataset_id,
            request,
            result,
            parent_record_id=parent_record_id,
            note=note,
            provenance=provenance,
            accept=True,
        )

    def record(self, record_id: int) -> EOSFitRecord:
        """Read one immutable fit record."""
        path = f"session/records/{int(record_id):06d}"
        if path not in self._h5:
            raise KeyError(f"EOS archive has no fit record {record_id}")
        return read_fit_record(self._h5[path])

    def records(
        self,
        slot: str | EOSResultSlot | None = None,
    ) -> tuple[EOSFitRecord, ...]:
        """Return all fit records, optionally filtered by result slot."""
        resolved = None if slot is None else EOSResultSlot.parse(slot)
        values = tuple(self.record(record_id) for record_id in self.record_ids)
        if resolved is None:
            return values
        return tuple(record for record in values if record.slot == resolved)

    def accept(
        self,
        record_id: int,
        *,
        note: str | None = None,
    ) -> EOSSlotState:
        """Select a successful fit record as the current accepted result."""
        self._require_writable()
        record = self.record(record_id)
        if not record.successful:
            raise ValueError("only successful EOS fit records can be accepted")
        slot = record.slot
        state = self.slot_state(slot)
        previous = state.accepted_record_id
        accepted_root = self._accepted_parent(slot, create=True)
        temporary_name = f"_pending_{slot.target}"
        if temporary_name in accepted_root:
            del accepted_root[temporary_name]
        pending = accepted_root.create_group(temporary_name)
        try:
            write_accepted_result(pending, record)
            if slot.target in accepted_root:
                del accepted_root[slot.target]
            accepted_root.move(temporary_name, slot.target)
            updated = EOSSlotState(
                slot=slot,
                status=EOSSlotStatus.ACCEPTED,
                accepted_record_id=record_id,
                last_record_id=state.last_record_id,
                attempted_record_ids=state.attempted_record_ids,
            )
            write_slot_state(self._slot_group(slot, create=True), updated)
            self._append_event(
                EOSStateEventType.RECORD_ACCEPTED,
                record_id=record_id,
                slot=slot,
                note=note,
                metadata={"previous_accepted_record_id": previous},
            )
            self._h5.flush()
        except Exception:
            if temporary_name in accepted_root:
                del accepted_root[temporary_name]
            raise
        return updated

    def reject(self, record_id: int, *, note: str | None = None) -> EOSStateEvent:
        """Append a rejection decision without deleting the fit record.

        Raises
        ------
        ValueError
            If the record is currently accepted. Acceptance must first be
            revoked explicitly with :meth:`unaccept` so the archive history is
            unambiguous.
        """
        record = self.record(record_id)
        state = self.slot_state(record.slot)
        if state.accepted_record_id == record_id:
            raise ValueError(
                "cannot reject the current accepted EOS record; revoke acceptance first"
            )
        return self._append_event(
            EOSStateEventType.RECORD_REJECTED,
            record_id=record_id,
            slot=record.slot,
            note=note,
        )

    def unaccept(
        self,
        slot: str | EOSResultSlot,
        *,
        note: str | None = None,
    ) -> EOSSlotState:
        """Revoke the current accepted result for one slot.

        The immutable fit record and its earlier acceptance event are retained.
        The materialized accepted-result view is removed and a
        ``record_unaccepted`` event is appended.

        Parameters
        ----------
        slot : str or EOSResultSlot
            Scientific slot whose current acceptance is revoked.
        note : str or None, optional
            Scientific reason for the revocation.

        Returns
        -------
        EOSSlotState
            Updated slot state with no accepted record.

        Raises
        ------
        ValueError
            If the slot has no current accepted record.
        """
        self._require_writable()
        resolved = EOSResultSlot.parse(slot)
        state = self.slot_state(resolved)
        record_id = state.accepted_record_id
        if record_id is None:
            raise ValueError(f"EOS result slot {resolved.key!r} has no accepted record")

        accepted_root = self._accepted_parent(resolved, create=True)
        if resolved.target in accepted_root:
            del accepted_root[resolved.target]
        updated = EOSSlotState(
            slot=resolved,
            status=(
                EOSSlotStatus.ATTEMPTED
                if state.attempted_record_ids
                else EOSSlotStatus.NOT_PROCESSED
            ),
            accepted_record_id=None,
            last_record_id=state.last_record_id,
            attempted_record_ids=state.attempted_record_ids,
        )
        write_slot_state(self._slot_group(resolved, create=True), updated)
        self._append_event(
            EOSStateEventType.RECORD_UNACCEPTED,
            record_id=record_id,
            slot=resolved,
            note=note,
        )
        self._h5.flush()
        return updated

    def mark_candidate(
        self,
        record_id: int,
        *,
        note: str | None = None,
    ) -> EOSStateEvent:
        """Bookmark a fit record without changing the accepted result."""
        record = self.record(record_id)
        return self._append_event(
            EOSStateEventType.RECORD_CANDIDATE,
            record_id=record_id,
            slot=record.slot,
            note=note,
        )

    def mark_used_as_initial_guess(
        self,
        record_id: int,
        *,
        child_record_id: int | None = None,
        note: str | None = None,
    ) -> EOSStateEvent:
        """Record explicit reuse of one result as a later initial guess."""
        record = self.record(record_id)
        metadata: dict[str, Any] = {}
        if child_record_id is not None:
            child = self.record(child_record_id)
            if child.parent_record_id != record_id:
                raise ValueError("child record does not identify the source as parent")
            metadata["child_record_id"] = child_record_id
        return self._append_event(
            EOSStateEventType.RECORD_USED_AS_INITIAL_GUESS,
            record_id=record_id,
            slot=record.slot,
            note=note,
            metadata=metadata,
        )

    def add_note(self, record_id: int, note: str) -> EOSStateEvent:
        """Append a note event without mutating the original fit record."""
        if not str(note).strip():
            raise ValueError("EOS archive note cannot be empty")
        record = self.record(record_id)
        return self._append_event(
            EOSStateEventType.NOTE_ADDED,
            record_id=record_id,
            slot=record.slot,
            note=note,
        )

    def record_no_operation(
        self,
        slot: str | EOSResultSlot,
        *,
        note: str | None = None,
    ) -> EOSStateEvent:
        """Record that a registered property was deliberately left unfitted."""
        resolved = EOSResultSlot.parse(slot)
        self.slot_state(resolved)
        return self._append_event(
            EOSStateEventType.NO_OPERATION,
            slot=resolved,
            note=note,
        )

    def events(self) -> tuple[EOSStateEvent, ...]:
        """Return all append-only archive state events."""
        root = self._h5["session/state_events"]
        return tuple(
            read_state_event(root[name]) for name in sorted_numeric_children(root)
        )

    def accepted(self, slot: str | EOSResultSlot) -> EOSFitRecord | None:
        """Return the current accepted fit record for one slot."""
        state = self.slot_state(slot)
        if state.accepted_record_id is None:
            return None
        return self.record(state.accepted_record_id)

    def accepted_result(self, slot: str | EOSResultSlot) -> EOSFitResult | None:
        """Return only the current accepted numerical result for one slot."""
        record = self.accepted(slot)
        return None if record is None else record.result

    def write_batch_manifest(self, payload: dict[str, Any]) -> None:
        """Persist one normalized declarative batch plan.

        Parameters
        ----------
        payload : dict
            Serialization-ready plan mapping. The manifest is write-once for
            one newly created archive.

        Raises
        ------
        OSError
            If the archive is read-only.
        ValueError
            If a batch manifest is already present.
        """
        self._require_writable()
        session = self._h5["session"]
        if "batch_plan" in session:
            raise ValueError("EOS archive already contains a batch manifest")
        text = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        session.create_dataset(
            "batch_plan",
            data=text,
            dtype=h5py.string_dtype(encoding="utf-8"),
        )
        self._h5.flush()

    def batch_manifest(self) -> dict[str, Any] | None:
        """Return the archived declarative batch plan, when present."""
        session = self._h5["session"]
        if "batch_plan" not in session:
            return None
        value = session["batch_plan"][()]
        if isinstance(value, bytes):
            value = value.decode("utf-8")
        return dict(json.loads(str(value)))

    def summary(self) -> dict[str, Any]:
        """Return a compact frontend-neutral archive summary."""
        return {
            "path": str(self.path),
            "schema_version": str(self._h5["metadata"].attrs["schema_version"]),
            "datasets": list(self.dataset_ids),
            "records": list(self.record_ids),
            "events": list(self.event_ids),
            "slots": {state.slot.key: state.as_dict() for state in self.slots()},
        }

    def size_info(
        self,
        *,
        warning_threshold_mib: float | None = 100.0,
    ) -> EOSArchiveSizeInfo:
        """Return current physical archive size and warning status.

        Parameters
        ----------
        warning_threshold_mib : float or None, optional
            Advisory threshold in mebibytes. ``None`` disables the warning
            comparison. This method never prevents archive growth.

        Returns
        -------
        EOSArchiveSizeInfo
            Current file size and threshold comparison.

        Raises
        ------
        ValueError
            If the threshold is not positive.
        """
        self.flush()
        threshold_bytes: int | None
        if warning_threshold_mib is None:
            threshold_bytes = None
        else:
            threshold = float(warning_threshold_mib)
            if threshold <= 0.0:
                raise ValueError("EOS archive warning threshold must be positive")
            threshold_bytes = int(round(threshold * 1024.0**2))
            if threshold_bytes <= 0:
                threshold_bytes = 1
        return EOSArchiveSizeInfo(
            size_bytes=self.path.stat().st_size,
            warning_threshold_bytes=threshold_bytes,
        )

    def inspect_record(self, record_id: int) -> EOSRecordInspection:
        """Return a derived scientific inspection view for one record.

        Parameters
        ----------
        record_id : int
            Immutable fit-record identifier.

        Returns
        -------
        EOSRecordInspection
            Record, effective disposition, lineage, events, and notes.
        """
        record = self.record(record_id)
        records = self.records(record.slot)
        events = self.events()
        return _build_record_inspection(
            record,
            state=self.slot_state(record.slot),
            slot_records=records,
            archive_events=events,
        )

    def inspect_slot(self, slot: str | EOSResultSlot) -> EOSSlotInspection:
        """Return all attempts and current state for one result slot."""
        state = self.slot_state(slot)
        records = self.records(state.slot)
        events = self.events()
        return EOSSlotInspection(
            state=state,
            records=tuple(
                _build_record_inspection(
                    record,
                    state=state,
                    slot_records=records,
                    archive_events=events,
                )
                for record in records
            ),
        )

    def inspect(
        self,
        *,
        warning_threshold_mib: float | None = 100.0,
    ) -> EOSArchiveInspection:
        """Return a complete frontend-neutral archive inspection snapshot."""
        states = self.slots()
        records = self.records()
        events = self.events()
        records_by_slot = {
            state.slot: tuple(record for record in records if record.slot == state.slot)
            for state in states
        }
        slots = tuple(
            EOSSlotInspection(
                state=state,
                records=tuple(
                    _build_record_inspection(
                        record,
                        state=state,
                        slot_records=records_by_slot[state.slot],
                        archive_events=events,
                    )
                    for record in records_by_slot[state.slot]
                ),
            )
            for state in states
        )
        return EOSArchiveInspection(
            path=self.path,
            schema_version=str(self._h5["metadata"].attrs["schema_version"]),
            dataset_ids=self.dataset_ids,
            slots=slots,
            event_count=len(events),
            size=self.size_info(warning_threshold_mib=warning_threshold_mib),
        )

    def _append_event(
        self,
        event_type: EOSStateEventType,
        *,
        record_id: int | None = None,
        slot: EOSResultSlot | None = None,
        note: str | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> EOSStateEvent:
        self._require_writable()
        session = self._h5["session"]
        event_id = int(session.attrs["next_event_id"])
        event = EOSStateEvent(
            event_id=event_id,
            event_type=event_type,
            record_id=record_id,
            slot=slot,
            note=note,
            metadata=dict(metadata or {}),
        )
        root = self._h5["session/state_events"]
        name = f"{event_id:06d}"
        group = root.create_group(name)
        try:
            write_state_event(group, event)
            session.attrs["next_event_id"] = event_id + 1
            self._h5.flush()
        except Exception:
            if name in root:
                del root[name]
            raise
        return event

    def _slot_group(self, slot: EOSResultSlot, *, create: bool) -> h5py.Group:
        root = self._h5["session/current"]
        if create:
            domain = root.require_group(slot.domain.value)
            return domain.require_group(slot.target)
        return root[slot.domain.value][slot.target]

    def _accepted_parent(self, slot: EOSResultSlot, *, create: bool) -> h5py.Group:
        root = self._h5["results/accepted"]
        if create:
            return root.require_group(slot.domain.value)
        return root[slot.domain.value]

    def _require_writable(self) -> None:
        if not self.writable:
            raise OSError("EOS archive is open read-only")


def infer_result_slots(dataset: EOSDataset) -> tuple[EOSResultSlot, ...]:
    """Infer scientifically possible result slots from available coordinates.

    Constant control coordinates do not create fitting slots. For example, a
    constant pressure in an isobaric V--T dataset is preserved as metadata but
    does not create P--V slots.
    """
    slots: list[EOSResultSlot] = []
    targets = tuple(name for name in ("volume", "a", "b", "c") if dataset.has(name))
    if dataset.has("pressure"):
        try:
            pressure_variable = dataset.coordinate_profile("pressure").is_variable
        except ValueError:
            pressure_variable = False
        if pressure_variable:
            for target in targets:
                if dataset.coordinate_profile(target).is_variable:
                    slots.append(EOSResultSlot.parse(f"pv/{target}"))
    if dataset.has("temperature"):
        try:
            temperature_variable = dataset.coordinate_profile("temperature").is_variable
        except ValueError:
            temperature_variable = False
        if temperature_variable:
            for target in targets:
                if dataset.coordinate_profile(target).is_variable:
                    slots.append(EOSResultSlot.parse(f"vt/{target}"))
    if all(dataset.has(name) for name in ("pressure", "volume", "temperature")):
        profiles = {
            name: dataset.coordinate_profile(name)
            for name in ("pressure", "volume", "temperature")
        }
        if all(profile.is_variable for profile in profiles.values()):
            slots.append(EOSResultSlot.parse("pvt/volume"))
    if dataset.has("energy") and dataset.has("volume"):
        if (
            dataset.coordinate_profile("energy").is_variable
            and dataset.coordinate_profile("volume").is_variable
        ):
            slots.append(EOSResultSlot.parse("ev/energy"))
    unique = {slot.key: slot for slot in slots}
    return tuple(unique[key] for key in sorted(unique))


def _build_record_inspection(
    record: EOSFitRecord,
    *,
    state: EOSSlotState,
    slot_records: tuple[EOSFitRecord, ...],
    archive_events: tuple[EOSStateEvent, ...],
) -> EOSRecordInspection:
    """Build one record view from an already loaded archive snapshot."""
    events = tuple(
        event for event in archive_events if event.record_id == record.record_id
    )
    child_ids = tuple(
        item.record_id
        for item in slot_records
        if item.parent_record_id == record.record_id
    )
    notes = tuple(
        value
        for value in (
            record.note,
            *(event.note for event in events),
        )
        if value is not None and str(value).strip()
    )
    used_as_initial_guess = any(
        event.event_type is EOSStateEventType.RECORD_USED_AS_INITIAL_GUESS
        for event in events
    )
    return EOSRecordInspection(
        record=record,
        disposition=_record_disposition(record.record_id, state, events),
        events=events,
        child_record_ids=child_ids,
        notes=notes,
        used_as_initial_guess=used_as_initial_guess,
    )


def _record_disposition(
    record_id: int,
    state: EOSSlotState,
    events: tuple[EOSStateEvent, ...],
) -> EOSRecordDisposition:
    """Derive one record's effective disposition from state and events."""
    if state.accepted_record_id == record_id:
        return EOSRecordDisposition.ACCEPTED
    decisions = {
        EOSStateEventType.RECORD_ACCEPTED,
        EOSStateEventType.RECORD_UNACCEPTED,
        EOSStateEventType.RECORD_REJECTED,
        EOSStateEventType.RECORD_CANDIDATE,
    }
    direct = tuple(event for event in events if event.event_type in decisions)
    last = None if not direct else direct[-1].event_type
    if last is EOSStateEventType.RECORD_REJECTED:
        return EOSRecordDisposition.REJECTED
    if last is EOSStateEventType.RECORD_CANDIDATE:
        return EOSRecordDisposition.CANDIDATE
    if last is EOSStateEventType.RECORD_UNACCEPTED:
        return EOSRecordDisposition.WITHDRAWN
    if any(event.event_type is EOSStateEventType.RECORD_ACCEPTED for event in direct):
        return EOSRecordDisposition.SUPERSEDED
    return EOSRecordDisposition.UNDECIDED


__all__ = ["EOSArchive", "infer_result_slots"]
