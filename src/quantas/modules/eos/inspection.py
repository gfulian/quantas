# -*- coding: utf-8 -*-

"""Passive archive-inspection contracts for EOS sessions."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any

from .history import EOSFitRecord, EOSSlotState, EOSStateEvent


class EOSRecordDisposition(str, Enum):
    """Effective scientific disposition of one immutable fit record.

    Attributes
    ----------
    UNDECIDED
        The fit was stored but no explicit scientific decision was recorded.
    CANDIDATE
        The record is bookmarked for later comparison.
    REJECTED
        The record was explicitly rejected.
    ACCEPTED
        The record is the current accepted result for its slot.
    SUPERSEDED
        The record was accepted previously but another record is current.
    WITHDRAWN
        Acceptance was explicitly revoked and no later decision replaced it.
    """

    UNDECIDED = "undecided"
    CANDIDATE = "candidate"
    REJECTED = "rejected"
    ACCEPTED = "accepted"
    SUPERSEDED = "superseded"
    WITHDRAWN = "withdrawn"


@dataclass(frozen=True, slots=True)
class EOSArchiveSizeInfo:
    """Physical size and warning state of one EOS archive.

    Parameters
    ----------
    size_bytes : int
        Current on-disk file size.
    warning_threshold_bytes : int or None
        Configured warning threshold. ``None`` disables size warnings.
    """

    size_bytes: int
    warning_threshold_bytes: int | None = None

    def __post_init__(self) -> None:
        """Validate byte counts."""
        if int(self.size_bytes) < 0:
            raise ValueError("EOS archive size cannot be negative")
        object.__setattr__(self, "size_bytes", int(self.size_bytes))
        if self.warning_threshold_bytes is not None:
            threshold = int(self.warning_threshold_bytes)
            if threshold <= 0:
                raise ValueError("EOS archive warning threshold must be positive")
            object.__setattr__(self, "warning_threshold_bytes", threshold)

    @property
    def size_mib(self) -> float:
        """Return the file size in mebibytes."""
        return self.size_bytes / (1024.0**2)

    @property
    def warning_threshold_mib(self) -> float | None:
        """Return the configured warning threshold in mebibytes."""
        if self.warning_threshold_bytes is None:
            return None
        return self.warning_threshold_bytes / (1024.0**2)

    @property
    def exceeds_warning_threshold(self) -> bool:
        """Return whether the archive reached the configured warning size."""
        return (
            self.warning_threshold_bytes is not None
            and self.size_bytes >= self.warning_threshold_bytes
        )

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready size description."""
        return {
            "size_bytes": self.size_bytes,
            "size_mib": self.size_mib,
            "warning_threshold_bytes": self.warning_threshold_bytes,
            "warning_threshold_mib": self.warning_threshold_mib,
            "exceeds_warning_threshold": self.exceeds_warning_threshold,
        }


@dataclass(frozen=True, slots=True)
class EOSRecordInspection:
    """Derived inspection view for one immutable EOS fit record.

    Parameters
    ----------
    record : EOSFitRecord
        Complete immutable scientific record.
    disposition : EOSRecordDisposition
        Effective current decision derived from append-only state events.
    events : tuple of EOSStateEvent
        Events directly associated with the record.
    child_record_ids : tuple of int
        Later records that identify this record as their parent.
    notes : tuple of str
        Creation and event notes in chronological order.
    used_as_initial_guess : bool
        Whether explicit reuse as an initial-guess source was recorded.
    """

    record: EOSFitRecord
    disposition: EOSRecordDisposition
    events: tuple[EOSStateEvent, ...] = ()
    child_record_ids: tuple[int, ...] = ()
    notes: tuple[str, ...] = ()
    used_as_initial_guess: bool = False

    @property
    def record_id(self) -> int:
        """Return the immutable record identifier."""
        return self.record.record_id

    @property
    def slot_key(self) -> str:
        """Return the stable result-slot key."""
        return self.record.slot.key

    @property
    def model_tag(self) -> str:
        """Return the normalized model tag used by the fit."""
        return str(getattr(self.record.request.model, "tag", self.record.request.model))

    @property
    def successful(self) -> bool:
        """Return whether the numerical fit succeeded."""
        return self.record.successful

    @property
    def is_current_accepted(self) -> bool:
        """Return whether this is the current accepted result."""
        return self.disposition is EOSRecordDisposition.ACCEPTED

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready inspection view."""
        return {
            "record_id": self.record_id,
            "dataset_id": self.record.dataset_id,
            "slot": self.record.slot.as_dict(),
            "model_tag": self.model_tag,
            "successful": self.successful,
            "fit_status": self.record.result.fit.status.value,
            "disposition": self.disposition.value,
            "parent_record_id": self.record.parent_record_id,
            "child_record_ids": list(self.child_record_ids),
            "used_as_initial_guess": self.used_as_initial_guess,
            "notes": list(self.notes),
            "events": [event.as_dict() for event in self.events],
            "created_at": self.record.created_at.isoformat(),
        }


@dataclass(frozen=True, slots=True)
class EOSSlotInspection:
    """Derived archive view for one scientific result slot.

    Parameters
    ----------
    state : EOSSlotState
        Compact current state materialized in the archive.
    records : tuple of EOSRecordInspection
        All attempts for the slot in creation order.
    """

    state: EOSSlotState
    records: tuple[EOSRecordInspection, ...] = ()

    @property
    def accepted(self) -> EOSRecordInspection | None:
        """Return the current accepted record view, when present."""
        accepted_id = self.state.accepted_record_id
        return next(
            (item for item in self.records if item.record_id == accepted_id),
            None,
        )

    @property
    def last(self) -> EOSRecordInspection | None:
        """Return the most recent attempt, when present."""
        last_id = self.state.last_record_id
        return next((item for item in self.records if item.record_id == last_id), None)

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready slot inspection."""
        return {
            "state": self.state.as_dict(),
            "records": [record.as_dict() for record in self.records],
        }


@dataclass(frozen=True, slots=True)
class EOSArchiveInspection:
    """Frontend-neutral inspection snapshot of one EOS archive.

    Parameters
    ----------
    path : Path
        Native HDF5 archive path.
    schema_version : str
        EOS archive schema version.
    dataset_ids : tuple of int
        Embedded datasets in creation order.
    slots : tuple of EOSSlotInspection
        Registered result slots and their record histories.
    event_count : int
        Number of append-only state events.
    size : EOSArchiveSizeInfo
        Current physical archive size.
    """

    path: Path
    schema_version: str
    dataset_ids: tuple[int, ...]
    slots: tuple[EOSSlotInspection, ...]
    event_count: int
    size: EOSArchiveSizeInfo

    @property
    def record_count(self) -> int:
        """Return the total number of stored fit attempts."""
        return sum(len(slot.records) for slot in self.slots)

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready archive inspection."""
        return {
            "path": str(self.path),
            "schema_version": self.schema_version,
            "dataset_ids": list(self.dataset_ids),
            "record_count": self.record_count,
            "event_count": self.event_count,
            "size": self.size.as_dict(),
            "slots": [slot.as_dict() for slot in self.slots],
        }


__all__ = [
    "EOSArchiveInspection",
    "EOSArchiveSizeInfo",
    "EOSRecordDisposition",
    "EOSRecordInspection",
    "EOSSlotInspection",
]
