# -*- coding: utf-8 -*-

"""Passive history and state contracts for persistent EOS archives."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Any

from .models import EOSFitDomain, EOSFitRequest, EOSFitResult


class EOSSlotStatus(str, Enum):
    """Current processing state of one scientific result slot.

    Attributes
    ----------
    NOT_PROCESSED
        No fit attempt has been recorded for the slot.
    ATTEMPTED
        One or more attempts exist, but no result is currently accepted.
    ACCEPTED
        A successful record is selected as the current scientific result.
    """

    NOT_PROCESSED = "not_processed"
    ATTEMPTED = "attempted_no_accepted_result"
    ACCEPTED = "accepted"


class EOSStateEventType(str, Enum):
    """Append-only state changes stored in an EOS archive."""

    RECORD_CREATED = "record_created"
    RECORD_ACCEPTED = "record_accepted"
    RECORD_UNACCEPTED = "record_unaccepted"
    RECORD_REJECTED = "record_rejected"
    RECORD_CANDIDATE = "record_candidate"
    RECORD_USED_AS_INITIAL_GUESS = "record_used_as_initial_guess"
    NOTE_ADDED = "note_added"
    NO_OPERATION = "no_operation"


@dataclass(frozen=True, slots=True)
class EOSResultSlot:
    """Identify one accepted-result location by scientific domain and target.

    Parameters
    ----------
    domain : EOSFitDomain or str
        Scientific fit domain.
    target : str
        Fitted property, for example ``"volume"`` or ``"c"``.
    """

    domain: EOSFitDomain
    target: str

    def __post_init__(self) -> None:
        """Normalize and validate the slot identifier."""
        object.__setattr__(self, "domain", EOSFitDomain(self.domain))
        target = str(self.target).strip().lower()
        if not target or "/" in target:
            raise ValueError("EOS result-slot target must be a non-empty path segment")
        object.__setattr__(self, "target", target)

    @property
    def key(self) -> str:
        """Return the stable ``domain/target`` archive key."""
        return f"{self.domain.value}/{self.target}"

    @classmethod
    def parse(cls, value: str | EOSResultSlot) -> EOSResultSlot:
        """Return a canonical slot from a key or existing slot.

        Parameters
        ----------
        value : str or EOSResultSlot
            Slot key such as ``"pv/volume"``.

        Returns
        -------
        EOSResultSlot
            Canonical slot object.
        """
        if isinstance(value, cls):
            return value
        parts = str(value).strip().lower().split("/")
        if len(parts) != 2:
            raise ValueError("EOS result slot must use the form 'domain/target'")
        return cls(EOSFitDomain(parts[0]), parts[1])

    @classmethod
    def from_request(cls, request: EOSFitRequest) -> EOSResultSlot:
        """Construct the slot addressed by a fit request."""
        return cls(request.domain, request.target)

    def as_dict(self) -> dict[str, str]:
        """Return a serialization-ready slot mapping."""
        return {"domain": self.domain.value, "target": self.target, "key": self.key}


@dataclass(frozen=True, slots=True)
class EOSFitRecord:
    """Immutable scientific record for one executed EOS fit.

    Parameters
    ----------
    record_id : int
        Monotonic archive-global record identifier.
    dataset_id : int
        Dataset used by the request.
    request : EOSFitRequest
        Complete normalized fit request.
    result : EOSFitResult
        Complete numerical result, including failures.
    parent_record_id : int or None, optional
        Earlier record used as the scientific parent or initial-guess source.
    created_at : datetime, optional
        UTC-aware record creation time.
    note : str or None, optional
        Optional creation note. Later notes are append-only state events.
    provenance : dict, optional
        Passive caller and workflow provenance.

    Raises
    ------
    ValueError
        If identifiers are invalid or request and result are inconsistent.
    """

    record_id: int
    dataset_id: int
    request: EOSFitRequest
    result: EOSFitResult
    parent_record_id: int | None = None
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    note: str | None = None
    provenance: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate immutable record identity and scientific consistency."""
        if int(self.record_id) <= 0 or int(self.dataset_id) <= 0:
            raise ValueError("EOS record and dataset identifiers must be positive")
        object.__setattr__(self, "record_id", int(self.record_id))
        object.__setattr__(self, "dataset_id", int(self.dataset_id))
        if self.parent_record_id is not None:
            parent = int(self.parent_record_id)
            if parent <= 0 or parent >= self.record_id:
                raise ValueError("parent_record_id must identify an earlier record")
            object.__setattr__(self, "parent_record_id", parent)
        if self.result.request.as_dict() != self.request.as_dict():
            raise ValueError("EOS fit record request must match result.request")
        timestamp = self.created_at
        if timestamp.tzinfo is None:
            timestamp = timestamp.replace(tzinfo=timezone.utc)
        object.__setattr__(self, "created_at", timestamp)
        object.__setattr__(self, "note", None if self.note is None else str(self.note))
        object.__setattr__(self, "provenance", dict(self.provenance))

    @property
    def slot(self) -> EOSResultSlot:
        """Return the scientific result slot addressed by this record."""
        return EOSResultSlot.from_request(self.request)

    @property
    def successful(self) -> bool:
        """Return whether the stored numerical result is usable."""
        return bool(self.result.fit.success)

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready complete record."""
        return {
            "record_id": self.record_id,
            "dataset_id": self.dataset_id,
            "parent_record_id": self.parent_record_id,
            "created_at": self.created_at.isoformat(),
            "slot": self.slot.as_dict(),
            "note": self.note,
            "provenance": dict(self.provenance),
            "request": self.request.as_dict(),
            "result": self.result.as_dict(),
        }


@dataclass(frozen=True, slots=True)
class EOSStateEvent:
    """Append-only decision or annotation associated with an EOS archive."""

    event_id: int
    event_type: EOSStateEventType
    record_id: int | None = None
    slot: EOSResultSlot | None = None
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    note: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize event identifiers, enums, and passive metadata."""
        if int(self.event_id) <= 0:
            raise ValueError("EOS state-event identifier must be positive")
        object.__setattr__(self, "event_id", int(self.event_id))
        object.__setattr__(self, "event_type", EOSStateEventType(self.event_type))
        if self.record_id is not None:
            record_id = int(self.record_id)
            if record_id <= 0:
                raise ValueError("EOS state-event record_id must be positive")
            object.__setattr__(self, "record_id", record_id)
        if self.slot is not None:
            object.__setattr__(self, "slot", EOSResultSlot.parse(self.slot))
        timestamp = self.created_at
        if timestamp.tzinfo is None:
            timestamp = timestamp.replace(tzinfo=timezone.utc)
        object.__setattr__(self, "created_at", timestamp)
        object.__setattr__(self, "note", None if self.note is None else str(self.note))
        object.__setattr__(self, "metadata", dict(self.metadata))

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready state event."""
        return {
            "event_id": self.event_id,
            "event_type": self.event_type.value,
            "record_id": self.record_id,
            "slot": None if self.slot is None else self.slot.as_dict(),
            "created_at": self.created_at.isoformat(),
            "note": self.note,
            "metadata": dict(self.metadata),
        }


@dataclass(frozen=True, slots=True)
class EOSSlotState:
    """Current compact state for one result slot.

    Parameters
    ----------
    slot : EOSResultSlot
        Scientific domain and target.
    status : EOSSlotStatus
        Current processing state.
    accepted_record_id, last_record_id : int or None
        Current accepted result and most recent attempt.
    attempted_record_ids : tuple of int, optional
        All attempts addressing the slot in creation order.
    """

    slot: EOSResultSlot
    status: EOSSlotStatus = EOSSlotStatus.NOT_PROCESSED
    accepted_record_id: int | None = None
    last_record_id: int | None = None
    attempted_record_ids: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        """Normalize current state and validate identifier consistency."""
        object.__setattr__(self, "slot", EOSResultSlot.parse(self.slot))
        object.__setattr__(self, "status", EOSSlotStatus(self.status))
        attempts = tuple(int(value) for value in self.attempted_record_ids)
        if any(value <= 0 for value in attempts):
            raise ValueError("attempted EOS record identifiers must be positive")
        object.__setattr__(self, "attempted_record_ids", attempts)
        for name in ("accepted_record_id", "last_record_id"):
            value = getattr(self, name)
            if value is not None:
                normalized = int(value)
                if normalized <= 0:
                    raise ValueError(f"{name} must be positive")
                object.__setattr__(self, name, normalized)
        if self.status is EOSSlotStatus.ACCEPTED and self.accepted_record_id is None:
            raise ValueError("accepted slot state requires accepted_record_id")
        if (
            self.accepted_record_id is not None
            and self.accepted_record_id not in attempts
        ):
            raise ValueError("accepted record must appear among attempted records")
        if self.last_record_id is not None and self.last_record_id not in attempts:
            raise ValueError("last record must appear among attempted records")

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-ready compact slot state."""
        return {
            "slot": self.slot.as_dict(),
            "status": self.status.value,
            "accepted_record_id": self.accepted_record_id,
            "last_record_id": self.last_record_id,
            "attempted_record_ids": list(self.attempted_record_ids),
        }


__all__ = [
    "EOSFitRecord",
    "EOSResultSlot",
    "EOSSlotState",
    "EOSSlotStatus",
    "EOSStateEvent",
    "EOSStateEventType",
]
