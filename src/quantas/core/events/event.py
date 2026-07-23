# -*- coding: utf-8 -*-

"""Structured event records emitted by Quantas workflows.

Events report progress, diagnostics, warnings, and results without coupling a
calculator to a terminal, notebook, or graphical interface.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Any

import numpy as np


class EventLevel(Enum):
    """
    Enumeration of the event levels supported by Quantas.

    Attributes
    ----------
    DEBUG
        Event used to report detailed information useful during development.
    INFO
        Event used to report standard information during a workflow.
    WARNING
        Event used to report a non-critical problem.
    ERROR
        Event used to report an error.
    PROGRESS
        Event used to report the progress of a workflow.
    RESULT
        Event used to report that a result, or part of it, is available.
    """

    DEBUG = "debug"
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    PROGRESS = "progress"
    RESULT = "result"


@dataclass(slots=True, frozen=True)
class Event:
    """
    Message emitted by a Quantas workflow.

    Parameters
    ----------
    message : str
        Textual message associated with the event.
    level : EventLevel, optional
        Type or severity level of the event.
    progress : float or None, optional
        Progress value associated with the event. When provided, it must be
        between 0 and 1.
    data : dict, optional
        Optional payload associated with the event. This can be used to store
        additional information, such as elapsed time, current step, or partial
        results.
    timestamp : datetime, optional
        Date and time at which the event was created.
    """

    message: str
    level: EventLevel = EventLevel.INFO
    progress: float | None = None
    data: dict[str, Any] = field(default_factory=dict)
    timestamp: datetime = field(default_factory=lambda: datetime.now(timezone.utc))

    def __post_init__(self) -> None:
        """
        Validate the event after initialization.

        Raises
        ------
        ValueError
            If the progress value is not between 0 and 1.
        """
        if self.progress is not None and not 0.0 <= self.progress <= 1.0:
            raise ValueError("progress must be between 0.0 and 1.0")


@dataclass(slots=True, frozen=True)
class EventRecord:
    """Serializable record of a Quantas workflow event.

    Parameters
    ----------
    message : str
        Textual event message.
    level : str
        Event level stored as its stable string value.
    progress : float or None, optional
        Numerical progress between zero and one.
    data : dict, optional
        Lightweight structured payload. Numerical arrays and active objects are
        represented by summaries to avoid duplicating scientific result data.
    timestamp : datetime, optional
        Time at which the original event was emitted.
    """

    message: str
    level: str = EventLevel.INFO.value
    progress: float | None = None
    data: dict[str, Any] = field(default_factory=dict)
    timestamp: datetime = field(default_factory=lambda: datetime.now(timezone.utc))

    @classmethod
    def from_event(cls, event: Event) -> "EventRecord":
        """Create a persistent record from an operational event.

        Parameters
        ----------
        event : Event
            Operational event emitted by a calculator.

        Returns
        -------
        EventRecord
            Frontend-neutral, serializable event record.
        """
        return cls(
            message=event.message,
            level=event.level.value,
            progress=event.progress,
            data=_summarize_event_value(event.data),
            timestamp=event.timestamp,
        )


def _summarize_event_value(value: Any) -> Any:
    """Return a lightweight serializable representation of event data."""
    from dataclasses import is_dataclass
    from enum import Enum

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, datetime):
        return value.isoformat()
    if isinstance(value, np.ndarray):
        return {
            "type": "ndarray",
            "shape": list(value.shape),
            "dtype": str(value.dtype),
        }
    if isinstance(value, dict):
        return {str(key): _summarize_event_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        if len(value) <= 32:
            return [_summarize_event_value(item) for item in value]
        return {"type": type(value).__name__, "length": len(value)}
    if is_dataclass(value):
        return {"type": type(value).__name__}
    return {"type": type(value).__name__}
