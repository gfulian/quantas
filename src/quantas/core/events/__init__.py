# -*- coding: utf-8 -*-

"""Frontend-neutral events and observers for Quantas workflows.

Calculators emit structured events; observers decide whether to display, store,
forward, or ignore them.
"""

from .event import Event, EventLevel, EventRecord
from .observer import CallbackObserver, ListObserver, NullObserver, Observer

__all__ = [
    "Event",
    "EventLevel",
    "EventRecord",
    "Observer",
    "NullObserver",
    "ListObserver",
    "CallbackObserver",
]
