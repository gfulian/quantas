# -*- coding: utf-8 -*-

"""Observers that consume frontend-neutral Quantas events.

The implementations support callbacks, in-memory collection, and no-op handling
without introducing user-interface dependencies into scientific workflows.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Protocol

from .event import Event


class Observer(Protocol):
    """
    Protocol for Quantas event observers.

    Any callable object accepting an :class:`Event` instance can be used as an
    observer. This makes it possible to connect calculators to command-line
    loggers, graphical widgets, notebooks, lists, or custom user callbacks.
    """

    def __call__(self, event: Event) -> None:
        """
        Receive an event emitted by a workflow.

        Parameters
        ----------
        event : Event
            Event emitted by a Quantas calculator.
        """
        ...


class NullObserver:
    """
    Observer that ignores all received events.

    This class is used as the default observer when no logging, reporting, or
    graphical update is required.
    """

    def __call__(self, event: Event) -> None:
        """
        Receive and ignore an event.

        Parameters
        ----------
        event : Event
            Event emitted by a Quantas calculator.
        """
        return None


@dataclass
class ListObserver:
    """
    Observer that stores all received events in a list.

    This observer is mainly useful for testing, debugging, and workflows where
    the event log has to be saved after the calculation.

    Attributes
    ----------
    events : list of Event
        List containing the received events.
    """

    events: list[Event] = field(default_factory=list)

    def __call__(self, event: Event) -> None:
        """
        Store an event in the internal list.

        Parameters
        ----------
        event : Event
            Event emitted by a Quantas calculator.
        """
        self.events.append(event)


@dataclass
class CallbackObserver:
    """
    Observer that forwards events to a user-defined callback.

    Parameters
    ----------
    callback : callable
        Function or callable object that receives an :class:`Event` instance.
    """

    callback: Callable[[Event], None]

    def __call__(self, event: Event) -> None:
        """
        Forward an event to the stored callback.

        Parameters
        ----------
        event : Event
            Event emitted by a Quantas calculator.
        """
        self.callback(event)
