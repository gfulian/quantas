# -*- coding: utf-8 -*-

"""Common frontend-neutral contracts exposed by the Quantas public API."""

from __future__ import annotations

from typing import TypeVar

from quantas.core.events import (
    CallbackObserver,
    Event,
    EventLevel,
    EventRecord,
    ListObserver,
    NullObserver,
    Observer,
)
from quantas.models import InputData, ReportTable, ResultData
from quantas.models.phonons import PhononInputData

from .plotting import PlotCollection

PayloadT = TypeVar("PayloadT")


def get_result_payload(
    result: ResultData,
    *,
    module: str,
    key: str,
    expected_type: type[PayloadT],
) -> PayloadT:
    """Return and type-check one module-specific result payload.

    Parameters
    ----------
    result : ResultData
        Complete Quantas result envelope.
    module : str
        Expected stable module identifier.
    key : str
        Payload key inside ``result.results``.
    expected_type : type
        Expected payload class.

    Returns
    -------
    PayloadT
        Validated module-specific payload.

    Raises
    ------
    TypeError
        If ``result`` is not a :class:`ResultData` object.
    ValueError
        If module metadata, payload key, or payload type is invalid.
    """
    if not isinstance(result, ResultData):
        raise TypeError("result must be a ResultData object")
    if result.metadata.module != module:
        raise ValueError(
            f"expected a {module!r} result, got {result.metadata.module!r}"
        )
    payload = result.results.get(key)
    if not isinstance(payload, expected_type):
        raise ValueError(
            f"result does not contain a valid {expected_type.__name__} payload"
        )
    return payload


def _public_dir(names: list[str] | tuple[str, ...]) -> list[str]:
    """Return a stable autocomplete view for one public namespace."""
    return sorted(names)


def __dir__() -> list[str]:
    """Return only supported common public names."""
    return _public_dir(__all__)


__all__ = [
    "CallbackObserver",
    "Event",
    "EventLevel",
    "EventRecord",
    "InputData",
    "ListObserver",
    "NullObserver",
    "Observer",
    "PhononInputData",
    "PlotCollection",
    "ReportTable",
    "ResultData",
    "get_result_payload",
]
