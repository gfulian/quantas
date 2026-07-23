# -*- coding: utf-8 -*-

"""Workflow event serialization for native Quantas HDF5 files."""

from __future__ import annotations

import h5py

from quantas.core.events import EventRecord
from quantas.io.hdf5.attrs import (
    datetime_to_string,
    decode_text,
    numeric_sort_key,
    parse_datetime,
)
from quantas.io.hdf5.datasets import read_mapping, write_mapping


def write_events(h5: h5py.File, events: list[EventRecord]) -> h5py.Group:
    """Write persistent workflow event records.

    Parameters
    ----------
    h5 : h5py.File
        Open destination file.
    events : list of EventRecord
        Events to serialize.

    Returns
    -------
    h5py.Group
        Created events group.
    """
    group = h5.create_group("events")
    group.attrs["count"] = len(events)
    for index, event in enumerate(events):
        item = group.create_group(str(index))
        item.attrs["message"] = event.message
        item.attrs["level"] = event.level
        item.attrs["timestamp"] = datetime_to_string(event.timestamp)
        if event.progress is not None:
            item.attrs["progress"] = float(event.progress)
        data_group = item.create_group("data")
        write_mapping(data_group, event.data)
    return group


def read_events(h5: h5py.File) -> list[EventRecord]:
    """Read persistent workflow event records.

    Parameters
    ----------
    h5 : h5py.File
        Open source file.

    Returns
    -------
    list of EventRecord
        Reconstructed event records.
    """
    if "events" not in h5:
        return []
    records: list[EventRecord] = []
    for key in sorted(h5["events"], key=numeric_sort_key):
        item = h5["events"][key]
        records.append(
            EventRecord(
                message=decode_text(item.attrs.get("message", "")),
                level=decode_text(item.attrs.get("level", "info")),
                progress=(
                    float(item.attrs["progress"]) if "progress" in item.attrs else None
                ),
                data=read_mapping(item["data"]) if "data" in item else {},
                timestamp=parse_datetime(item.attrs.get("timestamp")),
            )
        )
    return records
