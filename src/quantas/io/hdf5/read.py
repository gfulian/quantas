# -*- coding: utf-8 -*-

"""Public reading helpers for native Quantas HDF5 files."""

from __future__ import annotations

from quantas.io.hdf5.attrs import (
    decode_optional_text,
    decode_scalar,
    decode_text,
    decode_value,
    numeric_sort_key,
    parse_datetime,
    require_attr,
    require_text_attr,
    require_unit,
)
from quantas.io.hdf5.datasets import (
    read_array_dataset,
    read_group_mapping,
    read_mapping,
    read_node,
    require_group,
)
from quantas.io.hdf5.envelope import (
    read_input_data,
    read_options,
    read_result_metadata,
    read_warnings,
)
from quantas.io.hdf5.events import read_events

__all__ = [
    "decode_optional_text",
    "decode_scalar",
    "decode_text",
    "decode_value",
    "numeric_sort_key",
    "parse_datetime",
    "read_array_dataset",
    "read_events",
    "read_group_mapping",
    "read_input_data",
    "read_mapping",
    "read_node",
    "read_options",
    "read_result_metadata",
    "read_warnings",
    "require_attr",
    "require_group",
    "require_text_attr",
    "require_unit",
]
