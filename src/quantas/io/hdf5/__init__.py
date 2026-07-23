# -*- coding: utf-8 -*-

"""Shared HDF5 helpers for native Quantas result files.

The package is split into small modules by responsibility:

``read``
    Public helpers for generic reading, decoding and validation.
``write``
    Public helpers for generic writing.
``envelope``
    Metadata, input, options, diagnostics and warning groups.
``events``
    Workflow event serialization.
``datasets`` and ``attrs``
    Low-level dataset and attribute mechanics.

Module-specific scientific payloads remain in the corresponding workflow
packages.  The re-exports below preserve the historical
``from quantas.io.hdf5 import ...`` API.
"""

from __future__ import annotations

from quantas.io.hdf5.attrs import (
    INTERNAL_KIND,
    datetime_to_string,
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
    write_array_dataset,
    write_mapping,
    write_numeric_attribute,
    write_sequence,
    write_value,
)
from quantas.io.hdf5.envelope import (
    read_input_data,
    read_options,
    read_precision_metadata,
    read_result_metadata,
    read_warnings,
    write_diagnostics,
    write_input_data,
    write_options,
    write_precision_metadata,
    write_result_metadata,
)
from quantas.io.hdf5.events import read_events, write_events

__all__ = [
    "INTERNAL_KIND",
    "datetime_to_string",
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
    "read_precision_metadata",
    "read_result_metadata",
    "read_warnings",
    "require_attr",
    "require_group",
    "require_text_attr",
    "require_unit",
    "write_array_dataset",
    "write_diagnostics",
    "write_events",
    "write_input_data",
    "write_mapping",
    "write_numeric_attribute",
    "write_options",
    "write_precision_metadata",
    "write_result_metadata",
    "write_sequence",
    "write_value",
]
