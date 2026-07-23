# -*- coding: utf-8 -*-

"""Shared HDF5 envelope for native Quantas result files."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import h5py
import numpy as np

from quantas.io.hdf5.attrs import (
    datetime_to_string,
    decode_optional_text,
    decode_text,
    parse_datetime,
)
from quantas.io.hdf5.datasets import read_mapping, read_node, write_mapping
from quantas.core.numerics import NumericPrecisionPolicy
from quantas.models.data import InputData, ResultData, ResultMetadata


def write_result_metadata(h5: h5py.File, result: ResultData) -> h5py.Group:
    """Write generic Quantas result metadata.

    Parameters
    ----------
    h5 : h5py.File
        Open destination file.
    result : ResultData
        Generic result container.

    Returns
    -------
    h5py.Group
        Created metadata group.
    """
    group = h5.create_group("metadata")
    metadata = result.metadata
    group.attrs["program"] = metadata.program
    group.attrs["module"] = metadata.module
    group.attrs["method"] = metadata.method
    group.attrs["version"] = metadata.version
    group.attrs["schema_version"] = metadata.schema_version
    group.attrs["created_by"] = metadata.created_by or "unknown"
    group.attrs["created_at"] = datetime_to_string(metadata.created_at)
    return group


def read_result_metadata(h5: h5py.File) -> ResultMetadata:
    """Read generic Quantas result metadata.

    Parameters
    ----------
    h5 : h5py.File
        Open source file.

    Returns
    -------
    ResultMetadata
        Reconstructed metadata.

    Raises
    ------
    ValueError
        If the metadata group or timestamp is invalid.
    """
    if "metadata" not in h5:
        raise ValueError("Missing metadata group in Quantas result file.")
    attrs = h5["metadata"].attrs
    return ResultMetadata(
        program=decode_text(attrs.get("program", "quantas")),
        module=decode_text(attrs.get("module", "unknown")),
        method=decode_text(attrs.get("method", "unknown")),
        version=decode_text(attrs.get("version", "unknown")),
        created_by=decode_optional_text(attrs.get("created_by")),
        created_at=parse_datetime(attrs.get("created_at")),
        schema_version=decode_text(attrs.get("schema_version", "unknown")),
    )


def write_precision_metadata(h5: h5py.File) -> h5py.Group:
    """Write the numerical precision policy into the metadata group.

    Parameters
    ----------
    h5 : h5py.File
        Open destination file containing the standard metadata group.
    Returns
    -------
    h5py.Group
        Created ``metadata/numerics`` group.
    """
    metadata = h5.require_group("metadata")
    group = metadata.create_group("numerics")
    for key, value in NumericPrecisionPolicy().as_metadata().items():
        group.attrs[key] = value
    return group


def read_precision_metadata(h5: h5py.File) -> dict[str, str | int]:
    """Read numerical precision metadata from a native Quantas file.

    Historical files are reported exactly as stored, including the short-lived
    experimental single-precision storage metadata.  Missing metadata uses the
    original Quantas double-precision default.
    """
    default = NumericPrecisionPolicy().as_metadata()
    path = "metadata/numerics"
    if path not in h5:
        return default
    attrs = h5[path].attrs
    return {
        "working": decode_text(attrs.get("working", default["working"])),
        "working_dtype": decode_text(
            attrs.get("working_dtype", default["working_dtype"])
        ),
        "working_bits": int(attrs.get("working_bits", default["working_bits"])),
        "storage": decode_text(attrs.get("storage", default["storage"])),
        "storage_dtype": decode_text(
            attrs.get("storage_dtype", default["storage_dtype"])
        ),
        "storage_bits": int(attrs.get("storage_bits", default["storage_bits"])),
    }


def write_input_data(h5: h5py.File, input_data: InputData | None) -> h5py.Group:
    """Write normalized workflow input data.

    Parameters
    ----------
    h5 : h5py.File
        Open destination file.
    input_data : InputData or None
        Generic normalized input container.

    Returns
    -------
    h5py.Group
        Created input group.
    """
    group = h5.create_group("input")
    if input_data is None:
        return group
    if input_data.source is not None:
        group.attrs["source"] = str(input_data.source)
    if input_data.raw is not None:
        group.create_dataset("raw", data=input_data.raw)
    data_group = group.create_group("data")
    write_mapping(data_group, input_data.data)
    return group


def read_input_data(
    h5: h5py.File, source: str | Path | None = None
) -> InputData | None:
    """Read normalized workflow input data.

    Parameters
    ----------
    h5 : h5py.File
        Open source file.
    source : str or Path or None, optional
        Fallback source description.

    Returns
    -------
    InputData or None
        Reconstructed input data, or ``None`` when no input was stored.
    """
    if "input" not in h5:
        return None
    group = h5["input"]
    stored_source = group.attrs.get("source", source)
    raw = read_node(group["raw"]) if "raw" in group else None
    data = read_mapping(group["data"]) if "data" in group else {}
    return InputData(
        source=None if stored_source is None else decode_text(stored_source),
        raw=None if raw is None else str(raw),
        data=data,
    )


def write_options(h5: h5py.File, options: Mapping[str, Any]) -> h5py.Group:
    """Write workflow options to the standard options group.

    Parameters
    ----------
    h5 : h5py.File
        Open destination file.
    options : Mapping[str, Any]
        Options to serialize.

    Returns
    -------
    h5py.Group
        Created options group.
    """
    group = h5.create_group("options")
    write_mapping(group, options)
    return group


def read_options(h5: h5py.File) -> dict[str, Any]:
    """Read workflow options from the standard options group.

    Parameters
    ----------
    h5 : h5py.File
        Open source file.

    Returns
    -------
    dict
        Decoded options mapping.
    """
    return read_mapping(h5["options"]) if "options" in h5 else {}


def write_diagnostics(
    h5: h5py.File,
    result: ResultData,
    report_text: str | None = None,
) -> h5py.Group:
    """Write generic warnings and optional frontend report text.

    Parameters
    ----------
    h5 : h5py.File
        Open destination file.
    result : ResultData
        Generic result container.
    report_text : str or None, optional
        Report text produced by a frontend renderer.

    Returns
    -------
    h5py.Group
        Created diagnostics group. Module-specific diagnostics can be added as
        children by module payload writers.
    """
    group = h5.create_group("diagnostics")
    if result.warnings:
        text_dtype = h5py.string_dtype(encoding="utf-8")
        group.create_dataset(
            "warnings",
            data=np.asarray(result.warnings, dtype=text_dtype),
        )
    if report_text is not None:
        group.create_dataset("report_text", data=report_text)
    return group


def read_warnings(h5: h5py.File) -> list[str]:
    """Read generic workflow warnings from the diagnostics group.

    Parameters
    ----------
    h5 : h5py.File
        Open source file.

    Returns
    -------
    list of str
        Warning messages.
    """
    if "diagnostics" not in h5 or "warnings" not in h5["diagnostics"]:
        return []
    values = np.atleast_1d(h5["diagnostics/warnings"][()])
    return [decode_text(value) for value in values]
