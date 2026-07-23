# -*- coding: utf-8 -*-

"""Dataset and recursive mapping helpers for native Quantas HDF5 files."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Mapping, cast

import h5py
import numpy as np
from numpy.typing import NDArray

from quantas.core.numerics import cast_floating_array, cast_floating_scalar
from quantas.io.hdf5.attrs import (
    INTERNAL_KIND,
    decode_scalar,
    decode_text,
    decode_value,
    numeric_sort_key,
)


def write_mapping(group: h5py.Group, values: Mapping[str, Any]) -> None:
    """Write a plain mapping recursively into a HDF5 group.

    Parameters
    ----------
    group : h5py.Group
        Destination group.
    values : Mapping[str, Any]
        Values to serialize.
    """
    for key, value in values.items():
        write_value(group, str(key), value)


def write_value(group: h5py.Group, key: str, value: Any) -> None:
    """Write one recursively serializable value.

    Parameters
    ----------
    group : h5py.Group
        Destination group.
    key : str
        Child name or attribute name.
    value : Any
        Value to serialize.
    """
    if is_dataclass(value) and not isinstance(value, type):
        value = asdict(cast(Any, value))
    if isinstance(value, Enum):
        value = value.value
    if isinstance(value, Path):
        value = str(value)
    if isinstance(value, datetime):
        value = value.isoformat()

    if value is None:
        item = group.create_group(key)
        item.attrs[INTERNAL_KIND] = "none"
        return
    if isinstance(value, Mapping):
        item = group.create_group(key)
        item.attrs[INTERNAL_KIND] = "mapping"
        write_mapping(item, value)
        return
    if isinstance(value, (list, tuple)):
        write_sequence(group, key, value)
        return
    if isinstance(value, np.ndarray):
        write_array_dataset(group, key, value)
        return
    if isinstance(value, str):
        group.attrs[key] = value
        return
    if np.isscalar(value):
        group.attrs[key] = cast_floating_scalar(value)
        return

    group.attrs[key] = str(value)


def write_numeric_attribute(
    group: h5py.Group | h5py.Dataset,
    key: str,
    value: Any,
) -> None:
    """Write one numeric HDF5 attribute using Quantas double precision."""
    group.attrs[key] = cast_floating_scalar(value)


def write_sequence(
    group: h5py.Group, key: str, value: list[Any] | tuple[Any, ...]
) -> None:
    """Write a sequence as a dataset when possible, otherwise recursively.

    Parameters
    ----------
    group : h5py.Group
        Destination group.
    key : str
        Child name.
    value : list or tuple
        Sequence to serialize.
    """
    if not value:
        item = group.create_group(key)
        item.attrs[INTERNAL_KIND] = "sequence"
        return
    try:
        array = np.asarray(value)
    except (TypeError, ValueError):
        array = np.asarray(value, dtype=object)
    if array.dtype.kind != "O":
        write_array_dataset(group, key, array)
        return
    item = group.create_group(key)
    item.attrs[INTERNAL_KIND] = "sequence"
    for index, entry in enumerate(value):
        write_value(item, str(index), entry)


def write_array_dataset(
    group: h5py.Group,
    key: str,
    value: Any,
    *,
    unit: str | None = None,
    description: str | None = None,
    compression: bool | None = None,
) -> h5py.Dataset:
    """Write an array-like value as a HDF5 dataset.

    Parameters
    ----------
    group : h5py.Group
        Destination group.
    key : str
        Dataset name.
    value : Any
        Array-like value.
    unit : str or None, optional
        Unit label stored as dataset attribute.
    description : str or None, optional
        Human-readable dataset description stored as attribute.
    compression : bool or None, optional
        If ``True``, use gzip chunked compression for non-scalar arrays. If
        ``None``, compression is enabled only for arrays with more than one
        element when ``unit`` or ``description`` is supplied by module payload
        writers that requested shared dataset options.

    Returns
    -------
    h5py.Dataset
        Created dataset.
    """
    array = np.asarray(value)
    array = cast_floating_array(array)
    if array.dtype.kind in {"U", "O"}:
        string_dtype = h5py.string_dtype(encoding="utf-8")
        array = np.asarray(array, dtype=string_dtype)

    options: dict[str, Any] = {}
    use_compression = bool(compression)
    if use_compression and array.ndim > 0 and array.size > 1:
        options.update(
            chunks=True,
            compression="gzip",
            compression_opts=4,
            shuffle=True,
        )
    dataset = group.create_dataset(key, data=array, **options)
    if unit is not None:
        dataset.attrs["unit"] = unit
    if description is not None:
        dataset.attrs["description"] = description
    return dataset


def read_mapping(group: h5py.Group) -> dict[str, Any]:
    """Read a recursively serialized mapping from a HDF5 group.

    Parameters
    ----------
    group : h5py.Group
        Source group.

    Returns
    -------
    dict
        Decoded mapping.
    """
    result = {
        key: decode_scalar(value)
        for key, value in group.attrs.items()
        if key != INTERNAL_KIND
    }
    for key, value in group.items():
        result[key] = read_node(value)
    return result


def read_group_mapping(group: h5py.Group) -> dict[str, Any]:
    """Read attributes and child nodes from a HDF5 group.

    Parameters
    ----------
    group : h5py.Group
        Source group.

    Returns
    -------
    dict
        Decoded mapping. Quantas ``none`` and ``sequence`` markers are honored.
    """
    kind = decode_text(group.attrs.get(INTERNAL_KIND, "mapping"))
    if kind == "none":
        return None  # type: ignore[return-value]
    if kind == "sequence":
        return {
            key: read_node(group[key]) for key in sorted(group, key=numeric_sort_key)
        }

    data: dict[str, Any] = {}
    for key, value in group.attrs.items():
        if key == INTERNAL_KIND:
            continue
        if key.endswith("_is_none") and bool(value):
            data[key.removesuffix("_is_none")] = None
        else:
            data[key] = decode_scalar(value)
    for key, node in group.items():
        data[key] = read_node(node)
    return data


def read_node(node: h5py.Dataset | h5py.Group) -> Any:
    """Read one recursively serialized HDF5 node.

    Parameters
    ----------
    node : h5py.Dataset or h5py.Group
        Node to decode.

    Returns
    -------
    Any
        Decoded value.
    """
    if isinstance(node, h5py.Dataset):
        value = node[()]
        decoded = decode_value(value)
        if isinstance(decoded, np.ndarray) and decoded.dtype.kind in {"S", "O", "U"}:
            return np.vectorize(decode_scalar, otypes=[object])(decoded)
        return decoded

    kind = decode_text(node.attrs.get(INTERNAL_KIND, "mapping"))
    if kind == "none":
        return None
    if kind == "sequence":
        return [read_node(node[key]) for key in sorted(node, key=numeric_sort_key)]
    return read_group_mapping(node)


def require_group(group: h5py.Group | h5py.File, key: str) -> h5py.Group:
    """Return a required child group.

    Parameters
    ----------
    group : h5py.Group or h5py.File
        Parent object.
    key : str
        Required child name.

    Returns
    -------
    h5py.Group
        Child group.

    Raises
    ------
    ValueError
        If the child is missing or is not a group.
    """
    if key not in group or not isinstance(group[key], h5py.Group):
        raise ValueError(f"Missing required HDF5 group: {group.name}/{key}")
    return group[key]


def read_array_dataset(
    group: h5py.Group,
    key: str,
    dtype: Any | None = None,
    shape: tuple[int, ...] | None = None,
    *,
    readonly: bool = False,
    required: bool = True,
) -> NDArray[Any] | None:
    """Read and optionally validate one HDF5 dataset.

    Parameters
    ----------
    group : h5py.Group
        Source group.
    key : str
        Dataset name.
    dtype : Any or None, optional
        NumPy dtype used for conversion.
    shape : tuple or None, optional
        Expected shape.
    readonly : bool, optional
        If ``True``, return a read-only NumPy copy.
    required : bool, optional
        If ``False``, return ``None`` when the dataset is absent.

    Returns
    -------
    ndarray or None
        Dataset content.

    Raises
    ------
    ValueError
        If a required dataset is missing or has an unexpected shape.
    """
    if key not in group or not isinstance(group[key], h5py.Dataset):
        if required:
            raise ValueError(f"Missing required HDF5 dataset: {group.name}/{key}")
        return None
    value = group[key][()]
    array = np.asarray(value if dtype is None else np.asarray(value, dtype=dtype))
    if shape is not None and array.shape != shape:
        raise ValueError(
            f"Dataset {group.name}/{key} has shape {array.shape}; expected {shape}."
        )
    if readonly:
        array = np.array(array, copy=True)
        array.setflags(write=False)
    return array
