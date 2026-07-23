# -*- coding: utf-8 -*-

"""Convert and validate scalar attributes in native Quantas HDF5 files.

The helpers are shared by all readers and writers and remain independent of
module-specific result schemas.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

import h5py
import numpy as np

INTERNAL_KIND = "__quantas_kind__"


def decode_scalar(value: Any) -> Any:
    """Decode byte and NumPy scalar values.

    Parameters
    ----------
    value : Any
        HDF5 attribute or scalar dataset value.

    Returns
    -------
    Any
        Python scalar, decoded text, or the original value when no conversion is
        needed.
    """
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        return value.item()
    return value


def decode_value(value: Any) -> Any:
    """Decode a generic HDF5 value.

    Parameters
    ----------
    value : Any
        Scalar or array value returned by h5py.

    Returns
    -------
    Any
        Decoded Python/NumPy value. Byte string arrays are converted to Unicode
        arrays; object arrays are decoded element-wise.
    """
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8")
    if isinstance(value, np.ndarray):
        if value.dtype.kind == "S":
            return np.char.decode(value, "utf-8")
        if value.dtype.kind == "O":
            flat = [decode_value(item) for item in value.ravel()]
            return np.asarray(flat, dtype=object).reshape(value.shape)
        return value
    if isinstance(value, np.generic):
        return value.item()
    return value


def decode_text(value: Any) -> str:
    """Decode one HDF5 scalar as text.

    Parameters
    ----------
    value : Any
        HDF5 scalar or attribute value.

    Returns
    -------
    str
        Decoded text representation.
    """
    return str(decode_scalar(value))


def decode_optional_text(value: Any) -> str | None:
    """Decode optional text and normalize schema placeholders.

    Parameters
    ----------
    value : Any
        HDF5 scalar, attribute value, or ``None``.

    Returns
    -------
    str or None
        Decoded text, or ``None`` for missing/placeholder values.
    """
    if value is None:
        return None
    text = decode_text(value)
    return None if text in {"", "none", "unknown"} else text


def datetime_to_string(value: datetime | str | None) -> str:
    """Return a UTC-aware ISO-8601 timestamp string.

    Parameters
    ----------
    value : datetime, str, or None
        Timestamp to serialize.

    Returns
    -------
    str
        ISO-8601 timestamp string.
    """
    if value is None:
        return datetime.now(timezone.utc).isoformat()
    if isinstance(value, datetime):
        return value.isoformat()
    return str(value)


def parse_datetime(value: Any) -> datetime:
    """Parse a stored ISO-8601 timestamp.

    Parameters
    ----------
    value : Any
        Stored timestamp value.

    Returns
    -------
    datetime
        Parsed timestamp.

    Raises
    ------
    ValueError
        If the timestamp text is not valid ISO-8601.
    """
    if value is None:
        return datetime.now(timezone.utc)
    text = decode_text(value)
    try:
        return datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError(f"Invalid Quantas result timestamp: {text}") from exc


def numeric_sort_key(value: str) -> tuple[int, str]:
    """Sort numerical HDF5 group names before arbitrary names.

    Parameters
    ----------
    value : str
        HDF5 child name.

    Returns
    -------
    tuple
        Numeric-aware sorting key.
    """
    try:
        return int(value), value
    except ValueError:
        return 2**31 - 1, value


def require_attr(group: h5py.Group | h5py.Dataset, key: str) -> Any:
    """Return a required HDF5 attribute.

    Parameters
    ----------
    group : h5py.Group or h5py.Dataset
        Object containing the attribute.
    key : str
        Attribute name.

    Returns
    -------
    Any
        Attribute value.

    Raises
    ------
    ValueError
        If the attribute is missing.
    """
    if key not in group.attrs:
        raise ValueError(f"Missing required HDF5 attribute: {group.name}:{key}")
    return group.attrs[key]


def require_text_attr(
    group: h5py.Group | h5py.Dataset,
    key: str,
    expected: str,
) -> None:
    """Validate one textual HDF5 attribute.

    Parameters
    ----------
    group : h5py.Group or h5py.Dataset
        Object containing the attribute.
    key : str
        Attribute name.
    expected : str
        Required decoded value.

    Raises
    ------
    ValueError
        If the attribute is missing or has a different value.
    """
    value = decode_text(require_attr(group, key))
    if value != expected:
        raise ValueError(
            f"Unsupported value for {group.name}:{key}: {value!r}; "
            f"expected {expected!r}."
        )


def require_unit(dataset: h5py.Dataset, expected: str) -> None:
    """Validate the physical unit attached to a dataset.

    Parameters
    ----------
    dataset : h5py.Dataset
        Dataset whose ``unit`` attribute is checked.
    expected : str
        Required unit label.

    Raises
    ------
    ValueError
        If the unit is missing or different.
    """
    require_text_attr(dataset, "unit", expected)
