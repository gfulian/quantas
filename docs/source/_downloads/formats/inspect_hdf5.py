#!/usr/bin/env python3
"""Inspect a native Quantas HDF5 file without loading large arrays."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import h5py
import numpy as np


def _decode(value: Any) -> Any:
    """Decode bytes and NumPy scalars for readable terminal output."""
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray) and value.dtype.kind == "S":
        return np.char.decode(value, "utf-8")
    return value


def _format_value(value: Any, *, limit: int = 120) -> str:
    """Return a compact representation suitable for attributes."""
    text = repr(_decode(value))
    return text if len(text) <= limit else text[: limit - 3] + "..."


def _print_attributes(obj: h5py.Group | h5py.Dataset, indent: str) -> None:
    """Print all attributes attached to one HDF5 object."""
    for key in sorted(obj.attrs):
        print(f"{indent}@{key} = {_format_value(obj.attrs[key])}")


def _visit(name: str, obj: h5py.Group | h5py.Dataset) -> None:
    """Print one group or dataset visited by h5py."""
    depth = name.count("/")
    indent = "  " * depth
    path = "/" + name
    if isinstance(obj, h5py.Group):
        print(f"{indent}[group]   {path}")
        _print_attributes(obj, indent + "  ")
        return

    unit = _decode(obj.attrs.get("unit", ""))
    description = _decode(obj.attrs.get("description", ""))
    storage = []
    if obj.compression:
        storage.append(f"compression={obj.compression}")
    if obj.chunks:
        storage.append(f"chunks={obj.chunks}")
    details = [f"shape={obj.shape}", f"dtype={obj.dtype}"]
    if unit:
        details.append(f"unit={unit}")
    if description:
        details.append(f"description={description}")
    details.extend(storage)
    print(f"{indent}[dataset] {path} ({', '.join(details)})")
    _print_attributes(obj, indent + "  ")


def inspect(path: Path) -> None:
    """Inspect one HDF5 file and print the complete object tree."""
    with h5py.File(path, "r") as handle:
        print(f"File: {path}")
        metadata = handle.get("metadata")
        if isinstance(metadata, h5py.Group):
            module = _decode(metadata.attrs.get("module", "unknown"))
            schema = _decode(metadata.attrs.get("schema_version", "unknown"))
            version = _decode(metadata.attrs.get("version", "unknown"))
            print(f"Quantas module: {module}")
            print(f"Envelope schema: {schema}")
            print(f"Writer version: {version}")
        else:
            print("Warning: /metadata is absent; this may not be a native Quantas file.")
        print("\nHDF5 tree")
        print("---------")
        _print_attributes(handle, "")
        handle.visititems(_visit)


def main() -> int:
    """Run the command-line inspector."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path, help="Native Quantas HDF5 file")
    args = parser.parse_args()
    if not args.file.is_file():
        parser.error(f"file does not exist: {args.file}")
    try:
        inspect(args.file)
    except BrokenPipeError:
        return 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
