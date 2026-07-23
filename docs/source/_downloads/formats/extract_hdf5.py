#!/usr/bin/env python3
"""Extract selected datasets from a native Quantas HDF5 file to NPZ."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np


def _decode(value: Any) -> Any:
    """Convert HDF5 scalar metadata to JSON-compatible Python values."""
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        if value.dtype.kind == "S":
            return np.char.decode(value, "utf-8").tolist()
        return value.tolist()
    return value




def _portable_array(value: Any) -> np.ndarray:
    """Return an NPZ-safe array without object-dtype string payloads."""
    array = np.asarray(value)
    if array.dtype.kind == "S":
        return np.char.decode(array, "utf-8")
    if array.dtype.kind != "O":
        return array

    decoded: list[str] = []
    for item in array.ravel():
        if isinstance(item, (bytes, np.bytes_)):
            decoded.append(item.decode("utf-8", errors="replace"))
        elif isinstance(item, str):
            decoded.append(item)
        else:
            raise TypeError(
                "object-dtype dataset contains non-text values and cannot be "
                "stored safely without pickle"
            )
    return np.asarray(decoded, dtype=np.str_).reshape(array.shape)

def _npz_key(path: str) -> str:
    """Map an absolute HDF5 path to a stable NPZ member name."""
    return path.strip("/").replace("/", "__") or "root"


def _dataset_paths(handle: h5py.File) -> list[str]:
    """Return all dataset paths in deterministic order."""
    paths: list[str] = []

    def collect(name: str, obj: h5py.Group | h5py.Dataset) -> None:
        if isinstance(obj, h5py.Dataset):
            paths.append("/" + name)

    handle.visititems(collect)
    return sorted(paths)


def extract(
    source: Path,
    destination: Path,
    requested: list[str],
    *,
    extract_all: bool,
) -> tuple[Path, Path]:
    """Extract arrays and write a JSON manifest beside the NPZ file."""
    arrays: dict[str, np.ndarray] = {}
    manifest: dict[str, Any] = {
        "source": str(source),
        "root_attributes": {},
        "groups": {},
        "datasets": {},
        "metadata": {},
    }

    with h5py.File(source, "r") as handle:
        manifest["root_attributes"] = {
            key: _decode(value) for key, value in handle.attrs.items()
        }

        def collect_group(name: str, node: h5py.Group | h5py.Dataset) -> None:
            if isinstance(node, h5py.Group):
                manifest["groups"]["/" + name] = {
                    "attributes": {
                        key: _decode(value) for key, value in node.attrs.items()
                    }
                }

        handle.visititems(collect_group)

        metadata = handle.get("metadata")
        if isinstance(metadata, h5py.Group):
            manifest["metadata"] = {
                key: _decode(value) for key, value in metadata.attrs.items()
            }

        paths = _dataset_paths(handle) if extract_all else requested
        if not paths:
            raise ValueError("no datasets requested; use --dataset or --all")

        for raw_path in paths:
            path = raw_path if raw_path.startswith("/") else "/" + raw_path
            if path not in handle:
                raise KeyError(f"dataset does not exist: {path}")
            node = handle[path]
            if not isinstance(node, h5py.Dataset):
                raise TypeError(f"requested path is not a dataset: {path}")

            value = _portable_array(node[()])
            key = _npz_key(path)
            arrays[key] = value
            manifest["datasets"][key] = {
                "path": path,
                "shape": list(value.shape),
                "dtype": str(value.dtype),
                "unit": _decode(node.attrs.get("unit")),
                "description": _decode(node.attrs.get("description")),
                "attributes": {
                    name: _decode(item) for name, item in node.attrs.items()
                },
            }

    destination = destination.with_suffix(".npz")
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(destination, **arrays)
    manifest_path = destination.with_suffix(".json")
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return destination, manifest_path


def main() -> int:
    """Run the dataset extractor."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path, help="Native Quantas HDF5 file")
    parser.add_argument(
        "--dataset",
        action="append",
        default=[],
        help="Absolute HDF5 dataset path; repeat for multiple arrays",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Extract every dataset (may require substantial memory)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("quantas_extract.npz"),
        help="Destination NPZ path",
    )
    args = parser.parse_args()
    if not args.file.is_file():
        parser.error(f"file does not exist: {args.file}")
    try:
        npz_path, manifest_path = extract(
            args.file,
            args.output,
            args.dataset,
            extract_all=args.all,
        )
    except (KeyError, TypeError, ValueError) as exc:
        parser.error(str(exc))
    print(f"Arrays:   {npz_path}")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
