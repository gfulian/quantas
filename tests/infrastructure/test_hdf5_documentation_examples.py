"""Execution tests for the standalone native-HDF5 documentation utilities."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
INSPECTOR = ROOT / "examples/io/inspect_hdf5.py"
EXTRACTOR = ROOT / "examples/io/extract_hdf5.py"


def _native_fixture(path: Path) -> None:
    """Write a small native-like HDF5 tree for utility tests."""
    with h5py.File(path, "w") as handle:
        metadata = handle.create_group("metadata")
        metadata.attrs["program"] = "quantas"
        metadata.attrs["version"] = "2.0.0b6"
        metadata.attrs["module"] = "qha"
        metadata.attrs["schema_version"] = "2.0"
        results = handle.create_group("results")
        temperature = results.create_dataset(
            "temperature", data=np.asarray([0.0, 300.0], dtype=np.float64)
        )
        temperature.attrs["unit"] = "K"
        temperature.attrs["description"] = "Temperature grid"
        volume = results.create_dataset(
            "equilibrium_volume",
            data=np.asarray([[19.0], [19.1]], dtype=np.float64),
            compression="gzip",
        )
        volume.attrs["unit"] = "angstrom^3"
        volume.attrs["description"] = "Equilibrium volume"
        results.create_dataset("valid_mask", data=np.asarray([[True], [False]]))
        results.create_dataset(
            "labels",
            data=np.asarray(["state-a", "state-b"], dtype=h5py.string_dtype("utf-8")),
        )


def test_inspector_reports_schema_shapes_and_units(tmp_path: Path) -> None:
    """The inspector traverses metadata without printing array contents."""
    source = tmp_path / "result.hdf5"
    _native_fixture(source)
    completed = subprocess.run(
        [sys.executable, str(INSPECTOR), str(source)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "Quantas module: qha" in completed.stdout
    assert "Envelope schema: 2.0" in completed.stdout
    assert "/results/equilibrium_volume" in completed.stdout
    assert "shape=(2, 1)" in completed.stdout
    assert "unit=angstrom^3" in completed.stdout


def test_extractor_writes_npz_and_metadata_manifest(tmp_path: Path) -> None:
    """Selected datasets retain source paths and unit metadata."""
    source = tmp_path / "result.hdf5"
    destination = tmp_path / "subset.npz"
    _native_fixture(source)
    subprocess.run(
        [
            sys.executable,
            str(EXTRACTOR),
            str(source),
            "--dataset",
            "/results/temperature",
            "--dataset",
            "/results/equilibrium_volume",
            "--output",
            str(destination),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    with np.load(destination, allow_pickle=False) as arrays:
        assert arrays["results__temperature"].dtype == np.float64
        assert arrays["results__equilibrium_volume"].shape == (2, 1)
    manifest = json.loads(destination.with_suffix(".json").read_text(encoding="utf-8"))
    item = manifest["datasets"]["results__equilibrium_volume"]
    assert item["path"] == "/results/equilibrium_volume"
    assert item["unit"] == "angstrom^3"
    assert item["shape"] == [2, 1]
    assert manifest["groups"]["/metadata"]["attributes"]["module"] == "qha"


def test_extract_all_writes_pickle_free_text(tmp_path: Path) -> None:
    """The documented --all mode remains portable for HDF5 text datasets."""
    source = tmp_path / "result.hdf5"
    destination = tmp_path / "complete.npz"
    _native_fixture(source)
    subprocess.run(
        [
            sys.executable,
            str(EXTRACTOR),
            str(source),
            "--all",
            "--output",
            str(destination),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    with np.load(destination, allow_pickle=False) as arrays:
        assert arrays["results__labels"].tolist() == ["state-a", "state-b"]
        assert arrays["results__valid_mask"].dtype == np.bool_
