"""Tests for the shared native HDF5 helper layer."""

from __future__ import annotations

from datetime import datetime, timezone

import h5py
import numpy as np

from quantas.io import hdf5
from quantas.io.hdf5.read import read_mapping as read_mapping_from_read_module
from quantas.io.hdf5.write import write_mapping as write_mapping_from_write_module
from quantas.models import ResultData, ResultMetadata


def test_hdf5_package_reexports_mapping_helpers(tmp_path) -> None:
    """The split package keeps the historical import surface available."""
    filename = tmp_path / "mapping.hdf5"
    payload = {
        "name": "example",
        "array": np.arange(3.0),
        "nested": {"flag": True, "none": None},
        "items": ["a", "b"],
    }

    with h5py.File(filename, "w") as handle:
        group = handle.create_group("payload")
        hdf5.write_mapping(group, payload)

    with h5py.File(filename, "r") as handle:
        decoded = hdf5.read_mapping(handle["payload"])

    assert decoded["name"] == "example"
    np.testing.assert_allclose(decoded["array"], np.arange(3.0))
    assert decoded["nested"]["flag"] is True
    assert decoded["nested"]["none"] is None
    assert decoded["items"].tolist() == ["a", "b"]


def test_hdf5_read_write_modules_share_recursive_mapping(tmp_path) -> None:
    """The explicit read/write modules operate on the same neutral schema."""
    filename = tmp_path / "mapping_explicit.hdf5"
    payload = {"temperatures": [300.0, 600.0], "label": "QHA"}

    with h5py.File(filename, "w") as handle:
        group = handle.create_group("payload")
        write_mapping_from_write_module(group, payload)

    with h5py.File(filename, "r") as handle:
        decoded = read_mapping_from_read_module(handle["payload"])

    np.testing.assert_allclose(decoded["temperatures"], [300.0, 600.0])
    assert decoded["label"] == "QHA"


def test_hdf5_envelope_round_trips_metadata_without_payload_knowledge(tmp_path) -> None:
    """Generic metadata are handled without module-specific payload logic."""
    filename = tmp_path / "metadata.hdf5"
    created = datetime(2026, 1, 2, 3, 4, 5, tzinfo=timezone.utc)
    result = ResultData(
        metadata=ResultMetadata(
            program="quantas",
            module="dummy",
            method="audit",
            version="2.0.0",
            created_by="pytest",
            created_at=created,
            schema_version="2.0",
        )
    )

    with h5py.File(filename, "w") as handle:
        hdf5.write_result_metadata(handle, result)

    with h5py.File(filename, "r") as handle:
        metadata = hdf5.read_result_metadata(handle)

    assert metadata.program == "quantas"
    assert metadata.module == "dummy"
    assert metadata.method == "audit"
    assert metadata.created_by == "pytest"
    assert metadata.created_at == created
