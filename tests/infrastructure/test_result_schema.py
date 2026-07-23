"""Tests for the shared Quantas 2.0 result schema and event log."""

from __future__ import annotations

from datetime import datetime, timezone

import h5py

from quantas.core.events import EventRecord
from quantas.io.hdf5 import (
    read_events,
    read_input_data,
    read_options,
    read_result_metadata,
    read_warnings,
    write_diagnostics,
    write_events,
    write_input_data,
    write_options,
    write_result_metadata,
)
from quantas.models import (
    InputData,
    RESULT_SCHEMA_VERSION,
    ResultData,
    ResultMetadata,
    validate_result_schema,
)


def test_native_result_schema_roundtrip(tmp_path) -> None:
    """Generic metadata, input, options, warnings and events round-trip."""
    timestamp = datetime(2026, 7, 1, 12, 0, tzinfo=timezone.utc)
    result = ResultData(
        metadata=ResultMetadata(
            module="test",
            method="schema",
            created_at=timestamp,
            created_by="pytest",
        ),
        input_data=InputData(
            source="input.yaml",
            raw="job: schema",
            data={"jobname": "Schema test", "count": 2},
        ),
        options={"debug": False, "mode": "test"},
        warnings=["example warning"],
        events=[
            EventRecord(
                message="Half complete",
                level="progress",
                progress=0.5,
                data={"current": 1, "total": 2},
                timestamp=timestamp,
            )
        ],
    )
    filename = tmp_path / "schema.hdf5"

    with h5py.File(filename, "w") as h5:
        write_result_metadata(h5, result)
        write_input_data(h5, result.input_data)
        write_options(h5, result.options)
        h5.create_group("results")
        write_diagnostics(h5, result, report_text="schema report")
        write_events(h5, result.events)
        validate_result_schema(h5, expected_module="test")

    with h5py.File(filename, "r") as h5:
        metadata = read_result_metadata(h5)
        input_data = read_input_data(h5)
        options = read_options(h5)
        warnings = read_warnings(h5)
        events = read_events(h5)

        assert h5["diagnostics/report_text"].asstr()[()] == "schema report"

    assert metadata.schema_version == RESULT_SCHEMA_VERSION
    assert metadata.created_at == timestamp
    assert metadata.created_by == "pytest"
    assert input_data is not None
    assert input_data.source == "input.yaml"
    assert input_data.raw == "job: schema"
    assert input_data.data == {"jobname": "Schema test", "count": 2}
    assert options == {"debug": False, "mode": "test"}
    assert warnings == ["example warning"]
    assert events == result.events


def test_schema_version_is_independent_from_package_version() -> None:
    """The persisted result schema has its own explicit version number."""
    metadata = ResultMetadata(module="test", method="schema")

    assert metadata.schema_version == "2.0"
    assert metadata.version != metadata.schema_version
