"""Characterization tests for the fixed Quantas numerical precision policy."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from quantas.core.numerics import (
    FLOAT_DTYPE,
    NumericPrecisionPolicy,
    cast_floating_array,
)
from quantas.io.hdf5 import (
    read_precision_metadata,
    write_array_dataset,
    write_numeric_attribute,
    write_precision_metadata,
    write_result_metadata,
)
from quantas.models import ReportTable, ResultData, ResultMetadata
from quantas.renderers.tables import (
    format_numeric,
    render_table,
    resolve_numeric_format,
)


def test_precision_policy_is_fixed_to_double_precision() -> None:
    policy = NumericPrecisionPolicy()
    assert policy.working == "double"
    assert policy.storage == "double"
    assert policy.working_dtype == "float64"
    assert policy.storage_bits == 64
    assert FLOAT_DTYPE == np.dtype("float64")


def test_precision_policy_rejects_runtime_precision_arguments() -> None:
    with pytest.raises(TypeError):
        NumericPrecisionPolicy(storage="single")


def test_floating_cast_uses_float64_and_preserves_integer_arrays() -> None:
    values = np.array([1.0, 2.0], dtype=np.float32)
    integers = np.array([1, 2], dtype=np.int64)
    assert cast_floating_array(values).dtype == np.dtype("float64")
    assert cast_floating_array(integers).dtype == np.dtype("int64")


def test_hdf5_writers_store_double_precision_floats(tmp_path) -> None:
    filename = tmp_path / "precision.hdf5"
    with h5py.File(filename, "w") as h5:
        group = h5.create_group("data")
        dataset = write_array_dataset(group, "values", np.array([1.0, 2.0]))
        write_numeric_attribute(group, "threshold", 1.0e-8)
        assert dataset.dtype == np.dtype("float64")
        assert np.asarray(group.attrs["threshold"]).dtype == np.dtype("float64")


def test_precision_metadata_is_fixed_and_round_trips(tmp_path) -> None:
    filename = tmp_path / "metadata.hdf5"
    result = ResultData(metadata=ResultMetadata(module="test", method="precision"))
    with h5py.File(filename, "w") as h5:
        write_result_metadata(h5, result)
        write_precision_metadata(h5)
    with h5py.File(filename, "r") as h5:
        restored = read_precision_metadata(h5)
        assert h5["metadata/numerics"].attrs["storage_dtype"] == "float64"
    assert restored == NumericPrecisionPolicy().as_metadata()


def test_display_formats_do_not_change_storage_precision() -> None:
    value = np.float64(1.2345679)
    assert resolve_numeric_format("energy_ha") == ".12E"
    assert format_numeric(value, "energy_ha") == "1.234567900000E+00"
    assert format_numeric(298.15, "temperature") == "298.15"
    assert format_numeric(123.456789, ".4f") == "123.4568"


def test_plain_text_tables_remain_free_of_ansi_and_box_drawing() -> None:
    table = ReportTable(
        title="Energy",
        columns=["Quantity", "Value"],
        rows=[["E", -123.456789]],
        metadata={"column_formats": [None, "energy_ha"]},
    )
    text = render_table(table)
    assert "-1.234567890000E+02" in text
    assert "\x1b[" not in text
    assert not {"┏", "┓", "┗", "┛", "━", "┃"}.intersection(text)
