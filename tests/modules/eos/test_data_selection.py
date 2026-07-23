"""Non-destructive EOS data selection, groups, and crystal-system metadata."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from quantas.core.math.fitting import OLSOptions
from quantas.modules.eos import (
    EOSArchive,
    EOSCrystalSystem,
    EOSFitOptions,
    EOSFitRequest,
    EOSFitter,
    EOSSpecError,
    parse_eos_crystal_system,
    read_eos_input,
    read_eos_spec,
    resolve_eos_spec,
)
from quantas.modules.eos.report import (
    eos_data_selection_table,
    eos_data_table,
    eos_input_summary_table,
)
from quantas.renderers.tables import render_table


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text.strip() + "\n", encoding="utf-8")
    return path


def _selection_dataset(tmp_path: Path) -> Path:
    return _write(
        tmp_path,
        "selection.dat",
        """
        TITLE Grouped quartz-like data
        SYSTEM tetragonal
        FORMAT P V SIGP SIGV GROUP USE
        0.0 100.0 0.01 0.02 1 1
        1.0  99.0 0.01 0.02 1 1 *
        2.0  98.0 0.01 0.02 2 0
        3.0  97.0 0.01 0.02 2 1
        4.0  96.0 0.01 0.02 2 1
        """,
    )


def test_reader_combines_use_and_trailing_marker_non_destructively(
    tmp_path: Path,
) -> None:
    dataset = read_eos_input(_selection_dataset(tmp_path))

    assert dataset.npoints == 5
    assert dataset.selected_npoints == 3
    assert dataset.excluded_npoints == 2
    np.testing.assert_array_equal(
        dataset.selection_mask(), [True, False, False, True, True]
    )
    np.testing.assert_array_equal(dataset.groups, [1, 1, 2, 2, 2])
    assert dataset.group_ids == (1, 2)
    assert dataset.group_summary() == (
        {"group": 1, "total": 2, "selected": 1, "excluded": 1},
        {"group": 2, "total": 3, "selected": 2, "excluded": 1},
    )
    assert dataset.crystal_system is EOSCrystalSystem.TETRAGONAL
    assert dataset.metadata["independent_cell_axes"] == ("a", "c")
    assert dataset.metadata["selection"]["excluded_by_marker"] == 1


def test_dataset_default_selection_is_used_by_api_fit(tmp_path: Path) -> None:
    dataset = read_eos_input(_selection_dataset(tmp_path))
    request = EOSFitRequest(
        model="BM2",
        options=EOSFitOptions(solver_options=OLSOptions()),
    )

    result = EOSFitter().fit(dataset, request)

    assert result.fit.success
    assert result.fit.n_points == 3
    assert result.request.mask is not None
    np.testing.assert_array_equal(result.request.mask, [True, False, False, True, True])
    assert result.request.metadata["selection"]["base"] == "default"


def test_spec_can_select_groups_and_rows_from_default_or_all(tmp_path: Path) -> None:
    dataset = read_eos_input(_selection_dataset(tmp_path))
    spec = _write(
        tmp_path,
        "selection.spec",
        """
        # QUANTAS EOS SPEC 1

        [job default-group-2]
        domain = pv
        targets = volume
        model = BM2
        groups = 2
        solver = ols

        [job all-group-1]
        domain = pv
        targets = volume
        model = BM2
        selection_base = all
        include_groups = 1
        exclude_rows = 2
        solver = ols
        accept = no
        """,
    )

    resolved = resolve_eos_spec(read_eos_spec(spec), dataset)
    first, second = (job.request for job in resolved.plan.jobs)

    np.testing.assert_array_equal(first.mask, [False, False, False, True, True])
    np.testing.assert_array_equal(second.mask, [True, False, False, False, False])
    assert first.metadata["selection"]["selected"] == 2
    assert second.metadata["selection"]["base"] == "all"


def test_selection_spec_rejects_missing_unknown_or_empty_groups(
    tmp_path: Path,
) -> None:
    plain = _write(
        tmp_path,
        "plain.dat",
        """
        FORMAT P V
        0 100
        1 99
        2 98
        """,
    )
    dataset = read_eos_input(plain)
    missing = _write(
        tmp_path,
        "missing.spec",
        """
        # QUANTAS EOS SPEC 1
        [job volume]
        domain = pv
        targets = volume
        groups = 1
        """,
    )
    with pytest.raises(EOSSpecError, match="requires a GROUP column"):
        resolve_eos_spec(read_eos_spec(missing), dataset)

    grouped = read_eos_input(_selection_dataset(tmp_path))
    unknown = _write(
        tmp_path,
        "unknown.spec",
        """
        # QUANTAS EOS SPEC 1
        [job volume]
        domain = pv
        targets = volume
        groups = 9
        """,
    )
    with pytest.raises(EOSSpecError, match="unknown data group"):
        resolve_eos_spec(read_eos_spec(unknown), grouped)

    empty = _write(
        tmp_path,
        "empty.spec",
        """
        # QUANTAS EOS SPEC 1
        [job volume]
        domain = pv
        targets = volume
        groups = 1
        exclude_rows = 1
        """,
    )
    with pytest.raises(EOSSpecError, match="excludes every observation"):
        resolve_eos_spec(read_eos_spec(empty), grouped)


def test_selection_and_groups_round_trip_through_hdf5(tmp_path: Path) -> None:
    dataset = read_eos_input(_selection_dataset(tmp_path))
    path = tmp_path / "selection.hdf5"

    with EOSArchive.create(path, dataset=dataset) as archive:
        restored = archive.dataset(1)

    np.testing.assert_array_equal(restored.default_mask, dataset.default_mask)
    np.testing.assert_array_equal(restored.groups, dataset.groups)
    assert restored.crystal_system is EOSCrystalSystem.TETRAGONAL
    assert restored.metadata == dataset.metadata


def test_schema_1_0_archive_without_selection_remains_readable(
    tmp_path: Path,
) -> None:
    dataset = read_eos_input(_selection_dataset(tmp_path))
    path = tmp_path / "selection.hdf5"

    with EOSArchive.create(path, dataset=dataset):
        pass
    with h5py.File(path, "r+") as h5:
        h5["metadata"].attrs["schema_version"] = "1.0"
        del h5["input/datasets/000001/selection"]

    with EOSArchive(path) as archive:
        restored = archive.dataset(1)
        assert archive.inspect().schema_version == "1.0"

    np.testing.assert_array_equal(
        restored.default_mask, np.ones(dataset.npoints, dtype=np.bool_)
    )
    assert restored.groups is None


def test_selection_reports_are_explicit_and_grouped(tmp_path: Path) -> None:
    dataset = read_eos_input(_selection_dataset(tmp_path))

    summary = render_table(eos_input_summary_table(dataset, "result.hdf5"))
    data = render_table(eos_data_table(dataset))
    selection = render_table(eos_data_selection_table(dataset))

    assert "Selected by default" in summary
    assert "tetragonal" in summary
    assert "Group" in data and "Use" in data
    assert "Group 1" in selection and "Group 2" in selection
    assert "non-destructive" in selection


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("triclinic", EOSCrystalSystem.TRICLINIC),
        ("rhombohedral", EOSCrystalSystem.TRIGONAL),
        ("CUBIC", EOSCrystalSystem.CUBIC),
    ],
)
def test_crystal_system_aliases_are_canonical(
    value: str, expected: EOSCrystalSystem
) -> None:
    assert parse_eos_crystal_system(value) is expected


def test_reader_rejects_invalid_system_use_and_group(tmp_path: Path) -> None:
    invalid_system = _write(
        tmp_path,
        "bad-system.dat",
        "SYSTEM quasicrystal\nFORMAT P V\n0 1",
    )
    with pytest.raises(ValueError, match="Unsupported EOS crystal system"):
        read_eos_input(invalid_system)

    invalid_use = _write(
        tmp_path,
        "bad-use.dat",
        "FORMAT P V USE\n0 1 maybe",
    )
    with pytest.raises(ValueError, match="USE value"):
        read_eos_input(invalid_use)

    invalid_group = _write(
        tmp_path,
        "bad-group.dat",
        "FORMAT P V GROUP\n0 1 1.5",
    )
    with pytest.raises(ValueError, match="GROUP value"):
        read_eos_input(invalid_group)
