"""Tests for QHA pressure-volume inspection utilities."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.modules.qha.api import inspect_qha_input
from quantas.modules.qha.inspect import PressureVolumePreview, pressure_volume_preview
from quantas.modules.qha.models import QHAInput, QHAOptions


def make_static_input() -> QHAInput:
    """Return a small energy-volume dataset with a clear minimum."""
    volume = np.asarray([69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0], dtype=float)
    energy = -200.0 + 0.0025 * (volume - 72.0) ** 2
    return QHAInput(jobname="preview", volume=volume, energy=energy)


def test_pressure_volume_preview_returns_polynomial_and_eos_estimates() -> None:
    """Pressure preview should contain both polynomial and EOS estimates."""
    qha_input = make_static_input()
    options = QHAOptions(
        energy_unit="eV",
        volume_unit="angstrom",
        pressure_unit="GPa",
        energy_degree=2,
        eos="BM",
    )

    preview = pressure_volume_preview(qha_input, options)

    assert isinstance(preview, PressureVolumePreview)
    assert preview.success
    assert preview.polynomial is not None
    assert preview.eos is not None
    assert preview.polynomial.success
    assert preview.eos.success
    assert preview.polynomial.pressure.shape == qha_input.volume.shape
    assert preview.eos.pressure.shape == qha_input.volume.shape
    assert preview.pressure_unit == "GPa"


def test_pressure_volume_preview_table_rows_are_neutral_records() -> None:
    """Preview rows should be frontend-neutral dictionaries."""
    qha_input = make_static_input()
    options = QHAOptions(
        energy_unit="eV", volume_unit="angstrom", pressure_unit="GPa", energy_degree=2
    )

    rows = pressure_volume_preview(qha_input, options, include_eos=False).table_rows()

    assert len(rows) == qha_input.nvol
    assert set(rows[0]) == {"volume", "energy", "pressure_polynomial", "pressure_eos"}
    assert rows[0]["pressure_eos"] is None
    assert rows[0]["pressure_polynomial"] is not None


def test_inspect_qha_input_accepts_qha_input_objects() -> None:
    """Public API should inspect normalized QHAInput objects."""
    qha_input = make_static_input()
    options = QHAOptions(
        energy_unit="eV", volume_unit="angstrom", pressure_unit="GPa", energy_degree=2
    )

    preview = inspect_qha_input(qha_input, options, include_eos=False)

    assert preview.success
    assert preview.eos is None
    assert preview.polynomial is not None


def test_pressure_volume_preview_raises_for_missing_static_data() -> None:
    """Inspection should reject inputs without static energy-volume data."""
    with pytest.raises(ValueError):
        pressure_volume_preview(QHAInput(volume=np.asarray([1.0, 2.0])))


def test_pressure_volume_preview_reports_failed_eos_without_failing_polynomial() -> (
    None
):
    """A failed EOS estimate should not invalidate a successful polynomial preview."""
    qha_input = make_static_input()
    options = QHAOptions(
        energy_unit="eV", volume_unit="angstrom", pressure_unit="GPa", energy_degree=2
    )

    preview = pressure_volume_preview(qha_input, options, eos="not-an-eos")

    assert preview.success
    assert preview.polynomial is not None
    assert preview.polynomial.success
    assert preview.eos is not None
    assert not preview.eos.success
    assert preview.warnings
    assert all(row["pressure_eos"] is None for row in preview.table_rows())
