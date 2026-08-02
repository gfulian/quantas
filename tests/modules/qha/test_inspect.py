"""Tests for QHA pressure-volume inspection utilities."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.api import qha as public_qha
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


def test_public_inspection_builders_report_and_plot_selected_fits() -> None:
    """Public builders should expose tables and sampled E(V) curves."""
    options = public_qha.Options(
        energy_unit="eV",
        volume_unit="angstrom",
        pressure_unit="GPa",
        energy_degree=2,
        eos="BM3",
    )
    preview = public_qha.inspect(make_static_input(), options=options)

    tables = public_qha.build_inspection_report(preview)
    plots = public_qha.build_inspection_plots(preview, sample_points=51)

    assert [table.title for table in tables] == [
        "Pressure-volume preview",
        "Pressure-volume fit diagnostics",
        "Pressure-volume fit parameters",
    ]
    assert len(plots.plots) == 1
    plot = plots.plots[0]
    assert plot.key == "qha_input_energy_volume"
    assert [series.key for series in plot.series] == [
        "observed",
        "polynomial_fit",
        "eos_fit",
    ]
    assert np.array_equal(plot.series[0].x, preview.volume)
    assert np.array_equal(plot.series[0].y, preview.energy)
    assert plot.series[1].x.size == 51
    assert plot.series[2].x.size == 51
    assert plot.metadata["volume_interval"] == (69.0, 75.0)
    assert plot.metadata["extrapolated"] is False


def test_inspection_plot_omits_failed_fit_and_preserves_warning() -> None:
    """An unusable EOS should not prevent plotting the data and polynomial."""
    preview = public_qha.inspect(
        make_static_input(),
        options=public_qha.Options(
            energy_unit="eV",
            volume_unit="angstrom",
            pressure_unit="GPa",
            energy_degree=2,
        ),
        eos="not-an-eos",
    )

    plots = public_qha.build_inspection_plots(preview)

    assert [series.key for series in plots.plots[0].series] == [
        "observed",
        "polynomial_fit",
    ]
    assert plots.warnings


def test_inspection_builders_validate_public_inputs() -> None:
    """Inspection presentation builders should reject invalid contracts."""
    with pytest.raises(TypeError, match="Preview"):
        public_qha.build_inspection_report(object())  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="Preview"):
        public_qha.build_inspection_plots(object())  # type: ignore[arg-type]

    preview = public_qha.inspect(
        make_static_input(),
        options=public_qha.Options(energy_degree=2),
        include_eos=False,
    )
    with pytest.raises(ValueError, match="at least 2"):
        public_qha.build_inspection_plots(preview, sample_points=1)
