"""Reports, neutral plots, tabular export, and public module contract."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from quantas.api import seismic as public_api
from quantas.core.geometry import Hemisphere
from quantas.core.physics.seismic import MODE_INDEX, SamplingLevel, WaveMode
from quantas.models import SphericalMapSpec
from quantas.modules.seismic import (
    MODULE_CONTRACT,
    SeismicInput,
    SeismicOptions,
    SeismicPlotOptions,
    build_seismic_plots,
    build_seismic_report,
    read_seismic_hdf5,
    run_seismic,
    write_seismic_csv,
    write_seismic_hdf5,
)


DATA = (
    Path(__file__).parents[2]
    / "physics"
    / "seismic"
    / "data"
    / "hydroxylapatite.dat"
)


def _enhancement_result():
    return run_seismic(
        DATA,
        options=SeismicOptions(
            ntheta=5,
            nphi=9,
            hemisphere=Hemisphere.UPPER,
            level=SamplingLevel.ENHANCEMENT,
            batch_size=10,
            track_polarization_axes=True,
        ),
    )


@pytest.mark.module
@pytest.mark.seismic
def test_report_contains_extrema_anisotropy_and_directional_diagnostics() -> None:
    result = _enhancement_result()
    tables = build_seismic_report(result)
    by_title = {table.title: table for table in tables}

    expected = {
        "Seismic calculation summary",
        "Selected numerical options",
        "Isotropic elastic reference properties",
        "Mechanical stability",
        "Hill-average isotropic velocity reference",
        "Acoustic-mode conventions",
        "Sampled phase-velocity extrema",
        "Derived phase-velocity diagnostics",
        "Sampled group-velocity extrema",
        "Power-flow angle statistics",
        "Enhancement and focusing diagnostics",
        "Sampling and numerical diagnostics",
        "Sampled acoustic-axis candidates",
    }
    assert expected <= set(by_title)

    payload = result.results["seismic"]
    phase = payload.field.phase
    index = MODE_INDEX[WaveMode.V_P]
    values = phase.phase_speeds[:, index]
    mask = phase.valid_mask[:, index]
    minimum = float(np.min(values[mask]))
    maximum = float(np.max(values[mask]))
    expected_anisotropy = 200.0 * (maximum - minimum) / (maximum + minimum)

    phase_table = by_title["Sampled phase-velocity extrema"]
    vp_row = next(row for row in phase_table.rows if row[0] == "V_P")
    assert float(vp_row[1]) == pytest.approx(minimum, rel=2.0e-9)
    assert float(vp_row[5]) == pytest.approx(maximum, rel=2.0e-9)
    assert float(vp_row[9]) == pytest.approx(expected_anisotropy, rel=2.0e-9)
    assert str(vp_row[2]).startswith("[")
    assert 0.0 <= float(vp_row[3]) <= 180.0
    assert 0.0 <= float(vp_row[4]) < 360.0

    derived = by_title["Derived phase-velocity diagnostics"]
    assert [row[0] for row in derived.rows] == [
        "V_S1 - V_S2",
        "Directional S-wave anisotropy",
        "V_P / V_S1",
        "V_P / V_S2",
    ]
    assert derived.rows[0][-1] == "not applicable"
    assert derived.rows[1][-1] == "not applicable"
    assert float(derived.rows[2][-1]) > 0.0

    group = by_title["Sampled group-velocity extrema"]
    assert all(str(row[5]).startswith("[") for row in group.rows)
    assert all(str(row[10]).startswith("[") for row in group.rows)

    enhancement = by_title["Enhancement and focusing diagnostics"]
    assert all(int(row[6]) > 0 for row in enhancement.rows)
    assert enhancement.metadata["enhancement_representation"] == "log10(A)"

    candidates = by_title["Sampled acoustic-axis candidates"]
    assert candidates.metadata["total_unique_candidates"] >= 1
    assert candidates.rows[0][0] == "V_S2-V_S1"
    assert candidates.rows[0][1] == "[0.000000, 0.000000, 1.000000]"
    assert candidates.rows[0][-1] == 9


@pytest.mark.module
@pytest.mark.seismic
def test_isotropic_report_does_not_assign_an_arbitrary_extremum_direction() -> None:
    bulk = 80.0
    shear = 40.0
    stiffness = np.zeros((6, 6), dtype=float)
    stiffness[:3, :3] = bulk - 2.0 * shear / 3.0
    np.fill_diagonal(stiffness[:3, :3], bulk + 4.0 * shear / 3.0)
    np.fill_diagonal(stiffness[3:, 3:], shear)
    result = run_seismic(
        SeismicInput(
            jobname="Isotropic solid",
            stiffness=stiffness,
            density=3000.0,
        ),
        SeismicOptions(
            ntheta=4,
            nphi=7,
            level=SamplingLevel.PHASE,
            track_polarization_axes=False,
        ),
    )

    phase_table = next(
        table
        for table in build_seismic_report(result)
        if table.title == "Sampled phase-velocity extrema"
    )
    for row in phase_table.rows:
        assert row[2] == "all sampled directions"
        assert row[6] == "all sampled directions"
        assert float(row[9]) == pytest.approx(0.0, abs=1.0e-10)
        assert float(row[10]) == pytest.approx(1.0, abs=1.0e-10)


@pytest.mark.module
@pytest.mark.seismic
def test_neutral_plot_specs_include_scalar_maps_and_polarizations() -> None:
    result = _enhancement_result()
    collection = build_seismic_plots(
        result,
        SeismicPlotOptions(include_polarizations=True),
    )

    assert len(collection.plots) == 16
    assert collection.warnings == []
    assert all(isinstance(spec, SphericalMapSpec) for spec in collection.plots)
    assert all(spec.values.shape == (5, 9) for spec in collection.plots)
    assert all(spec.hemisphere == "upper" for spec in collection.plots)
    assert all(
        spec.metadata["azimuth_endpoint_included"] is False for spec in collection.plots
    )

    shear = next(spec for spec in collection.plots if spec.key == "shear_anisotropy")
    assert shear.value_axis.unit == "%"
    assert len(shear.markers) == 2
    assert len(shear.axis_layers) == 1
    layer = shear.axis_layers[0]
    assert layer.key == "v_s1_polarization"
    assert layer.metadata["axial"] is True
    assert layer.axes.shape[1] == 3
    assert np.allclose(np.linalg.norm(layer.axes, axis=1), 1.0)


@pytest.mark.module
@pytest.mark.seismic
def test_plot_builder_warns_when_requested_quantities_are_unavailable() -> None:
    result = run_seismic(
        DATA,
        SeismicOptions(
            ntheta=3,
            nphi=5,
            level=SamplingLevel.PHASE,
            track_polarization_axes=False,
        ),
    )
    collection = build_seismic_plots(
        result,
        SeismicPlotOptions(
            properties=("phase_v_p", "group_v_p", "shear_anisotropy"),
            include_polarizations=True,
        ),
    )

    assert [spec.key for spec in collection.plots] == [
        "phase_v_p",
        "shear_anisotropy",
    ]
    assert len(collection.warnings) == 2
    assert "polarization" in collection.warnings[0].lower()
    assert "group_v_p" in collection.warnings[1]


@pytest.mark.module
@pytest.mark.seismic
def test_long_form_csv_export_has_named_mode_resolved_columns(tmp_path: Path) -> None:
    result = _enhancement_result()
    output = write_seismic_csv(result, tmp_path / "seismic_fields")

    assert output.suffix == ".csv"
    with output.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream))
    header = rows[4]
    assert "mode" in header
    assert "phase_speed_km_s" in header
    assert "group_speed_km_s" in header
    assert "power_flow_angle_degree" in header
    assert "log10_enhancement" in header
    assert "tracked_polarization_x" in header
    assert "enhancement" not in header

    payload = result.results["seismic"]
    data_rows = rows[5:]
    assert len(data_rows) == payload.field.n_points * 3
    assert data_rows[0][1] == "V_S2"
    assert data_rows[1][1] == "V_S1"
    assert data_rows[2][1] == "V_P"


@pytest.mark.module
@pytest.mark.seismic
def test_report_and_plot_specs_build_from_typed_hdf5_round_trip(
    tmp_path: Path,
) -> None:
    result = _enhancement_result()
    filename = write_seismic_hdf5(result, tmp_path / "seismic_result")
    restored = read_seismic_hdf5(filename)

    original_tables = build_seismic_report(result)
    restored_tables = build_seismic_report(restored)
    assert [table.title for table in restored_tables] == [
        table.title for table in original_tables
    ]
    original_phase = next(
        table
        for table in original_tables
        if table.title == "Sampled phase-velocity extrema"
    )
    restored_phase = next(
        table
        for table in restored_tables
        if table.title == "Sampled phase-velocity extrema"
    )
    assert restored_phase.rows == original_phase.rows

    plots = build_seismic_plots(restored)
    assert len(plots.plots) == 16
    assert plots.warnings == []


@pytest.mark.module
@pytest.mark.seismic
def test_public_contract_and_top_level_facade_are_complete() -> None:
    result = _enhancement_result()

    assert MODULE_CONTRACT.name == "seismic"
    assert MODULE_CONTRACT.result_key == "seismic"
    MODULE_CONTRACT.validate_result(result)
    assert MODULE_CONTRACT.build_report is build_seismic_report
    assert MODULE_CONTRACT.build_plots is build_seismic_plots

    assert public_api.run is not None
    assert public_api.build_report(result) == build_seismic_report(result)
    assert public_api.build_plots is not None
    assert public_api.build_summary is not None
    assert public_api.build_surfaces is not None
    assert public_api.write_csv is not None
