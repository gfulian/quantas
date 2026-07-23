"""Matplotlib rendering and Click boundary for the seismic workflow."""

from __future__ import annotations

from pathlib import Path

import gc

import matplotlib
import numpy as np
import pytest
from click.testing import CliRunner

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from quantas.cli.main import main
from quantas.core.geometry import Hemisphere
from quantas.core.physics.seismic import SamplingLevel, WaveMode
from quantas.models import SphericalMapSpec, SphericalSummarySpec, SurfacePlotSpec
from quantas.modules.seismic import (
    SeismicOptions,
    SeismicPlotOptions,
    SeismicSurfaceOptions,
    build_seismic_plots,
    build_seismic_summary,
    build_seismic_surfaces,
    read_seismic_hdf5,
    run_seismic,
    write_seismic_hdf5,
)
from quantas.renderers.plots import (
    MatplotlibOptions,
    render_plot,
    render_plot_collection,
)


@pytest.fixture(autouse=True)
def _close_matplotlib_figures():
    """Release Matplotlib figures after each renderer test."""
    yield
    plt.close("all")
    gc.collect()


DATA = (
    Path(__file__).parents[2]
    / "physics"
    / "seismic"
    / "data"
    / "hydroxylapatite.dat"
)


def _result(level: SamplingLevel = SamplingLevel.ENHANCEMENT):
    return run_seismic(
        DATA,
        SeismicOptions(
            ntheta=5,
            nphi=9,
            level=level,
            batch_size=12,
            track_polarization_axes=True,
        ),
    )


def test_spherical_renderer_draws_overlays(
    tmp_path: Path,
) -> None:
    collection = build_seismic_plots(
        _result(),
        SeismicPlotOptions(
            properties=("shear_anisotropy",),
            polarization_stride=2,
            include_polarizations=True,
        ),
    )
    rendered = render_plot_collection(
        collection,
        MatplotlibOptions(
            output_dir=tmp_path,
            axis_label_mode="cartesian",
            polarization_color="purple",
            polarization_line_width=2.0,
            polarization_scale=0.08,
            close=False,
        ),
    )

    artifact = rendered.artifacts[0]
    assert artifact.kind == "spherical_map"
    assert artifact.path == (tmp_path / "seismic_2d_shear_anisotropy_polarization.png")
    assert artifact.path.exists()
    assert isinstance(artifact.spec, SphericalMapSpec)
    map_axis = artifact.figure.axes[0]
    labels = {text.get_text() for text in map_axis.texts}
    assert {"+x", "−x", "+y", "−y", "+z"} <= labels
    assert any("min" in label for label in labels)
    assert any("max" in label for label in labels)
    assert len(map_axis.collections) >= 3


def test_summary_and_three_dimensional_surface_specs_are_structured() -> None:
    result = _result()
    summary = build_seismic_summary(result)
    assert isinstance(summary, SphericalSummarySpec)
    assert len(summary.maps) == 6
    assert [item.key for item in summary.maps[:3]] == [
        "phase_v_p",
        "phase_v_s1",
        "phase_v_s2",
    ]

    surfaces = build_seismic_surfaces(
        result,
        SeismicSurfaceOptions(
            modes=(WaveMode.V_P,),
            surface_types=("phase", "slowness", "group"),
            include_polarizations=True,
            polarization_stride=2,
        ),
    )
    assert surfaces.warnings == []
    assert len(surfaces.plots) == 3
    assert all(isinstance(spec, SurfacePlotSpec) for spec in surfaces.plots)
    assert all(len(spec.layers) == 2 for spec in surfaces.plots)
    assert all(len(spec.vector_layers) == 1 for spec in surfaces.plots)
    assert all(spec.metadata["geometry"] == "unit_sphere" for spec in surfaces.plots)
    assert [spec.metadata["surface_type"] for spec in surfaces.plots] == [
        "phase",
        "slowness",
        "group",
    ]
    phase = surfaces.plots[0]
    sampled = phase.layers[0]
    assert sampled.x.shape == (5, 10)
    assert np.nanmax(np.sqrt(sampled.x**2 + sampled.y**2 + sampled.z**2)) > 0.0


def test_summary_and_surface_render_with_crystal_frame_labels(tmp_path: Path) -> None:
    result = _result(SamplingLevel.GROUP)
    summary = build_seismic_summary(result)
    summary_artifact = render_plot(
        summary,
        MatplotlibOptions(output_dir=tmp_path, close=False),
    )
    assert summary_artifact.kind == "spherical_summary"
    assert summary_artifact.path == tmp_path / "seismic_summary.png"
    assert summary_artifact.path.exists()

    surfaces = build_seismic_surfaces(
        result,
        SeismicSurfaceOptions(
            modes=(WaveMode.V_S1,),
            surface_types=("group",),
            include_polarizations=True,
            polarization_stride=2,
        ),
    )
    artifact = render_plot_collection(
        surfaces,
        MatplotlibOptions(
            output_dir=tmp_path,
            axis_label_mode="crystal",
            polarization_color="green",
            close=False,
        ),
    ).artifacts[0]
    axis = artifact.figure.axes[0]
    labels = {text.get_text() for text in axis.texts}
    assert {"$[100]$", "$[010]$", "$[001]$"} <= labels
    assert axis.axison is False
    assert artifact.figure._suptitle is None
    assert artifact.path == (
        tmp_path / "seismic_3d_group_v_s1_unit_sphere_polarization.png"
    )
    assert artifact.path.exists()


def test_full_sphere_map_renders_upper_and_lower_tensor_axes(tmp_path: Path) -> None:
    result = run_seismic(
        DATA,
        SeismicOptions(
            ntheta=7,
            nphi=10,
            hemisphere=Hemisphere.FULL,
            level=SamplingLevel.PHASE,
            track_polarization_axes=False,
        ),
    )
    collection = build_seismic_plots(
        result,
        SeismicPlotOptions(properties=("phase_v_p",)),
    )
    artifact = render_plot_collection(
        collection,
        MatplotlibOptions(output_dir=tmp_path, close=False),
    ).artifacts[0]
    map_axes = artifact.figure.axes[:2]
    labels = [{text.get_text() for text in axis.texts} for axis in map_axes]
    assert "+z" in labels[0]
    assert "−z" in labels[1]


def test_seismic_run_command_writes_report_and_hdf5_without_figures(
    tmp_path: Path,
) -> None:
    output = tmp_path / "apatite.hdf5"
    report_path = tmp_path / "apatite.log"
    response = CliRunner().invoke(
        main,
        [
            "seismic",
            "run",
            str(DATA),
            "--ntheta",
            "4",
            "--nphi",
            "7",
            "--level",
            "phase",
            "--batch-size",
            "5",
            "--rotate-xyz",
            "0",
            "0",
            "30",
            "--quiet",
            "--output",
            str(output),
            "--report",
            str(report_path),
            "--verbosity",
            "extended",
        ],
    )

    assert response.exit_code == 0, response.output
    assert output.exists()
    assert report_path.exists()
    assert list(tmp_path.glob("*.png")) == []
    restored = read_seismic_hdf5(output)
    assert restored.options["ntheta"] == 4
    assert restored.options["nphi"] == 7
    assert restored.options["batch_size"] == 5
    assert restored.options["level"] == "phase"
    assert restored.options["rotation"]["kind"] == "xyz"
    np.testing.assert_allclose(
        restored.options["rotation"]["matrix"],
        restored.results["seismic"].metadata["tensor_frame"]["component_transform"],
    )
    report = report_path.read_text(encoding="utf-8")
    assert "Stiffness matrix before rotation" in report
    assert "Tensor component transformation" in report
    assert "Stiffness matrix after rotation" in report
    assert "Sampled phase-velocity extrema" in report


def test_plot_command_writes_requested_outputs(
    tmp_path: Path,
) -> None:
    hdf5 = write_seismic_hdf5(_result(), tmp_path / "result")
    base = tmp_path / "figure"
    response = CliRunner().invoke(
        main,
        [
            "seismic",
            "plot",
            str(hdf5),
            "--2d",
            "--property",
            "phase_v_p",
            "--3d",
            "--surface",
            "phase",
            "--mode",
            "v_p",
            "--summary",
            "--axis-labels",
            "crystal",
            "--polarization-color",
            "purple",
            "--polarization-linewidth",
            "2.0",
            "--output",
            str(base),
        ],
    )

    assert response.exit_code == 0, response.output
    assert "Generating figure: figure_seismic_2d_phase_v_p.png" in response.output
    assert (
        "Generating figure: figure_seismic_3d_phase_v_p_unit_sphere.png"
        in response.output
    )
    expected = {
        tmp_path / "figure_seismic_2d_phase_v_p.png",
        tmp_path / "figure_seismic_2d_phase_v_p_polarization.png",
        tmp_path / "figure_seismic_3d_phase_v_p_unit_sphere.png",
        tmp_path / "figure_seismic_3d_phase_v_p_unit_sphere_polarization.png",
        tmp_path / "figure_seismic_summary.png",
        tmp_path / "figure_seismic_summary_polarization.png",
    }
    assert expected <= set(tmp_path.glob("*.png"))
    assert all(str(path) in response.output for path in expected)


def test_seismic_plot_defaults_to_summary_and_export_writes_csv(tmp_path: Path) -> None:
    hdf5 = write_seismic_hdf5(_result(SamplingLevel.PHASE), tmp_path / "result")
    plot_response = CliRunner().invoke(
        main,
        ["seismic", "plot", str(hdf5), "--output", str(tmp_path / "default")],
    )
    assert plot_response.exit_code == 0, plot_response.output
    assert (tmp_path / "default_seismic_summary.png").exists()

    csv_file = tmp_path / "fields.csv"
    export_response = CliRunner().invoke(
        main,
        ["seismic", "export", str(hdf5), "--output", str(csv_file)],
    )
    assert export_response.exit_code == 0, export_response.output
    assert csv_file.exists()
    assert "phase_speed_km_s" in csv_file.read_text(encoding="utf-8")


def test_root_help_exposes_seismic_run_plot_and_export() -> None:
    root = CliRunner().invoke(main, ["--help"])
    assert root.exit_code == 0
    assert "seismic" in root.output

    group = CliRunner().invoke(main, ["seismic", "--help"])
    assert group.exit_code == 0
    assert "run" in group.output
    assert "plot" in group.output
    assert "export" in group.output

    plot_help = CliRunner().invoke(main, ["seismic", "plot", "--help"])
    assert plot_help.exit_code == 0
    for option in (
        "--2d",
        "--3d",
        "--summary",
        "--projection",
        "--cmap",
        "--axis-labels",
        "--polarization-color",
        "--polarization-linewidth",
        "--geometry",
    ):
        assert option in plot_help.output


def test_properties_use_labels_and_colormaps() -> None:
    collection = build_seismic_plots(_result())
    by_key = {spec.key: spec for spec in collection.plots}

    assert by_key["phase_v_p"].value_axis.label == (r"$V_P$ ($\mathrm{km\,s^{-1}}$)")
    assert by_key["phase_v_s1"].value_axis.label == (
        r"$V_{S1}$ ($\mathrm{km\,s^{-1}}$)"
    )
    assert by_key["phase_v_s2"].value_axis.label == (
        r"$V_{S2}$ ($\mathrm{km\,s^{-1}}$)"
    )
    assert by_key["power_flow_v_p"].value_axis.label == r"$\psi_P$ ($^\circ$)"
    assert by_key["log10_enhancement_v_p"].value_axis.label == (r"$\log_{10}(A_P)$")
    assert by_key["phase_v_p"].colormap == "seismic"
    assert by_key["power_flow_v_p"].colormap == "quantas_powerflow"
    assert by_key["log10_enhancement_v_p"].colormap == "quantas_enhancement"


def test_3d_default_builds_all_unit_sphere_properties() -> (
    None
):
    surfaces = build_seismic_surfaces(_result())

    assert len(surfaces.plots) == 16
    assert all(isinstance(spec, SurfacePlotSpec) for spec in surfaces.plots)
    assert all(spec.metadata["geometry"] == "unit_sphere" for spec in surfaces.plots)
    for spec in surfaces.plots:
        sampled = spec.layers[0]
        radius = np.sqrt(sampled.x**2 + sampled.y**2 + sampled.z**2)
        assert np.allclose(radius[np.isfinite(radius)], 1.0, atol=1.0e-12)


def test_physical_and_unit_sphere_geometry_are_selectable() -> None:
    result = _result(SamplingLevel.GROUP)
    unit = build_seismic_surfaces(
        result,
        SeismicSurfaceOptions(
            modes=(WaveMode.V_P,),
            surface_types=("phase", "slowness", "group"),
            geometry="unit_sphere",
        ),
    )
    physical = build_seismic_surfaces(
        result,
        SeismicSurfaceOptions(
            modes=(WaveMode.V_P,),
            surface_types=("phase", "slowness", "group"),
            geometry="physical",
        ),
    )

    for spec in unit.plots:
        sampled = spec.layers[0]
        radius = np.sqrt(sampled.x**2 + sampled.y**2 + sampled.z**2)
        assert np.allclose(radius[np.isfinite(radius)], 1.0, atol=1.0e-12)
    physical_ranges = []
    for spec in physical.plots:
        sampled = spec.layers[0]
        radius = np.sqrt(sampled.x**2 + sampled.y**2 + sampled.z**2)
        physical_ranges.append((float(np.nanmin(radius)), float(np.nanmax(radius))))
    assert all(
        not np.isclose(low, 1.0) or not np.isclose(high, 1.0)
        for low, high in physical_ranges
    )


def test_polarization_variants_have_distinct_neutral_filenames() -> None:
    result = _result()
    plain = build_seismic_plots(
        result,
        SeismicPlotOptions(
            properties=("phase_v_p",),
            include_polarizations=False,
        ),
    ).plots[0]
    polarized = build_seismic_plots(
        result,
        SeismicPlotOptions(
            properties=("phase_v_p",),
            include_polarizations=True,
            polarization_stride=2,
        ),
    ).plots[0]

    assert plain.filename_stem == "seismic_2d_phase_v_p"
    assert plain.axis_layers == []
    assert polarized.filename_stem == "seismic_2d_phase_v_p_polarization"
    assert len(polarized.axis_layers) == 1

    plain_surface = build_seismic_surfaces(
        result,
        SeismicSurfaceOptions(
            properties=("phase_v_p",),
            include_polarizations=False,
        ),
    ).plots[0]
    polarized_surface = build_seismic_surfaces(
        result,
        SeismicSurfaceOptions(
            properties=("phase_v_p",),
            include_polarizations=True,
            polarization_stride=2,
        ),
    ).plots[0]
    assert plain_surface.filename_stem == "seismic_3d_phase_v_p_unit_sphere"
    assert plain_surface.vector_layers == []
    assert polarized_surface.filename_stem == (
        "seismic_3d_phase_v_p_unit_sphere_polarization"
    )
    assert len(polarized_surface.vector_layers) == 1


def test_derived_ratio_has_no_mode_polarization_or_physical_radius() -> None:
    result = _result()
    key = "phase_v_p_over_v_s1"

    map_spec = build_seismic_plots(
        result,
        SeismicPlotOptions(
            properties=(key,),
            include_polarizations=True,
            polarization_stride=2,
        ),
    ).plots[0]
    assert map_spec.axis_layers == []
    assert map_spec.filename_stem == f"seismic_2d_{key}"

    surface = build_seismic_surfaces(
        result,
        SeismicSurfaceOptions(
            properties=(key,),
            geometry="physical",
            include_polarizations=True,
            polarization_stride=2,
        ),
    ).plots[0]
    sampled = surface.layers[0]
    radius = np.sqrt(sampled.x**2 + sampled.y**2 + sampled.z**2)
    assert np.allclose(radius[np.isfinite(radius)], 1.0, atol=1.0e-12)
    assert surface.vector_layers == []
    assert surface.metadata["geometry"] == "unit_sphere"
    assert surface.metadata["mode"] is None


def test_3d_renderer_uses_full_surface_mesh(tmp_path: Path) -> None:
    result = run_seismic(
        DATA,
        SeismicOptions(
            ntheta=61,
            nphi=101,
            level=SamplingLevel.PHASE,
            track_polarization_axes=False,
        ),
    )
    spec = build_seismic_surfaces(
        result,
        SeismicSurfaceOptions(
            properties=("phase_v_p",),
            geometry="physical",
            include_polarizations=False,
        ),
    ).plots[0]
    artifact = render_plot(
        spec,
        MatplotlibOptions(output_dir=tmp_path, close=False),
    )

    surface_axis = artifact.figure.axes[0]
    first_surface = surface_axis.collections[0]
    assert first_surface._facecolor3d.shape[0] >= 6000
    xlim = surface_axis.get_xlim3d()
    assert xlim[1] > abs(xlim[0])
    assert surface_axis.axison is False


def test_cli_warns_when_physical_3d_polarizations_are_requested(
    tmp_path: Path,
) -> None:
    hdf5 = write_seismic_hdf5(_result(SamplingLevel.PHASE), tmp_path / "physical")
    response = CliRunner().invoke(
        main,
        [
            "seismic",
            "plot",
            str(hdf5),
            "--3d",
            "--property",
            "phase_v_p",
            "--geometry",
            "physical",
            "--polarization",
            "--dpi",
            "40",
            "--output",
            str(tmp_path / "physical_plot"),
        ],
    )

    assert response.exit_code == 0, response.output
    assert (
        "Generating figure: physical_plot_seismic_3d_phase_v_p_physical.png"
        in response.output
    )
    assert (
        "3D polarization overlays are supported only on unit-sphere" in response.output
    )
    assert (tmp_path / "physical_plot_seismic_3d_phase_v_p_physical.png").exists()
    assert not (
        tmp_path / "physical_plot_seismic_3d_phase_v_p_physical_polarization.png"
    ).exists()


def test_physical_3d_surfaces_disable_polarization_overlays() -> None:
    result = _result(SamplingLevel.GROUP)
    surfaces = build_seismic_surfaces(
        result,
        SeismicSurfaceOptions(
            modes=(WaveMode.V_P,),
            surface_types=("phase", "group"),
            geometry="physical",
            include_polarizations=True,
            polarization_stride=2,
        ),
    )

    assert len(surfaces.warnings) == 1
    assert "unit-sphere" in surfaces.warnings[0]
    assert "Physical-geometry" in surfaces.warnings[0]
    assert len(surfaces.plots) == 2
    assert all(isinstance(spec, SurfacePlotSpec) for spec in surfaces.plots)
    assert all(spec.metadata["geometry"] == "physical" for spec in surfaces.plots)
    assert all(spec.vector_layers == [] for spec in surfaces.plots)
    assert all("polarization" not in spec.filename_stem for spec in surfaces.plots)


def test_3d_unit_sphere_polarizations_are_clipped_to_visible_side(
    tmp_path: Path,
) -> None:
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    result = _result(SamplingLevel.GROUP)
    spec = build_seismic_surfaces(
        result,
        SeismicSurfaceOptions(
            modes=(WaveMode.V_P,),
            surface_types=("phase",),
            geometry="unit_sphere",
            include_polarizations=True,
            polarization_stride=1,
        ),
    ).plots[0]
    assert isinstance(spec, SurfacePlotSpec)
    original_count = spec.vector_layers[0].origins.shape[0]

    artifact = render_plot(
        spec,
        MatplotlibOptions(output_dir=tmp_path, close=False),
    )
    axis = artifact.figure.axes[0]
    line_collections = [
        item for item in axis.collections if isinstance(item, Line3DCollection)
    ]
    assert line_collections
    rendered_count = sum(len(item._segments3d) for item in line_collections)
    assert 0 < rendered_count < original_count


def test_seismic_figures_do_not_render_titles(tmp_path: Path) -> None:
    result = _result()
    map_spec = build_seismic_plots(
        result, SeismicPlotOptions(properties=("phase_v_p",))
    ).plots[0]
    map_artifact = render_plot(
        map_spec, MatplotlibOptions(output_dir=tmp_path, close=False)
    )
    assert map_artifact.figure._suptitle is None
    assert map_artifact.figure.axes[0].get_title() == ""

    summary_artifact = render_plot(
        build_seismic_summary(result),
        MatplotlibOptions(output_dir=tmp_path, close=False),
    )
    assert summary_artifact.figure._suptitle is None
    assert all(axis.get_title() == "" for axis in summary_artifact.figure.axes)


def test_seismic_plot_3d_without_property_selection_renders_all_fields(
    tmp_path: Path,
) -> None:
    hdf5 = write_seismic_hdf5(_result(), tmp_path / "all_fields")
    response = CliRunner().invoke(
        main,
        [
            "seismic",
            "plot",
            str(hdf5),
            "--3d",
            "--geometry",
            "unit-sphere",
            "--no-polarization",
            "--dpi",
            "40",
            "--output",
            str(tmp_path / "all3d"),
        ],
    )

    assert response.exit_code == 0, response.output
    figures = set(tmp_path.glob("all3d_seismic_3d_*_unit_sphere.png"))
    assert len(figures) == 16
    assert tmp_path / "all3d_seismic_3d_phase_v_p_unit_sphere.png" in figures
    assert (
        tmp_path / "all3d_seismic_3d_log10_enhancement_v_s2_unit_sphere.png" in figures
    )
