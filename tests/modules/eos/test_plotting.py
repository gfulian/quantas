from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=True)

from click.testing import CliRunner
import numpy as np

from quantas.cli.main import main
from quantas.core.math.fitting import EffectiveVarianceOptions, OLSOptions
from quantas.core.physics.eos import PVTEOS, PVTModel
from quantas.models import LinePlotSpec
from quantas.modules.eos import (
    EOSArchive,
    EOSDataset,
    EOSFitOptions,
    EOSFitRequest,
    EOSFitter,
    EOSPlotOptions,
    EOSPlotter,
    ParameterConstraint,
    read_eos_input,
)
from quantas.renderers.plots import MatplotlibOptions, render_plot_collection


DATA = Path(__file__).with_name("data")


def _quartz_archive(tmp_path: Path) -> Path:
    dataset = read_eos_input(DATA / "PV_quartz.dat")
    request = EOSFitRequest(
        model="BM3",
        domain="pv",
        target="volume",
        options=EOSFitOptions(
            solver_options=EffectiveVarianceOptions(max_iterations=50)
        ),
    )
    result = EOSFitter().fit(dataset, request)
    path = tmp_path / "quartz.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        archive.store_fit(1, request, result)
    return path


def _pvt_archive(tmp_path: Path) -> Path:
    model = PVTModel("BM3", "linear", "berman:quadratic")
    pressure_parameters = {"K0": 160.0, "KP": 4.2, "KPP": -0.02, "V0": 100.0}
    thermal_parameters = {
        "V0": 100.0,
        "temperature_ref": 300.0,
        "alpha0": 3.0e-5,
        "alpha1": 1.0e-8,
    }
    coupling_parameters = {"dK0_dT": -0.02}
    temperature = np.repeat([300.0, 600.0, 900.0], 6)
    pressure = np.tile(np.linspace(0.0, 10.0, 6), 3)
    volume = PVTEOS().volume(
        model,
        pressure_parameters,
        thermal_parameters,
        coupling_parameters,
        pressure,
        temperature,
    )
    dataset = EOSDataset(
        jobname="Synthetic PVT",
        columns={
            "pressure": pressure,
            "temperature": temperature,
            "volume": volume,
            "sigma_pressure": np.full(pressure.size, 0.01),
            "sigma_temperature": np.full(pressure.size, 1.0),
            "sigma_volume": np.full(pressure.size, 0.005),
        },
        units={
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3",
            "sigma_pressure": "GPa",
            "sigma_temperature": "K",
            "sigma_volume": "angstrom^3",
        },
    )
    constraints = (
        ParameterConstraint.free("K0", 158.0),
        ParameterConstraint.free("KP", 4.0),
        ParameterConstraint.free("V0", 99.8),
        ParameterConstraint.fixed("temperature_ref", 300.0),
        ParameterConstraint.free("alpha0", 2.8e-5),
        ParameterConstraint.free("alpha1", 0.8e-8),
        ParameterConstraint.free("dK0_dT", -0.018),
    )
    request = EOSFitRequest(
        model=model,
        domain="pvt",
        constraints=constraints,
        options=EOSFitOptions(solver_options=OLSOptions(max_iterations=5000)),
    )
    result = EOSFitter().fit(dataset, request)
    path = tmp_path / "pvt.hdf5"
    with EOSArchive.create(path, dataset=dataset) as archive:
        archive.store_fit(1, request, result)
    return path


def test_pv_plotter_builds_fit_residual_and_normalized_specs(tmp_path: Path) -> None:
    plotter = EOSPlotter.from_archive(_quartz_archive(tmp_path))

    collection = plotter.build(options=EOSPlotOptions(point_size=4.0, curve_width=2.2))

    assert plotter.available_plot_types() == (
        "fit",
        "residuals",
        "standardized-residuals",
        "normalized-pressure",
    )
    assert all(isinstance(spec, LinePlotSpec) for spec in collection.plots)
    by_key = {spec.key: spec for spec in collection.plots}
    assert by_key["fit"].x_axis.label == "Volume (Å³)"
    assert by_key["fit"].y_axis.label == "Pressure (GPa)"
    assert by_key["fit"].series[0].style.line_width == 2.2
    data = next(
        series for series in by_key["fit"].series if series.key.startswith("data")
    )
    assert data.style.marker_size == 4.0
    assert data.x_error is not None
    assert data.y_error is not None
    assert "$f_E$" in by_key["normalized_pressure"].x_axis.label
    assert "$F_E$" in by_key["normalized_pressure"].y_axis.label


def test_pvt_plotter_builds_coverage_isotherms_isobars_and_residuals(
    tmp_path: Path,
) -> None:
    plotter = EOSPlotter.from_archive(_pvt_archive(tmp_path))

    collection = plotter.build(
        options=EOSPlotOptions(isotherms=(300.0, 900.0), isobars=(0.0, 10.0))
    )

    keys = {spec.key for spec in collection.plots}
    assert {
        "coverage",
        "isotherms",
        "isobars",
        "residuals_vs_pressure",
        "residuals_vs_temperature",
    } <= keys
    isotherms = next(spec for spec in collection.plots if spec.key == "isotherms")
    assert isotherms.x_axis.label == "Volume (Å³)"
    assert isotherms.y_axis.label == "Pressure (GPa)"
    assert [series.label for series in isotherms.series[:2]] == ["300 K", "900 K"]
    isobars = next(spec for spec in collection.plots if spec.key == "isobars")
    assert [series.label for series in isobars.series[:2]] == ["0 GPa", "10 GPa"]


def test_matplotlib_renderer_draws_eos_errorbars_and_writes_files(
    tmp_path: Path,
) -> None:
    collection = EOSPlotter.from_archive(_quartz_archive(tmp_path)).build(
        ("fit", "normalized-pressure"),
        options=EOSPlotOptions(show_uncertainties=True),
    )
    output = tmp_path / "plots"

    rendered = render_plot_collection(
        collection,
        MatplotlibOptions(output_dir=output, dpi=90, close=True),
    )

    assert len(rendered.artifacts) == 2
    assert all(
        artifact.path is not None and artifact.path.is_file()
        for artifact in rendered.artifacts
    )


def test_eos_plot_title_figure_size_and_tight_layout_are_configurable(
    tmp_path: Path,
) -> None:
    collection = EOSPlotter.from_archive(_quartz_archive(tmp_path)).build(
        ("fit",),
        options=EOSPlotOptions(show_title=False),
    )

    rendered = render_plot_collection(
        collection,
        MatplotlibOptions(
            figure_size=(7.2, 5.1),
            tight_layout=True,
            close=False,
        ),
    )

    artifact = rendered.artifacts[0]
    assert artifact.spec.title == ""
    assert artifact.figure.axes[0].get_title() == ""
    assert np.allclose(artifact.figure.get_size_inches(), [7.2, 5.1])


def test_normalized_pressure_hides_zero_pressure_point_by_default(
    tmp_path: Path,
) -> None:
    plotter = EOSPlotter.from_archive(_quartz_archive(tmp_path))

    hidden = plotter.build(("normalized-pressure",)).plots[0]
    shown = plotter.build(
        ("normalized-pressure",),
        options=EOSPlotOptions(show_zero_pressure_point=True),
    ).plots[0]

    hidden_data = np.concatenate(
        [series.y for series in hidden.series if series.key.startswith("data")]
    )
    shown_data = np.concatenate(
        [series.y for series in shown.series if series.key.startswith("data")]
    )
    assert not np.any(np.isclose(hidden_data, 0.0, atol=1.0e-14))
    assert np.any(np.isclose(shown_data, 0.0, atol=1.0e-14))
    assert hidden.metadata["zero_pressure_points_hidden"] == 1
    assert shown.metadata["zero_pressure_points_shown"] is True


def test_matplotlib_typography_is_configurable(tmp_path: Path) -> None:
    collection = EOSPlotter.from_archive(_quartz_archive(tmp_path)).build(("fit",))

    rendered = render_plot_collection(
        collection,
        MatplotlibOptions(
            axis_label_font_size=16.0,
            legend_font_size=14.0,
            title_font_size=17.0,
            tick_label_font_size=12.0,
            tight_layout=True,
            close=False,
        ),
    )

    axis = rendered.artifacts[0].figure.axes[0]
    assert axis.xaxis.label.get_fontsize() == 16.0
    assert axis.yaxis.label.get_fontsize() == 16.0
    assert axis.title.get_fontsize() == 17.0
    assert all(label.get_fontsize() == 12.0 for label in axis.get_xticklabels())
    legend = axis.get_legend()
    assert legend is not None
    assert all(text.get_fontsize() == 14.0 for text in legend.get_texts())


def test_eos_plot_cli_writes_selected_figures_and_documents_options(
    tmp_path: Path,
) -> None:
    archive = _quartz_archive(tmp_path)
    output = tmp_path / "cli_plots"

    result = CliRunner().invoke(
        main,
        [
            "eos",
            "plot",
            str(archive),
            "--plot",
            "fit",
            "--plot",
            "normalized-pressure",
            "--point-size",
            "4",
            "--curve-width",
            "2.1",
            "--dpi",
            "90",
            "--no-title",
            "--figure-width",
            "7.2",
            "--figure-height",
            "5.1",
            "--axis-label-font-size",
            "16",
            "--legend-font-size",
            "14",
            "--title-font-size",
            "16",
            "--tick-label-font-size",
            "12",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "EOS plot written to" in result.output
    assert len(list(output.glob("*.png"))) == 2

    help_result = CliRunner().invoke(main, ["eos", "plot", "--help"])
    assert help_result.exit_code == 0
    for heading in (
        "Record selection:",
        "Plot selection:",
        "Data presentation:",
        "Style:",
        "P-V-T curves:",
        "Output and reporting:",
    ):
        assert heading in help_result.output
    assert "--uncertainties" in help_result.output
    assert "--curve-width" in help_result.output
    assert "--title" in help_result.output
    assert "--figure-width" in help_result.output
    assert "--figure-height" in help_result.output
    assert "--axis-label-font-size" in help_result.output
    assert "--legend-font-size" in help_result.output
    assert "--title-font-size" in help_result.output
    assert "--tick-label-font-size" in help_result.output
    assert "--zero-pressure-point" in help_result.output
    assert "--dpi" in help_result.output


def test_vt_fit_curve_does_not_extend_below_zero_kelvin(tmp_path: Path) -> None:
    """V-T plot padding must respect the physical lower bound at zero kelvin."""
    dataset = read_eos_input(DATA / "rutile.dat")
    request = EOSFitRequest(
        model="berman:quadratic",
        domain="vt",
        options=EOSFitOptions(
            solver_options=EffectiveVarianceOptions(
                max_iterations=30,
                inner_max_iterations=5000,
            )
        ),
    )
    result = EOSFitter().fit(dataset, request)
    archive_path = tmp_path / "rutile.hdf5"
    with EOSArchive.create(archive_path, dataset=dataset) as archive:
        archive.store_fit(1, request, result)

    fit = EOSPlotter.from_archive(archive_path).build(("fit",)).plots[0]
    curve = fit.series[0]

    assert np.all(curve.x >= 0.0)
    assert curve.x[0] == np.finfo(np.float64).tiny
