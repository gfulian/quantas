"""Tests for neutral thermoelastic plot foundations and Matplotlib rendering."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=True)

import numpy as np

from quantas.core.physics.elasticity import (
    cold_finite_strain_component,
    reconstruct_stiffness_from_components,
)
from quantas.core.physics.eos import EnergyEOS
from quantas.core.physics.units import energy_to_pressure
from quantas.models import (
    ColoredPathSeries,
    ContourPlotSpec,
    PanelPlotSpec,
    PlotBand,
    PlotCollection,
    PlotMask,
    PlotSpan,
    ResultData,
    ResultMetadata,
    SecondaryAxis,
)
from quantas.models.structures import CrystalStructure, SymmetryMetadata
from quantas.modules.qha.models import QHAResult
from quantas.modules.thermoelasticity import (
    ElasticVolumePoint,
    ElasticVolumeSeries,
    ThermoelasticContext,
    ThermoelasticDepthProfile,
    ThermoelasticInput,
    ThermoelasticOptions,
    build_thermoelastic_domain_plot,
    build_thermoelastic_fit_plots,
    build_thermoelastic_profile_plots,
    build_thermoelastic_pt_plots,
    analyze_thermoelastic_result,
    run_thermoelastic_context,
)
from quantas.modules.thermoelasticity.plot import (
    ThermoelasticProfilePlotOptions,
    component_style,
    resolve_components,
)
from quantas.renderers.plots import MatplotlibOptions, render_plot_collection


def _plot_result() -> ResultData:
    volumes = np.linspace(90.0, 110.0, 9)
    v0 = 100.0
    k0_gpa = 120.0
    kp = 4.4
    pressure_factor = float(energy_to_pressure(1.0, "Ha", "A", "GPa"))
    energies = EnergyEOS().evaluate(
        "BM3",
        volumes,
        [-200.0, k0_gpa / pressure_factor, kp, v0],
    )
    component_parameters = {
        "C11": (220.0, 6.0),
        "C12": (85.0, 3.5),
        "C44": (72.0, 1.8),
    }
    points: list[ElasticVolumePoint] = []
    mass = 4.0e-25
    for index, (volume, energy) in enumerate(zip(volumes, energies, strict=True)):
        components = {
            label: float(
                cold_finite_strain_component(
                    volume,
                    reference_volume=v0,
                    bulk_modulus=k0_gpa,
                    bulk_modulus_derivative=kp,
                    reference_component=c0,
                    component_pressure_derivative=cp,
                    wallace_delta={"C11": -3.0, "C12": -1.0, "C44": -1.0}[label],
                    order=3,
                )
            )
            for label, (c0, cp) in component_parameters.items()
        }
        matrix = reconstruct_stiffness_from_components(components, "cubic")
        a = volume ** (1.0 / 3.0)
        points.append(
            ElasticVolumePoint(
                source=f"point-{index}.out",
                pressure=float(index - 4),
                stress_pressure=float(index - 4),
                volume=float(volume),
                density=mass / (volume * 1.0e-30),
                energy=float(energy),
                stiffness=matrix,
                lattice=np.diag([a, a, a]),
            )
        )
    reference = CrystalStructure(
        lattice=points[4].lattice,
        fractional_positions=np.asarray([[0.0, 0.0, 0.0]]),
        atomic_numbers=np.asarray([12]),
    )
    series = ElasticVolumeSeries(
        points=tuple(points),
        reference_structure=reference,
        symmetry=SymmetryMetadata(
            space_group_number=221,
            international_symbol="Pm-3m",
            point_group="m-3m",
        ),
        elastic_symmetry="cubic",
        reference_index=4,
    )
    temperature = np.asarray([300.0, 500.0, 700.0])
    pressure = np.asarray([0.0, 5.0, 10.0])
    equilibrium_volume = np.asarray(
        [
            [100.5, 97.0, 94.0],
            [102.0, 98.5, 95.5],
            [103.5, 100.0, 97.0],
        ],
        dtype=np.float64,
    )
    qha = QHAResult(
        jobname="synthetic",
        temperature=temperature,
        pressure=pressure,
        volume=volumes,
        static_energy=energies,
        equilibrium_volume=equilibrium_volume,
        uncertainties={"sigma_VT": np.full_like(equilibrium_volume, 0.02)},
    )
    qha_result = ResultData(
        metadata=ResultMetadata(module="qha", method="quasi-harmonic"),
        options={
            "temperature_unit": "K",
            "pressure_unit": "GPa",
            "volume_unit": "A",
            "energy_unit": "Ha",
        },
        results={"qha": qha},
    )
    context = ThermoelasticContext(
        input_data=ThermoelasticInput(
            jobname="synthetic QSA",
            elastic_series=series,
        ),
        qha_result_data=qha_result,
        qha=qha,
        extrapolation_mask=np.zeros_like(equilibrium_volume, dtype=np.bool_),
    )
    profile = ThermoelasticDepthProfile.linear(
        name="lithosphere",
        depth_min=0.0,
        depth_max=100.0,
        npoints=6,
        pressure_gradient=0.05,
        temperature_at_depth_min=300.0,
        temperature_gradient=2.0,
    )
    calibrated = run_thermoelastic_context(
        context,
        options=ThermoelasticOptions(),
    )
    return analyze_thermoelastic_result(calibrated, profiles=[profile])


def test_reusable_plot_primitives_are_present_in_profile_spec() -> None:
    result = _plot_result()
    collection = build_thermoelastic_profile_plots(
        result,
        profile_name="lithosphere",
        components=["C11", "C12"],
    )
    spec = collection.plots[0]

    assert spec.invert_y_axis
    assert spec.colored_paths
    assert isinstance(spec.colored_paths[0], ColoredPathSeries)
    assert not spec.bands
    assert spec.secondary_axes
    assert isinstance(spec.secondary_axes[0], SecondaryAxis)


def test_profile_pressure_axis_uses_two_decimal_labels() -> None:
    collection = build_thermoelastic_profile_plots(
        _plot_result(),
        profile_name="lithosphere",
        components=["C11"],
    )
    pressure_axis = collection.plots[0].secondary_axes[0]

    assert pressure_axis.key == "pressure"
    assert all(label.count(".") == 1 for label in pressure_axis.labels)
    assert all(len(label.rsplit(".", 1)[1]) == 2 for label in pressure_axis.labels)


def test_auto_uncertainty_bands_single_component() -> None:
    result = _plot_result()
    single = build_thermoelastic_profile_plots(
        result,
        profile_name="lithosphere",
        components=["C11"],
    ).plots[0]
    multiple = build_thermoelastic_profile_plots(
        result,
        profile_name="lithosphere",
        components=["C11", "C12"],
    ).plots[0]

    assert single.bands
    assert isinstance(single.bands[0], PlotBand)
    assert not multiple.bands


def test_colored_profile_markers_have_black_edges() -> None:
    collection = build_thermoelastic_profile_plots(
        _plot_result(),
        profile_name="lithosphere",
        components=["C11"],
    )
    path = collection.plots[0].colored_paths[0]

    assert path.style.marker_edge_color == "black"
    assert path.style.marker_edge_width == 0.7


def test_component_styles_are_stable_and_groups_resolve() -> None:
    payload = _plot_result().results["thermoelasticity"]

    assert component_style("C11") == component_style("11")
    assert component_style("C66").marker != component_style("C11").marker
    assert component_style("C45").marker is None
    assert component_style("C45").marker not in {"$1$", "$2$", "$3$", "$4$"}
    assert resolve_components(payload, group="independent") == (
        "C11",
        "C12",
        "C44",
    )
    assert resolve_components(payload, group="normal") == ("C11", "C22", "C33")


def test_fit_builder_creates_fit_and_residual_panel() -> None:
    collection = build_thermoelastic_fit_plots(
        _plot_result(),
        components=["C11"],
    )
    spec = collection.plots[0]

    assert isinstance(spec, PanelPlotSpec)
    assert len(spec.panels) == 2
    assert spec.panels[0].bands
    assert spec.panels[0].series[0].label == "Observed"
    assert spec.panels[1].key == "C11_residuals"


def test_pt_builder_maps_values_and_extrapolation() -> None:
    result = _plot_result()
    payload = result.results["thermoelasticity"]
    payload.extrapolation_mask[0, -1] = True
    payload.qha_extrapolation_mask[-1, 0] = True
    collection = build_thermoelastic_pt_plots(
        result,
        components=["C11", "C12"],
    )
    spec = collection.plots[0]

    assert isinstance(spec, PanelPlotSpec)
    assert all(isinstance(panel, ContourPlotSpec) for panel in spec.panels)
    assert all(panel.masks for panel in spec.panels)
    assert isinstance(spec.panels[0].masks[0], PlotMask)


def test_relative_profile_uncertainty_is_zero_at_reference() -> None:
    collection = build_thermoelastic_profile_plots(
        _plot_result(),
        profile_name="lithosphere",
        components=["C11"],
        options=ThermoelasticProfilePlotOptions(
            mode="relative",
            color_by="component",
            uncertainty="band",
        ),
    )
    spec = collection.plots[0]
    band = spec.bands[0]

    assert band.lower[0] == band.upper[0] == 0.0
    assert spec.metadata["reference_depth_km"] == 0.0


def test_domain_builder_overlays_depth_colored_profile() -> None:
    collection = build_thermoelastic_domain_plot(_plot_result())
    spec = collection.plots[0]

    assert isinstance(spec, ContourPlotSpec)
    assert spec.colored_paths
    assert spec.colored_paths[0].value_axis.key == "depth"
    assert spec.colored_paths[0].style.colormap == "plasma"
    assert spec.colored_paths[0].style.marker_edge_color == "black"


def test_matplotlib_renders_all_thermoelastic_plot_families(tmp_path: Path) -> None:
    result = _plot_result()
    collections: list[PlotCollection] = [
        build_thermoelastic_fit_plots(result, components=["C11"]),
        build_thermoelastic_pt_plots(result, components=["C11"]),
        build_thermoelastic_profile_plots(
            result,
            profile_name="lithosphere",
            components=["C11", "C12"],
        ),
        build_thermoelastic_domain_plot(result),
    ]
    artifacts = []
    for collection in collections:
        rendered = render_plot_collection(
            collection,
            MatplotlibOptions(output_dir=tmp_path, close=True),
        )
        artifacts.extend(rendered.artifacts)

    assert len(artifacts) == 4
    assert all(
        artifact.path is not None and artifact.path.exists() for artifact in artifacts
    )


def test_generic_renderer_handles_masks_spans_and_secondary_axes(
    tmp_path: Path,
) -> None:
    result = _plot_result()
    profile_collection = build_thermoelastic_profile_plots(
        result,
        profile_name="lithosphere",
        components=["C11"],
    )
    profile_spec = profile_collection.plots[0]
    profile_spec.spans.append(
        PlotSpan(
            key="test_span",
            label="Diagnostic interval",
            axis="y",
            start=20.0,
            end=40.0,
            hatch="///",
        )
    )
    pt_collection = build_thermoelastic_pt_plots(result, components=["C11"])
    contour = pt_collection.plots[0]
    assert isinstance(contour, ContourPlotSpec)
    contour.masks.append(
        PlotMask(
            key="test_mask",
            label="Diagnostic mask",
            x=contour.x,
            y=contour.y,
            mask=np.eye(contour.y.size, contour.x.size, dtype=np.bool_),
        )
    )

    first = render_plot_collection(
        PlotCollection(plots=[profile_spec]),
        MatplotlibOptions(output_dir=tmp_path, close=True),
    )
    second = render_plot_collection(
        PlotCollection(plots=[contour]),
        MatplotlibOptions(output_dir=tmp_path, close=True),
    )

    assert first.artifacts[0].path is not None
    assert second.artifacts[0].path is not None
