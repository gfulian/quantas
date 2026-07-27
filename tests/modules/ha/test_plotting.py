from __future__ import annotations

import matplotlib

matplotlib.use("Agg", force=True)

import numpy as np
import pytest

from quantas.models import ContourPlotSpec, LinePlotSpec, PlotCollection
from quantas.modules.ha.models import HAResult
from quantas.modules.ha.plot import (
    HAPlotOptions,
    build_ha_plot_collection,
    build_thermodynamic_contour_spec,
    build_thermodynamic_plot_spec,
)
from quantas.renderers.plots import (
    MatplotlibOptions,
    render_plot,
    render_plot_collection,
)


@pytest.fixture
def ha_result() -> HAResult:
    temperature = np.array([0.0, 100.0, 200.0], dtype=float)
    volume = np.array([18.0, 19.0], dtype=float)
    free_energy = np.array(
        [
            [-10.0, -9.8],
            [-9.5, -9.3],
            [-9.0, -8.8],
        ],
        dtype=float,
    )
    entropy = np.array([0.0, 0.1, 0.2], dtype=float)
    return HAResult(
        jobname="MgO",
        temperature=temperature,
        volume=volume,
        free_energy=free_energy,
        entropy=entropy,
    )


def test_build_thermodynamic_plot_spec_is_frontend_neutral(ha_result: HAResult) -> None:
    spec = build_thermodynamic_plot_spec(ha_result, "F")

    assert isinstance(spec, LinePlotSpec)
    assert spec.key == "free_energy"
    assert spec.x_axis.label == "Temperature (K)"
    assert len(spec.series) == 2
    np.testing.assert_allclose(spec.series[0].y, [-10.0, -9.5, -9.0])


def test_rendered_thermodynamic_property_returns_figure(ha_result: HAResult) -> None:
    artifact = render_plot(build_thermodynamic_plot_spec(ha_result, "F"))
    figure = artifact.figure

    assert figure.axes
    assert figure.axes[0].get_xlabel() == "Temperature (K)"
    assert len(figure.axes[0].lines) == 2


def test_rendered_thermodynamic_property_saves_file(
    tmp_path,
    ha_result: HAResult,
) -> None:
    spec = build_thermodynamic_plot_spec(ha_result, "free_energy")
    artifact = render_plot(
        spec,
        MatplotlibOptions(
            output_dir=tmp_path,
            filename_prefix="ha_",
            close=True,
        ),
    )

    assert artifact.path == tmp_path / "ha_free_energy.png"
    assert artifact.path.exists()
    assert artifact.path.stat().st_size > 0
    assert artifact.figure.axes[0].get_title() == "Helmholtz free energy"


def test_builder_accepts_one_dimensional_property(ha_result: HAResult) -> None:
    spec = build_thermodynamic_plot_spec(ha_result, "entropy")

    assert len(spec.series) == 1


def test_builder_rejects_unknown_property(ha_result: HAResult) -> None:
    with pytest.raises(KeyError):
        build_thermodynamic_plot_spec(ha_result, "not_a_property")


def test_builder_requires_temperature() -> None:
    result = HAResult(free_energy=np.array([[1.0]]))

    with pytest.raises(ValueError, match="temperatures"):
        build_thermodynamic_plot_spec(result, "F")


def test_builder_rejects_incompatible_shape(ha_result: HAResult) -> None:
    ha_result.free_energy = np.ones((2, 2), dtype=float)

    with pytest.raises(ValueError, match="temperature dimension"):
        build_thermodynamic_plot_spec(ha_result, "F")


def test_ha_collection_renders_requested_properties(
    tmp_path, ha_result: HAResult
) -> None:
    collection = build_ha_plot_collection(ha_result, properties=["F", "S"])

    assert isinstance(collection, PlotCollection)
    rendered = render_plot_collection(
        collection,
        MatplotlibOptions(
            output_dir=tmp_path,
            filename_prefix="ha_",
            close=True,
        ),
    )

    assert {artifact.path.name for artifact in rendered.artifacts} == {
        "ha_free_energy.png",
        "ha_entropy.png",
    }


def test_zero_point_plot_broadcasts_temperature_independent_volume_series(ha_result: HAResult) -> None:
    ha_result.zero_point_energy = np.array([[0.10, 0.20]], dtype=np.float64)

    spec = build_thermodynamic_plot_spec(ha_result, "Uzp")

    assert len(spec.series) == 2
    np.testing.assert_allclose(spec.series[0].y, [0.10, 0.10, 0.10])
    np.testing.assert_allclose(spec.series[1].y, [0.20, 0.20, 0.20])


def test_temperature_sections_select_exact_native_volumes(ha_result: HAResult) -> None:
    spec = build_thermodynamic_plot_spec(
        ha_result,
        "F",
        options=HAPlotOptions(
            curve_axis="temperature",
            selected_volumes=(19.0,),
        ),
    )

    assert spec.x_axis.key == "temperature"
    assert len(spec.series) == 1
    np.testing.assert_allclose(spec.series[0].x, [0.0, 100.0, 200.0])
    np.testing.assert_allclose(spec.series[0].y, [-9.8, -9.3, -8.8])
    assert spec.series[0].metadata["volume"] == 19.0


def test_volume_sections_select_exact_native_temperatures(ha_result: HAResult) -> None:
    spec = build_thermodynamic_plot_spec(
        ha_result,
        "F",
        options=HAPlotOptions(
            curve_axis="volume",
            selected_temperatures=(100.0,),
        ),
    )

    assert spec.x_axis.key == "volume"
    assert spec.filename_stem == "free_energy_vs_volume"
    assert len(spec.series) == 1
    np.testing.assert_allclose(spec.series[0].x, [18.0, 19.0])
    np.testing.assert_allclose(spec.series[0].y, [-9.5, -9.3])
    assert spec.series[0].metadata["temperature_native"] == 100.0


def test_volume_temperature_contour_preserves_native_grid(ha_result: HAResult) -> None:
    spec = build_thermodynamic_contour_spec(ha_result, "F")

    assert isinstance(spec, ContourPlotSpec)
    assert spec.x_axis.key == "temperature"
    assert spec.y_axis.key == "volume"
    np.testing.assert_allclose(spec.x, [0.0, 100.0, 200.0])
    np.testing.assert_allclose(spec.y, [18.0, 19.0])
    np.testing.assert_allclose(
        spec.z,
        [[-10.0, -9.5, -9.0], [-9.8, -9.3, -8.8]],
    )


def test_ha_sections_reject_coordinates_absent_from_native_grid(
    ha_result: HAResult,
) -> None:
    with pytest.raises(ValueError, match="not present in the native grid"):
        build_thermodynamic_plot_spec(
            ha_result,
            "F",
            options=HAPlotOptions(
                curve_axis="volume",
                selected_temperatures=(150.0,),
            ),
        )


def test_ha_collection_can_build_line_and_contour_without_mutation(
    ha_result: HAResult,
) -> None:
    original = np.array(ha_result.free_energy, copy=True)
    collection = build_ha_plot_collection(
        ha_result,
        properties=("F",),
        options=HAPlotOptions(include_contours=True),
    )

    assert [type(item).__name__ for item in collection.plots] == [
        "LinePlotSpec",
        "ContourPlotSpec",
    ]
    np.testing.assert_array_equal(ha_result.free_energy, original)


def test_volume_sections_require_multiple_sampled_volumes(
    ha_result: HAResult,
) -> None:
    ha_result.volume = np.array([18.0], dtype=np.float64)
    ha_result.free_energy = np.asarray(ha_result.free_energy)[:, :1]

    with pytest.raises(ValueError, match="at least two sampled volumes"):
        build_thermodynamic_plot_spec(
            ha_result,
            "F",
            options=HAPlotOptions(curve_axis="volume"),
        )


def test_ha_sections_reject_duplicate_native_temperatures(
    ha_result: HAResult,
) -> None:
    ha_result.temperature = np.array([0.0, 100.0, 100.0], dtype=np.float64)

    with pytest.raises(ValueError, match="duplicate coordinates"):
        build_thermodynamic_plot_spec(ha_result, "F")
