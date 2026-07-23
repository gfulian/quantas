from __future__ import annotations

import matplotlib

matplotlib.use("Agg", force=True)

import numpy as np
import pytest

from quantas.models import LinePlotSpec, PlotCollection
from quantas.modules.ha.models import HAResult
from quantas.modules.ha.plot import (
    build_ha_plot_collection,
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
