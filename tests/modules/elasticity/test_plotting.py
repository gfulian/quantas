"""Tests for neutral elasticity polar-plot builders."""

from __future__ import annotations

import numpy as np

from quantas.models import PlotCollection, PolarPlotSpec
from quantas.modules.elasticity.models import ElasticityResult
from quantas.modules.elasticity.plot import (
    build_elasticity_plot_collection,
    build_elasticity_2d_plot_spec,
)
from quantas.renderers.plots import MatplotlibOptions, render_plot_collection


def _young_result() -> ElasticityResult:
    angle = np.linspace(0.0, 2.0 * np.pi, 9)
    values = 100.0 + 10.0 * np.cos(2.0 * angle)
    return ElasticityResult(
        properties_2d={
            "xy": {"phi": angle, "young_modulus": values},
            "xz": {"theta": angle, "young_modulus": values + 5.0},
            "yz": {"theta": angle, "young_modulus": values + 10.0},
        }
    )


def test_builder_creates_three_frontend_neutral_panels() -> None:
    spec = build_elasticity_2d_plot_spec(_young_result(), "young")
    assert isinstance(spec, PolarPlotSpec)
    assert spec.metadata["module"] == "elasticity"
    assert [panel.key for panel in spec.panels] == ["xy", "xz", "yz"]
    assert all(panel.angle_unit == "degree" for panel in spec.panels)
    np.testing.assert_allclose(
        spec.panels[0].series[0].x,
        np.linspace(0.0, 360.0, 9),
    )


def test_collection_renders_and_saves_polar_plot(tmp_path) -> None:
    collection = build_elasticity_plot_collection(_young_result())
    assert isinstance(collection, PlotCollection)
    rendered = render_plot_collection(
        collection,
        MatplotlibOptions(
            output_dir=tmp_path,
            filename_prefix="MgO_",
            close=True,
        ),
    )
    assert rendered.artifacts[0].path == tmp_path / "MgO_2d_young.png"
    assert rendered.artifacts[0].path.exists()
    assert len(rendered.artifacts[0].figure.axes) == 3


def test_phase_velocity_is_not_an_elasticity_plot_property() -> None:
    collection = build_elasticity_plot_collection(
        ElasticityResult(
            properties_2d={
                plane: {
                    "theta": np.zeros(3),
                    "phi": np.zeros(3),
                    "phase_velocity": np.ones((3, 3)),
                }
                for plane in ("xy", "xz", "yz")
            }
        )
    )
    assert collection.plots == []
