"""Shared Matplotlib figure-preset contracts."""

from __future__ import annotations

import matplotlib
import pytest

matplotlib.use("Agg")

from quantas.models import (
    LinePlotSpec,
    PlotAxis,
    PlotCollection,
    PlotSeries,
    PlotSeriesStyle,
)
from quantas.renderers.plots import (
    FIGURE_PRESET_NAMES,
    MatplotlibOptions,
    figure_preset,
    render_plot_collection,
)

pytestmark = [pytest.mark.elasticity, pytest.mark.plotting]


def test_named_presets_define_resolution_and_monochrome_policy() -> None:
    assert FIGURE_PRESET_NAMES == ("screen", "publication", "monochrome")
    assert figure_preset("SCREEN").dpi == 150
    assert figure_preset("publication").dpi == 300
    assert figure_preset("monochrome").monochrome is True


def test_explicit_renderer_options_override_preset_defaults() -> None:
    options = MatplotlibOptions.from_preset(
        "publication",
        dpi=450,
        title_font_size=15.0,
        savefig_kwargs={"transparent": True},
    )
    assert options.preset == "publication"
    assert options.dpi == 450
    assert options.title_font_size == 15.0
    assert options.savefig_kwargs == {"bbox_inches": "tight", "transparent": True}


def test_invalid_preset_is_rejected() -> None:
    with pytest.raises(ValueError, match="unknown figure preset"):
        MatplotlibOptions.from_preset("poster")
    with pytest.raises(ValueError, match="unknown figure preset"):
        MatplotlibOptions(preset="poster")  # type: ignore[arg-type]


def test_monochrome_preset_converts_lines_and_markers() -> None:
    collection = PlotCollection(
        plots=[
            LinePlotSpec(
                key="preset-check",
                title="Preset check",
                filename_stem="preset-check",
                x_axis=PlotAxis(key="x", label="x"),
                y_axis=PlotAxis(key="y", label="y"),
                series=[
                    PlotSeries(
                        key="first",
                        x=[0.0, 1.0],
                        y=[0.0, 1.0],
                        label="first",
                        style=PlotSeriesStyle(color="tab:red", marker="o"),
                    ),
                    PlotSeries(
                        key="second",
                        x=[0.0, 1.0],
                        y=[1.0, 0.0],
                        label="second",
                        style=PlotSeriesStyle(color="tab:blue", marker="s"),
                    ),
                ],
            )
        ]
    )
    rendered = render_plot_collection(
        collection,
        MatplotlibOptions.from_preset("monochrome", close=False),
    )
    figure = rendered.artifacts[0].figure
    lines = figure.axes[0].lines
    assert [line.get_color() for line in lines] == ["black", "black"]
    assert [line.get_markerfacecolor() for line in lines] == ["white", "white"]
    assert lines[0].get_linestyle() != lines[1].get_linestyle()
    import matplotlib.pyplot as plt

    plt.close(figure)
