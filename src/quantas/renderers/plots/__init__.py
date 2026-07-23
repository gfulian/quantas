# -*- coding: utf-8 -*-

"""Plot renderers for neutral Quantas plot specifications."""

from .presets import FIGURE_PRESET_NAMES, FigurePreset, FigurePresetSpec, figure_preset
from .matplotlib import (
    MatplotlibOptions,
    PlotArtifact,
    PlotRenderResult,
    render_plot,
    render_plot_collection,
    validate_colormap_name,
)

__all__ = [
    "FIGURE_PRESET_NAMES",
    "FigurePreset",
    "FigurePresetSpec",
    "figure_preset",
    "MatplotlibOptions",
    "PlotArtifact",
    "PlotRenderResult",
    "render_plot",
    "render_plot_collection",
    "validate_colormap_name",
]
