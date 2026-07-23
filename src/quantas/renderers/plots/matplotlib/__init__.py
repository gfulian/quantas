# -*- coding: utf-8 -*-

"""Matplotlib renderer for frontend-neutral Quantas plot specifications.

The backend is split by plot family so that line/contour/polar plots,
spherical maps, three-dimensional surfaces, colors, and file output can evolve
independently while preserving the public :func:`render_plot` and
:func:`render_plot_collection` API.
"""

from __future__ import annotations

from .colors import validate_colormap_name
from .dispatcher import render_plot, render_plot_collection
from .options import MatplotlibOptions, PlotArtifact, PlotRenderResult

__all__ = [
    "MatplotlibOptions",
    "PlotArtifact",
    "PlotRenderResult",
    "render_plot",
    "render_plot_collection",
    "validate_colormap_name",
]
