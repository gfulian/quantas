# -*- coding: utf-8 -*-

"""Public dispatcher for the Matplotlib plotting backend."""

from __future__ import annotations

from itertools import cycle
from typing import Any, Protocol

import numpy as np

from quantas.models.plot import (
    ContourPlotSpec,
    LinePlotSpec,
    PanelPlotSpec,
    PlotCollection,
    PlotSpec,
    PolarPlotSpec,
    SphericalMapSpec,
    SphericalSummarySpec,
    SurfacePlotSpec,
)

from .basic import (
    _render_contour_plot,
    _render_line_plot,
    _render_panel_plot,
    _render_polar_plot,
)
from .common import pyplot
from .io import save_figure
from .options import MatplotlibOptions, PlotArtifact, PlotRenderResult
from .spherical import _render_spherical_map, _render_spherical_summary
from .surface import _render_surface_plot


class _FigureLike(Protocol):
    """Minimal Matplotlib figure surface required by the dispatcher."""

    axes: list[Any]
    legends: list[Any]


_MONO_LINE_STYLES = ("-", "--", "-.", ":")


def _grayscale(colors: Any) -> np.ndarray:
    """Return RGBA colors converted to luminance-preserving grayscale."""
    values = np.asarray(colors, dtype=float)
    if values.size == 0:
        return values
    converted = values.copy()
    luminance = (
        0.2126 * converted[..., 0]
        + 0.7152 * converted[..., 1]
        + 0.0722 * converted[..., 2]
    )
    converted[..., 0] = luminance
    converted[..., 1] = luminance
    converted[..., 2] = luminance
    return converted


def _apply_monochrome(figure: _FigureLike) -> None:
    """Convert one rendered figure to a grayscale-safe presentation."""
    for axis in figure.axes:
        styles = cycle(_MONO_LINE_STYLES)
        for line in axis.lines:
            line.set_color("black")
            if line.get_linestyle() not in {"", "None", "none"}:
                line.set_linestyle(next(styles))
            if line.get_marker() not in {None, "", "None", "none", " "}:
                line.set_markerfacecolor("white")
                line.set_markeredgecolor("black")
        for collection in axis.collections:
            if getattr(collection, "get_array", lambda: None)() is not None:
                collection.set_cmap("Greys")
            get_face = getattr(collection, "get_facecolors", None)
            set_face = getattr(collection, "set_facecolors", None)
            if get_face is not None and set_face is not None:
                face = get_face()
                if np.asarray(face).size:
                    set_face(_grayscale(face))
            get_edge = getattr(collection, "get_edgecolors", None)
            set_edge = getattr(collection, "set_edgecolors", None)
            if get_edge is not None and set_edge is not None:
                edge = get_edge()
                if np.asarray(edge).size:
                    set_edge(_grayscale(edge))
        for image in axis.images:
            image.set_cmap("Greys")
        for patch in axis.patches:
            try:
                patch.set_facecolor(_grayscale([patch.get_facecolor()])[0])
                patch.set_edgecolor(_grayscale([patch.get_edgecolor()])[0])
            except (TypeError, ValueError):
                continue
        for text in axis.texts:
            text.set_color("black")
        axis.xaxis.label.set_color("black")
        axis.yaxis.label.set_color("black")
        axis.title.set_color("black")
        if hasattr(axis, "zaxis"):
            axis.zaxis.label.set_color("black")
        for spine in axis.spines.values():
            spine.set_color("black")


def _apply_typography(figure: _FigureLike, options: MatplotlibOptions) -> None:
    """Apply optional renderer-level font sizes to one figure.

    Parameters
    ----------
    figure : matplotlib.figure.Figure
        Figure whose axes, legends, and titles are updated.
    options : MatplotlibOptions
        Renderer options containing optional font-size overrides.

    Returns
    -------
    None
        The figure is modified in place.
    """
    for axis in figure.axes:
        if options.axis_label_font_size is not None:
            axis.xaxis.label.set_fontsize(options.axis_label_font_size)
            axis.yaxis.label.set_fontsize(options.axis_label_font_size)
            if hasattr(axis, "zaxis"):
                axis.zaxis.label.set_fontsize(options.axis_label_font_size)
        if options.title_font_size is not None:
            axis.title.set_fontsize(options.title_font_size)
        if options.tick_label_font_size is not None:
            axis.tick_params(
                axis="both",
                which="both",
                labelsize=options.tick_label_font_size,
            )
            if hasattr(axis, "zaxis"):
                axis.tick_params(
                    axis="z",
                    which="both",
                    labelsize=options.tick_label_font_size,
                )
        legend = axis.get_legend()
        if legend is not None and options.legend_font_size is not None:
            for text in legend.get_texts():
                text.set_fontsize(options.legend_font_size)
            title = legend.get_title()
            if title is not None:
                title.set_fontsize(options.legend_font_size)

    if options.title_font_size is not None:
        suptitle = getattr(figure, "_suptitle", None)
        if suptitle is not None:
            suptitle.set_fontsize(options.title_font_size)
    if options.legend_font_size is not None:
        for legend in figure.legends:
            for text in legend.get_texts():
                text.set_fontsize(options.legend_font_size)
            title = legend.get_title()
            if title is not None:
                title.set_fontsize(options.legend_font_size)


def render_plot(
    spec: PlotSpec,
    options: MatplotlibOptions | None = None,
) -> PlotArtifact:
    """Render one neutral plot specification with Matplotlib.

    Parameters
    ----------
    spec : PlotSpec
        Neutral plot specification.
    options : MatplotlibOptions or None, optional
        Rendering and file-output options.

    Returns
    -------
    PlotArtifact
        Rendered figure and optional output path.

    Raises
    ------
    TypeError
        If the specification type is unsupported.
    """
    opts = options or MatplotlibOptions()
    plt = pyplot()

    if isinstance(spec, LinePlotSpec):
        figure = _render_line_plot(plt, spec)
    elif isinstance(spec, ContourPlotSpec):
        figure = _render_contour_plot(plt, spec)
    elif isinstance(spec, PolarPlotSpec):
        figure = _render_polar_plot(plt, spec)
    elif isinstance(spec, PanelPlotSpec):
        figure = _render_panel_plot(plt, spec)
    elif isinstance(spec, SurfacePlotSpec):
        figure = _render_surface_plot(plt, spec, opts)
    elif isinstance(spec, SphericalMapSpec):
        figure = _render_spherical_map(plt, spec, opts)
    elif isinstance(spec, SphericalSummarySpec):
        figure = _render_spherical_summary(plt, spec, opts)
    else:
        raise TypeError(f"unsupported plot specification: {type(spec).__name__}")

    if opts.figure_size is not None:
        figure.set_size_inches(*opts.figure_size, forward=True)
    _apply_typography(figure, opts)
    if opts.monochrome:
        _apply_monochrome(figure)
    if opts.tight_layout:
        figure.tight_layout()

    path = save_figure(figure, spec, opts)
    if opts.show:
        plt.show()
    if opts.close:
        plt.close(figure)
    return PlotArtifact(spec=spec, figure=figure, path=path)


def render_plot_collection(
    collection: PlotCollection,
    options: MatplotlibOptions | None = None,
) -> PlotRenderResult:
    """Render all specifications in a neutral plot collection.

    Parameters
    ----------
    collection : PlotCollection
        Neutral specifications and preparation warnings.
    options : MatplotlibOptions or None, optional
        Rendering and file-output options.

    Returns
    -------
    PlotRenderResult
        Concrete Matplotlib artifacts and propagated warnings.
    """
    opts = options or MatplotlibOptions()
    return PlotRenderResult(
        artifacts=[render_plot(spec, opts) for spec in collection.plots],
        warnings=list(collection.warnings),
    )


__all__ = ["render_plot", "render_plot_collection"]
