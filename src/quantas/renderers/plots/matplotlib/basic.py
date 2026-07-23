# -*- coding: utf-8 -*-

"""Cartesian line, contour, panel, and polar renderers for Matplotlib."""

from __future__ import annotations

from math import ceil
from typing import Any

import numpy as np

from quantas.models.plot import (
    ColoredPathSeries,
    ContourPlotSpec,
    LinePlotSpec,
    PanelPlotSpec,
    PlotBand,
    PlotMask,
    PlotSeries,
    PlotSpan,
    PolarPlotSpec,
    ScalarBackground,
    SecondaryAxis,
)

from .colors import resolve_colormap
from .common import LINE_STYLES, apply_x_limits, apply_y_limits


def _render_line_plot(plt: Any, spec: LinePlotSpec) -> Any:
    """Render a Cartesian line specification.

    Parameters
    ----------
    plt : module
        Imported :mod:`matplotlib.pyplot` module.
    spec : LinePlotSpec
        Neutral line-plot specification.

    Returns
    -------
    matplotlib.figure.Figure
        Rendered Matplotlib figure.
    """
    figure, axis = plt.subplots()
    _draw_line_plot(plt, figure, axis, spec)
    figure.tight_layout()
    return figure


def _draw_line_plot(plt: Any, figure: Any, axis: Any, spec: LinePlotSpec) -> None:
    """Draw one line specification into an existing Matplotlib axis."""
    for span in spec.spans:
        _draw_span(axis, span)
    for band in spec.bands:
        _draw_band(axis, band)
    for series in spec.series:
        _draw_series(axis, series)

    colorbar_pad = _line_colorbar_pad(spec)
    rendered_colorbars: set[str] = set()
    for path in spec.colored_paths:
        key = path.value_axis.key
        show_colorbar = path.style.show_colorbar and key not in rendered_colorbars
        _draw_colored_path(
            plt,
            figure,
            axis,
            path,
            show_colorbar=show_colorbar,
            colorbar_pad=colorbar_pad,
        )
        if show_colorbar:
            rendered_colorbars.add(key)

    axis.set_xlabel(spec.x_axis.label)
    axis.set_ylabel(spec.y_axis.label)
    axis.set_title(spec.title)
    apply_x_limits(axis, spec.x_axis.limits)
    apply_y_limits(axis, spec.y_axis.limits)
    if spec.backgrounds:
        x_limits = axis.get_xlim()
        y_limits = axis.get_ylim()
        for background in spec.backgrounds:
            _draw_scalar_background(
                plt,
                figure,
                axis,
                background,
                colorbar_pad=colorbar_pad,
            )
        axis.set_xlim(x_limits)
        axis.set_ylim(y_limits)
    if spec.invert_x_axis:
        axis.invert_xaxis()
    if spec.invert_y_axis:
        axis.invert_yaxis()
    if spec.grid:
        axis.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    for secondary in spec.secondary_axes:
        _draw_secondary_axis(axis, secondary)
    if spec.show_legend and (
        spec.series or spec.bands or spec.colored_paths or spec.spans
    ):
        handles, labels = axis.get_legend_handles_labels()
        if labels:
            axis.legend(
                handles,
                labels,
                title=spec.legend_title,
                fontsize="small",
                ncol=max(1, spec.legend_columns),
            )


def _line_colorbar_pad(spec: LinePlotSpec) -> float:
    """Return colorbar spacing that clears prepared right-side axes."""
    has_right_axis = any(
        secondary.orientation == "y" and secondary.location == "right"
        for secondary in spec.secondary_axes
    )
    return 0.30 if has_right_axis else 0.04


def _draw_series(axis: Any, series: PlotSeries) -> None:
    """Draw one ordinary Cartesian series."""
    kwargs: dict[str, Any] = {
        "label": series.label,
        "linestyle": LINE_STYLES[series.style.line_style],
        "linewidth": series.style.line_width,
        "alpha": series.style.alpha,
    }
    if series.style.color is not None:
        kwargs["color"] = series.style.color
    if series.style.marker is not None:
        kwargs["marker"] = series.style.marker
    if series.style.marker_size is not None:
        kwargs["markersize"] = series.style.marker_size
    if series.style.marker_edge_color is not None:
        kwargs["markeredgecolor"] = series.style.marker_edge_color
    if series.style.marker_edge_width is not None:
        kwargs["markeredgewidth"] = series.style.marker_edge_width
    if series.x_error is not None or series.y_error is not None:
        axis.errorbar(
            series.x,
            series.y,
            xerr=series.x_error,
            yerr=series.y_error,
            capsize=series.style.errorbar_capsize,
            elinewidth=series.style.errorbar_line_width,
            **kwargs,
        )
    else:
        axis.plot(series.x, series.y, **kwargs)


def _draw_band(axis: Any, band: PlotBand) -> None:
    """Draw one vertical or horizontal confidence band."""
    kwargs: dict[str, Any] = {
        "label": band.label,
        "alpha": band.style.alpha,
        "linewidth": band.style.edge_width,
        "linestyle": LINE_STYLES[band.style.line_style],
    }
    if band.style.color is not None:
        kwargs["facecolor"] = band.style.color
    if band.style.edge_color is not None:
        kwargs["edgecolor"] = band.style.edge_color
    if band.orientation == "vertical":
        axis.fill_between(
            band.coordinates,
            band.lower,
            band.upper,
            **kwargs,
        )
    else:
        axis.fill_betweenx(
            band.coordinates,
            band.lower,
            band.upper,
            **kwargs,
        )


def _draw_span(axis: Any, span: PlotSpan) -> None:
    """Draw one highlighted interval behind line data."""
    kwargs: dict[str, Any] = {
        "label": span.label,
        "alpha": span.alpha,
        "hatch": span.hatch,
        "zorder": 0,
    }
    if span.color is not None:
        kwargs["facecolor"] = span.color
    if span.axis == "x":
        axis.axvspan(span.start, span.end, **kwargs)
    else:
        axis.axhspan(span.start, span.end, **kwargs)


def _draw_colored_path(
    plt: Any,
    figure: Any,
    axis: Any,
    path: ColoredPathSeries,
    *,
    show_colorbar: bool,
    colorbar_pad: float = 0.04,
) -> None:
    """Draw one scalar-colored Cartesian path."""
    x = np.asarray(path.x, dtype=np.float64)
    y = np.asarray(path.y, dtype=np.float64)
    values = np.asarray(path.values, dtype=np.float64)
    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        return

    cmap = resolve_colormap(plt, path.style.colormap)
    if path.style.value_limits is None:
        lower = float(np.min(finite_values))
        upper = float(np.max(finite_values))
        if np.isclose(lower, upper):
            pad = max(abs(lower), 1.0) * 1.0e-12
            lower -= pad
            upper += pad
    else:
        lower, upper = path.style.value_limits
    norm = plt.matplotlib.colors.Normalize(vmin=lower, vmax=upper)

    mappable: Any
    if x.size >= 2:
        points = np.column_stack((x, y))
        segments = np.stack((points[:-1], points[1:]), axis=1)
        segment_values = 0.5 * (values[:-1] + values[1:])
        line = plt.matplotlib.collections.LineCollection(
            segments,
            cmap=cmap,
            norm=norm,
            linewidths=path.style.line_width,
            linestyles=LINE_STYLES[path.style.line_style],
            alpha=path.style.alpha,
        )
        line.set_array(segment_values)
        axis.add_collection(line)
        axis.update_datalim(points)
        axis.autoscale_view()
        mappable = line
    else:
        mappable = plt.matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    if path.style.marker is not None or path.style.marker_size is not None:
        marker = path.style.marker or "o"
        size = (path.style.marker_size or 4.0) ** 2
        scatter = axis.scatter(
            x,
            y,
            c=values,
            cmap=cmap,
            norm=norm,
            marker=marker,
            s=size,
            alpha=path.style.alpha,
            edgecolors=path.style.marker_edge_color,
            linewidths=path.style.marker_edge_width,
            zorder=3,
        )
        mappable = scatter

    if path.label:
        proxy_color = cmap(norm(float(np.nanmedian(finite_values))))
        axis.plot(
            [],
            [],
            label=path.label,
            color=proxy_color,
            linestyle=LINE_STYLES[path.style.line_style],
            linewidth=path.style.line_width,
            marker=path.style.marker,
            markersize=path.style.marker_size,
            markeredgecolor=path.style.marker_edge_color,
            markeredgewidth=path.style.marker_edge_width,
        )
    if show_colorbar:
        colorbar = figure.colorbar(mappable, ax=axis, pad=colorbar_pad)
        colorbar.set_label(path.value_axis.label)


def _draw_secondary_axis(axis: Any, secondary: SecondaryAxis) -> None:
    """Draw one explicit secondary tick mapping."""
    if secondary.orientation == "y":
        target = axis.twinx()
        target.set_ylim(axis.get_ylim())
        target.set_yticks(secondary.positions)
        target.set_yticklabels(secondary.labels)
        target.set_ylabel(secondary.label)
        if secondary.location == "left":
            target.yaxis.set_label_position("left")
            target.yaxis.tick_left()
            target.spines["left"].set_position(("outward", 55))
            target.spines["right"].set_visible(False)
        else:
            target.yaxis.set_label_position("right")
            target.yaxis.tick_right()
            target.spines["right"].set_position(("outward", 28))
            target.set_ylabel(secondary.label, labelpad=8)
    else:
        target = axis.twiny()
        target.set_xlim(axis.get_xlim())
        target.set_xticks(secondary.positions)
        target.set_xticklabels(secondary.labels)
        target.set_xlabel(secondary.label)
        if secondary.location == "bottom":
            target.xaxis.set_label_position("bottom")
            target.xaxis.tick_bottom()
            target.spines["bottom"].set_position(("outward", 45))
            target.spines["top"].set_visible(False)


def _draw_scalar_background(
    plt: Any,
    figure: Any,
    axis: Any,
    background: ScalarBackground,
    *,
    colorbar_pad: float = 0.04,
) -> None:
    """Draw a scalar field varying along one Cartesian axis."""
    coordinates = np.asarray(background.coordinates, dtype=np.float64)
    values = np.asarray(background.values, dtype=np.float64)
    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        return
    if background.value_limits is None:
        lower = float(np.min(finite_values))
        upper = float(np.max(finite_values))
        if np.isclose(lower, upper):
            pad = max(abs(lower), 1.0) * 1.0e-12
            lower -= pad
            upper += pad
    else:
        lower, upper = background.value_limits
    cmap = resolve_colormap(plt, background.colormap)
    norm = plt.matplotlib.colors.Normalize(vmin=lower, vmax=upper)
    edges = _coordinate_edges(coordinates)

    if background.axis == "y":
        xmin, xmax = axis.get_xlim()
        mesh = axis.pcolormesh(
            np.asarray([xmin, xmax], dtype=np.float64),
            edges,
            values[:, np.newaxis],
            cmap=cmap,
            norm=norm,
            shading="flat",
            alpha=background.alpha,
            zorder=-10,
        )
    else:
        ymin, ymax = axis.get_ylim()
        mesh = axis.pcolormesh(
            edges,
            np.asarray([ymin, ymax], dtype=np.float64),
            values[np.newaxis, :],
            cmap=cmap,
            norm=norm,
            shading="flat",
            alpha=background.alpha,
            zorder=-10,
        )
    if background.show_colorbar:
        colorbar = figure.colorbar(mesh, ax=axis, pad=colorbar_pad)
        colorbar.set_label(background.value_axis.label)


def _coordinate_edges(coordinates: np.ndarray) -> np.ndarray:
    """Return cell edges for monotonic one-dimensional center coordinates."""
    if coordinates.size == 1:
        return np.asarray(
            [coordinates[0] - 0.5, coordinates[0] + 0.5],
            dtype=np.float64,
        )
    midpoint = 0.5 * (coordinates[:-1] + coordinates[1:])
    first = coordinates[0] - (midpoint[0] - coordinates[0])
    last = coordinates[-1] + (coordinates[-1] - midpoint[-1])
    return np.concatenate(([first], midpoint, [last])).astype(np.float64)


def _render_contour_plot(plt: Any, spec: ContourPlotSpec) -> Any:
    """Render a Cartesian filled-contour specification.

    Parameters
    ----------
    plt : module
        Imported :mod:`matplotlib.pyplot` module.
    spec : ContourPlotSpec
        Neutral contour specification.

    Returns
    -------
    matplotlib.figure.Figure
        Rendered Matplotlib figure.
    """
    figure, axis = plt.subplots()
    _draw_contour_plot(plt, figure, axis, spec)
    figure.tight_layout()
    return figure


def _draw_contour_plot(
    plt: Any,
    figure: Any,
    axis: Any,
    spec: ContourPlotSpec,
) -> None:
    """Draw one contour specification into an existing Matplotlib axis."""
    xgrid, ygrid = np.meshgrid(spec.x, spec.y)
    filled_levels = (
        spec.levels if spec.mode == "discrete" else max(128, spec.levels * 8)
    )
    kwargs: dict[str, Any] = {"cmap": resolve_colormap(plt, spec.colormap)}
    if spec.value_limits is not None:
        kwargs["vmin"], kwargs["vmax"] = spec.value_limits
    if spec.center is not None:
        finite = np.asarray(spec.z, dtype=np.float64)
        finite = finite[np.isfinite(finite)]
        if finite.size:
            lower = (
                float(np.min(finite))
                if spec.value_limits is None
                else spec.value_limits[0]
            )
            upper = (
                float(np.max(finite))
                if spec.value_limits is None
                else spec.value_limits[1]
            )
            if lower < spec.center < upper:
                kwargs.pop("vmin", None)
                kwargs.pop("vmax", None)
                kwargs["norm"] = plt.matplotlib.colors.TwoSlopeNorm(
                    vmin=lower,
                    vcenter=spec.center,
                    vmax=upper,
                )
    filled = axis.contourf(
        xgrid,
        ygrid,
        spec.z,
        levels=filled_levels,
        **kwargs,
    )
    if spec.isolines:
        lines = axis.contour(
            xgrid,
            ygrid,
            spec.z,
            levels=spec.levels,
            colors="black",
            linewidths=0.55,
            alpha=0.75,
        )
        if spec.isoline_labels:
            axis.clabel(lines, inline=True, fontsize=8)
    colorbar = figure.colorbar(filled, ax=axis)
    colorbar.set_label(spec.value_axis.label)

    mask_handles: list[Any] = []
    for mask in spec.masks:
        handle = _draw_mask(plt, axis, mask)
        if handle is not None:
            mask_handles.append(handle)
    for series in spec.series:
        _draw_series(axis, series)
    rendered_colorbars: set[str] = set()
    for path in spec.colored_paths:
        key = path.value_axis.key
        show_colorbar = path.style.show_colorbar and key not in rendered_colorbars
        _draw_colored_path(plt, figure, axis, path, show_colorbar=show_colorbar)
        if show_colorbar:
            rendered_colorbars.add(key)

    axis.set_xlabel(spec.x_axis.label)
    axis.set_ylabel(spec.y_axis.label)
    axis.set_title(spec.title)
    apply_x_limits(axis, spec.x_axis.limits)
    apply_y_limits(axis, spec.y_axis.limits)
    handles, labels = axis.get_legend_handles_labels()
    if mask_handles:
        handles.extend(mask_handles)
        labels.extend(handle.get_label() for handle in mask_handles)
    if labels:
        axis.legend(handles, labels, fontsize="small")


def _draw_mask(plt: Any, axis: Any, mask: PlotMask) -> Any | None:
    """Draw one hatched boolean mask and return a legend proxy."""
    if not np.any(mask.mask):
        return None
    xgrid, ygrid = np.meshgrid(mask.x, mask.y)
    colors = [mask.color or "none"]
    axis.contourf(
        xgrid,
        ygrid,
        mask.mask.astype(np.float64),
        levels=[0.5, 1.5],
        colors=colors,
        hatches=[mask.hatch],
        alpha=mask.alpha,
    )
    return plt.matplotlib.patches.Patch(
        facecolor=mask.color or "none",
        edgecolor="0.35",
        hatch=mask.hatch,
        label=mask.label,
        alpha=max(mask.alpha, 0.8 if mask.color == "none" else mask.alpha),
    )


def _render_panel_plot(plt: Any, spec: PanelPlotSpec) -> Any:
    """Render a Cartesian multi-panel specification."""
    panel_count = len(spec.panels)
    columns = min(spec.columns, panel_count)
    rows = int(ceil(panel_count / columns))
    figure, axes = plt.subplots(
        rows,
        columns,
        squeeze=False,
        sharex=spec.share_x,
        sharey=spec.share_y,
        figsize=(6.0 * columns, 4.8 * rows),
    )
    flat_axes = axes.ravel()
    for axis, panel in zip(flat_axes, spec.panels, strict=False):
        if isinstance(panel, LinePlotSpec):
            _draw_line_plot(plt, figure, axis, panel)
        else:
            _draw_contour_plot(plt, figure, axis, panel)
    for axis in flat_axes[panel_count:]:
        axis.set_visible(False)
    if spec.title:
        figure.suptitle(spec.title)
    figure.tight_layout()
    return figure


def _render_polar_plot(plt: Any, spec: PolarPlotSpec) -> Any:
    """Render a multi-panel polar specification.

    Parameters
    ----------
    plt : module
        Imported :mod:`matplotlib.pyplot` module.
    spec : PolarPlotSpec
        Neutral polar specification.

    Returns
    -------
    matplotlib.figure.Figure
        Rendered Matplotlib figure.

    Raises
    ------
    ValueError
        If the specification contains no polar panels.
    """
    panel_count = len(spec.panels)
    if panel_count == 0:
        raise ValueError("a polar plot requires at least one panel")
    figure, axes = plt.subplots(
        1,
        panel_count,
        figsize=(6 * panel_count, 6),
        subplot_kw={"projection": "polar"},
        squeeze=False,
    )
    flat_axes = axes.ravel()
    legend_handles: list[Any] = []
    legend_labels: list[str] = []

    for axis, panel in zip(flat_axes, spec.panels, strict=True):
        axis.set_theta_zero_location(panel.theta_zero_location)
        axis.set_theta_direction(panel.theta_direction)
        axis.grid(panel.grid, linestyle="dotted")
        for series in panel.series:
            angle = np.asarray(series.x, dtype=np.float64)
            if panel.angle_unit == "degree":
                angle = np.radians(angle)
            kwargs: dict[str, Any] = {
                "label": series.label,
                "linestyle": LINE_STYLES[series.style.line_style],
                "linewidth": series.style.line_width,
            }
            if series.style.color is not None:
                kwargs["color"] = series.style.color
            if series.style.marker is not None:
                kwargs["marker"] = series.style.marker
            if series.style.marker_size is not None:
                kwargs["markersize"] = series.style.marker_size
            if series.style.marker_edge_color is not None:
                kwargs["markeredgecolor"] = series.style.marker_edge_color
            if series.style.marker_edge_width is not None:
                kwargs["markeredgewidth"] = series.style.marker_edge_width
            line = axis.plot(angle, series.y, **kwargs)[0]
            if series.label not in legend_labels:
                legend_handles.append(line)
                legend_labels.append(series.label)
        if panel.radial_limit is not None:
            axis.set_rmax(panel.radial_limit)
        axis.text(
            0.2,
            -0.05,
            panel.title,
            ha="center",
            va="center",
            transform=axis.transAxes,
            fontsize=12,
        )

    if spec.show_legend and legend_handles:
        figure.legend(
            legend_handles,
            legend_labels,
            fontsize=12,
            ncol=max(1, spec.legend_columns),
        )
    figure.tight_layout()
    return figure


__all__ = [
    "_render_contour_plot",
    "_render_line_plot",
    "_render_panel_plot",
    "_render_polar_plot",
]
