# -*- coding: utf-8 -*-

"""CLI command for thermoelastic finite-strain fit plots."""

from __future__ import annotations

from pathlib import Path
from typing import cast

import click

from quantas.cli.grouped_options import GroupedCommand, grouped_option
from quantas.cli.thermoelastic_plot_common import (
    component_options,
    make_style,
    read_archive,
    render_collection,
    render_options,
    show_components,
    style_options,
)
from quantas.api.thermoelasticity import (
    ComponentGroup as ThermoelasticComponentGroup,
    FitPlotOptions as ThermoelasticFitPlotOptions,
    build_fit_plots as build_thermoelastic_fit_plots,
)


@click.command(name="fit", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@component_options
@grouped_option(
    "--residuals/--no-residuals",
    group="Fit diagnostics",
    default=True,
    show_default=True,
    help="Include a residual panel beneath each fitted component.",
)
@grouped_option(
    "--confidence-band/--no-confidence-band",
    group="Fit diagnostics",
    default=True,
    show_default=True,
    help="Display the propagated confidence band around the fitted curve.",
)
@grouped_option(
    "--confidence",
    group="Fit diagnostics",
    type=click.FloatRange(min=0.0, max=1.0, min_open=True, max_open=True),
    default=0.95,
    show_default=True,
    help="Central confidence probability used by --confidence-band.",
)
@grouped_option(
    "--symmetry-spread/--no-symmetry-spread",
    group="Fit diagnostics",
    default=False,
    show_default=True,
    help="Show the numerical spread among symmetry-equivalent observed entries.",
)
@grouped_option(
    "--curve-points",
    group="Fit diagnostics",
    hidden=True,
    type=click.IntRange(min=20),
    default=300,
    show_default=True,
    help="Number of samples used for each smooth fitted curve.",
)
@style_options
@render_options(default_width=7.0, default_height=6.5)
def fit_plot(
    archive: Path,
    components: tuple[str, ...],
    component_group: str,
    list_components: bool,
    residuals: bool,
    confidence_band: bool,
    confidence: float,
    symmetry_spread: bool,
    curve_points: int,
    preset: str,
    line_width: float | None,
    marker_size: float | None,
    marker_edge_color: str,
    marker_edge_width: float,
    errorbar_width: float,
    errorbar_capsize: float,
    grid: bool | None,
    title: bool | None,
    output_dir: Path | None,
    image_format: str,
    dpi: int | None,
    transparent: bool,
    show: bool,
    figure_width: float,
    figure_height: float,
    axis_label_font_size: float | None,
    legend_font_size: float | None,
    title_font_size: float | None,
    tick_label_font_size: float | None,
) -> None:
    """Plot observed and fitted independent elastic components versus volume."""
    result_data, payload = read_archive(archive)
    if list_components:
        show_components(payload)
        return
    style = make_style(
        preset=preset,
        line_width=line_width,
        marker_size=marker_size,
        marker_edge_color=marker_edge_color,
        marker_edge_width=marker_edge_width,
        errorbar_width=errorbar_width,
        errorbar_capsize=errorbar_capsize,
        grid=grid,
        title=title,
    )
    collection = build_thermoelastic_fit_plots(
        result_data,
        components=components or None,
        component_group=cast(ThermoelasticComponentGroup, component_group.lower()),
        options=ThermoelasticFitPlotOptions(
            style=style,
            residuals=residuals,
            confidence=confidence if confidence_band else None,
            curve_points=curve_points,
            show_symmetry_spread=symmetry_spread,
        ),
    )
    render_collection(
        archive,
        collection,
        family="fit",
        output_dir=output_dir,
        image_format=image_format,
        preset=preset,
        dpi=dpi,
        transparent=transparent,
        show=show,
        figure_size=(figure_width, figure_height),
        axis_label_font_size=axis_label_font_size,
        legend_font_size=legend_font_size,
        title_font_size=title_font_size,
        tick_label_font_size=tick_label_font_size,
    )


__all__ = ["fit_plot"]
