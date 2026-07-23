# -*- coding: utf-8 -*-

"""CLI command for thermoelastic pressure-temperature contour plots."""

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
    PTPlotOptions as ThermoelasticPTPlotOptions,
    PTQuantity as ThermoelasticPTQuantity,
    PlotLayout as ThermoelasticPlotLayout,
    TensorCondition as ThermoelasticTensorCondition,
    build_pt_plots as build_thermoelastic_pt_plots,
)
from quantas.renderers.plots.matplotlib import validate_colormap_name


@click.command(name="pt", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@component_options
@grouped_option(
    "--tensor-condition",
    group="Mapped quantity",
    type=click.Choice(["isothermal", "adiabatic"], case_sensitive=False),
    default="isothermal",
    show_default=True,
    help="Plot the isothermal C^T or adiabatic C^S tensor field.",
)
@grouped_option(
    "--quantity",
    group="Mapped quantity",
    type=click.Choice(
        ["value", "uncertainty", "relative-uncertainty"],
        case_sensitive=False,
    ),
    default="value",
    show_default=True,
    help="Quantity represented by the contour colors.",
)
@grouped_option(
    "--layout",
    group="Mapped quantity",
    type=click.Choice(["auto", "facets", "separate"], case_sensitive=False),
    default="auto",
    show_default=True,
    help="Multi-component arrangement.",
)
@grouped_option(
    "--colormap",
    group="Contours",
    default=None,
    help="Matplotlib colormap. Default: automatic sequential/diverging choice.",
)
@grouped_option(
    "--levels",
    group="Contours",
    type=click.IntRange(min=2),
    default=12,
    show_default=True,
    help="Number of main contour levels.",
)
@grouped_option(
    "--isolines/--no-isolines",
    group="Contours",
    default=True,
    show_default=True,
    help="Overlay contour lines.",
)
@grouped_option(
    "--isoline-labels/--no-isoline-labels",
    group="Contours",
    default=True,
    show_default=True,
    help="Label contour lines when isolines are enabled.",
)
@grouped_option(
    "--extrapolation/--no-extrapolation",
    group="Domain diagnostics",
    default=True,
    show_default=True,
    help="Overlay QHA-coordinate and elastic-volume extrapolation masks.",
)
@grouped_option(
    "--panel-columns",
    group="Mapped quantity",
    hidden=True,
    type=click.IntRange(min=1),
    default=2,
    show_default=True,
    help="Preferred number of columns in faceted figures.",
)
@style_options
@render_options(default_width=7.0, default_height=5.5)
def pt_plot(
    archive: Path,
    components: tuple[str, ...],
    component_group: str,
    list_components: bool,
    tensor_condition: str,
    quantity: str,
    layout: str,
    colormap: str | None,
    levels: int,
    isolines: bool,
    isoline_labels: bool,
    extrapolation: bool,
    panel_columns: int,
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
    """Plot thermoelastic stiffness or uncertainty on the archived P-T grid."""
    result_data, payload = read_archive(archive)
    if list_components:
        show_components(payload)
        return
    if colormap is not None:
        validate_colormap_name(colormap)
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
    collection = build_thermoelastic_pt_plots(
        result_data,
        components=components or None,
        component_group=cast(ThermoelasticComponentGroup, component_group.lower()),
        options=ThermoelasticPTPlotOptions(
            style=style,
            tensor_condition=cast(
                ThermoelasticTensorCondition, tensor_condition.lower()
            ),
            quantity=cast(ThermoelasticPTQuantity, quantity.lower()),
            layout=cast(ThermoelasticPlotLayout, layout.lower()),
            colormap=colormap,
            levels=levels,
            isolines=isolines,
            isoline_labels=isoline_labels,
            show_extrapolation=extrapolation,
            panel_columns=panel_columns,
        ),
    )
    render_collection(
        archive,
        collection,
        family="pt",
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


__all__ = ["pt_plot"]
