# -*- coding: utf-8 -*-

"""CLI command comparing isothermal and adiabatic thermoelastic stiffness."""

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
    ComparePlotOptions as ThermoelasticComparePlotOptions,
    ComponentGroup as ThermoelasticComponentGroup,
    PlotLayout as ThermoelasticPlotLayout,
    build_compare_plots as build_thermoelastic_compare_plots,
)


@click.command(name="compare", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@component_options
@grouped_option(
    "--pressure",
    group="Fixed thermodynamic coordinate",
    type=float,
    default=None,
    help="Hold pressure fixed in GPa and vary temperature.",
)
@grouped_option(
    "--temperature",
    group="Fixed thermodynamic coordinate",
    type=float,
    default=None,
    help="Hold temperature fixed in K and vary pressure.",
)
@grouped_option(
    "--layout",
    group="Component arrangement",
    type=click.Choice(["auto", "facets", "separate"], case_sensitive=False),
    default="auto",
    show_default=True,
    help=(
        "Arrange multiple components automatically, in facets, or as "
        "separate figures."
    ),
)
@style_options
@render_options(default_width=7.0, default_height=6.5)
def compare_plot(
    archive: Path,
    components: tuple[str, ...],
    component_group: str,
    list_components: bool,
    pressure: float | None,
    temperature: float | None,
    layout: str,
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
    """Compare C^T and C^S at fixed pressure or fixed temperature."""
    result_data, payload = read_archive(archive)
    if list_components:
        show_components(payload)
        return
    if (pressure is None) == (temperature is None):
        raise click.UsageError("set exactly one of --pressure or --temperature")
    collection = build_thermoelastic_compare_plots(
        result_data,
        components=components or None,
        component_group=cast(ThermoelasticComponentGroup, component_group.lower()),
        options=ThermoelasticComparePlotOptions(
            style=make_style(
                preset=preset,
                line_width=line_width,
                marker_size=marker_size,
                marker_edge_color=marker_edge_color,
                marker_edge_width=marker_edge_width,
                errorbar_width=errorbar_width,
                errorbar_capsize=errorbar_capsize,
                grid=grid,
                title=title,
            ),
            fixed_pressure=pressure,
            fixed_temperature=temperature,
            layout=cast(ThermoelasticPlotLayout, layout.lower()),
        ),
    )
    render_collection(
        archive,
        collection,
        family="compare",
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


__all__ = ["compare_plot"]
