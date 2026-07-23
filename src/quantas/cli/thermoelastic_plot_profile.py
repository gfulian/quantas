# -*- coding: utf-8 -*-

"""CLI command for thermoelastic geological-profile plots."""

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
    show_profiles,
    style_options,
)
from quantas.api.thermoelasticity import (
    ComponentGroup as ThermoelasticComponentGroup,
    PlotLayout as ThermoelasticPlotLayout,
    ProfileBackground as ThermoelasticProfileBackground,
    ProfileColor as ThermoelasticProfileColor,
    ProfileMode as ThermoelasticProfileMode,
    ProfilePlotOptions as ThermoelasticProfilePlotOptions,
    TensorCondition as ThermoelasticTensorCondition,
    UncertaintyMode as ThermoelasticUncertaintyMode,
    build_profile_plots as build_thermoelastic_profile_plots,
)
from quantas.renderers.plots.matplotlib import validate_colormap_name


@click.command(name="profile", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@grouped_option(
    "--profile-name",
    group="Profile selection",
    default=None,
    metavar="NAME",
    help="Archived profile name. Optional only when exactly one profile exists.",
)
@grouped_option(
    "--list-profiles",
    group="Profile selection",
    is_flag=True,
    default=False,
    help="List archived profiles and exit.",
)
@component_options
@grouped_option(
    "--tensor-condition",
    group="Profile coordinates",
    type=click.Choice(["isothermal", "adiabatic"], case_sensitive=False),
    default="isothermal",
    show_default=True,
    help="Plot the isothermal C^T or adiabatic C^S tensor profile.",
)
@grouped_option(
    "--mode",
    group="Profile coordinates",
    type=click.Choice(["absolute", "relative"], case_sensitive=False),
    default="absolute",
    show_default=True,
    help="Plot absolute stiffness or percentage change from a reference depth.",
)
@grouped_option(
    "--reference-depth",
    group="Profile coordinates",
    type=float,
    default=None,
    help="Reference depth in km for relative profiles. Default: first point.",
)
@grouped_option(
    "--layout",
    group="Profile coordinates",
    type=click.Choice(["auto", "overlay", "facets", "separate"], case_sensitive=False),
    default="auto",
    show_default=True,
    help="Multi-component arrangement.",
)
@grouped_option(
    "--max-overlay-components",
    group="Profile coordinates",
    hidden=True,
    type=click.IntRange(min=1),
    default=6,
    show_default=True,
    help="Maximum number of overlaid components selected by layout=auto.",
)
@grouped_option(
    "--panel-columns",
    group="Profile coordinates",
    hidden=True,
    type=click.IntRange(min=1),
    default=2,
    show_default=True,
    help="Preferred number of columns in faceted figures.",
)
@grouped_option(
    "--uncertainty",
    group="Uncertainty presentation",
    type=click.Choice(["auto", "none", "band", "bars"], case_sensitive=False),
    default="auto",
    show_default=True,
    help=(
        "One-sigma uncertainty style. Auto uses a band for one component and "
        "omits uncertainty when several Cij are requested together."
    ),
)
@grouped_option(
    "--color-by",
    group="Temperature encoding",
    type=click.Choice(["component", "temperature", "none"], case_sensitive=False),
    default="temperature",
    show_default=True,
    help="Variable controlling profile-line colors.",
)
@grouped_option(
    "--background",
    group="Temperature encoding",
    type=click.Choice(["none", "temperature"], case_sensitive=False),
    default="none",
    show_default=True,
    help="Optional scalar background behind the depth profile.",
)
@grouped_option(
    "--temperature-colormap",
    group="Temperature encoding",
    default="plasma",
    show_default=True,
    help="Colormap used for temperature-colored paths or backgrounds.",
)
@grouped_option(
    "--pressure-axis/--no-pressure-axis",
    group="Axes and annotations",
    default=True,
    show_default=True,
    help="Show a right-hand pressure axis when P(depth) is monotonic.",
)
@grouped_option(
    "--extrapolation/--no-extrapolation",
    group="Axes and annotations",
    default=True,
    show_default=True,
    help="Mark depth intervals requiring QHA or elastic extrapolation.",
)
@style_options
@render_options(default_width=8.0, default_height=7.0)
def profile_plot(
    archive: Path,
    profile_name: str | None,
    list_profiles: bool,
    components: tuple[str, ...],
    component_group: str,
    list_components: bool,
    tensor_condition: str,
    mode: str,
    reference_depth: float | None,
    layout: str,
    max_overlay_components: int,
    panel_columns: int,
    uncertainty: str,
    color_by: str,
    background: str,
    temperature_colormap: str,
    pressure_axis: bool,
    extrapolation: bool,
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
    """Plot elastic components along one archived geological profile."""
    result_data, payload = read_archive(archive)
    if list_profiles:
        show_profiles(payload)
        return
    if list_components:
        show_components(payload)
        return
    validate_colormap_name(temperature_colormap)
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
    collection = build_thermoelastic_profile_plots(
        result_data,
        profile_name=profile_name,
        components=components or None,
        component_group=cast(ThermoelasticComponentGroup, component_group.lower()),
        options=ThermoelasticProfilePlotOptions(
            style=style,
            tensor_condition=cast(
                ThermoelasticTensorCondition, tensor_condition.lower()
            ),
            mode=cast(ThermoelasticProfileMode, mode.lower()),
            layout=cast(ThermoelasticPlotLayout, layout.lower()),
            uncertainty=cast(ThermoelasticUncertaintyMode, uncertainty.lower()),
            color_by=cast(ThermoelasticProfileColor, color_by.lower()),
            background=cast(ThermoelasticProfileBackground, background.lower()),
            reference_depth=reference_depth,
            show_pressure_axis=pressure_axis,
            show_extrapolation=extrapolation,
            max_overlay_components=max_overlay_components,
            panel_columns=panel_columns,
            temperature_colormap=temperature_colormap,
        ),
    )
    render_collection(
        archive,
        collection,
        family="profile",
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


__all__ = ["profile_plot"]
