# -*- coding: utf-8 -*-

"""CLI command for thermoelastic pressure-temperature domain plots."""

from __future__ import annotations

from pathlib import Path

import click

from quantas.cli.grouped_options import GroupedCommand, grouped_option
from quantas.cli.thermoelastic_plot_common import (
    make_style,
    read_archive,
    render_collection,
    render_options,
    show_profiles,
    style_options,
)
from quantas.api.thermoelasticity import (
    DomainPlotOptions as ThermoelasticDomainPlotOptions,
    build_domain_plot as build_thermoelastic_domain_plot,
)
from quantas.renderers.plots.matplotlib import validate_colormap_name


@click.command(name="domain", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@grouped_option(
    "--profile-name",
    "profile_names",
    group="Profile selection",
    multiple=True,
    metavar="NAME",
    help="Archived profile to overlay. Repeat as needed. Default: all profiles.",
)
@grouped_option(
    "--list-profiles",
    group="Profile selection",
    is_flag=True,
    default=False,
    help="List archived profiles and exit.",
)
@grouped_option(
    "--profiles/--no-profiles",
    group="Profile selection",
    default=True,
    show_default=True,
    help="Overlay archived geological paths.",
)
@grouped_option(
    "--color-profiles-by-depth/--uniform-profile-color",
    group="Geological-path presentation",
    default=True,
    show_default=True,
    help="Color geological paths by depth or use ordinary line colors.",
)
@grouped_option(
    "--profile-colormap",
    group="Geological-path presentation",
    default="plasma",
    show_default=True,
    help="Depth colormap for geological paths, independent of the volume map.",
)
@grouped_option(
    "--colormap",
    group="Equilibrium-volume field",
    default="viridis",
    show_default=True,
    help="Colormap used for the equilibrium-volume background.",
)
@grouped_option(
    "--levels",
    group="Equilibrium-volume field",
    type=click.IntRange(min=2),
    default=12,
    show_default=True,
    help="Number of main volume contour levels.",
)
@grouped_option(
    "--extrapolation/--no-extrapolation",
    group="Domain diagnostics",
    default=True,
    show_default=True,
    help="Overlay QHA-coordinate and elastic-volume extrapolation masks.",
)
@style_options
@render_options(default_width=7.0, default_height=5.8)
def domain_plot(
    archive: Path,
    profile_names: tuple[str, ...],
    list_profiles: bool,
    profiles: bool,
    color_profiles_by_depth: bool,
    profile_colormap: str,
    colormap: str,
    levels: int,
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
    """Plot QHA and elastic validity domains with geological P-T paths."""
    result_data, payload = read_archive(archive)
    if list_profiles:
        show_profiles(payload)
        return
    validate_colormap_name(colormap)
    validate_colormap_name(profile_colormap)
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
    collection = build_thermoelastic_domain_plot(
        result_data,
        profile_names=profile_names or None,
        options=ThermoelasticDomainPlotOptions(
            style=style,
            colormap=colormap,
            levels=levels,
            show_profiles=profiles,
            color_profiles_by_depth=color_profiles_by_depth,
            profile_colormap=profile_colormap,
            show_extrapolation=extrapolation,
        ),
    )
    render_collection(
        archive,
        collection,
        family="domain",
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


__all__ = ["domain_plot"]
