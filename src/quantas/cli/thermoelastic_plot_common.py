# -*- coding: utf-8 -*-

"""Shared Click options and rendering helpers for thermoelastic plots."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any, TypeVar, cast

import click

from quantas.cli.contracts import OUTPUT_GROUP, figure_preset_option
from quantas.cli.grouped_options import grouped_option
from quantas.cli.messages import quantas_finish, quantas_title
from quantas.cli.output import CLIOutput
from quantas.core.events import EventLevel
from quantas.models import PlotCollection, ReportTable, ResultData
from quantas.api.thermoelasticity import (
    PlotPreset as ThermoelasticPlotPreset,
    PlotStyleOptions as ThermoelasticPlotStyleOptions,
    read_result as read_thermoelastic_hdf5,
    resolve_components,
)
from quantas.api.thermoelasticity import Result as ThermoelasticResult
from quantas.renderers.plots import MatplotlibOptions, render_plot_collection

_F = TypeVar("_F", bound=Callable[..., Any])
_COMPONENT_GROUPS = (
    "independent",
    "normal",
    "shear",
    "coupling",
    "offdiagonal",
    "all",
)


def component_options(function: _F) -> _F:
    """Attach shared elastic-component selection options."""
    decorators = (
        grouped_option(
            "--component",
            "components",
            group="Component selection",
            multiple=True,
            metavar="CIJ",
            help=(
                "Elastic component to plot, for example C11 or C33. Repeat as "
                "needed; explicit components override --component-group."
            ),
        ),
        grouped_option(
            "--component-group",
            group="Component selection",
            type=click.Choice(_COMPONENT_GROUPS, case_sensitive=False),
            default="independent",
            show_default=True,
            help="Semantic component group used when --component is absent.",
        ),
        grouped_option(
            "--list-components",
            group="Component selection",
            is_flag=True,
            default=False,
            help="List components available in the archive and exit.",
        ),
    )
    for decorator in reversed(decorators):
        function = decorator(function)
    return function


def style_options(function: _F) -> _F:
    """Attach common frontend-neutral plot style controls."""
    decorators = (
        figure_preset_option(parameter_name="preset", group="Style preset"),
        grouped_option(
            "--line-width",
            group="Lines and markers",
            hidden=True,
            type=click.FloatRange(min=0.0, min_open=True),
            default=None,
            help="Component line width. Default is supplied by the style preset.",
        ),
        grouped_option(
            "--marker-size",
            group="Lines and markers",
            hidden=True,
            type=click.FloatRange(min=0.0, min_open=True),
            default=None,
            help="Marker size in typographic points. Default: style preset.",
        ),
        grouped_option(
            "--marker-edge-color",
            group="Lines and markers",
            hidden=True,
            default="black",
            show_default=True,
            help=(
                "Marker border color. Use 'none' for borderless markers; colored "
                "thermoelastic plots use black borders by default."
            ),
        ),
        grouped_option(
            "--marker-edge-width",
            group="Lines and markers",
            hidden=True,
            type=click.FloatRange(min=0.0),
            default=0.7,
            show_default=True,
            help="Marker border width in typographic points.",
        ),
        grouped_option(
            "--errorbar-width",
            group="Lines and markers",
            hidden=True,
            type=click.FloatRange(min=0.0, min_open=True),
            default=1.0,
            show_default=True,
            help="Error-bar stroke width.",
        ),
        grouped_option(
            "--errorbar-capsize",
            group="Lines and markers",
            hidden=True,
            type=click.FloatRange(min=0.0),
            default=2.5,
            show_default=True,
            help="Error-bar cap size in typographic points.",
        ),
        grouped_option(
            "--grid/--no-grid",
            group="Axes and annotations",
            default=None,
            help="Override the style-preset Cartesian-grid choice.",
        ),
        grouped_option(
            "--title/--no-title",
            group="Axes and annotations",
            default=None,
            help="Override the style-preset title choice.",
        ),
    )
    for decorator in reversed(decorators):
        function = decorator(function)
    return function


def render_options(
    *,
    default_width: float,
    default_height: float,
) -> Callable[[_F], _F]:
    """Return a decorator attaching shared Matplotlib output controls."""

    def decorate(function: _F) -> _F:
        decorators = (
            grouped_option(
                "-o",
                "--output",
                "output_dir",
                group=OUTPUT_GROUP,
                type=click.Path(file_okay=False, path_type=Path),
                default=None,
                help="Output directory. Default: archive base name + '_plots'.",
            ),
            grouped_option(
                "--format",
                "image_format",
                group=OUTPUT_GROUP,
                type=click.Choice(["png", "pdf", "svg"], case_sensitive=False),
                default="png",
                show_default=True,
                help="Image format.",
            ),
            grouped_option(
                "--dpi",
                group=OUTPUT_GROUP,
                type=click.IntRange(min=1),
                default=None,
                help="Override the raster resolution supplied by the figure preset.",
            ),
            grouped_option(
                "--transparent/--opaque",
                group=OUTPUT_GROUP,
                default=False,
                show_default=True,
                help="Save figures with a transparent or opaque background.",
            ),
            grouped_option(
                "--show",
                group=OUTPUT_GROUP,
                is_flag=True,
                default=False,
                help="Show figures interactively after generation.",
            ),
            grouped_option(
                "--figure-width",
                group="Figure geometry and typography",
                hidden=True,
                type=click.FloatRange(min=0.0, min_open=True),
                default=default_width,
                show_default=True,
                help="Figure width in inches.",
            ),
            grouped_option(
                "--figure-height",
                group="Figure geometry and typography",
                hidden=True,
                type=click.FloatRange(min=0.0, min_open=True),
                default=default_height,
                show_default=True,
                help="Figure height in inches.",
            ),
            grouped_option(
                "--axis-label-font-size",
                group="Figure geometry and typography",
                hidden=True,
                type=click.FloatRange(min=0.0, min_open=True),
                default=None,
                help="Override the value supplied by the figure preset.",
            ),
            grouped_option(
                "--legend-font-size",
                group="Figure geometry and typography",
                hidden=True,
                type=click.FloatRange(min=0.0, min_open=True),
                default=None,
                help="Override the value supplied by the figure preset.",
            ),
            grouped_option(
                "--title-font-size",
                group="Figure geometry and typography",
                hidden=True,
                type=click.FloatRange(min=0.0, min_open=True),
                default=None,
                help="Override the value supplied by the figure preset.",
            ),
            grouped_option(
                "--tick-label-font-size",
                group="Figure geometry and typography",
                hidden=True,
                type=click.FloatRange(min=0.0, min_open=True),
                default=None,
                help="Override the value supplied by the figure preset.",
            ),
        )
        for decorator in reversed(decorators):
            function = decorator(function)
        return function

    return decorate


def read_archive(archive: Path) -> tuple[ResultData, ThermoelasticResult]:
    """Read one thermoelastic HDF5 archive and return its passive payload."""
    try:
        result_data = read_thermoelastic_hdf5(archive)
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc
    payload = result_data.results.get("thermoelasticity")
    if not isinstance(payload, ThermoelasticResult):
        raise click.ClickException(
            "archive does not contain a thermoelasticity result payload"
        )
    return result_data, payload


def make_style(
    *,
    preset: str,
    line_width: float | None,
    marker_size: float | None,
    marker_edge_color: str,
    marker_edge_width: float,
    errorbar_width: float,
    errorbar_capsize: float,
    grid: bool | None,
    title: bool | None,
) -> ThermoelasticPlotStyleOptions:
    """Build validated common plot style options from CLI values."""
    normalized_preset = "analysis" if preset.lower() == "screen" else preset.lower()
    return ThermoelasticPlotStyleOptions(
        preset=cast(ThermoelasticPlotPreset, normalized_preset),
        show_title=title,
        grid=grid,
        line_width=line_width,
        marker_size=marker_size,
        marker_edge_color=marker_edge_color,
        marker_edge_width=marker_edge_width,
        errorbar_width=errorbar_width,
        errorbar_capsize=errorbar_capsize,
    )


def render_collection(
    archive: Path,
    collection: PlotCollection,
    *,
    family: str,
    output_dir: Path | None,
    image_format: str,
    preset: str,
    dpi: int | None,
    transparent: bool,
    show: bool,
    figure_size: tuple[float, float],
    axis_label_font_size: float | None,
    legend_font_size: float | None,
    title_font_size: float | None,
    tick_label_font_size: float | None,
) -> None:
    """Render one prepared plot collection and report generated artifacts."""
    destination = output_dir or archive.with_name(f"{archive.stem}_plots")
    try:
        rendered = render_plot_collection(
            collection,
            MatplotlibOptions.from_preset(
                preset,
                output_dir=destination,
                filename_prefix="",
                image_format=image_format.lower(),
                dpi=dpi,
                figure_size=figure_size,
                tight_layout=True,
                axis_label_font_size=axis_label_font_size,
                legend_font_size=legend_font_size,
                title_font_size=title_font_size,
                tick_label_font_size=tick_label_font_size,
                show=show,
                close=not show,
                savefig_kwargs={"transparent": transparent},
            ),
        )
    except click.ClickException:
        raise
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc

    output = CLIOutput(show_progress=False)
    try:
        output.message(quantas_title(), bold=True)
        for warning in rendered.warnings:
            output.message(str(warning), level=EventLevel.WARNING)
        for artifact in rendered.artifacts:
            if artifact.path is not None:
                output.message(
                    f"Thermoelastic {family} plot written to: {artifact.path}"
                )
        output.message(quantas_finish())
        output.save()
    finally:
        output.close()


def show_components(result: ThermoelasticResult) -> None:
    """Print all non-zero and independent elastic components in one archive."""
    available = resolve_components(result, group="all")
    independent = set(result.independent_labels)
    rows = [[label, "yes" if label in independent else "no"] for label in available]
    output = CLIOutput(show_progress=False)
    try:
        output.table(
            ReportTable(
                "Available thermoelastic components",
                ["Component", "Independent fit"],
                rows,
            )
        )
        output.save()
    finally:
        output.close()


def show_profiles(result: ThermoelasticResult) -> None:
    """Print archived geological profile coverage."""
    rows = [
        [
            name,
            profile.depth.size,
            float(profile.depth[0]),
            float(profile.depth[-1]),
            float(profile.pressure.min()),
            float(profile.pressure.max()),
            float(profile.temperature.min()),
            float(profile.temperature.max()),
        ]
        for name, profile in sorted(result.profiles.items())
    ]
    output = CLIOutput(show_progress=False)
    try:
        output.table(
            ReportTable(
                "Archived thermoelastic profiles",
                [
                    "Profile",
                    "Points",
                    "z min (km)",
                    "z max (km)",
                    "P min (GPa)",
                    "P max (GPa)",
                    "T min (K)",
                    "T max (K)",
                ],
                rows,
                metadata={
                    "notes": [
                        "No profiles are stored in this archive."
                        if not rows
                        else "Use --profile-name to select one path."
                    ]
                },
            )
        )
        output.save()
    finally:
        output.close()


__all__ = [
    "component_options",
    "make_style",
    "read_archive",
    "render_collection",
    "render_options",
    "show_components",
    "show_profiles",
    "style_options",
]
