# -*- coding: utf-8 -*-

"""EOS post-fit diagnostic, plotting, and calculation commands."""

from __future__ import annotations

from pathlib import Path

import click
from quantas.cli.contracts import OUTPUT_GROUP, figure_preset_option
from quantas.cli.eos_common import GroupedCommand, grouped_option
from quantas.cli.eos_helpers import (
    _combine_pvt_coordinates,
    _combined_grid,
)
from quantas.cli.messages import quantas_finish, quantas_title
from quantas.core.events import EventLevel
from quantas.cli.output import CLIOutput
from quantas.api.eos import (
    PlotOptions as EOSPlotOptions,
    available_plot_types,
    build_plots as build_eos_plots,
    calculate as calculate_eos,
    diagnose as diagnose_eos,
    calculation_summary_table,
    calculation_table,
    diagnostic_summary_table,
    diagnostic_table,
    record_domain,
    write_calculation_csv as write_eos_calculation_csv,
    write_diagnostic_csv as write_eos_diagnostic_csv,
)
from quantas.renderers.plots import MatplotlibOptions, render_plot_collection


@click.command(name="diagnose", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@grouped_option(
    "--slot",
    group="Record selection",
    default=None,
    metavar="DOMAIN/TARGET",
    help="Accepted result slot. Optional only when exactly one accepted result exists.",
)
@grouped_option(
    "--record-id",
    group="Record selection",
    type=click.IntRange(min=1),
    default=None,
    help="Inspect an explicit immutable record instead of the accepted result.",
)
@grouped_option(
    "--no-normalized-pressure",
    group="Diagnostics",
    is_flag=True,
    default=False,
    help="Suppress finite-strain and normalized-pressure columns.",
)
@grouped_option(
    "-o",
    "--output",
    group=OUTPUT_GROUP,
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Optional CSV destination for all diagnostic rows.",
)
@grouped_option(
    "-f",
    "--force",
    group=OUTPUT_GROUP,
    is_flag=True,
    default=False,
    help="Replace an existing CSV file. The command never prompts.",
)
def diagnose(
    archive: Path,
    slot: str | None,
    record_id: int | None,
    no_normalized_pressure: bool,
    output: Path | None,
    force: bool,
) -> None:
    """Inspect residuals and finite-strain diagnostics from an EOS archive."""
    terminal = CLIOutput(show_progress=False)
    try:
        terminal.message(quantas_title(), bold=True)
        diagnostic = diagnose_eos(
            archive,
            slot=slot,
            record_id=record_id,
            include_normalized_pressure=not no_normalized_pressure,
        )
        terminal.tables(
            (
                diagnostic_summary_table(diagnostic),
                diagnostic_table(diagnostic),
            )
        )
        if output is not None:
            destination = write_eos_diagnostic_csv(diagnostic, output, overwrite=force)
            terminal.message(f"EOS diagnostic CSV written to: {destination}")
        terminal.message(quantas_finish())
        terminal.save()
    except FileExistsError as exc:
        raise click.UsageError(
            f"output file exists: {exc}; use --force to replace it"
        ) from exc
    except click.ClickException:
        raise
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc
    finally:
        terminal.close()


_EOS_PLOT_DESCRIPTIONS = {
    "fit": "Observed data and fitted P-V or V-T curve.",
    "residuals": "Physical residuals against pressure or temperature.",
    "standardized-residuals": "Dimensionless standardized residuals.",
    "normalized-pressure": "Finite-strain normalized-pressure diagnostic.",
    "coverage": "Pressure-temperature coverage of a P-V-T dataset.",
    "isotherms": "Calculated P-V isotherms with P-V-T observations.",
    "isobars": "Calculated V-T isobars with P-V-T observations.",
}


@click.command(name="plot", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@grouped_option(
    "--slot",
    group="Record selection",
    default=None,
    metavar="DOMAIN/TARGET",
    help="Accepted result slot. Optional only when exactly one accepted result exists.",
)
@grouped_option(
    "--record-id",
    group="Record selection",
    type=click.IntRange(min=1),
    default=None,
    help="Plot an explicit immutable record instead of the accepted result.",
)
@grouped_option(
    "--plot",
    "plot_types",
    group="Plot selection",
    type=click.Choice(
        [
            "all",
            "fit",
            "residuals",
            "standardized-residuals",
            "normalized-pressure",
            "coverage",
            "isotherms",
            "isobars",
        ],
        case_sensitive=False,
    ),
    multiple=True,
    help="Plot type. Repeat for several types. Default: all available plots.",
)
@grouped_option(
    "--list-plots",
    group="Plot selection",
    is_flag=True,
    default=False,
    help="List plot types available for the selected record and exit.",
)
@figure_preset_option(group="Style")
@grouped_option(
    "--uncertainties/--no-uncertainties",
    group="Data presentation",
    default=True,
    show_default=True,
    help="Draw input one-sigma error bars when available.",
)
@grouped_option(
    "--excluded/--no-excluded",
    group="Data presentation",
    default=True,
    show_default=True,
    help="Display observations excluded from the fit.",
)
@grouped_option(
    "--group-data/--no-group-data",
    group="Data presentation",
    default=True,
    show_default=True,
    help="Use separate data series for input groups.",
)
@grouped_option(
    "--point-size",
    group="Style",
    type=click.FloatRange(min=0.0, min_open=True),
    default=5.0,
    show_default=True,
    help="Observation marker size in typographic points.",
)
@grouped_option(
    "--curve-width",
    group="Style",
    type=click.FloatRange(min=0.0, min_open=True),
    default=1.8,
    show_default=True,
    help="Fitted and calculated curve width.",
)
@grouped_option(
    "--errorbar-width",
    group="Style",
    type=click.FloatRange(min=0.0, min_open=True),
    default=1.0,
    show_default=True,
    help="Error-bar stroke width.",
)
@grouped_option(
    "--errorbar-capsize",
    group="Style",
    type=click.FloatRange(min=0.0),
    default=2.5,
    show_default=True,
    help="Error-bar cap size in typographic points.",
)
@grouped_option(
    "--curve-points",
    group="Style",
    type=click.IntRange(min=20),
    default=300,
    show_default=True,
    help="Number of samples used for smooth calculated curves.",
)
@grouped_option(
    "--curve-color",
    group="Style",
    default="black",
    show_default=True,
    help="Matplotlib-compatible color used for the primary fitted curve.",
)
@grouped_option(
    "--excluded-color",
    group="Style",
    default="0.55",
    show_default=True,
    help="Matplotlib-compatible color used for excluded observations.",
)
@grouped_option(
    "--grid/--no-grid",
    group="Style",
    default=True,
    show_default=True,
    help="Display Cartesian grids.",
)
@grouped_option(
    "--legend/--no-legend",
    group="Style",
    default=True,
    show_default=True,
    help="Display legends when a plot contains several series.",
)
@grouped_option(
    "--title/--no-title",
    group="Style",
    default=True,
    show_default=True,
    help="Display figure titles. Disable for publication figures using captions.",
)
@grouped_option(
    "--zero-pressure-point/--no-zero-pressure-point",
    group="Data presentation",
    default=False,
    show_default=True,
    help=(
        "Display finite zero-pressure observations in normalized-pressure plots. "
        "They are hidden by default to focus the diagnostic strain range."
    ),
)
@grouped_option(
    "--axis-label-font-size",
    group="Style",
    type=click.FloatRange(min=0.0, min_open=True),
    default=None,
    help="Axis-label font size. Default: figure preset.",
)
@grouped_option(
    "--legend-font-size",
    group="Style",
    type=click.FloatRange(min=0.0, min_open=True),
    default=None,
    help="Legend font size. Default: figure preset.",
)
@grouped_option(
    "--title-font-size",
    group="Style",
    type=click.FloatRange(min=0.0, min_open=True),
    default=None,
    help="Title font size. Default: figure preset.",
)
@grouped_option(
    "--tick-label-font-size",
    group="Style",
    type=click.FloatRange(min=0.0, min_open=True),
    default=None,
    help="Font size for axis tick labels. Default: Matplotlib setting.",
)
@grouped_option(
    "--figure-width",
    group="Style",
    type=click.FloatRange(min=0.0, min_open=True),
    default=None,
    help="Figure width in inches. Default: plot-specific preset geometry.",
)
@grouped_option(
    "--figure-height",
    group="Style",
    type=click.FloatRange(min=0.0, min_open=True),
    default=None,
    help="Figure height in inches. Default: plot-specific preset geometry.",
)
@grouped_option(
    "--isotherm",
    "isotherms",
    group="P-V-T curves",
    type=float,
    multiple=True,
    help="Temperature in kelvin for a P-V-T isotherm. Repeat as needed.",
)
@grouped_option(
    "--isobar",
    "isobars",
    group="P-V-T curves",
    type=float,
    multiple=True,
    help="Pressure for a P-V-T isobar. Repeat as needed.",
)
@grouped_option(
    "-o",
    "--output",
    "output_dir",
    group=OUTPUT_GROUP,
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Output directory. Default: archive base name + '_plots'.",
)
@grouped_option(
    "--format",
    "image_format",
    group=OUTPUT_GROUP,
    type=click.Choice(["png", "pdf", "svg"], case_sensitive=False),
    default="png",
    show_default=True,
    help="Output image format.",
)
@grouped_option(
    "--dpi",
    group=OUTPUT_GROUP,
    type=click.IntRange(min=1),
    default=None,
    help="Override the raster resolution supplied by the figure preset.",
)
@grouped_option(
    "--show",
    group=OUTPUT_GROUP,
    is_flag=True,
    default=False,
    help="Show figures on screen after generation.",
)
def plot(
    archive: Path,
    slot: str | None,
    record_id: int | None,
    plot_types: tuple[str, ...],
    list_plots: bool,
    figure_preset: str,
    uncertainties: bool,
    excluded: bool,
    group_data: bool,
    point_size: float,
    curve_width: float,
    errorbar_width: float,
    errorbar_capsize: float,
    curve_points: int,
    curve_color: str,
    excluded_color: str,
    grid: bool,
    legend: bool,
    title: bool,
    zero_pressure_point: bool,
    axis_label_font_size: float | None,
    legend_font_size: float | None,
    title_font_size: float | None,
    tick_label_font_size: float | None,
    figure_width: float | None,
    figure_height: float | None,
    isotherms: tuple[float, ...],
    isobars: tuple[float, ...],
    output_dir: Path | None,
    image_format: str,
    dpi: int | None,
    show: bool,
) -> None:
    """Generate EOS plots from an immutable HDF5 fit record."""
    if (figure_width is None) != (figure_height is None):
        raise click.UsageError(
            "--figure-width and --figure-height must be supplied together"
        )
    figure_size = (
        None
        if figure_width is None or figure_height is None
        else (figure_width, figure_height)
    )
    try:
        available = available_plot_types(
            archive, slot=slot, record_id=record_id
        )
        if list_plots:
            for kind in available:
                click.echo(f"{kind:28s} {_EOS_PLOT_DESCRIPTIONS[kind]}")
            return
        if output_dir is None:
            stem = archive.with_suffix("")
            output_dir = stem.with_name(f"{stem.name}_plots")
        preparation = EOSPlotOptions(
            show_uncertainties=uncertainties,
            show_excluded=excluded,
            group_data=group_data,
            point_size=point_size,
            curve_width=curve_width,
            errorbar_width=errorbar_width,
            errorbar_capsize=errorbar_capsize,
            curve_points=curve_points,
            curve_color=curve_color,
            excluded_color=excluded_color,
            show_legend=legend,
            show_title=title,
            show_zero_pressure_point=zero_pressure_point,
            grid=grid,
            isotherms=isotherms,
            isobars=isobars,
        )
        collection = build_eos_plots(
            archive,
            plot_types or None,
            slot=slot,
            record_id=record_id,
            options=preparation,
        )
        render_options = MatplotlibOptions.from_preset(
            figure_preset,
            output_dir=output_dir,
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
        )
        result = render_plot_collection(collection, render_options)
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc
    terminal = CLIOutput(show_progress=False)
    try:
        terminal.message(quantas_title(), bold=True)
        for warning in result.warnings:
            terminal.message(str(warning), level=EventLevel.WARNING)
        for artifact in result.artifacts:
            if artifact.path is not None:
                terminal.message(f"EOS plot written to: {artifact.path}")
        terminal.message(quantas_finish())
        terminal.save()
    finally:
        terminal.close()


@click.command(name="calculate", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@grouped_option(
    "--slot",
    group="Record selection",
    default=None,
    metavar="DOMAIN/TARGET",
    help="Accepted result slot. Optional only when exactly one accepted result exists.",
)
@grouped_option(
    "--record-id",
    group="Record selection",
    type=click.IntRange(min=1),
    default=None,
    help="Evaluate an explicit immutable record instead of the accepted result.",
)
@grouped_option(
    "--pressure",
    group="State coordinates",
    type=float,
    multiple=True,
    help="Pressure value. Repeat for several values.",
)
@grouped_option(
    "--pressure-range",
    group="State coordinates",
    default=None,
    metavar="START:STOP:STEP",
    help="Inclusive pressure grid.",
)
@grouped_option(
    "--coordinate",
    group="State coordinates",
    type=float,
    multiple=True,
    help="Volume or fitted-axis value. Repeat for several values.",
)
@grouped_option(
    "--coordinate-range",
    group="State coordinates",
    default=None,
    metavar="START:STOP:STEP",
    help="Inclusive volume or fitted-axis grid.",
)
@grouped_option(
    "--temperature",
    group="State coordinates",
    type=float,
    multiple=True,
    help="Temperature value in kelvin. Repeat for several values.",
)
@grouped_option(
    "--temperature-range",
    group="State coordinates",
    default=None,
    metavar="START:STOP:STEP",
    help="Inclusive temperature grid in kelvin.",
)
@grouped_option(
    "--pairwise",
    group="State coordinates",
    is_flag=True,
    default=False,
    help="Pair P/V and T values element-wise instead of constructing a Cartesian grid.",
)
@grouped_option(
    "--no-uncertainty",
    group="Uncertainty propagation",
    is_flag=True,
    default=False,
    help="Do not propagate the fitted parameter covariance.",
)
@grouped_option(
    "--relative-step",
    group="Uncertainty propagation",
    type=click.FloatRange(min=0.0, min_open=True),
    default=1.0e-5,
    show_default=True,
    help="Relative finite-difference step for parameter-covariance propagation.",
)
@grouped_option(
    "-o",
    "--output",
    group=OUTPUT_GROUP,
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Optional CSV destination for all calculated states.",
)
@grouped_option(
    "-f",
    "--force",
    group=OUTPUT_GROUP,
    is_flag=True,
    default=False,
    help="Replace an existing CSV file. The command never prompts.",
)
def calculate(
    archive: Path,
    slot: str | None,
    record_id: int | None,
    pressure: tuple[float, ...],
    pressure_range: str | None,
    coordinate: tuple[float, ...],
    coordinate_range: str | None,
    temperature: tuple[float, ...],
    temperature_range: str | None,
    pairwise: bool,
    no_uncertainty: bool,
    relative_step: float,
    output: Path | None,
    force: bool,
) -> None:
    """Calculate fitted EOS properties at selected states or regular grids."""
    terminal = CLIOutput(show_progress=False)
    try:
        terminal.message(quantas_title(), bold=True)
        domain = record_domain(archive, slot=slot, record_id=record_id)
        p_values = _combined_grid(pressure, pressure_range, "pressure")
        x_values = _combined_grid(coordinate, coordinate_range, "coordinate")
        t_values = _combined_grid(temperature, temperature_range, "temperature")
        if p_values is not None and x_values is not None:
            raise click.UsageError("use pressure or coordinate values, not both")
        if domain.value == "pvt":
            if t_values is None:
                raise click.UsageError("P-V-T calculation requires temperature values")
            primary = p_values if p_values is not None else x_values
            if primary is None:
                raise click.UsageError(
                    "P-V-T calculation requires pressure or coordinate values"
                )
            primary, t_values = _combine_pvt_coordinates(
                primary, t_values, pairwise=pairwise
            )
            p_values = primary if p_values is not None else None
            x_values = primary if x_values is not None else None
        calculation = calculate_eos(
            archive,
            slot=slot,
            record_id=record_id,
            pressure=p_values,
            volume=x_values,
            temperature=t_values,
            propagate_uncertainty=not no_uncertainty,
            relative_step=relative_step,
        )
        terminal.tables(
            (
                calculation_summary_table(calculation),
                calculation_table(calculation),
            )
        )
        if output is not None:
            destination = write_eos_calculation_csv(
                calculation, output, overwrite=force
            )
            terminal.message(f"EOS calculation CSV written to: {destination}")
        terminal.message(quantas_finish())
        terminal.save()
    except FileExistsError as exc:
        raise click.UsageError(
            f"output file exists: {exc}; use --force to replace it"
        ) from exc
    except click.ClickException:
        raise
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc
    finally:
        terminal.close()


__all__ = ["calculate", "diagnose", "plot"]
