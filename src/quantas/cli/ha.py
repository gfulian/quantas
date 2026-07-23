# -*- coding: utf-8 -*-

"""
Command-line interface for harmonic-approximation workflows.

The commands defined here are thin frontend wrappers around the public HA API.
They do not implement numerical logic directly.
"""

from __future__ import annotations

from pathlib import Path

import click

from quantas.cli.reference_help import apply_reference_help

from quantas.cli.contracts import (
    DOMAIN_GROUP,
    PLOTTING_GROUP,
    UNITS_GROUP,
    VALIDATION_GROUP,
    default_hdf5_path,
    default_report_path,
    figure_preset_option,
    force_option,
    output_option,
    parse_verbosity,
    progress_option,
    quiet_option,
    report_option,
    verbosity_option,
)
from quantas.cli.grouped_options import GroupedCommand, grouped_option
from quantas.cli.messages import (
    confirm,
    echo,
    echo_error,
    echo_highlight,
    echo_warning,
    quantas_error,
    quantas_finish,
    quantas_title,
)
from quantas.api.ha import (
    Options as HAOptions,
    Result as HAResult,
    build_plots as build_ha_plots,
    read_result as read_ha_hdf5,
    run as run_ha,
    write_result as write_ha_hdf5,
)
from quantas.cli.ha_observer import HATextObserver
from quantas.cli.phonon_input import phonon_inpgen
from quantas.modules.ha.io.export import HATableExport
from quantas.references import (
    module_citation_keys,
    render_citation_notice,
)
from quantas.renderers.plots import MatplotlibOptions, render_plot_collection


@click.group(name="ha")
def ha() -> None:
    """Harmonic-approximation thermodynamic analysis."""
    return


ha.add_command(phonon_inpgen)


@ha.command(name="run", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@grouped_option(
    "-T",
    "--temperature",
    group=DOMAIN_GROUP,
    nargs=3,
    type=float,
    metavar="MIN MAX STEP",
    default=(298.15, 298.15, 1.0),
    show_default="298.15 298.15 1.0",
    help="Temperature range.",
)
@grouped_option(
    "--eunit",
    group=UNITS_GROUP,
    type=click.Choice(["Ha", "eV", "Ry"], case_sensitive=True),
    default="Ha",
    show_default=True,
    help="Measurement unit for energy values.",
)
@grouped_option(
    "--vunit",
    group=UNITS_GROUP,
    type=click.Choice(["A", "bohr"], case_sensitive=True),
    default="A",
    show_default=True,
    help="Measurement unit for volume values.",
)
@grouped_option(
    "--funit",
    group=UNITS_GROUP,
    type=click.Choice(["cm-1", "cm^-1", "THz", "Hz"], case_sensitive=True),
    default="cm-1",
    show_default="cm^-1",
    help="Measurement unit for phonon frequency values.",
)
@grouped_option(
    "--tunit",
    group=UNITS_GROUP,
    type=click.Choice(["K", "C"], case_sensitive=True),
    default="K",
    show_default=True,
    help="Measurement unit for temperature values.",
)
@grouped_option(
    "-B",
    "--benchmark",
    group=VALIDATION_GROUP,
    is_flag=True,
    default=False,
    help="Render backend timing events during the run.",
)
@grouped_option(
    "-p",
    "--plot",
    group=PLOTTING_GROUP,
    is_flag=True,
    default=False,
    help="Create static Matplotlib plots after the calculation.",
)
@grouped_option(
    "--plot-property",
    "plot_property",
    group=PLOTTING_GROUP,
    default=None,
    help=(
        "HA property to plot after the run, for example F, Cv, S, Utot, "
        "or all. Default: standard compact set."
    ),
)
@grouped_option(
    "--plot-unit",
    group=PLOTTING_GROUP,
    default=None,
    help="Energy unit used for plotted energy-like values. Defaults to stored result units.",
)
@figure_preset_option(
    option_name="--plot-preset",
    parameter_name="plot_preset",
    group=PLOTTING_GROUP,
)
@grouped_option(
    "--dpi",
    group=PLOTTING_GROUP,
    type=click.IntRange(min=1),
    default=None,
    help="Override the raster resolution supplied by the figure preset.",
)
@output_option(
    help="Output HDF5 file. Default: input file base name + '_HA.hdf5'.",
)
@force_option()
@report_option()
@verbosity_option()
@quiet_option()
@progress_option()
def run(
    filename: Path,
    output: Path | None,
    temperature: tuple[float, float, float],
    eunit: str,
    vunit: str,
    funit: str,
    tunit: str,
    benchmark: bool,
    plot: bool,
    plot_property: str | None,
    plot_unit: str | None,
    plot_preset: str,
    dpi: int | None,
    force: bool,
    report: Path | None,
    verbosity: str,
    quiet: bool,
    progress: bool,
) -> None:
    """
    Run a harmonic-approximation calculation from a Quantas YAML input file.
    """
    destination = default_hdf5_path(filename, output, suffix="_HA")
    report = default_report_path(filename, report)
    report_verbosity = parse_verbosity(verbosity)

    echo_highlight(quantas_title(), silent=quiet)

    options = HAOptions(
        temperature_min=temperature[0],
        temperature_max=temperature[1],
        temperature_step=temperature[2],
        energy_unit=eunit,
        volume_unit=vunit,
        frequency_unit="cm^-1" if funit == "cm-1" else funit,
        temperature_unit=tunit,
    )
    observer = HATextObserver(
        report_file=report,
        silent=quiet,
        show_progress=progress,
        verbosity=report_verbosity,
        include_timing=benchmark,
    )

    try:
        result = run_ha(filename, options=options, observer=observer)
    except Exception as exc:
        observer.close()
        echo_error(quantas_error(), bold=True)
        echo_error(str(exc))
        raise click.Abort() from exc

    observer.output.text_block(render_citation_notice(module_citation_keys("ha")))
    observer.save()

    overwrite = not destination.exists() or force
    if destination.exists() and not force:
        overwrite = confirm(
            f"Output file {destination} exists. Would you like to overwrite it?",
            default=False,
        )
        if not overwrite:
            echo("Results not saved", silent=quiet)

    if overwrite:
        write_ha_hdf5(result, destination, report_text=observer.text())

    if plot:
        basename = destination.with_suffix("")
        collection = build_ha_plots(
            result,
            properties=plot_property,
            unit=plot_unit,
        )
        rendered = render_plot_collection(
            collection,
            MatplotlibOptions.from_preset(
                plot_preset,
                output_dir=basename.parent,
                filename_prefix=f"{basename.name}_",
                dpi=dpi,
                close=True,
            ),
        )
        for warning in rendered.warnings:
            echo_warning(str(warning))
        for artifact in rendered.artifacts:
            if artifact.path is not None:
                echo(f"Plot saved to: {artifact.path}")

    echo_highlight(quantas_finish(), silent=quiet)


@ha.command(name="export", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "-o",
    "--output",
    "outfile",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Output table file. Default: input file base name + '.dat'.",
)
@click.option(
    "--property",
    "property_name",
    default="F",
    show_default=True,
    help="HA property to export, for example F, Cv, S, Utot, or free_energy.",
)
@click.option(
    "--unit",
    default=None,
    help="Output unit for energy-like data, for example Ha, eV, Ry, or kJ/mol.",
)
@click.option(
    "--ask-unit",
    is_flag=True,
    default=False,
    help="Prompt for the output unit during export.",
)
def export(
    filename: Path,
    outfile: Path | None,
    property_name: str,
    unit: str | None,
    ask_unit: bool,
) -> None:
    """
    Export a selected HA property from a Quantas HDF5 result file.
    """
    if outfile is None:
        outfile = filename.with_suffix(".dat")

    try:
        result = read_ha_hdf5(filename)
        ha_result = result.results["ha"]
        if ask_unit:
            unit = click.prompt(
                "Output unit",
                default=unit or _stored_property_unit(ha_result, property_name),
                show_default=True,
            )
        HATableExport().export(
            result,
            outfile,
            property_name=property_name,
            unit=unit,
        )
    except Exception as exc:
        echo_error(str(exc))
        raise click.Abort() from exc

    echo(f"Table saved to: {outfile}")


@ha.command(name="plot", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "-o",
    "--output",
    "outbase",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Output figure basename. Default: input file without extension.",
)
@click.option(
    "--property",
    "property_name",
    default=None,
    help=(
        "HA property to plot, for example F, Cv, S, Utot, or all. "
        "Default: standard compact set."
    ),
)
@figure_preset_option()
@click.option(
    "--format",
    "image_format",
    type=click.Choice(["png", "pdf", "svg"], case_sensitive=False),
    default="png",
    show_default=True,
    help="Output image format.",
)
@click.option(
    "--dpi",
    type=click.IntRange(min=1),
    default=None,
    help="Override the raster resolution supplied by the figure preset.",
)
@click.option(
    "--show",
    is_flag=True,
    default=False,
    help="Show figures on screen after generation.",
)
@click.option(
    "--unit",
    default=None,
    help="Energy unit used for plotted energy-like values. Defaults to stored HDF5 units.",
)
def plot(
    filename: Path,
    outbase: Path | None,
    property_name: str | None,
    figure_preset: str,
    image_format: str,
    dpi: int | None,
    show: bool,
    unit: str | None,
) -> None:
    """
    Generate static Matplotlib plots from a Quantas HA HDF5 result file.
    """
    if outbase is None:
        outbase = filename.with_suffix("")

    try:
        result = read_ha_hdf5(filename)
        collection = build_ha_plots(
            result,
            properties=property_name,
            unit=unit,
        )
        rendered = render_plot_collection(
            collection,
            MatplotlibOptions.from_preset(
                figure_preset,
                output_dir=outbase.parent,
                filename_prefix=f"{outbase.name}_",
                image_format=image_format.lower(),
                dpi=dpi,
                show=show,
                close=not show,
            ),
        )
        plot_files = [
            artifact.path
            for artifact in rendered.artifacts
            if artifact.path is not None
        ]
    except Exception as exc:
        echo_error(str(exc))
        raise click.Abort() from exc

    for plot_file in plot_files:
        echo(f"Plot saved to: {plot_file}")


def _stored_property_unit(result: HAResult, property_name: str) -> str:
    """
    Return the stored unit label for a selected HA property.

    Parameters
    ----------
    result : HAResult
        Harmonic-approximation result object.
    property_name : str
        Requested property name or historical key.

    Returns
    -------
    str
        Stored unit label, used as the default CLI prompt value.
    """
    units = result.metadata.get("units", {})
    key = property_name.strip()
    if key in {"S", "entropy"}:
        return str(units.get("entropy", "unknown"))
    if key in {"Cv", "isochoric_heat_capacity"}:
        return str(units.get("heat_capacity", "unknown"))
    if key in {"V", "volume"}:
        return str(units.get("volume", "unknown"))
    if key in {"T", "temperature"}:
        return str(units.get("temperature", "K"))
    return str(units.get("energy", "unknown"))

apply_reference_help(ha, ('ha',))
