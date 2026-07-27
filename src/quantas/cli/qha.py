# -*- coding: utf-8 -*-

"""Command-line interface for quasi-harmonic approximation workflows.

The commands defined here are thin frontend wrappers around the public QHA API.
They parse command-line options, attach text observers, and delegate numerical
work to the library layer.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import click

from quantas.cli.reference_help import apply_reference_help

from quantas.cli.contracts import (
    DOMAIN_GROUP,
    NUMERICAL_GROUP,
    OUTPUT_GROUP,
    SCIENTIFIC_GROUP,
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
from quantas.api.common import EventLevel
from quantas.api.qha import (
    CurveAxis as QHACurveAxis,
    FitFailurePolicy as QHAFitFailurePolicy,
    Minimization as QHAMinimization,
    ModeContinuity as QHAModeContinuity,
    Options as QHAOptions,
    PlotOptions as QHAPlotOptions,
    PolynomialDerivativeMethod as QHAPolynomialDerivativeMethod,
    Scheme as QHAScheme,
    ThermalExpansionMethod as QHAThermalExpansionMethod,
    available_energy_eos,
    build_plots as build_qha_plots,
    inspect as inspect_qha_input,
    list_plot_properties as list_available_plot_properties,
    read_input as read_qha_input,
    read_result as read_qha_hdf5,
    run as run_qha,
    write_result as write_qha_hdf5,
    write_table as write_qha_table,
)
from quantas.cli.output import CLIOutput
from quantas.cli.qha_observer import QHATextObserver
from quantas.cli.phonon_input import phonon_inpgen
from quantas.cli.qha_render import preview_report_tables
from quantas.renderers.plots import MatplotlibOptions, render_plot_collection
from quantas.references import module_citation_keys, render_citation_notice

_QHA_ENERGY_EOS_CHOICES = available_energy_eos()


@click.group(name="qha")
def qha() -> None:
    """Quasi-harmonic approximation thermodynamic analysis."""
    return


qha.add_command(phonon_inpgen)


@qha.command(name="inspect", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@grouped_option(
    "--eos",
    group=SCIENTIFIC_GROUP,
    type=click.Choice(_QHA_ENERGY_EOS_CHOICES, case_sensitive=True),
    default="BM3",
    show_default=True,
    help="Equation-of-state family and order used for the preview fit.",
)
@grouped_option(
    "--degree",
    group=NUMERICAL_GROUP,
    type=int,
    default=3,
    show_default=True,
    help="Polynomial degree used for the preview polynomial fit.",
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
    "--punit",
    group=UNITS_GROUP,
    type=click.Choice(["GPa", "kbar"], case_sensitive=True),
    default="GPa",
    show_default=True,
    help="Measurement unit for pressure values.",
)
@grouped_option(
    "--no-polynomial",
    group=VALIDATION_GROUP,
    is_flag=True,
    default=False,
    help="Do not include the polynomial pressure estimate.",
)
@grouped_option(
    "--no-eos",
    group=VALIDATION_GROUP,
    is_flag=True,
    default=False,
    help="Do not include the EOS pressure estimate.",
)
@grouped_option(
    "-q",
    "--quiet",
    "silent",
    group=OUTPUT_GROUP,
    is_flag=True,
    default=False,
    help="Do not print output on screen.",
)
@grouped_option(
    "-r",
    "--report",
    group=OUTPUT_GROUP,
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Text report file. Default: do not write a report file.",
)
def inspect(
    filename: Path,
    eos: str,
    degree: int,
    eunit: str,
    vunit: str,
    punit: str,
    no_polynomial: bool,
    no_eos: bool,
    silent: bool,
    report: Path | None,
) -> None:
    """Inspect static energy-volume data before a full QHA run."""
    options = QHAOptions(
        energy_unit=eunit,
        volume_unit=vunit,
        pressure_unit=punit,
        eos=eos,
        energy_degree=degree,
    )

    try:
        preview = inspect_qha_input(
            filename,
            options=options,
            include_polynomial=not no_polynomial,
            include_eos=not no_eos,
            polynomial_degree=degree,
            eos=eos,
        )
    except Exception as exc:
        echo_error(quantas_error(), bold=True)
        echo_error(str(exc))
        raise click.Abort() from exc

    output = CLIOutput(report_file=report, silent=silent)
    output.tables(preview_report_tables(preview, include_diagnostics=True))
    for warning in preview.warnings:
        output.message(warning, level=EventLevel.WARNING)
    output.save()


@qha.command(name="run", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@grouped_option(
    "-S",
    "--scheme",
    group=SCIENTIFIC_GROUP,
    type=click.Choice(["freq", "td"], case_sensitive=True),
    default="freq",
    show_default=True,
    help="QHA interpolation scheme.",
)
@grouped_option(
    "--mode-continuity",
    group=SCIENTIFIC_GROUP,
    type=click.Choice(
        ["verified", "assumed", "unknown", "unreliable"],
        case_sensitive=True,
    ),
    default=None,
    help="Override the phonon-mode continuity status stored in the input file.",
)
@grouped_option(
    "-N",
    "--minimization",
    group=SCIENTIFIC_GROUP,
    type=click.Choice(["poly", "eos"], case_sensitive=True),
    default="poly",
    show_default=True,
    help="Volume minimization method.",
)
@grouped_option(
    "-E",
    "--eos",
    group=SCIENTIFIC_GROUP,
    type=click.Choice(_QHA_ENERGY_EOS_CHOICES, case_sensitive=True),
    default="BM3",
    show_default=True,
    help="Equation-of-state family and order used for EOS minimization.",
)
@grouped_option(
    "--gruneisen/--no-gruneisen",
    "calculate_gruneisen",
    group=SCIENTIFIC_GROUP,
    default=True,
    show_default=True,
    help="Calculate the thermodynamic Gruneisen parameter.",
)
@grouped_option(
    "--mode-gruneisen/--no-mode-gruneisen",
    "calculate_mode_gruneisen",
    group=SCIENTIFIC_GROUP,
    default=True,
    show_default=True,
    help="Calculate mode-resolved Gruneisen parameters for the frequency scheme.",
)
@grouped_option(
    "--thermal-expansion",
    "thermal_expansion_method",
    group=SCIENTIFIC_GROUP,
    type=click.Choice(
        ["mixed_derivative", "mode_gruneisen", "numerical"],
        case_sensitive=True,
    ),
    default="mixed_derivative",
    show_default=True,
    help=(
        "Method used to calculate volumetric thermal expansion. "
        "mode_gruneisen is available only with --scheme=freq."
    ),
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
    "-P",
    "--pressure",
    group=DOMAIN_GROUP,
    nargs=3,
    type=float,
    metavar="MIN MAX STEP",
    default=(0.0, 0.0, 1.0),
    show_default="0.0 0.0 1.0",
    help="Pressure range.",
)
@grouped_option(
    "-D",
    "--energy-degree",
    group=NUMERICAL_GROUP,
    type=click.IntRange(min=2),
    default=3,
    show_default=True,
    help=(
        "Polynomial degree for static energy, free energy, and "
        "thermodynamic-property fits."
    ),
)
@grouped_option(
    "-F",
    "--frequency-degree",
    group=NUMERICAL_GROUP,
    type=click.IntRange(min=1),
    default=3,
    show_default=True,
    help="Polynomial degree for mode-resolved frequency fits.",
)
@grouped_option(
    "--poly-derivatives",
    "polynomial_derivative_method",
    group=NUMERICAL_GROUP,
    type=click.Choice(["local_grid", "analytic"], case_sensitive=True),
    default="local_grid",
    show_default=True,
    help="Method used to calculate KT and Kp after polynomial minimization.",
)
@grouped_option(
    "--poly-grid-points",
    "polynomial_grid_points",
    group=NUMERICAL_GROUP,
    type=click.IntRange(min=3),
    default=5,
    show_default=True,
    help="Odd number of volumes in the local polynomial derivative grid.",
)
@grouped_option(
    "--poly-grid-separation",
    "polynomial_grid_separation",
    group=NUMERICAL_GROUP,
    type=click.FloatRange(min=0.0, min_open=True),
    default=0.05,
    show_default=True,
    help="Adjacent local-grid spacing as a percentage of equilibrium volume.",
)
@grouped_option(
    "--gruneisen-min-cv-fraction",
    group=VALIDATION_GROUP,
    type=click.FloatRange(min=0.0),
    default=1.0e-2,
    show_default=True,
    help="Minimum Cv/Dulong-Petit fraction for the macroscopic Gruneisen ratio.",
)
@grouped_option(
    "--max-failures",
    group=VALIDATION_GROUP,
    type=click.IntRange(min=1),
    default=5,
    show_default=True,
    help="Maximum number of consecutive failed fits before stopping.",
)
@grouped_option(
    "--failure-policy",
    group=VALIDATION_GROUP,
    type=click.Choice(["continue", "stop", "raise"], case_sensitive=True),
    default="stop",
    show_default=True,
    help="Policy applied when local fits fail.",
)
@grouped_option(
    "--eunit", group=UNITS_GROUP, default="Ha", show_default=True, help="Energy unit."
)
@grouped_option(
    "--vunit", group=UNITS_GROUP, default="A", show_default=True, help="Volume unit."
)
@grouped_option(
    "--funit",
    group=UNITS_GROUP,
    default="cm^-1",
    show_default=True,
    help="Frequency unit.",
)
@grouped_option(
    "--tunit",
    group=UNITS_GROUP,
    default="K",
    show_default=True,
    help="Temperature unit.",
)
@grouped_option(
    "--punit",
    group=UNITS_GROUP,
    default="GPa",
    show_default=True,
    help="Pressure unit.",
)
@output_option(
    help="Output HDF5 file. Default: input file base name + '_QHA.hdf5'.",
)
@force_option()
@report_option()
@verbosity_option()
@quiet_option()
@progress_option()
def run(
    filename: Path,
    scheme: str,
    mode_continuity: str | None,
    minimization: str,
    eos: str,
    calculate_gruneisen: bool,
    calculate_mode_gruneisen: bool,
    thermal_expansion_method: str,
    temperature: tuple[float, float, float],
    pressure: tuple[float, float, float],
    energy_degree: int,
    frequency_degree: int,
    polynomial_derivative_method: str,
    polynomial_grid_points: int,
    polynomial_grid_separation: float,
    gruneisen_min_cv_fraction: float,
    max_failures: int,
    failure_policy: str,
    eunit: str,
    vunit: str,
    funit: str,
    tunit: str,
    punit: str,
    output: Path | None,
    force: bool,
    report: Path | None,
    verbosity: str,
    quiet: bool,
    progress: bool,
) -> None:
    """Run a QHA calculation from a Quantas YAML input file."""
    destination = default_hdf5_path(filename, output, suffix="_QHA")
    report = default_report_path(filename, report)
    report_verbosity = parse_verbosity(verbosity)

    echo_highlight(quantas_title(), silent=quiet)

    options = QHAOptions(
        temperature_min=temperature[0],
        temperature_max=temperature[1],
        temperature_step=temperature[2],
        pressure_min=pressure[0],
        pressure_max=pressure[1],
        pressure_step=pressure[2],
        scheme=cast(QHAScheme, scheme),
        minimization=cast(QHAMinimization, minimization),
        eos=eos,
        energy_degree=energy_degree,
        free_energy_degree=energy_degree,
        frequency_degree=frequency_degree,
        polynomial_derivative_method=cast(
            QHAPolynomialDerivativeMethod, polynomial_derivative_method
        ),
        polynomial_grid_points=polynomial_grid_points,
        polynomial_grid_separation=polynomial_grid_separation,
        calculate_gruneisen=calculate_gruneisen,
        calculate_mode_gruneisen=calculate_mode_gruneisen,
        thermal_expansion_method=cast(
            QHAThermalExpansionMethod, thermal_expansion_method
        ),
        gruneisen_min_cv_fraction=gruneisen_min_cv_fraction,
        energy_unit=eunit,
        volume_unit=vunit,
        frequency_unit="cm^-1" if funit == "cm-1" else funit,
        temperature_unit=tunit,
        pressure_unit=punit,
        debug=report_verbosity.includes_debug,
        max_consecutive_failures=max_failures,
        fit_failure_policy=cast(QHAFitFailurePolicy, failure_policy),
    )
    observer = QHATextObserver(
        report_file=report,
        silent=quiet,
        show_progress=progress,
        verbosity=report_verbosity,
    )

    try:
        calculation_input = filename
        if mode_continuity is not None:
            input_data = read_qha_input(filename)
            continuity = cast(QHAModeContinuity, mode_continuity)
            input_data.mode_continuity = continuity
            input_data.metadata["mode_continuity"] = continuity
            result = run_qha(input_data, options=options, observer=observer)
        else:
            result = run_qha(calculation_input, options=options, observer=observer)
    except Exception as exc:
        observer.close()
        echo_error(quantas_error(), bold=True)
        echo_error(str(exc))
        raise click.Abort() from exc

    observer.output.text_block(render_citation_notice(module_citation_keys("qha")))
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
        write_qha_hdf5(result, destination, report_text=observer.text())

    echo_highlight(quantas_finish(), silent=quiet)


@qha.command(name="plot", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "-p",
    "--property",
    "property_names",
    multiple=True,
    help=(
        "QHA property to plot, for example VT, KT, KS, alphaV, Cp, Cv, "
        "or heat_capacities. Can be used more than once. Default: "
        "standard properties."
    ),
)
@click.option(
    "--axis",
    "curve_axis",
    type=click.Choice(["temperature", "pressure"], case_sensitive=False),
    default="temperature",
    show_default=True,
    help="Independent variable used for one-dimensional line sections.",
)
@click.option(
    "--pressure",
    "selected_pressures",
    type=float,
    multiple=True,
    help=(
        "Exact native pressure included in temperature sections. May be "
        "repeated; default: all pressures."
    ),
)
@click.option(
    "--temperature",
    "selected_temperatures",
    type=float,
    multiple=True,
    help=(
        "Exact native temperature included in pressure sections. May be "
        "repeated; default: all temperatures."
    ),
)
@click.option(
    "--2d",
    "include_2d",
    is_flag=True,
    default=False,
    help=(
        "Generate filled pressure-temperature contour maps when enough "
        "data are available."
    ),
)
@click.option(
    "--cmap",
    type=click.Choice(
        ["viridis", "plasma", "inferno", "magma", "cividis", "turbo"],
        case_sensitive=False,
    ),
    default="viridis",
    show_default=True,
    help="Colormap used for two-dimensional contour maps.",
)
@click.option(
    "--contour-mode",
    type=click.Choice(["discrete", "smooth"], case_sensitive=False),
    default="smooth",
    show_default=True,
    help="Filled contour rendering mode.",
)
@click.option(
    "--levels",
    type=int,
    default=12,
    show_default=True,
    help="Number of contour levels or isolines.",
)
@click.option(
    "--isolines/--no-isolines",
    default=True,
    show_default=True,
    help="Draw contour lines on two-dimensional maps.",
)
@click.option(
    "--isoline-labels/--no-isoline-labels",
    default=True,
    show_default=True,
    help="Label contour lines with isovalues.",
)
@click.option(
    "--temperature-unit",
    type=click.Choice(
        ["K", "C", "Celsius", "Fahrenheit", "Rankine"], case_sensitive=False
    ),
    default=None,
    help="Display temperature unit used only for plotting. Default: HDF5 result unit.",
)
@click.option(
    "--pressure-unit",
    type=click.Choice(
        ["Pa", "kPa", "MPa", "GPa", "bar", "kbar", "Mbar"], case_sensitive=False
    ),
    default=None,
    help="Display pressure unit used only for plotting. Default: HDF5 result unit.",
)
@click.option(
    "--energy-unit",
    type=click.Choice(["J/mol", "kJ/mol", "Ha", "eV", "Ry"], case_sensitive=False),
    default=None,
    help=(
        "Display energy unit used only for plotting energy, entropy and "
        "heat-capacity quantities."
    ),
)
@click.option(
    "--dulong-petit/--no-dulong-petit",
    default=True,
    show_default=True,
    help=(
        "Draw the Dulong-Petit limit on Cv and combined heat-capacity "
        "plots when atom-count metadata are available."
    ),
)
@click.option(
    "-o",
    "--output",
    "output_dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Output directory for plot files. Default: input file base name + '_plots'.",
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
    "--list-properties",
    is_flag=True,
    default=False,
    help="List available properties and exit without generating figures.",
)
def plot(
    filename: Path,
    property_names: tuple[str, ...],
    curve_axis: str,
    selected_pressures: tuple[float, ...],
    selected_temperatures: tuple[float, ...],
    include_2d: bool,
    cmap: str,
    contour_mode: str,
    levels: int,
    isolines: bool,
    isoline_labels: bool,
    temperature_unit: str | None,
    pressure_unit: str | None,
    energy_unit: str | None,
    dulong_petit: bool,
    output_dir: Path | None,
    figure_preset: str,
    image_format: str,
    dpi: int | None,
    show: bool,
    list_properties: bool,
) -> None:
    """Generate QHA plots from a HDF5 result file."""
    try:
        result = read_qha_hdf5(filename)
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc

    if list_properties:
        for key, attribute, description in list_available_plot_properties(result):
            echo(f"{key:18s} {attribute:30s} {description}")
        return

    if output_dir is None:
        stem = filename.with_suffix("")
        output_dir = stem.with_name(f"{stem.name}_plots")

    plot_options = QHAPlotOptions(
        include_contours=include_2d,
        cmap=cmap.lower(),
        contour_mode=contour_mode.lower(),
        levels=levels,
        isolines=isolines,
        isoline_labels=isoline_labels,
        temperature_unit=temperature_unit,
        pressure_unit=pressure_unit,
        energy_unit=energy_unit,
        include_dulong_petit=dulong_petit,
        curve_axis=cast(QHACurveAxis, curve_axis.lower()),
        selected_pressures=selected_pressures or None,
        selected_temperatures=selected_temperatures or None,
    )
    render_options = MatplotlibOptions.from_preset(
        figure_preset,
        output_dir=output_dir,
        image_format=image_format.lower(),
        dpi=dpi,
        show=show,
        close=not show,
    )
    try:
        collection = build_qha_plots(
            result,
            properties=property_names or None,
            options=plot_options,
        )
        plot_result = render_plot_collection(collection, render_options)
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc

    for warning in plot_result.warnings:
        echo_warning(str(warning))
    for artifact in plot_result.artifacts:
        if artifact.path is not None:
            echo(f"Written {artifact.path}")


@qha.command(name="export", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "-o",
    "--output",
    "outfile",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Output table file. Default: input file base name + '_table.dat'.",
)
@click.option(
    "--property",
    "property_name",
    default=None,
    help=(
        "Single QHA property to export, for example VT, KT, G, or "
        "equilibrium_volume. Use 'structure' to export equilibrium lattice "
        "parameters and anisotropic thermal expansion only. Default: export "
        "all pressure-temperature properties."
    ),
)
@click.option(
    "--format",
    "table_format",
    type=click.Choice(["txt", "csv"], case_sensitive=False),
    default="txt",
    show_default=True,
    help="Output table format.",
)
@click.option(
    "--no-uncertainty",
    is_flag=True,
    default=False,
    help="Do not include uncertainty columns, even when available.",
)
def export(
    filename: Path,
    outfile: Path | None,
    property_name: str | None,
    table_format: str,
    no_uncertainty: bool,
) -> None:
    """Export QHA pressure-temperature tables from a HDF5 result file."""
    suffix = ".csv" if table_format.lower() == "csv" else ".dat"
    if outfile is None:
        stem = filename.with_suffix("")
        tag = property_name if property_name is not None else "table"
        outfile = stem.with_name(f"{stem.name}_{tag}").with_suffix(suffix)

    try:
        result = read_qha_hdf5(filename)
        outfile = write_qha_table(
            result,
            outfile,
            property_name=property_name,
            include_uncertainty=not no_uncertainty,
            file_format=cast(Literal["txt", "csv"], table_format.lower()),
        )
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc

    echo(f"Written {outfile}")

apply_reference_help(qha, ("qha",))
