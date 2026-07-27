# -*- coding: utf-8 -*-

"""Click adapters for the Quantas elasticity workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import click

from quantas.cli.reference_help import apply_reference_help

from quantas.cli.contracts import (
    DOMAIN_GROUP,
    PLOTTING_GROUP,
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
from quantas.cli.elasticity_observer import ElasticityTextObserver
from quantas.cli.grouped_options import GroupedCommand, grouped_option
from quantas.cli.tensor_rotation import (
    resolve_tensor_rotation,
    tensor_rotation_options,
)
from quantas.cli.messages import (
    confirm,
    echo,
    echo_error,
    echo_highlight,
    echo_warning,
    quantas_error,
    quantas_finish,
    quantas_title,
    prompt,
)
from quantas.api.elasticity import (
    Input as ElasticityInput,
    InputInterface as ElasticityInputInterface,
    Options as ElasticityOptions,
    PlotProperty as ElasticityPlotProperty,
    SurfaceGeometry as ElasticitySurfaceGeometry,
    SurfaceOptions as ElasticitySurfaceOptions,
    build_2d_plots as build_elasticity_2d_plots,
    build_3d_plots as build_elasticity_3d_plots,
    create_input as create_elasticity_input,
    read_result as read_elasticity_hdf5,
    run as run_elasticity,
    write_result as write_elasticity_hdf5,
    write_table as write_elasticity_table,
)
from quantas.api.plotting import PlotCollection
from quantas.renderers.plots import (
    MatplotlibOptions,
    render_plot_collection,
    validate_colormap_name,
)
from quantas.references import module_citation_keys, render_citation_notice


_PROPERTY_CHOICE = click.Choice(
    ["young", "compressibility", "shear", "poisson"],
    case_sensitive=False,
)
_GEOMETRY_CHOICE = click.Choice(
    ["physical", "unit-sphere"],
    case_sensitive=False,
)


@click.group(name="elasticity")
def elasticity() -> None:
    """Second-order elastic-tensor analysis."""


@elasticity.command(name="run", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@grouped_option(
    "--2d",
    "calculate_2d",
    group=DOMAIN_GROUP,
    is_flag=True,
    default=False,
    help="Calculate and store properties on the principal Cartesian planes.",
)
@grouped_option(
    "--3d",
    "calculate_3d",
    group=DOMAIN_GROUP,
    is_flag=True,
    default=False,
    help="Calculate and store three-dimensional directional data.",
)
@grouped_option(
    "--ntheta",
    group=DOMAIN_GROUP,
    type=click.IntRange(min=2),
    default=None,
    help="Polar samples: default 361 for 2D and 61 for 3D.",
)
@grouped_option(
    "--nphi",
    group=DOMAIN_GROUP,
    type=click.IntRange(min=3),
    default=None,
    help="Azimuthal samples for 3D. Default: 2 * ntheta - 1.",
)
@grouped_option(
    "--property",
    "properties",
    group=DOMAIN_GROUP,
    type=_PROPERTY_CHOICE,
    multiple=True,
    help="3D property to store or plot. Repeat to select multiple properties.",
)
@grouped_option(
    "-p",
    "--plot/--no-plot",
    group=PLOTTING_GROUP,
    default=False,
    show_default=True,
    help="Render requested 2D/3D figures after the calculation.",
)
@grouped_option(
    "--geometry",
    group=PLOTTING_GROUP,
    type=_GEOMETRY_CHOICE,
    default="physical",
    show_default=True,
    help="Three-dimensional radial representation.",
)
@grouped_option(
    "--color-mode",
    group=PLOTTING_GROUP,
    type=click.Choice(["solid", "property"]),
    default="property",
    show_default=True,
    help="Use monochromatic surfaces or a property colormap.",
)
@grouped_option(
    "--cmap",
    group=PLOTTING_GROUP,
    default="viridis",
    show_default=True,
    help="Colormap name.",
)
@grouped_option(
    "--mesh/--no-mesh",
    "show_mesh",
    group=PLOTTING_GROUP,
    default=False,
    show_default=True,
    help="Draw mesh edges on 3D surfaces.",
)
@grouped_option(
    "--mesh-color",
    group=PLOTTING_GROUP,
    default="black",
    show_default=True,
    help="Mesh edge color for 3D surfaces.",
)
@grouped_option(
    "--mesh-linewidth",
    "mesh_line_width",
    group=PLOTTING_GROUP,
    type=click.FloatRange(min=0.0),
    default=0.5,
    show_default=True,
    help="Mesh edge line width for 3D surfaces.",
)
@grouped_option(
    "--title/--no-title",
    "show_title",
    group=PLOTTING_GROUP,
    default=False,
    show_default=True,
    help="Display a title above 3D figures.",
)
@figure_preset_option(
    option_name="--plot-preset",
    group=PLOTTING_GROUP,
)
@grouped_option(
    "-F",
    "--format",
    "image_format",
    group=PLOTTING_GROUP,
    type=click.Choice(["png", "pdf", "svg"]),
    default="png",
    show_default=True,
    help="Static figure format.",
)
@grouped_option(
    "--dpi",
    group=PLOTTING_GROUP,
    type=click.IntRange(min=1),
    default=None,
    help="Override the raster resolution supplied by the figure preset.",
)
@output_option(
    help="Output HDF5 file. Default: input base name + '_elasticity.hdf5'.",
)
@force_option()
@report_option()
@verbosity_option()
@quiet_option()
@progress_option()
@tensor_rotation_options
def run(
    filename: Path,
    calculate_2d: bool,
    calculate_3d: bool,
    ntheta: int | None,
    nphi: int | None,
    properties: tuple[str, ...],
    plot: bool,
    geometry: str,
    color_mode: str,
    cmap: str,
    show_mesh: bool,
    mesh_color: str,
    mesh_line_width: float,
    show_title: bool,
    figure_preset: str,
    image_format: str,
    dpi: int | None,
    output: Path | None,
    force: bool,
    report: Path | None,
    verbosity: str,
    quiet: bool,
    progress: bool,
    rotate_xyz: tuple[float, float, float] | None,
    rotation_matrix: Path | None,
) -> None:
    """Run an elasticity calculation from a Quantas text input file."""
    if plot and not (calculate_2d or calculate_3d):
        raise click.ClickException("The --plot option requires --2d or --3d.")
    if nphi is not None and not calculate_3d:
        raise click.ClickException("The --nphi option requires --3d.")
    if plot and calculate_3d and color_mode == "property":
        _validate_colormap(cmap)

    destination = default_hdf5_path(filename, output, suffix="_elasticity")
    report = default_report_path(filename, report)
    report_verbosity = parse_verbosity(verbosity)

    ntheta_2d = 361 if ntheta is None else ntheta
    ntheta_3d = 61 if ntheta is None else ntheta
    nphi_3d = 2 * ntheta_3d - 1 if nphi is None else nphi
    selected = _selected_properties(properties)
    surface_options = (
        ElasticitySurfaceOptions(
            ntheta=ntheta_3d,
            nphi=nphi_3d,
            properties=selected,
        )
        if calculate_3d
        else None
    )

    rotation = resolve_tensor_rotation(rotate_xyz, rotation_matrix)

    echo_highlight(quantas_title(), silent=quiet)
    options = ElasticityOptions(
        calculate_2d=calculate_2d,
        ntheta=ntheta_2d,
        calculate_3d=calculate_3d,
        surface_options=surface_options,
        rotation=rotation,
    )
    observer = ElasticityTextObserver(
        report_file=report,
        silent=quiet,
        show_progress=progress,
        verbosity=report_verbosity,
    )

    try:
        result = run_elasticity(filename, options=options, observer=observer)
    except Exception as exc:
        observer.close()
        echo_error(quantas_error(), bold=True)
        echo_error(str(exc))
        raise click.Abort() from exc

    observer.output.text_block(
        render_citation_notice(module_citation_keys("elasticity"))
    )
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
        write_elasticity_hdf5(result, destination, report_text=observer.text())

    if plot and calculate_2d:
        _render_collection(
            build_elasticity_2d_plots(result, selected),
            destination,
            figure_preset=figure_preset,
            dpi=dpi,
            image_format=image_format,
            show_title=show_title,
        )
    if plot and calculate_3d:
        _render_collection(
            build_elasticity_3d_plots(
                result,
                properties=selected,
                geometry=cast(ElasticitySurfaceGeometry, geometry.replace("-", "_")),
                color_mode=cast(Literal["solid", "property"], color_mode),
                colormap=cmap,
                show_mesh=show_mesh,
                mesh_color=mesh_color,
                mesh_line_width=mesh_line_width,
            ),
            destination,
            figure_preset=figure_preset,
            dpi=dpi,
            image_format=image_format,
            show_title=show_title,
        )

    echo_highlight(quantas_finish(), silent=quiet)


@elasticity.command(name="export", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "-o",
    "--output",
    "outfile",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Output table file. Default: input base name + '.dat'.",
)
def export(filename: Path, outfile: Path | None) -> None:
    """Export 2D elasticity data from a Quantas HDF5 result file."""
    outfile = filename.with_suffix(".dat") if outfile is None else outfile
    try:
        result = read_elasticity_hdf5(filename)
        outfile = write_elasticity_table(result, outfile)
    except Exception as exc:
        echo_error(quantas_error(), bold=True)
        echo_error(str(exc))
        raise click.Abort() from exc
    echo(f"Elasticity 2D data exported to: {outfile}")


@elasticity.command(name="inpgen", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "-o",
    "--output",
    "outfile",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Output input file. Default: input base name + '_elasticity_input.dat'.",
)
@click.option(
    "-i",
    "--interface",
    type=click.Choice(["crystal", "vasp"], case_sensitive=True),
    default="crystal",
    show_default=True,
    help="Interface used to read the external output file.",
)
def inpgen(filename: Path, outfile: Path | None, interface: str) -> None:
    """Generate a Quantas elasticity input from CRYSTAL or VASP output."""
    if outfile is None:
        stem = filename.with_suffix("")
        outfile = stem.with_name(stem.name + "_elasticity_input").with_suffix(".dat")
    else:
        outfile = outfile.with_suffix(".dat")
    echo_highlight(quantas_title())
    try:
        jobname = prompt(
            "\nPlease enter a short description for the input file",
            default="Unknown",
            show_default=False,
        )
        outfile = create_elasticity_input(
            filename,
            outfile,
            interface=cast(ElasticityInputInterface, interface),
            jobname=jobname,
        )
    except Exception as exc:
        echo_error(quantas_error(), bold=True)
        echo_error(str(exc))
        raise click.Abort() from exc
    echo(f"Elasticity input file written to: {outfile}")
    echo_highlight(render_citation_notice(module_citation_keys("elasticity")))
    echo_highlight(quantas_finish())


@elasticity.command(name="plot", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "-o",
    "--output",
    "outfile",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Base name for output figures. Default: input file base name.",
)
@click.option("--2d", "plot_2d", is_flag=True, help="Render 2D principal-plane plots.")
@click.option("--3d", "plot_3d", is_flag=True, help="Render 3D elasticity surfaces.")
@click.option(
    "--ntheta",
    type=click.IntRange(min=2),
    default=None,
    help="Override the number of polar samples used to rebuild plotted surfaces.",
)
@click.option(
    "--nphi",
    type=click.IntRange(min=3),
    default=None,
    help="Override the number of azimuthal samples used for 3D surfaces.",
)
@click.option(
    "--property",
    "properties",
    type=_PROPERTY_CHOICE,
    multiple=True,
    help="Elastic property to plot. Repeat to select multiple properties.",
)
@click.option(
    "--geometry",
    type=_GEOMETRY_CHOICE,
    default="physical",
    show_default=True,
    help="Render physical-radius surfaces or values on a unit sphere.",
)
@click.option(
    "--color-mode",
    type=click.Choice(["solid", "property"]),
    default="property",
    show_default=True,
    help="Color 3D surfaces uniformly or by the plotted property.",
)
@click.option(
    "--cmap",
    default="viridis",
    show_default=True,
    help="Matplotlib colormap used when --color-mode=property.",
)
@click.option(
    "--mesh/--no-mesh",
    "show_mesh",
    default=False,
    show_default=True,
    help="Draw mesh edges on 3D surfaces.",
)
@click.option(
    "--mesh-color",
    default="black",
    show_default=True,
    help="Mesh edge color for 3D surfaces.",
)
@click.option(
    "--mesh-linewidth",
    "mesh_line_width",
    type=click.FloatRange(min=0.0),
    default=0.5,
    show_default=True,
    help="Mesh edge line width for 3D surfaces.",
)
@click.option(
    "--title/--no-title",
    "show_title",
    default=False,
    show_default=True,
    help="Show or hide figure titles.",
)
@figure_preset_option()
@click.option(
    "--format",
    "image_format",
    type=click.Choice(["png", "pdf", "svg"]),
    default="png",
    show_default=True,
    help="Static figure format.",
)
@click.option(
    "--dpi",
    type=click.IntRange(min=1),
    default=None,
    help="Override the raster resolution supplied by the figure preset.",
)
def plot(
    filename: Path,
    outfile: Path | None,
    plot_2d: bool,
    plot_3d: bool,
    ntheta: int | None,
    nphi: int | None,
    properties: tuple[str, ...],
    geometry: str,
    color_mode: str,
    cmap: str,
    show_mesh: bool,
    mesh_color: str,
    mesh_line_width: float,
    show_title: bool,
    figure_preset: str,
    image_format: str,
    dpi: int | None,
) -> None:
    """Create 2D and/or 3D elasticity plots from a native HDF5 result."""
    if not (plot_2d or plot_3d):
        raise click.ClickException("Select at least one of --2d or --3d.")
    if nphi is not None and not plot_3d:
        raise click.ClickException("The --nphi option requires --3d.")
    if plot_3d and color_mode == "property":
        _validate_colormap(cmap)

    basename = filename.with_suffix("") if outfile is None else outfile.with_suffix("")
    selected = _selected_properties(properties)
    try:
        result_data = read_elasticity_hdf5(filename)
        payload = result_data.results["elasticity"]
        if plot_2d:
            if ntheta is not None or not payload.has_2d_data():
                if payload.stiffness is None:
                    raise click.ClickException("The result has no stiffness matrix.")
                result_data_2d = run_elasticity(
                    ElasticityInput(
                        jobname=payload.jobname,
                        stiffness=payload.stiffness,
                    ),
                    options=ElasticityOptions(
                        calculate_2d=True,
                        ntheta=361 if ntheta is None else ntheta,
                    ),
                )
            else:
                result_data_2d = result_data
            _render_collection(
                build_elasticity_2d_plots(result_data_2d, selected),
                basename,
                figure_preset=figure_preset,
                dpi=dpi,
                image_format=image_format,
                prefix_name=basename.name,
                show_title=show_title,
            )
        if plot_3d:
            if ntheta is None and nphi is None and payload.has_3d_data():
                surface_options = None
            else:
                ntheta_3d = 61 if ntheta is None else ntheta
                nphi_3d = 2 * ntheta_3d - 1 if nphi is None else nphi
                surface_options = ElasticitySurfaceOptions(
                    ntheta=ntheta_3d,
                    nphi=nphi_3d,
                    properties=selected,
                )
            _render_collection(
                build_elasticity_3d_plots(
                    result_data,
                    surface_options,
                    properties=selected,
                    geometry=cast(
                        ElasticitySurfaceGeometry, geometry.replace("-", "_")
                    ),
                    color_mode=cast(Literal["solid", "property"], color_mode),
                    colormap=cmap,
                    show_mesh=show_mesh,
                    mesh_color=mesh_color,
                    mesh_line_width=mesh_line_width,
                ),
                basename,
                figure_preset=figure_preset,
                dpi=dpi,
                image_format=image_format,
                prefix_name=basename.name,
                show_title=show_title,
            )
    except click.ClickException:
        raise
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc


def _selected_properties(
    properties: tuple[str, ...],
) -> tuple[ElasticityPlotProperty, ...]:
    """Return unique requested property names or the complete default set."""
    selected = properties or ("young", "compressibility", "shear", "poisson")
    return tuple(
        cast(ElasticityPlotProperty, name.lower()) for name in dict.fromkeys(selected)
    )


def _render_collection(
    collection: PlotCollection,
    output_base: Path,
    *,
    figure_preset: str,
    dpi: int | None,
    image_format: str,
    prefix_name: str | None = None,
    show_title: bool = False,
) -> None:
    """Render and report one neutral plot collection."""
    rendered = render_plot_collection(
        collection,
        MatplotlibOptions.from_preset(
            figure_preset,
            output_dir=output_base.parent,
            filename_prefix=f"{prefix_name or output_base.stem}_",
            image_format=image_format,
            dpi=dpi,
            close=True,
            show_title=show_title,
        ),
    )
    for warning in rendered.warnings:
        echo_warning(str(warning))
    for artifact in rendered.artifacts:
        if artifact.path is not None:
            echo(f"Plot saved to: {artifact.path}")


def _validate_colormap(name: str) -> None:
    """Raise a clean Click error for an unknown Matplotlib colormap."""
    try:
        validate_colormap_name(name)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

apply_reference_help(elasticity, ('elasticity',))
