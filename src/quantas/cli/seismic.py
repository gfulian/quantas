# -*- coding: utf-8 -*-

"""Click commands for sampled seismic-wave calculations and plots."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import click

from quantas.cli.reference_help import apply_reference_help

from quantas.cli.contracts import (
    ADVANCED_GROUP,
    DOMAIN_GROUP,
    NUMERICAL_GROUP,
    SCIENTIFIC_GROUP,
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
from quantas.cli.general import MostSimilarCommandGroup
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
from quantas.cli.output import CLIOutput
from quantas.cli.seismic_observer import SeismicProgressObserver
from quantas.cli.tensor_rotation import (
    resolve_tensor_rotation,
    tensor_rotation_options,
)
from quantas.core.geometry import Hemisphere
from quantas.core.physics.seismic import SamplingLevel, WaveMode
from quantas.api.plotting import (
    PlotCollection,
    SphericalMapSpec,
    SphericalProjection,
    SphericalSummarySpec,
    SurfacePlotSpec,
)
from quantas.api.seismic import (
    Options as SeismicOptions,
    PlotOptions as SeismicPlotOptions,
    SurfaceGeometry,
    SurfaceOptions as SeismicSurfaceOptions,
    SurfaceType,
    build_plots as build_seismic_plots,
    build_report as build_seismic_report,
    build_summary as build_seismic_summary,
    build_surfaces as build_seismic_surfaces,
    read_result as read_seismic_hdf5,
    run as run_seismic,
    write_csv as write_seismic_csv,
    write_result as write_seismic_hdf5,
)
from quantas.renderers.plots import (
    MatplotlibOptions,
    render_plot,
)
from quantas.references import module_citation_keys, render_citation_notice


_PROPERTY_KEYS = (
    "phase_v_p",
    "phase_v_s1",
    "phase_v_s2",
    "shear_anisotropy",
    "shear_splitting",
    "phase_v_p_over_v_s1",
    "phase_v_p_over_v_s2",
    "group_v_p",
    "group_v_s1",
    "group_v_s2",
    "power_flow_v_p",
    "power_flow_v_s1",
    "power_flow_v_s2",
    "log10_enhancement_v_p",
    "log10_enhancement_v_s1",
    "log10_enhancement_v_s2",
)


@click.group(cls=MostSimilarCommandGroup)
def seismic() -> None:
    """Calculate, inspect, export, and plot acoustic-wave properties."""


@seismic.command(name="run", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@grouped_option(
    "--hemisphere",
    group=SCIENTIFIC_GROUP,
    type=click.Choice([item.value for item in Hemisphere]),
    default=Hemisphere.UPPER.value,
    show_default=True,
    help="Sample the selected spherical hemisphere.",
)
@grouped_option(
    "-L",
    "--level",
    group=SCIENTIFIC_GROUP,
    type=click.Choice([item.value for item in SamplingLevel]),
    default=SamplingLevel.ENHANCEMENT.value,
    show_default=True,
    help="Highest acoustic-property level calculated and persisted.",
)
@grouped_option(
    "--track-polarizations/--no-track-polarizations",
    group=SCIENTIFIC_GROUP,
    default=True,
    show_default=True,
    help="Track shear-wave polarization axes across the sampled sphere.",
)
@grouped_option(
    "--ntheta",
    group=DOMAIN_GROUP,
    type=click.IntRange(min=2),
    default=91,
    show_default=True,
    help="Number of polar samples.",
)
@grouped_option(
    "--nphi",
    group=DOMAIN_GROUP,
    type=click.IntRange(min=3),
    default=181,
    show_default=True,
    help="Number of azimuthal samples.",
)
@grouped_option(
    "--batch-size",
    group=NUMERICAL_GROUP,
    type=click.IntRange(min=1),
    default=512,
    show_default=True,
    help="Directions solved per vectorized numerical batch.",
)
@grouped_option(
    "--eigenvalue-rtol",
    group=ADVANCED_GROUP,
    hidden=True,
    type=float,
    default=1.0e-10,
    show_default=True,
    help="Relative tolerance for small negative Christoffel eigenvalues.",
)
@grouped_option(
    "--eigenvalue-atol",
    group=ADVANCED_GROUP,
    hidden=True,
    type=float,
    default=1.0e-12,
    show_default=True,
    help="Absolute tolerance for small negative Christoffel eigenvalues.",
)
@grouped_option(
    "--degeneracy-rtol",
    group=ADVANCED_GROUP,
    hidden=True,
    type=float,
    default=1.0e-8,
    show_default=True,
    help="Relative tolerance used to identify acoustic-mode degeneracies.",
)
@grouped_option(
    "--degeneracy-atol",
    group=ADVANCED_GROUP,
    hidden=True,
    type=float,
    default=1.0e-10,
    show_default=True,
    help="Absolute tolerance used to identify acoustic-mode degeneracies.",
)
@grouped_option(
    "--pseudoinverse-rcond",
    group=ADVANCED_GROUP,
    hidden=True,
    type=float,
    default=1.0e-10,
    show_default=True,
    help="Relative singular-value cutoff used in analytical eigenvalue Hessians.",
)
@grouped_option(
    "--caustic-rtol",
    group=ADVANCED_GROUP,
    hidden=True,
    type=float,
    default=1.0e-10,
    show_default=True,
    help="Relative tolerance used to identify possible caustic candidates.",
)
@grouped_option(
    "--caustic-atol",
    group=ADVANCED_GROUP,
    hidden=True,
    type=float,
    default=1.0e-12,
    show_default=True,
    help="Absolute tolerance used to identify possible caustic candidates.",
)
@output_option(
    help="HDF5 result file. Default: input base name + '_seismic.hdf5'.",
)
@force_option()
@report_option()
@verbosity_option()
@quiet_option()
@progress_option()
@tensor_rotation_options
def run(
    filename: Path,
    hemisphere: str,
    level: str,
    track_polarizations: bool,
    ntheta: int,
    nphi: int,
    batch_size: int,
    eigenvalue_rtol: float,
    eigenvalue_atol: float,
    degeneracy_rtol: float,
    degeneracy_atol: float,
    pseudoinverse_rcond: float,
    caustic_rtol: float,
    caustic_atol: float,
    output: Path | None,
    force: bool,
    report: Path | None,
    verbosity: str,
    quiet: bool,
    progress: bool,
    rotate_xyz: tuple[float, float, float] | None,
    rotation_matrix: Path | None,
) -> None:
    """Run SEISMIC, print a detailed report, and save an HDF5 result."""
    destination = default_hdf5_path(filename, output, suffix="_seismic")
    report = default_report_path(filename, report)
    report_verbosity = parse_verbosity(verbosity)

    rotation = resolve_tensor_rotation(rotate_xyz, rotation_matrix)

    if not quiet:
        echo_highlight(quantas_title())
    options = SeismicOptions(
        ntheta=ntheta,
        nphi=nphi,
        hemisphere=Hemisphere(hemisphere),
        level=SamplingLevel(level),
        batch_size=batch_size,
        track_polarization_axes=track_polarizations,
        eigenvalue_rtol=eigenvalue_rtol,
        eigenvalue_atol=eigenvalue_atol,
        degeneracy_rtol=degeneracy_rtol,
        degeneracy_atol=degeneracy_atol,
        pseudoinverse_rcond=pseudoinverse_rcond,
        caustic_rtol=caustic_rtol,
        caustic_atol=caustic_atol,
        rotation=rotation,
    )
    output_service = CLIOutput(
        report_file=report,
        silent=quiet,
        show_progress=progress,
    )
    observer = SeismicProgressObserver(
        silent=quiet,
        show_progress=progress,
        output=output_service,
    )
    try:
        result = run_seismic(filename, options=options, observer=observer)
        output_service.tables(
            build_seismic_report(result, level=report_verbosity.value)
        )
        output_service.text_block(
            render_citation_notice(module_citation_keys("seismic"))
        )
        output_service.save()
        report_text = output_service.text()
    except Exception as exc:
        observer.close()
        echo_error(quantas_error(), bold=True)
        echo_error(str(exc))
        raise click.Abort() from exc

    overwrite = not destination.exists() or force
    if destination.exists() and not force:
        overwrite = confirm(
            f"Output file {destination} exists. Would you like to overwrite it?",
            default=False,
        )
    if overwrite:
        write_seismic_hdf5(result, destination, report_text=report_text)
        if not quiet:
            echo(f"SEISMIC result written to: {destination}")
            echo(f"SEISMIC report written to: {report}")
    elif not quiet:
        echo("Results not saved")

    if not quiet:
        echo_highlight(quantas_finish())


@seismic.command(name="plot", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "-o",
    "--output",
    "outfile",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Base name for output figures. Default: HDF5 base name.",
)
@click.option("--2d", "plot_2d", is_flag=True, help="Render spherical 2D maps.")
@click.option("--3d", "plot_3d", is_flag=True, help="Render 3D acoustic surfaces.")
@click.option("--summary", is_flag=True, help="Render the six-panel seismic summary.")
@click.option(
    "--property",
    "properties",
    type=click.Choice(_PROPERTY_KEYS),
    multiple=True,
    help="Scalar property for 2D or 3D output. Repeat to select multiple properties.",
)
@click.option(
    "--surface",
    "surface_types",
    type=click.Choice(["phase", "slowness", "group"]),
    multiple=True,
    help="3D surface geometry. Repeat to select multiple types.",
)
@click.option(
    "--mode",
    "modes",
    type=click.Choice(
        [mode.value for mode in (WaveMode.V_P, WaveMode.V_S1, WaveMode.V_S2)]
    ),
    multiple=True,
    help="Acoustic mode used for 3D surfaces. Repeat to select multiple modes.",
)
@click.option(
    "--projection",
    type=click.Choice(["equal_area", "stereographic"]),
    default="equal_area",
    show_default=True,
    help="Spherical projection used for two-dimensional maps.",
)
@click.option(
    "--cmap",
    default=None,
    help="Override the property-dependent default with a Matplotlib colormap.",
)
@click.option(
    "--levels",
    type=click.IntRange(min=2),
    default=16,
    show_default=True,
    help="Number of filled-contour levels used for two-dimensional maps.",
)
@click.option(
    "--extrema/--no-extrema",
    default=True,
    show_default=True,
    help="Show or hide minimum and maximum markers.",
)
@click.option(
    "--polarization/--no-polarization",
    default=True,
    show_default=True,
    help=(
        "Also render separate polarization-overlay variants whenever a "
        "tracked polarization field is available."
    ),
)
@click.option(
    "--polarization-stride", type=click.IntRange(min=1), default=8, show_default=True,
    help="Sampling stride between plotted polarization axes.",
)
@click.option(
    "--polarization-color",
    default="black",
    show_default=True,
    help="Color of plotted polarization axes.",
)
@click.option(
    "--polarization-linewidth",
    type=click.FloatRange(min=0.01),
    default=1.2,
    show_default=True,
    help="Line width of plotted polarization axes.",
)
@click.option(
    "--polarization-scale",
    type=click.FloatRange(min=0.001),
    default=0.07,
    show_default=True,
    help="Length scale of plotted polarization axes.",
)
@click.option(
    "--geometry",
    type=click.Choice(["unit-sphere", "physical"]),
    default="unit-sphere",
    show_default=True,
    help="Use a unit sphere or the natural physical acoustic geometry in 3D.",
)
@click.option(
    "--axis-labels",
    type=click.Choice(["cartesian", "crystal"]),
    default="cartesian",
    show_default=True,
    help="Use x/y/z labels or [100]/[010]/[001] tensor-frame labels.",
)
@click.option(
    "--complete-surface/--sampled-domain",
    default=True,
    show_default=True,
    help="Complete hemispherical 3D results by antipodal symmetry.",
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
@click.option("--show", is_flag=True, help="Display figures after rendering.")
def plot(
    filename: Path,
    outfile: Path | None,
    plot_2d: bool,
    plot_3d: bool,
    summary: bool,
    properties: tuple[str, ...],
    surface_types: tuple[str, ...],
    modes: tuple[str, ...],
    projection: str,
    cmap: str,
    levels: int,
    extrema: bool,
    polarization: bool,
    polarization_stride: int,
    polarization_color: str,
    polarization_linewidth: float,
    polarization_scale: float,
    axis_labels: str,
    geometry: str,
    complete_surface: bool,
    figure_preset: str,
    image_format: str,
    dpi: int | None,
    show: bool,
) -> None:
    """Render 2D maps, 3D surfaces, or a compact summary from HDF5."""
    if not (plot_2d or plot_3d or summary):
        summary = True
    base = filename.with_suffix("") if outfile is None else outfile.with_suffix("")
    render_options = MatplotlibOptions.from_preset(
        figure_preset,
        output_dir=base.parent,
        filename_prefix=f"{base.name}_",
        image_format=image_format,
        dpi=dpi,
        show=show,
        close=not show,
        axis_label_mode=cast(Literal["cartesian", "crystal"], axis_labels),
        polarization_color=polarization_color,
        polarization_line_width=polarization_linewidth,
        polarization_scale=polarization_scale,
        annotate_extrema=extrema,
        savefig_kwargs={"bbox_inches": "tight"},
    )
    try:
        result = read_seismic_hdf5(filename)
        base_map_options = SeismicPlotOptions(
            properties=properties or None,
            projection=cast(SphericalProjection, projection),
            colormap=cmap,
            levels=levels,
            polarization_stride=polarization_stride,
            include_extrema_markers=extrema,
            include_polarizations=False,
        )
        polarized_map_options = SeismicPlotOptions(
            properties=properties or None,
            projection=cast(SphericalProjection, projection),
            colormap=cmap,
            levels=levels,
            polarization_stride=polarization_stride,
            include_extrema_markers=extrema,
            include_polarizations=True,
        )
        collections: list[PlotCollection] = []
        if plot_2d:
            collections.append(build_seismic_plots(result, base_map_options))
            if polarization:
                polarized_maps = build_seismic_plots(result, polarized_map_options)
                polarized_maps.plots = [
                    spec
                    for spec in polarized_maps.plots
                    if isinstance(spec, SphericalMapSpec) and spec.axis_layers
                ]
                if polarized_maps.plots or polarized_maps.warnings:
                    collections.append(polarized_maps)
        if plot_3d:
            selected_modes = tuple(WaveMode(mode) for mode in modes) or (
                WaveMode.V_P,
                WaveMode.V_S1,
                WaveMode.V_S2,
            )
            selected_surfaces = tuple(cast(SurfaceType, item) for item in surface_types)
            surface_properties: tuple[str, ...] | None
            if properties:
                surface_properties = tuple(properties)
            elif selected_surfaces:
                surface_properties = ()
            else:
                surface_properties = None
            collections.append(
                build_seismic_surfaces(
                    result,
                    SeismicSurfaceOptions(
                        properties=surface_properties,
                        modes=selected_modes,
                        surface_types=selected_surfaces,
                        geometry=cast(SurfaceGeometry, geometry.replace("-", "_")),
                        colormap=cmap,
                        complete_antipodal_surface=complete_surface,
                        include_polarizations=False,
                        polarization_stride=polarization_stride,
                        polarization_color=polarization_color,
                        polarization_line_width=polarization_linewidth,
                        polarization_scale=polarization_scale,
                    ),
                )
            )
            if polarization:
                polarized_surfaces = build_seismic_surfaces(
                    result,
                    SeismicSurfaceOptions(
                        properties=surface_properties,
                        modes=selected_modes,
                        surface_types=selected_surfaces,
                        geometry=cast(SurfaceGeometry, geometry.replace("-", "_")),
                        colormap=cmap,
                        complete_antipodal_surface=complete_surface,
                        include_polarizations=True,
                        polarization_stride=polarization_stride,
                        polarization_color=polarization_color,
                        polarization_line_width=polarization_linewidth,
                        polarization_scale=polarization_scale,
                    ),
                )
                polarized_surfaces.plots = [
                    spec
                    for spec in polarized_surfaces.plots
                    if isinstance(spec, SurfacePlotSpec) and spec.vector_layers
                ]
                if polarized_surfaces.plots or polarized_surfaces.warnings:
                    collections.append(polarized_surfaces)
        paths: list[Path] = []
        for collection in collections:
            _render_cli_collection(collection, render_options, paths)
        if summary:
            base_summary = build_seismic_summary(result, base_map_options)
            _render_cli_spec(base_summary, render_options, paths)
            if polarization:
                polarized_summary = build_seismic_summary(
                    result,
                    polarized_map_options,
                )
                if isinstance(polarized_summary, SphericalSummarySpec) and any(
                    map_spec.axis_layers for map_spec in polarized_summary.maps
                ):
                    _render_cli_spec(polarized_summary, render_options, paths)
    except Exception as exc:
        echo_error(quantas_error(), bold=True)
        echo_error(str(exc))
        raise click.Abort() from exc

    for path in paths:
        echo(f"Figure written to: {path}")


def _render_cli_collection(
    collection: PlotCollection,
    options: MatplotlibOptions,
    paths: list[Path],
) -> None:
    """Render one collection while reporting each figure before generation."""
    for warning in collection.warnings:
        echo_warning(warning)
    for spec in collection.plots:
        if isinstance(spec, (SphericalMapSpec, SphericalSummarySpec, SurfacePlotSpec)):
            _render_cli_spec(spec, options, paths)
        else:
            raise TypeError(
                f"Unsupported seismic plot specification: {type(spec).__name__}"
            )


def _render_cli_spec(
    spec: SphericalMapSpec | SphericalSummarySpec | SurfacePlotSpec,
    options: MatplotlibOptions,
    paths: list[Path],
) -> None:
    """Render one seismic plot and report the target filename to the user."""
    suffix = options.image_format.lower().lstrip(".")
    prefix = options.filename_prefix
    echo(f"Generating figure: {prefix}{spec.filename_stem}.{suffix}")
    artifact = render_plot(spec, options)
    if artifact.path is not None:
        paths.append(artifact.path)


@seismic.command(name="export", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "-o",
    "--output",
    "outfile",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Output CSV file. Default: input base name + '.csv'.",
)
def export(filename: Path, outfile: Path | None) -> None:
    """Export sampled acoustic fields from a Quantas HDF5 result."""
    output = filename.with_suffix(".csv") if outfile is None else outfile
    try:
        result = read_seismic_hdf5(filename)
        output = write_seismic_csv(result, output)
    except Exception as exc:
        echo_error(quantas_error(), bold=True)
        echo_error(str(exc))
        raise click.Abort() from exc
    echo(f"SEISMIC data exported to: {output}")

apply_reference_help(seismic, ('seismic',))
