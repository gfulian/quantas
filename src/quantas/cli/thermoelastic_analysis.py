# -*- coding: utf-8 -*-

"""Simplified thermoelastic point, grid, and geological-profile analysis CLI."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import click
import numpy as np

from quantas.cli.contracts import (
    DOMAIN_GROUP,
    OUTPUT_GROUP,
    PLOTTING_GROUP,
    SCIENTIFIC_GROUP,
    VALIDATION_GROUP,
)
from quantas.cli.grouped_options import GroupedCommand, grouped_option
from quantas.cli.messages import quantas_finish, quantas_title
from quantas.cli.output import CLIOutput
from quantas.cli.thermoelastic_common import (
    approve_output_replacement,
    thermoelastic_payload,
)
from quantas.cli.thermoelastic_plot_common import render_collection
from quantas.models import ReportTable
from quantas.api.thermoelasticity import (
    DepthProfile as ThermoelasticDepthProfile,
    ProfilePlotOptions as ThermoelasticProfilePlotOptions,
    PlotStyleOptions as ThermoelasticPlotStyleOptions,
    analyze_profiles as analyze_thermoelastic_profiles,
    analyze_grid as analyze_thermoelastic_result,
    build_profile_plots as build_thermoelastic_profile_plots,
    build_profile_preset as build_thermoelastic_profile_preset,
    point_table as thermoelastic_point_table,
    profile_presets as thermoelastic_profile_presets,
    read_depth_profile as read_thermoelastic_depth_profile,
    read_profile_spec as read_thermoelastic_profile_spec,
    read_result as read_thermoelastic_hdf5,
    regular_grid,
    write_grid_table as write_thermoelastic_grid_table,
    write_profile_table as write_thermoelastic_profile_table,
    write_result as write_thermoelastic_hdf5,
    write_state_input as write_thermoelastic_state_input,
)

_PROFILE_PRESET_NAMES = tuple(preset.name for preset in thermoelastic_profile_presets())


@click.group(name="analysis")
def analysis() -> None:
    """Evaluate a calibrated thermoelastic model at points, grids, or profiles."""


@analysis.command(name="point", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("pressure", type=float)
@click.argument("temperature", type=float)
@grouped_option(
    "--tensor-condition",
    group=SCIENTIFIC_GROUP,
    type=click.Choice(["isothermal", "adiabatic"], case_sensitive=False),
    default="adiabatic",
    show_default=True,
    help="Stiffness condition shown and written to an optional reusable input.",
)
@grouped_option(
    "--extrapolation",
    group=VALIDATION_GROUP,
    type=click.Choice(["fail", "warn", "allow"], case_sensitive=False),
    default="fail",
    show_default=True,
    help="Policy for states outside QHA or elastic-volume support.",
)
@grouped_option(
    "-o",
    "--output",
    group=OUTPUT_GROUP,
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Optional shared Elasticity/SEISMIC text input containing the selected "
        "stiffness tensor and density."
    ),
)
@grouped_option(
    "-f",
    "--force",
    group=OUTPUT_GROUP,
    is_flag=True,
    default=False,
    help="Replace OUTPUT without prompting.",
)
def point_analysis(
    archive: Path,
    pressure: float,
    temperature: float,
    tensor_condition: str,
    extrapolation: str,
    output: Path | None,
    force: bool,
) -> None:
    """Evaluate one state and optionally write an Elasticity/SEISMIC input."""
    result_data = read_thermoelastic_hdf5(archive)
    analyzed = analyze_thermoelastic_result(
        result_data,
        pressure=np.asarray([pressure], dtype=np.float64),
        temperature=np.asarray([temperature], dtype=np.float64),
        extrapolation_policy=extrapolation.lower(),
    )
    payload = thermoelastic_payload(analyzed)
    terminal = CLIOutput(show_progress=False)
    try:
        terminal.message(quantas_title(), bold=True)
        terminal.table(
            thermoelastic_point_table(
                payload,
                tensor_condition=cast(
                    Literal["isothermal", "adiabatic"],
                    tensor_condition.lower(),
                ),
            )
        )
        if output is not None:
            destination = output.with_suffix(output.suffix or ".dat")
            if approve_output_replacement(destination, force):
                write_thermoelastic_state_input(
                    payload,
                    destination,
                    tensor_condition=cast(
                        Literal["isothermal", "adiabatic"],
                        tensor_condition.lower(),
                    ),
                    overwrite=True,
                )
                terminal.message(f"Elasticity/SEISMIC input written to: {destination}")
        terminal.message(quantas_finish())
        terminal.save()
    finally:
        terminal.close()


@analysis.command(name="grid", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@grouped_option(
    "-P",
    "--pressure",
    group=DOMAIN_GROUP,
    nargs=3,
    type=float,
    metavar="MIN MAX STEP",
    required=True,
    help="Pressure grid in GPa.",
)
@grouped_option(
    "-T",
    "--temperature",
    group=DOMAIN_GROUP,
    nargs=3,
    type=float,
    metavar="MIN MAX STEP",
    required=True,
    help="Temperature grid in K.",
)
@grouped_option(
    "--extrapolation",
    group=VALIDATION_GROUP,
    type=click.Choice(["fail", "warn", "allow"], case_sensitive=False),
    default="warn",
    show_default=True,
    help="Policy outside QHA or elastic-volume support.",
)
@grouped_option(
    "-o",
    "--output",
    "outfile",
    group=OUTPUT_GROUP,
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Analysis HDF5. Default: fit archive base name + '_GRID.hdf5'.",
)
@grouped_option(
    "--table-output",
    group=OUTPUT_GROUP,
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Optional wide CSV/text table with one row per P-T state.",
)
@grouped_option(
    "--table-format",
    group=OUTPUT_GROUP,
    type=click.Choice(["csv", "text"], case_sensitive=False),
    default="csv",
    show_default=True,
    help="Format of the optional wide table output.",
)
@grouped_option(
    "--tensor-condition",
    group=SCIENTIFIC_GROUP,
    type=click.Choice(["isothermal", "adiabatic", "both"], case_sensitive=False),
    default="both",
    show_default=True,
    help="Select isothermal, adiabatic, or both stiffness tensors.",
)
@grouped_option(
    "--uncertainties/--no-uncertainties",
    group=VALIDATION_GROUP,
    default=True,
    show_default=True,
    help="Include or omit one-sigma tensor uncertainties in the table.",
)
@grouped_option(
    "-f",
    "--force",
    group=OUTPUT_GROUP,
    is_flag=True,
    default=False,
    help="Replace existing outputs without prompting.",
)
def grid_analysis(
    archive: Path,
    pressure: tuple[float, float, float],
    temperature: tuple[float, float, float],
    extrapolation: str,
    outfile: Path | None,
    table_output: Path | None,
    table_format: str,
    tensor_condition: str,
    uncertainties: bool,
    force: bool,
) -> None:
    """Evaluate a rectangular P-T grid and save an analysis HDF5 archive."""
    destination = (
        archive.with_name(f"{archive.stem}_GRID.hdf5")
        if outfile is None
        else outfile.with_suffix(".hdf5")
    )
    if not approve_output_replacement(destination, force):
        click.echo("Results not saved")
        return
    source = read_thermoelastic_hdf5(archive)
    analyzed = analyze_thermoelastic_result(
        source,
        pressure=regular_grid(*pressure),
        temperature=regular_grid(*temperature),
        extrapolation_policy=extrapolation.lower(),
    )
    write_thermoelastic_hdf5(analyzed, destination)
    payload = thermoelastic_payload(analyzed)
    terminal = CLIOutput(show_progress=False)
    try:
        terminal.message(f"Thermoelastic grid archive written to: {destination}")
        if table_output is not None:
            table_path = _normalized_table_path(table_output, table_format)
            if approve_output_replacement(table_path, force):
                write_thermoelastic_grid_table(
                    payload,
                    table_path,
                    tensor_condition=cast(
                        Literal["isothermal", "adiabatic", "both"],
                        tensor_condition.lower(),
                    ),
                    file_format=cast(Literal["csv", "text"], table_format.lower()),
                    include_uncertainties=uncertainties,
                    overwrite=True,
                )
                terminal.message(f"Wide grid table written to: {table_path}")
        terminal.message(
            f"Next: quantas thermoelasticity plot pt {destination} --component C11"
        )
        terminal.save()
    finally:
        terminal.close()


@analysis.command(name="profile", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@grouped_option(
    "--profile",
    "profile_file",
    group="Profile definition",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Complete tabulated depth/P/T profile.",
)
@grouped_option(
    "--profile-spec",
    group="Profile definition",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="YAML profile composing independent pressure and temperature models.",
)
@grouped_option(
    "--preset",
    "profile_preset",
    group="Profile definition",
    type=click.Choice(_PROFILE_PRESET_NAMES, case_sensitive=False),
    default=None,
    help="Built-in geological profile preset.",
)
@grouped_option(
    "--list-presets",
    group="Profile definition",
    is_flag=True,
    default=False,
    help="List built-in profiles and exit.",
)
@grouped_option(
    "--linear-profile",
    group="Synthetic profile",
    default=None,
    metavar="NAME",
    help="Create a controlled linear profile for testing.",
)
@grouped_option(
    "--depth",
    group="Synthetic profile",
    nargs=3,
    type=float,
    metavar="MIN MAX STEP",
    default=None,
    help="Synthetic-profile depth grid in km as MIN MAX STEP.",
)
@grouped_option(
    "--pressure-at-depth-min",
    group="Synthetic profile",
    type=float,
    default=0.0,
    show_default=True,
    help="Pressure in GPa at the minimum synthetic-profile depth.",
)
@grouped_option(
    "--pressure-gradient",
    group="Synthetic profile",
    type=float,
    default=0.03,
    show_default=True,
    help="GPa km^-1.",
)
@grouped_option(
    "--temperature-at-depth-min",
    group="Synthetic profile",
    type=float,
    default=298.15,
    show_default=True,
    help="Temperature in K at the minimum synthetic-profile depth.",
)
@grouped_option(
    "--temperature-gradient",
    group="Synthetic profile",
    type=float,
    default=0.5,
    show_default=True,
    help="K km^-1.",
)
@grouped_option(
    "--extrapolation",
    group=VALIDATION_GROUP,
    type=click.Choice(["fail", "warn", "allow"], case_sensitive=False),
    default="warn",
    show_default=True,
    help="Policy for states outside QHA or elastic-volume support.",
)
@grouped_option(
    "-o",
    "--output",
    "outfile",
    group=OUTPUT_GROUP,
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Profile HDF5. Default is derived from fit archive and profile name.",
)
@grouped_option(
    "--table-output",
    group=OUTPUT_GROUP,
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Wide profile table. Default: HDF5 base name + '.csv' or '.dat'.",
)
@grouped_option(
    "--table-format",
    group=OUTPUT_GROUP,
    type=click.Choice(["csv", "text"], case_sensitive=False),
    default="csv",
    show_default=True,
    help="Format of the optional wide table output.",
)
@grouped_option(
    "--tensor-condition",
    group=SCIENTIFIC_GROUP,
    type=click.Choice(["isothermal", "adiabatic", "both"], case_sensitive=False),
    default="both",
    show_default=True,
    help="Select isothermal, adiabatic, or both stiffness tensors.",
)
@grouped_option(
    "--uncertainties/--no-uncertainties",
    group=VALIDATION_GROUP,
    default=True,
    show_default=True,
    help="Include or omit one-sigma tensor uncertainties in the table.",
)
@grouped_option(
    "--plot/--no-plot",
    group=PLOTTING_GROUP,
    default=False,
    show_default=True,
    help="Generate a standard diagnostic profile plot after evaluation.",
)
@grouped_option(
    "-f",
    "--force",
    group=OUTPUT_GROUP,
    is_flag=True,
    default=False,
    help="Replace existing outputs without prompting.",
)
def profile_analysis(
    archive: Path,
    profile_file: Path | None,
    profile_spec: Path | None,
    profile_preset: str | None,
    list_presets: bool,
    linear_profile: str | None,
    depth: tuple[float, float, float] | None,
    pressure_at_depth_min: float,
    pressure_gradient: float,
    temperature_at_depth_min: float,
    temperature_gradient: float,
    extrapolation: str,
    outfile: Path | None,
    table_output: Path | None,
    table_format: str,
    tensor_condition: str,
    uncertainties: bool,
    plot: bool,
    force: bool,
) -> None:
    """Evaluate one geological profile, save HDF5, and write a wide table."""
    if list_presets:
        _show_presets()
        return
    profile = _resolve_one_profile(
        profile_file=profile_file,
        profile_spec=profile_spec,
        profile_preset=profile_preset,
        linear_profile=linear_profile,
        depth=depth,
        pressure_at_depth_min=pressure_at_depth_min,
        pressure_gradient=pressure_gradient,
        temperature_at_depth_min=temperature_at_depth_min,
        temperature_gradient=temperature_gradient,
    )
    destination = (
        archive.with_name(f"{archive.stem}_{profile.name}.hdf5")
        if outfile is None
        else outfile.with_suffix(".hdf5")
    )
    table_path = (
        destination.with_suffix(".csv" if table_format.lower() == "csv" else ".dat")
        if table_output is None
        else _normalized_table_path(table_output, table_format)
    )
    for path in (destination, table_path):
        if not approve_output_replacement(path, force):
            click.echo("Results not saved")
            return
    source = read_thermoelastic_hdf5(archive)
    analyzed = analyze_thermoelastic_profiles(
        source,
        (profile,),
        extrapolation_policy=extrapolation.lower(),
    )
    write_thermoelastic_hdf5(analyzed, destination)
    payload = thermoelastic_payload(analyzed)
    write_thermoelastic_profile_table(
        payload,
        table_path,
        profile_name=profile.name,
        tensor_condition=cast(
            Literal["isothermal", "adiabatic", "both"],
            tensor_condition.lower(),
        ),
        file_format=cast(Literal["csv", "text"], table_format.lower()),
        include_uncertainties=uncertainties,
        overwrite=True,
    )
    terminal = CLIOutput(show_progress=False)
    try:
        terminal.message(f"Thermoelastic profile archive written to: {destination}")
        terminal.message(f"Wide profile table written to: {table_path}")
        if plot:
            collection = build_thermoelastic_profile_plots(
                analyzed,
                profile_name=profile.name,
                options=ThermoelasticProfilePlotOptions(
                    style=ThermoelasticPlotStyleOptions(preset="analysis"),
                    tensor_condition=(
                        "adiabatic"
                        if payload.profiles[profile.name].stiffness_adiabatic
                        is not None
                        else "isothermal"
                    ),
                    layout="facets",
                ),
            )
            render_collection(
                destination,
                collection,
                family="profile",
                output_dir=destination.with_name(destination.stem + "_plots"),
                image_format="png",
                preset="screen",
                dpi=150,
                transparent=False,
                show=False,
                figure_size=(7.0, 7.0),
                axis_label_font_size=14.0,
                legend_font_size=11.0,
                title_font_size=14.0,
                tick_label_font_size=11.0,
            )
        terminal.message(
            "The CSV/text table is the recommended interface for custom figures."
        )
        terminal.save()
    finally:
        terminal.close()


def _resolve_one_profile(
    *,
    profile_file: Path | None,
    profile_spec: Path | None,
    profile_preset: str | None,
    linear_profile: str | None,
    depth: tuple[float, float, float] | None,
    pressure_at_depth_min: float,
    pressure_gradient: float,
    temperature_at_depth_min: float,
    temperature_gradient: float,
) -> ThermoelasticDepthProfile:
    selected = sum(
        value is not None
        for value in (profile_file, profile_spec, profile_preset, linear_profile)
    )
    if selected != 1:
        raise click.UsageError(
            "choose exactly one profile source: --profile, --profile-spec, "
            "--preset, or --linear-profile"
        )
    if profile_file is not None:
        return read_thermoelastic_depth_profile(profile_file)
    if profile_spec is not None:
        return read_thermoelastic_profile_spec(profile_spec)
    if profile_preset is not None:
        return build_thermoelastic_profile_preset(profile_preset)
    if depth is None:
        raise click.UsageError("--linear-profile requires --depth MIN MAX STEP")
    values = regular_grid(*depth)
    offset = values - values[0]
    return ThermoelasticDepthProfile(
        name=cast(str, linear_profile),
        depth=values,
        pressure=pressure_at_depth_min + pressure_gradient * offset,
        temperature=temperature_at_depth_min + temperature_gradient * offset,
        metadata={
            "kind": "linear_gradients",
            "pressure_gradient_GPa_per_km": pressure_gradient,
            "temperature_gradient_K_per_km": temperature_gradient,
        },
    )


def _show_presets() -> None:
    rows = [
        [
            preset.name,
            f"{preset.depth_min:g}-{preset.depth_max:g}",
            f"{preset.pressure_model} + {preset.temperature_model}",
        ]
        for preset in thermoelastic_profile_presets()
    ]
    output = CLIOutput(show_progress=False)
    try:
        output.table(
            ReportTable(
                "Built-in thermoelastic profiles",
                ["Preset", "Depth (km)", "Models"],
                rows,
            )
        )
        output.save()
    finally:
        output.close()


def _normalized_table_path(path: Path, file_format: str) -> Path:
    suffix = ".csv" if file_format.lower() == "csv" else ".dat"
    return (
        path.with_suffix(suffix)
        if path.suffix.lower() not in {".csv", ".dat"}
        else path
    )


__all__ = ["analysis"]
