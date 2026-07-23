# -*- coding: utf-8 -*-

"""Complete tensor-export command for thermoelasticity."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import click
import numpy as np

from quantas.cli.contracts import OUTPUT_GROUP
from quantas.cli.grouped_options import GroupedCommand, grouped_option
from quantas.cli.output import CLIOutput
from quantas.cli.thermoelastic_common import (
    require_output_replacement,
    thermoelastic_payload,
)
from quantas.models import ReportTable
from quantas.api.thermoelasticity import (
    analyze_grid as analyze_thermoelastic_result,
    grid_info_table as thermoelastic_grid_info_table,
    read_result as read_thermoelastic_hdf5,
    regular_grid,
    write_tensor_export as write_thermoelastic_tensor_export,
)
from quantas.api.thermoelasticity import Result as ThermoelasticResult


@click.command(name="export", cls=GroupedCommand)
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@grouped_option(
    "--list-grid",
    group="Archive inspection",
    is_flag=True,
    default=False,
    help="Show available P-T ranges, steps, and archived profiles.",
)
@grouped_option(
    "--all",
    "export_all",
    group="Tensor selection",
    is_flag=True,
    default=False,
    help="Export the complete reconstructed grid and every archived profile.",
)
@grouped_option(
    "--profiles",
    "profiles_only",
    group="Tensor selection",
    is_flag=True,
    default=False,
    help="Export every archived depth profile, without the rectangular P-T grid.",
)
@grouped_option(
    "--point",
    group="Tensor selection",
    nargs=2,
    type=float,
    metavar="P_GPA T_K",
    default=None,
    help="Export only one reconstructed pressure-temperature state.",
)
@grouped_option(
    "-P",
    "--pressure",
    group="Tensor selection",
    nargs=3,
    type=float,
    metavar="MIN MAX STEP",
    default=None,
    help="Pressure range in GPa for a grid-only export.",
)
@grouped_option(
    "-T",
    "--temperature",
    group="Tensor selection",
    nargs=3,
    type=float,
    metavar="MIN MAX STEP",
    default=None,
    help="Temperature range in K for a grid-only export.",
)
@grouped_option(
    "--profile-name",
    "profile_names",
    group="Tensor selection",
    multiple=True,
    help="Export only the named archived profile. Repeat as needed.",
)
@grouped_option(
    "--tensor-condition",
    group="Tensor selection",
    type=click.Choice(["isothermal", "adiabatic", "both"], case_sensitive=False),
    default="isothermal",
    show_default=True,
    help="Export C^T, C^S, or both tensor conditions.",
)
@grouped_option(
    "-o",
    "--output",
    group=OUTPUT_GROUP,
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Export destination; a selection-specific name is used by default.",
)
@grouped_option(
    "--format",
    "file_format",
    group=OUTPUT_GROUP,
    type=click.Choice(["txt", "csv"], case_sensitive=False),
    default="txt",
    show_default=True,
    help="Tensor export format.",
)
@grouped_option(
    "--uncertainties/--no-uncertainties",
    group=OUTPUT_GROUP,
    default=True,
    show_default=True,
    help="Include one-sigma tensor uncertainties.",
)
@grouped_option(
    "-f",
    "--force",
    group=OUTPUT_GROUP,
    is_flag=True,
    default=False,
    help="Replace an existing export file.",
)
@grouped_option(
    "--extrapolation",
    "extrapolation_policy",
    group="Validation policy",
    type=click.Choice(["fail", "warn", "allow"], case_sensitive=False),
    default="warn",
    show_default=True,
    help="Policy for custom states outside the archived QHA or elastic range.",
)
def export(
    archive: Path,
    list_grid: bool,
    export_all: bool,
    profiles_only: bool,
    point: tuple[float, float] | None,
    pressure: tuple[float, float, float] | None,
    temperature: tuple[float, float, float] | None,
    profile_names: tuple[str, ...],
    tensor_condition: str,
    output: Path | None,
    file_format: str,
    uncertainties: bool,
    force: bool,
    extrapolation_policy: str,
) -> None:
    """Inspect coverage and export an explicit tensor selection."""
    terminal = CLIOutput(show_progress=False)
    try:
        result_data = read_thermoelastic_hdf5(archive)
        payload = thermoelastic_payload(result_data)
        if list_grid:
            terminal.table(thermoelastic_grid_info_table(payload))
            if payload.profiles:
                terminal.table(_archived_profiles_table(payload))

        mode = _resolve_export_mode(
            export_all=export_all,
            profiles_only=profiles_only,
            point=point,
            pressure=pressure,
            temperature=temperature,
            profile_names=profile_names,
        )
        if mode is None:
            if list_grid:
                terminal.save()
                return
            raise click.UsageError(
                "choose one export selection: --all, --point, -P/-T, "
                "--profiles, or --profile-name"
            )

        include_grid = mode in {"all", "point", "grid"}
        include_profiles = mode in {"all", "profiles"}
        if mode in {"point", "grid"}:
            target_pressure, target_temperature = _resolve_pt_selection(
                point=point,
                pressure=pressure,
                temperature=temperature,
                archive=archive,
            )
            result_data = analyze_thermoelastic_result(
                result_data,
                temperature=target_temperature,
                pressure=target_pressure,
                extrapolation_policy=extrapolation_policy,
            )
            payload = thermoelastic_payload(result_data)
            payload.profiles = {}
        elif mode == "all" and payload.stiffness_isothermal is None:
            archived_profiles = dict(payload.profiles)
            result_data = analyze_thermoelastic_result(
                result_data,
                temperature=payload.temperature,
                pressure=payload.pressure,
                extrapolation_policy=extrapolation_policy,
            )
            payload = thermoelastic_payload(result_data)
            payload.profiles = archived_profiles
        elif mode == "profiles":
            if not payload.profiles:
                raise click.UsageError("archive does not contain depth profiles")
            if profile_names:
                missing = sorted(set(profile_names) - set(payload.profiles))
                if missing:
                    raise click.UsageError(
                        "profile(s) not present in archive: " + ", ".join(missing)
                    )
                payload.profiles = {
                    name: payload.profiles[name] for name in profile_names
                }

        suffix = ".csv" if file_format.lower() == "csv" else ".dat"
        destination = (
            archive.with_name(f"{archive.stem}_{mode}_tensors{suffix}")
            if output is None
            else output
        )
        require_output_replacement(destination, force)
        written = write_thermoelastic_tensor_export(
            result_data,
            destination,
            file_format=file_format,
            tensor_condition=cast(
                Literal["isothermal", "adiabatic", "both"],
                tensor_condition.lower(),
            ),
            include_uncertainties=uncertainties,
            include_grid=include_grid,
            include_profiles=include_profiles,
            overwrite=True,
        )
        terminal.message(f"Thermoelastic tensor export written to: {written}")
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


def _resolve_pt_selection(
    *,
    point: tuple[float, float] | None,
    pressure: tuple[float, float, float] | None,
    temperature: tuple[float, float, float] | None,
    archive: Path,
) -> tuple[np.ndarray, np.ndarray]:
    if point is not None and (pressure is not None or temperature is not None):
        raise click.UsageError(
            "--point cannot be combined with pressure or temperature ranges"
        )
    source = thermoelastic_payload(read_thermoelastic_hdf5(archive))
    if point is not None:
        return (
            np.asarray([point[0]], dtype=np.float64),
            np.asarray([point[1]], dtype=np.float64),
        )
    pressure_values = (
        source.pressure
        if pressure is None
        else regular_grid(pressure[0], pressure[1], pressure[2])
    )
    temperature_values = (
        source.temperature
        if temperature is None
        else regular_grid(temperature[0], temperature[1], temperature[2])
    )
    return np.asarray(pressure_values), np.asarray(temperature_values)


def _resolve_export_mode(
    *,
    export_all: bool,
    profiles_only: bool,
    point: tuple[float, float] | None,
    pressure: tuple[float, float, float] | None,
    temperature: tuple[float, float, float] | None,
    profile_names: tuple[str, ...],
) -> Literal["all", "point", "grid", "profiles"] | None:
    """Validate mutually exclusive export selectors and return the mode."""
    custom_grid = pressure is not None or temperature is not None
    selected = sum(
        (
            int(export_all),
            int(profiles_only),
            int(point is not None),
            int(custom_grid),
            int(bool(profile_names)),
        )
    )
    if selected == 0:
        return None
    if selected > 1:
        raise click.UsageError(
            "export selectors are mutually exclusive: choose only one of "
            "--all, --point, -P/-T, --profiles, or --profile-name"
        )
    if export_all:
        return "all"
    if profiles_only or profile_names:
        return "profiles"
    if point is not None:
        return "point"
    return "grid"


def _archived_profiles_table(result: ThermoelasticResult) -> ReportTable:
    """Return a compact list of depth profiles stored in one archive."""
    rows = [
        [
            name,
            profile.depth.size,
            float(profile.depth[0]),
            float(profile.depth[-1]),
            float(np.min(profile.pressure)),
            float(np.max(profile.pressure)),
            float(np.min(profile.temperature)),
            float(np.max(profile.temperature)),
        ]
        for name, profile in result.profiles.items()
    ]
    return ReportTable(
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
    )


__all__ = ["export"]
