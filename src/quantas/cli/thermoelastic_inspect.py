# -*- coding: utf-8 -*-

"""User-oriented inspection of thermoelastic HDF5 archives."""

from __future__ import annotations

from pathlib import Path

import click
import numpy as np

from quantas.cli.output import CLIOutput
from quantas.models import ReportTable
from quantas.api.thermoelasticity import read_result as read_thermoelastic_hdf5
from quantas.api.thermoelasticity import Result as ThermoelasticResult


@click.command(name="inspect")
@click.argument("archive", type=click.Path(exists=True, dir_okay=False, path_type=Path))
def inspect_archive(archive: Path) -> None:
    """Explain archive stage, scientific coverage, and useful next commands."""
    result_data = read_thermoelastic_hdf5(archive)
    payload = result_data.results.get("thermoelasticity")
    if not isinstance(payload, ThermoelasticResult):
        raise click.ClickException("archive is not a thermoelastic result")
    stage = _stage(payload)
    available = ["isothermal"]
    if payload.stiffness_adiabatic is not None or any(
        profile.stiffness_adiabatic is not None for profile in payload.profiles.values()
    ):
        available.append("adiabatic")
    elif (
        payload.isochoric_heat_capacity_cell is not None
        and payload.thermal_expansion_tensor is not None
    ):
        available.append("adiabatic on analysis")
    quality_levels: dict[str, int] = {}
    for component in payload.component_fits.values():
        level = "unclassified" if component.quality is None else component.quality.level
        quality_levels[level] = quality_levels.get(level, 0) + 1
    stable, unstable, indeterminate = _stability_counts(payload)
    output = CLIOutput(show_progress=False)
    try:
        output.table(
            ReportTable(
                "Thermoelastic archive summary",
                ["Property", "Value"],
                [
                    ["Archive", str(archive)],
                    ["Stage", stage],
                    ["Job", payload.jobname],
                    ["Reference EOS", payload.reference_eos.eos],
                    ["Independent Cij", " ".join(payload.independent_labels)],
                    ["Tensor conditions", ", ".join(available)],
                    ["Frame status", _frame_status(payload)],
                    ["P support (GPa)", _range(payload.pressure)],
                    ["T support (K)", _range(payload.temperature)],
                    ["Profiles", ", ".join(payload.profiles) or "none"],
                    ["Fit quality", _mapping_text(quality_levels)],
                    [
                        "Stable / unstable / indeterminate",
                        f"{stable} / {unstable} / {indeterminate}",
                    ],
                ],
            )
        )
        if payload.profiles:
            rows = []
            for name, profile in payload.profiles.items():
                rows.append(
                    [
                        name,
                        profile.depth.size,
                        f"{profile.depth[0]:g}-{profile.depth[-1]:g}",
                        f"{profile.pressure.min():g}-{profile.pressure.max():g}",
                        f"{profile.temperature.min():g}-{profile.temperature.max():g}",
                        int(np.count_nonzero(profile.qha_extrapolation_mask)),
                        int(np.count_nonzero(profile.elastic_extrapolation_mask)),
                    ]
                )
            output.table(
                ReportTable(
                    "Archived profiles",
                    [
                        "Profile",
                        "Points",
                        "Depth (km)",
                        "P (GPa)",
                        "T (K)",
                        "QHA ext.",
                        "Elastic ext.",
                    ],
                    rows,
                )
            )
        output.table(_next_steps(archive, stage, payload))
        output.save()
    finally:
        output.close()


def _frame_status(payload: ThermoelasticResult) -> str:
    """Return the archived tensor-frame normalization status."""
    value = payload.metadata.get("frame_normalization")
    if isinstance(value, dict):
        status = value.get("status")
        if isinstance(status, str) and status:
            return status
    legacy = payload.metadata.get("frame_status")
    return legacy if isinstance(legacy, str) and legacy else "not recorded"


def _stage(payload: ThermoelasticResult) -> str:
    has_grid = payload.stiffness_isothermal is not None
    has_profiles = bool(payload.profiles)
    if has_grid and has_profiles:
        return "combined grid and profile analysis"
    if has_grid:
        return "pressure-temperature grid analysis"
    if has_profiles:
        return "geological profile analysis"
    return "model calibration"


def _stability_counts(payload: ThermoelasticResult) -> tuple[int, int, int]:
    stable = unstable = indeterminate = 0
    fields = []
    if payload.stability is not None:
        fields.append(payload.stability)
    fields.extend(
        profile.stability
        for profile in payload.profiles.values()
        if profile.stability is not None
    )
    for field in fields:
        stable += int(np.count_nonzero(field.stable_mask))
        unstable += int(np.count_nonzero(field.unstable_mask))
        indeterminate += int(np.count_nonzero(field.indeterminate_mask))
    return stable, unstable, indeterminate


def _next_steps(
    archive: Path,
    stage: str,
    payload: ThermoelasticResult,
) -> ReportTable:
    rows: list[list[str]] = []
    if stage == "model calibration":
        rows.extend(
            [
                ["Inspect fit figures", f"quantas thermoelasticity plot fit {archive}"],
                [
                    "Evaluate one state",
                    f"quantas thermoelasticity analysis point {archive} 0 300",
                ],
                [
                    "Evaluate a P-T grid",
                    "quantas thermoelasticity analysis grid "
                    f"{archive} -P 0 10 1 -T 300 1500 100",
                ],
                [
                    "Evaluate a profile",
                    "quantas thermoelasticity analysis profile "
                    f"{archive} --preset mantle-katsura-2022",
                ],
            ]
        )
    elif stage == "pressure-temperature grid analysis":
        rows.extend(
            [
                [
                    "Plot maps",
                    f"quantas thermoelasticity plot pt {archive} --component C11",
                ],
                [
                    "Compare C^T and C^S",
                    "quantas thermoelasticity plot compare "
                    f"{archive} --component C11 --pressure 0",
                ],
                [
                    "Export wide table",
                    f"quantas thermoelasticity table grid {archive} -o grid.csv",
                ],
            ]
        )
    else:
        name = next(iter(payload.profiles), "PROFILE")
        rows.extend(
            [
                [
                    "Plot profile",
                    "quantas thermoelasticity plot profile "
                    f"{archive} --profile-name {name}",
                ],
                [
                    "Export wide table",
                    "quantas thermoelasticity table profile "
                    f"{archive} --profile-name {name} -o profile.csv",
                ],
            ]
        )
    return ReportTable("Common next commands", ["Purpose", "Command"], rows)


def _range(values: np.ndarray) -> str:
    array = np.asarray(values, dtype=np.float64)
    return "none" if array.size == 0 else f"{array.min():g}-{array.max():g}"


def _mapping_text(values: dict[str, int]) -> str:
    return (
        ", ".join(f"{key}={value}" for key, value in sorted(values.items())) or "none"
    )


__all__ = ["inspect_archive"]
