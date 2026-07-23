# -*- coding: utf-8 -*-

"""Native HDF5 export for sampled seismic-wave results."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from quantas.io.hdf5 import (
    write_diagnostics,
    write_events,
    write_input_data,
    write_options,
    write_precision_metadata,
    write_result_metadata,
)
from quantas.io.path import ensure_suffix
from quantas.models import BasicExport, BasicHDF5Export, ResultData
from quantas.core.physics.seismic import MODE_INDEX, MODE_ORDER, MODE_SYMBOLS
from quantas.modules.seismic.io.hdf5_payload import (
    write_seismic_diagnostics,
    write_seismic_payload,
)
from quantas.modules.seismic.models import SeismicResult


class SeismicHDF5Export(BasicHDF5Export):
    """Export a complete seismic calculation to native Quantas HDF5."""

    def export(
        self,
        result: ResultData,
        filename: str | Path,
        report_text: str | None = None,
    ) -> None:
        """Write a structured seismic result.

        Parameters
        ----------
        result : ResultData
            Generic result containing a ``seismic`` payload.
        filename : str or Path
            Destination file.
        report_text : str or None, optional
            Optional frontend-generated report text.

        Raises
        ------
        ValueError
            If the generic result does not contain a seismic payload.
        """
        payload = result.results.get("seismic")
        if result.metadata.module != "seismic" or not isinstance(
            payload, SeismicResult
        ):
            raise ValueError("ResultData does not contain a valid seismic result.")

        output = ensure_suffix(filename, ".hdf5")
        with h5py.File(output, "w") as h5:
            write_result_metadata(h5, result)
            write_precision_metadata(h5)
            write_input_data(h5, result.input_data)
            write_options(h5, result.options)
            write_seismic_payload(h5, payload)
            diagnostics = write_diagnostics(h5, result, report_text=report_text)
            write_seismic_diagnostics(diagnostics, payload)
            write_events(h5, result.events)
        self.completed = True
        self.error = None


class SeismicTableExport(BasicExport):
    """Export sampled seismic fields as a named long-form CSV table."""

    def export(self, result: ResultData, filename: str | Path) -> None:
        """Write one row for every sampled direction and acoustic mode.

        Parameters
        ----------
        result : ResultData
            Generic result containing a ``seismic`` payload.
        filename : str or Path
            Destination CSV file.

        Raises
        ------
        ValueError
            If the result does not contain a seismic payload.
        """
        payload = result.results.get("seismic")
        if result.metadata.module != "seismic" or not isinstance(
            payload, SeismicResult
        ):
            raise ValueError("ResultData does not contain a valid seismic result.")

        output = ensure_suffix(filename, ".csv")
        columns = _seismic_table_columns(payload)
        with output.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(["# Quantas sampled seismic-wave export"])
            writer.writerow(["# Job", payload.jobname])
            writer.writerow(["# Density / kg m^-3", f"{payload.density:.12g}"])
            writer.writerow(["# Mode order", "V_S2, V_S1, V_P"])
            writer.writerow(columns)
            for row in _iter_seismic_table_rows(payload):
                writer.writerow(row)
        self.completed = True
        self.error = None

    def export_from_hdf5(
        self,
        input_file: str | Path,
        output_file: str | Path,
    ) -> None:
        """Export a native seismic HDF5 result as long-form CSV.

        Parameters
        ----------
        input_file : str or Path
            Native Quantas seismic HDF5 result.
        output_file : str or Path
            Destination CSV file.
        """
        from quantas.modules.seismic.io.reader import read_seismic_hdf5

        self.export(read_seismic_hdf5(input_file), output_file)


def _seismic_table_columns(result: SeismicResult) -> list[str]:
    """Return named columns available for one seismic result."""
    columns = [
        "point_index",
        "mode",
        "theta_rad",
        "theta_degree",
        "phi_rad",
        "phi_degree",
        "wave_normal_x",
        "wave_normal_y",
        "wave_normal_z",
        "eigenvalue_km2_s2",
        "phase_speed_km_s",
        "polarization_x",
        "polarization_y",
        "polarization_z",
        "eigenvalue_gap_km2_s2",
        "relative_eigenvalue_gap",
        "phase_valid",
        "phase_clamped",
        "phase_degenerate",
    ]
    if result.field.group is not None:
        columns.extend(
            [
                "group_velocity_x_km_s",
                "group_velocity_y_km_s",
                "group_velocity_z_km_s",
                "group_speed_km_s",
                "ray_direction_x",
                "ray_direction_y",
                "ray_direction_z",
                "power_flow_angle_rad",
                "power_flow_angle_degree",
                "group_valid",
                "group_resolved",
            ]
        )
    if result.field.enhancement is not None:
        columns.extend(
            [
                "area_factor",
                "log10_enhancement",
                "enhancement_valid",
                "enhancement_resolved",
                "enhancement_finite",
                "caustic_candidate",
            ]
        )
    if result.field.tracking is not None:
        columns.extend(
            [
                "tracked_polarization_x",
                "tracked_polarization_y",
                "tracked_polarization_z",
                "tracking_resolved",
            ]
        )
    return columns


def _iter_seismic_table_rows(result: SeismicResult) -> Any:
    """Yield long-form rows for every direction and acoustic mode."""
    grid = result.grid
    phase = result.field.phase
    theta = grid.theta_grid.reshape(-1)
    phi = grid.phi_grid.reshape(-1)
    tracked_axes, tracked_resolved = _tracked_mode_polarizations(result)

    for point in range(phase.n_points):
        for mode in MODE_ORDER:
            index = MODE_INDEX[mode]
            row: list[Any] = [
                point,
                MODE_SYMBOLS[mode],
                theta[point],
                np.degrees(theta[point]),
                phi[point],
                np.degrees(phi[point]),
                *phase.directions[point],
                phase.eigenvalues[point, index],
                phase.phase_speeds[point, index],
                *phase.polarizations[point, index],
                phase.mode_eigenvalue_gaps[point, index],
                phase.mode_relative_eigenvalue_gaps[point, index],
                bool(phase.valid_mask[point, index]),
                bool(phase.clamped_mask[point, index]),
                bool(phase.degeneracy_mask[point, index]),
            ]
            if result.field.group is not None:
                group = result.field.group
                row.extend(
                    [
                        *group.group_velocities[point, index],
                        group.group_speeds[point, index],
                        *group.ray_directions[point, index],
                        group.power_flow_angles[point, index],
                        np.degrees(group.power_flow_angles[point, index]),
                        bool(group.valid_mask[point, index]),
                        bool(group.resolved_mask[point, index]),
                    ]
                )
            if result.field.enhancement is not None:
                enhancement = result.field.enhancement
                row.extend(
                    [
                        enhancement.area_factors[point, index],
                        enhancement.log10_enhancement[point, index],
                        bool(enhancement.valid_mask[point, index]),
                        bool(enhancement.resolved_mask[point, index]),
                        bool(enhancement.finite_mask[point, index]),
                        bool(enhancement.caustic_candidate_mask[point, index]),
                    ]
                )
            if tracked_axes is not None and tracked_resolved is not None:
                row.extend(
                    [
                        *tracked_axes[point, index],
                        bool(tracked_resolved[point, index]),
                    ]
                )
            yield row


def _tracked_mode_polarizations(
    result: SeismicResult,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    """Map continuity branches back to local acoustic-mode identifiers."""
    tracking = result.field.tracking
    if tracking is None:
        return None, None
    n_points = result.field.n_points
    axes = np.full((n_points, 3, 3), np.nan, dtype=float)
    resolved = np.zeros((n_points, 3), dtype=bool)
    for branch in range(3):
        mode_indices = tracking.branch_mode_indices[:, branch]
        valid = (mode_indices >= 0) & tracking.resolved_mask[:, branch]
        points = np.flatnonzero(valid)
        axes[points, mode_indices[points]] = tracking.polarizations[points, branch]
        resolved[points, mode_indices[points]] = True
    return axes, resolved
