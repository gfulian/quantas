# -*- coding: utf-8 -*-

"""HDF5 and tabular exporters for elasticity results."""

from __future__ import annotations

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
from quantas.modules.elasticity.io.reader import read_elasticity_hdf5
from quantas.modules.elasticity.io.hdf5_payload import write_elasticity_payload


class ElasticityHDF5Export(BasicHDF5Export):
    """Export elasticity calculations to native Quantas HDF5."""

    def export(
        self,
        result: ResultData,
        filename: str | Path,
        report_text: str | None = None,
    ) -> None:
        """Write a complete elasticity result using the shared schema.

        Parameters
        ----------
        result : ResultData
            Generic result containing an ``elasticity`` payload.
        filename : str or Path
            Destination HDF5 file.
        report_text : str or None, optional
            Optional frontend-generated report text.
        """
        filename = ensure_suffix(filename, ".hdf5")
        elasticity = result.results["elasticity"]

        with h5py.File(filename, "w") as h5:
            write_result_metadata(h5, result)
            write_precision_metadata(h5)
            write_input_data(h5, result.input_data)
            write_options(h5, result.options)
            write_elasticity_payload(h5, elasticity)
            write_diagnostics(h5, result, report_text=report_text)
            write_events(h5, result.events)
        self.completed = True


class ElasticityTableExport(BasicExport):
    """Export principal-plane elasticity data as a neutral text table."""

    def export(self, result: ResultData, filename: str | Path) -> None:
        """Export 2D directional data from a generic Quantas result."""
        filename = ensure_suffix(filename, ".dat")
        elasticity = result.results["elasticity"]

        with filename.open("w", encoding="utf-8") as stream:
            stream.write("# Quantas elasticity 2D directional data export\n")
            stream.write(f"# Job: {elasticity.jobname}\n")
            stream.write(f"# Crystal system: {elasticity.crystal_system}\n")
            stream.write("#\n")
            if not elasticity.properties_2d:
                stream.write("# No 2D directional data available.\n")
            else:
                for plane, plane_data in elasticity.properties_2d.items():
                    self._write_plane(stream, plane, plane_data)

        self.completed = True

    def export_from_hdf5(
        self,
        input_file: str | Path,
        output_file: str | Path,
    ) -> None:
        """Export 2D directional data from a Quantas elasticity HDF5 file."""
        output_file = ensure_suffix(output_file, ".dat")
        result = read_elasticity_hdf5(input_file)
        self.export(result, output_file)

    def _write_plane(self, stream: Any, plane: str, plane_data: dict[str, Any]) -> None:
        """Write one principal-plane data block."""
        theta = plane_data.get("theta")
        phi = plane_data.get("phi")
        if theta is None or phi is None:
            return

        stream.write(f"\n# Plane: {plane}\n")
        columns = ["theta_deg", "phi_deg"]
        arrays = [np.degrees(np.asarray(theta)), np.degrees(np.asarray(phi))]

        for key, value in plane_data.items():
            if key in {"theta", "phi"}:
                continue
            array = np.asarray(value)
            if array.ndim == 1:
                columns.append(key)
                arrays.append(array)
            elif array.ndim == 2:
                for index in range(array.shape[1]):
                    columns.append(f"{key}_{index + 1}")
                    arrays.append(array[:, index])

        stream.write("# " + " ".join(f"{name:>22s}" for name in columns) + "\n")
        data = np.column_stack(arrays)
        np.savetxt(stream, data, fmt="%22.10e")
