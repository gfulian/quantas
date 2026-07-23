# -*- coding: utf-8 -*-

"""Readers for Quantas elasticity input and HDF5 result files."""

from __future__ import annotations

from pathlib import Path
from typing import TypedDict

import h5py
import numpy as np

from quantas.core.physics.elasticity import validate_stiffness_matrix
from quantas.io.hdf5 import (
    read_events,
    read_input_data,
    read_options,
    read_result_metadata,
    read_warnings,
)
from quantas.models import BasicReader, ResultData, validate_result_schema
from quantas.modules.elasticity.io.hdf5_payload import read_elasticity_payload


class _ElasticityInputData(TypedDict):
    """Typed payload for the plain-text elasticity reader."""

    job: str
    stiffness: np.ndarray | None


class ElasticityInputFileReader(BasicReader[None]):
    """Read a Quantas second-order elasticity text input file."""

    def __init__(self, elasticity_input: str | Path | None = None) -> None:
        super().__init__()
        self._data: _ElasticityInputData = {
            "job": "Unknown",
            "stiffness": None,
        }
        if elasticity_input is not None:
            self.load(elasticity_input)

    @property
    def jobname(self) -> str:
        """Return the input job description."""
        return self._data["job"]

    @property
    def stiffness(self) -> np.ndarray:
        """Return the stiffness matrix in GPa."""
        stiffness = self._data["stiffness"]
        if stiffness is None:
            raise RuntimeError("elastic stiffness data are not available")
        return stiffness

    def load(self, filename: str | Path) -> None:
        """Parse a Quantas elasticity input file."""
        filename = Path(filename)
        self.completed = False
        self.error = None
        self._data = {
            "job": "Unknown",
            "stiffness": None,
        }
        with filename.open("r", encoding="utf-8") as stream:
            lines = stream.readlines()
        if len(lines) < 6:
            self.error = "File length is too short to contain the elastic matrix."
            return

        start = self._set_jobname(lines[0])
        matrix = self._read_stiffness_matrix(lines, start)
        if matrix is None:
            return
        self._data["stiffness"] = matrix.copy()
        if self.error is None:
            self.completed = True

    def _set_jobname(self, line: str) -> int:
        """Read an optional job name and return the matrix start line."""
        try:
            float(line.strip().split()[0])
        except (ValueError, IndexError):
            self._data["job"] = line.rstrip()
            return 1
        return 0

    def _read_stiffness_matrix(
        self,
        lines: list[str],
        start: int,
    ) -> np.ndarray | None:
        """Read a full or triangular ``6 x 6`` stiffness matrix."""
        try:
            matrix_lines = [lines[index + start] for index in range(6)]
        except IndexError:
            self.error = (
                "The elastic matrix does not have the expected (6, 6) shape. "
                "Please check the input file and retry."
            )
            return None

        try:
            rows = [
                [float(value) for value in matrix_line.split()]
                for matrix_line in matrix_lines
            ]
        except ValueError:
            self.error = "Elastic constants must be valid floating-point numbers."
            return None

        row_lengths = [len(row) for row in rows]
        matrix = self._assemble_matrix(rows, row_lengths)
        if matrix is None:
            self.error = (
                "The elastic matrix must be full, upper triangular, "
                "or lower triangular."
            )
            return None
        try:
            validate_stiffness_matrix(matrix, copy=False)
        except ValueError:
            self.error = "The elastic matrix must be symmetric or triangular."
            return None
        return matrix

    def _assemble_matrix(
        self,
        rows: list[list[float]],
        row_lengths: list[int],
    ) -> np.ndarray | None:
        """Assemble a full, upper-triangular, or lower-triangular matrix."""
        matrix = np.zeros((6, 6), dtype=float)
        if row_lengths == [6] * 6:
            matrix[:, :] = np.asarray(rows, dtype=float)
            return matrix
        if row_lengths == [6, 5, 4, 3, 2, 1]:
            for index, row in enumerate(rows):
                matrix[index, index:] = row
            return matrix + np.triu(matrix, 1).T
        if row_lengths == [1, 2, 3, 4, 5, 6]:
            for index, row in enumerate(rows):
                matrix[index, : index + 1] = row
            return matrix + np.tril(matrix, -1).T
        return None


class ElasticityHDF5Reader(BasicReader):
    """Read a complete Quantas elasticity result from native HDF5."""

    def load(self, filename: str | Path) -> ResultData:
        """Load and validate one elasticity HDF5 result file.

        Parameters
        ----------
        filename : str or Path
            Quantas elasticity HDF5 file.

        Returns
        -------
        ResultData
            Generic result containing an ``ElasticityResult`` under the
            ``"elasticity"`` key.

        Raises
        ------
        OSError
            If the file cannot be opened.
        ValueError
            If required metadata or result groups are missing.
        """
        path = Path(filename)
        self.completed = False
        self.error = None
        try:
            with h5py.File(path, "r") as h5:
                self._validate_file(h5)
                result = ResultData(
                    metadata=read_result_metadata(h5),
                    input_data=read_input_data(h5, path),
                    options=read_options(h5),
                    results={"elasticity": read_elasticity_payload(h5["results"])},
                    warnings=read_warnings(h5),
                    events=read_events(h5),
                )
        except Exception as exc:
            self.error = str(exc)
            raise
        self.completed = True
        return result

    def _validate_file(self, h5: h5py.File) -> None:
        """Validate program, module, and required HDF5 groups."""
        validate_result_schema(h5, expected_module="elasticity")


def read_elasticity_hdf5(filename: str | Path) -> ResultData:
    """Read a complete Quantas elasticity HDF5 result."""
    return ElasticityHDF5Reader().load(filename)
