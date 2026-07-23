# -*- coding: utf-8 -*-

"""Readers for seismic text inputs and native Quantas HDF5 results."""

from __future__ import annotations

from pathlib import Path
import h5py
import numpy as np
from numpy.typing import NDArray

from quantas.core.physics.elasticity import (
    StiffnessSymmetryCriterion,
    validate_stiffness_matrix,
)
from quantas.io.hdf5 import (
    read_events,
    read_input_data,
    read_options,
    read_result_metadata,
    read_warnings,
)
from quantas.models import BasicReader, InputData, ResultData, validate_result_schema
from quantas.modules.seismic.io.hdf5_payload import read_seismic_payload

_MODE_ORDER = "v_s2,v_s1,v_p"
_PAIR_ORDER = "v_s2-v_s1,v_s1-v_p"
_BRANCH_ORDER = "shear_a,shear_b,p"


class SeismicInputFileReader(BasicReader[InputData]):
    """Read stiffness and density data for a seismic calculation.

    The text format contains an optional job description, six lines defining a
    full or triangular stiffness matrix in GPa, and one density value in
    kg m^-3.
    """

    def __init__(self, seismic_input: str | Path | None = None) -> None:
        super().__init__()
        self._jobname = "Unknown"
        self._stiffness: NDArray[np.float64] | None = None
        self._density = 0.0
        self._source: Path | None = None
        self._raw: str | None = None
        if seismic_input is not None:
            self.load(seismic_input)

    @property
    def jobname(self) -> str:
        """Return the calculation description."""
        return self._jobname

    @property
    def stiffness(self) -> NDArray[np.float64] | None:
        """Return the stiffness matrix in GPa."""
        return None if self._stiffness is None else self._stiffness.copy()

    @property
    def density(self) -> float:
        """Return the material density in kg m^-3."""
        return self._density

    @property
    def raw(self) -> str | None:
        """Return the original text content."""
        return self._raw

    def load(self, filename: str | Path) -> InputData:
        """Parse a seismic text input file.

        Parameters
        ----------
        filename : str or Path
            Source text file.

        Returns
        -------
        InputData
            Generic parsed input representation. Reader properties expose the
            normalized seismic values.
        """
        self.completed = False
        self.error = None
        self._jobname = "Unknown"
        self._stiffness = None
        self._density = 0.0
        self._source = Path(filename)
        self._raw = self._source.read_text(encoding="utf-8")
        lines = [
            line.strip()
            for line in self._raw.splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        parsed = InputData(source=self._source, raw=self._raw, data={})
        if len(lines) < 7:
            self.error = (
                "The input must contain six stiffness rows and one density value."
            )
            return parsed

        start = self._read_jobname(lines[0])
        if len(lines) < start + 7:
            self.error = (
                "The input must contain six stiffness rows followed by density."
            )
            return parsed

        matrix = self._read_stiffness(lines[start : start + 6])
        if matrix is None:
            return parsed
        density = self._read_density(lines[start + 6])
        if density is None:
            return parsed

        self._stiffness = matrix
        self._density = density
        self.completed = True
        parsed.data.update(
            {
                "jobname": self._jobname,
                "stiffness": matrix.copy(),
                "density": density,
                "stiffness_unit": "GPa",
                "density_unit": "kg m^-3",
            }
        )
        return parsed

    def _read_jobname(self, first_line: str) -> int:
        """Read an optional job name and return the matrix start index."""
        try:
            float(first_line.split()[0])
        except (ValueError, IndexError):
            self._jobname = first_line
            return 1
        return 0

    def _read_stiffness(
        self,
        matrix_lines: list[str],
    ) -> NDArray[np.float64] | None:
        """Read a full, upper-triangular, or lower-triangular matrix."""
        try:
            rows = [[float(value) for value in line.split()] for line in matrix_lines]
        except ValueError:
            self.error = "Elastic constants must be valid floating-point numbers."
            return None
        if not all(np.all(np.isfinite(row)) for row in rows):
            self.error = "Elastic constants must contain finite values."
            return None

        lengths = [len(row) for row in rows]
        matrix = np.zeros((6, 6), dtype=float)
        if lengths == [6] * 6:
            matrix[:, :] = np.asarray(rows, dtype=float)
        elif lengths == [6, 5, 4, 3, 2, 1]:
            for index, row in enumerate(rows):
                matrix[index, index:] = row
            matrix += np.triu(matrix, 1).T
        elif lengths == [1, 2, 3, 4, 5, 6]:
            for index, row in enumerate(rows):
                matrix[index, : index + 1] = row
            matrix += np.tril(matrix, -1).T
        else:
            self.error = (
                "The stiffness matrix must be full, upper triangular, "
                "or lower triangular."
            )
            return None

        try:
            validate_stiffness_matrix(
                matrix,
                symmetry_tolerance=1.0e-8,
                symmetry_criterion=StiffnessSymmetryCriterion.ELEMENTWISE,
                copy=False,
            )
        except ValueError:
            self.error = "The stiffness matrix must be symmetric."
            return None
        return matrix

    def _read_density(self, line: str) -> float | None:
        """Read and validate a positive density value."""
        try:
            density = float(line.split()[0])
        except (ValueError, IndexError):
            self.error = "Density must be a valid floating-point number."
            return None
        if not np.isfinite(density) or density <= 0.0:
            self.error = "Density must be finite and positive."
            return None
        return density


class SeismicHDF5Reader(BasicReader[ResultData]):
    """Read a complete native Quantas seismic result."""

    def load(self, filename: str | Path) -> ResultData:
        """Load and validate one seismic HDF5 result file.

        Parameters
        ----------
        filename : str or Path
            Native Quantas seismic HDF5 file.

        Returns
        -------
        ResultData
            Generic result containing a :class:`SeismicResult` payload.

        Raises
        ------
        OSError
            If the file cannot be opened.
        ValueError
            If required metadata, units, or datasets are invalid.
        """
        path = Path(filename)
        self.completed = False
        self.error = None
        try:
            with h5py.File(path, "r") as h5:
                validate_result_schema(h5, expected_module="seismic")
                payload = read_seismic_payload(h5)
                result = ResultData(
                    metadata=read_result_metadata(h5),
                    input_data=read_input_data(h5, path),
                    options=read_options(h5),
                    results={"seismic": payload},
                    warnings=read_warnings(h5),
                    events=read_events(h5),
                )
        except Exception as exc:
            self.error = str(exc)
            raise
        self.completed = True
        return result


def read_seismic_hdf5(filename: str | Path) -> ResultData:
    """Read a complete native Quantas seismic HDF5 file."""
    return SeismicHDF5Reader().load(filename)
