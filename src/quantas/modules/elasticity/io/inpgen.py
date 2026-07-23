# -*- coding: utf-8 -*-

"""Generate Quantas elasticity inputs from supported electronic-structure codes."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from quantas.interfaces.crystal.elasticity import CrystalElasticityReader
from quantas.interfaces.vasp.elasticity import VASPElasticityReader

ElasticityReader = CrystalElasticityReader | VASPElasticityReader
ElasticityReaderType = type[CrystalElasticityReader] | type[VASPElasticityReader]


READERS: dict[str, ElasticityReaderType] = {
    "crystal": CrystalElasticityReader,
    "vasp": VASPElasticityReader,
}


class ElasticityInputCreator:
    """Create Quantas elasticity input files from supported external outputs.

    Parameters
    ----------
    interface : str
        Name of the external-code interface. Supported values are ``crystal``
        and ``vasp``.
    """

    def __init__(self, interface: str) -> None:
        """Initialize an input creator for one external-code interface."""
        if interface not in READERS:
            raise KeyError(f"Unsupported elasticity interface: {interface}")

        self.interface = interface
        self.reader: ElasticityReader | None = None
        self.completed = False
        self.error: str | None = None

    def read(self, filename: str | Path) -> tuple[bool, str | None]:
        """Read an external output file.

        Parameters
        ----------
        filename : str or Path
            Path to the external output file.

        Returns
        -------
        tuple
            Completion flag and error message.
        """
        reader_class = READERS[self.interface]
        reader = reader_class(filename)
        self.reader = reader

        self.completed = reader.completed
        self.error = reader.error

        return self.completed, self.error

    def write(self, filename: str | Path, jobname: str = "Unknown") -> None:
        """Write a Quantas elasticity input file.

        Parameters
        ----------
        filename : str or Path
            Output filename.
        jobname : str, optional
            Short description of the input file.

        Raises
        ------
        RuntimeError
            If no external output was successfully read.
        """
        if self.reader is None or not self.completed:
            raise RuntimeError("No completed elasticity input data are available.")

        filename = Path(filename)

        with filename.open("w", encoding="utf-8") as stream:
            stream.write(f"{jobname}\n")
            self._write_triangular_matrix(stream, self.reader.stiffness)

            if self.reader.density > 0.0:
                stream.write(f"{self.reader.density:.3f}\n")

    @staticmethod
    def _write_triangular_matrix(stream, matrix: np.ndarray) -> None:
        """Write the upper triangular part of the elastic matrix.

        Parameters
        ----------
        stream : file-like object
            Open text stream.
        matrix : ndarray
            Stiffness matrix in Voigt notation.
        """
        for i in range(6):
            values = [f"{matrix[i, j]:12.6f}" for j in range(i, 6)]
            indent = " " * (13 * i)
            stream.write(indent + " ".join(values) + "\n")

    @staticmethod
    def _write_full_matrix(stream, matrix: np.ndarray) -> None:
        """Write the full elastic matrix.

        Parameters
        ----------
        stream : file-like object
            Open text stream.
        matrix : ndarray
            Stiffness matrix in Voigt notation.
        """
        for i in range(6):
            values = [f"{matrix[i, j]:12.6f}" for j in range(6)]
            stream.write(" ".join(values) + "\n")
