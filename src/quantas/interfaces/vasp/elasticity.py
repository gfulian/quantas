# -*- coding: utf-8 -*-

"""VASP OUTCAR reader for second-order elastic constants."""

from __future__ import annotations

from pathlib import Path
from typing import TypedDict

import numpy as np

from quantas.models.reader import BasicReader

_CLAMPED_ELASTIC_MODULI = "SYMMETRIZED ELASTIC MODULI (kBar)"
_RELAXED_ELASTIC_MODULI = "TOTAL ELASTIC MODULI (kBar)"


class _ElasticityData(TypedDict):
    """Typed payload stored by elasticity interface readers."""

    stiffness: np.ndarray
    density: float


class VASPElasticityReader(BasicReader[None]):
    """Read a stiffness matrix from a VASP OUTCAR file."""

    def __init__(self, filename: str | Path | None = None) -> None:
        super().__init__()
        self._data: _ElasticityData = {
            "stiffness": np.zeros((6, 6), dtype=float),
            "density": 0.0,
        }
        if filename is not None:
            self.load(filename)

    @property
    def stiffness(self) -> np.ndarray:
        """Return the elastic stiffness matrix in Voigt notation, in GPa."""
        return self._data["stiffness"]

    @property
    def density(self) -> float:
        """Return the crystal density in kg m^-3 when available."""
        return float(self._data["density"])

    def load(self, filename: str | Path) -> None:
        """Read elastic constants from a VASP OUTCAR file.

        Parameters
        ----------
        filename : str or Path
            VASP OUTCAR text file.
        """
        path = Path(filename)
        self.completed = False
        self.error = None
        self._data = {
            "stiffness": np.zeros((6, 6), dtype=float),
            "density": 0.0,
        }

        if not self.is_elasticity_output(path):
            self.error = (
                f"'{path}' does not appear to be a VASP output containing "
                "elastic moduli."
            )
            return

        lines = path.read_text(encoding="utf-8").splitlines()
        start = self.elasticity_start_line(path)
        for row in range(6):
            values = lines[start + row].split()
            for column in range(6):
                try:
                    self._data["stiffness"][row, column] = (
                        float(values[column + 1]) / 10.0
                    )
                except (ValueError, IndexError):
                    self.error = (
                        "Problem collecting elastic stiffness element "
                        f"({row + 1}, {column + 1}) from the OUTCAR file."
                    )
                    return

        stiffness = self._data["stiffness"]
        stiffness[[3, 5]] = stiffness[[5, 3]]
        stiffness[:, [3, 5]] = stiffness[:, [5, 3]]
        self.completed = True

    def is_elasticity_output(self, filename: str | Path) -> bool:
        """Return whether a VASP OUTCAR contains an elastic-moduli section."""
        with Path(filename).open("r", encoding="utf-8") as stream:
            return any(_CLAMPED_ELASTIC_MODULI in line for line in stream)

    def elasticity_start_line(self, filename: str | Path) -> int:
        """Return the first data row, preferring relaxed-ion elastic moduli.

        Raises
        ------
        ValueError
            If no elastic-moduli section is found.
        """
        clamped_line: int | None = None
        relaxed_line: int | None = None
        with Path(filename).open("r", encoding="utf-8") as stream:
            for index, line in enumerate(stream):
                if _CLAMPED_ELASTIC_MODULI in line:
                    clamped_line = index + 3
                if _RELAXED_ELASTIC_MODULI in line:
                    relaxed_line = index + 3
        if relaxed_line is not None:
            return relaxed_line
        if clamped_line is not None:
            return clamped_line
        raise ValueError("Elastic moduli section not found in VASP OUTCAR.")
