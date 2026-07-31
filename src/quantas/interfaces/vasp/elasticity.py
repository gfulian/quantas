# -*- coding: utf-8 -*-

"""VASP OUTCAR reader for second-order elastic constants."""

from __future__ import annotations

import re
from pathlib import Path
from typing import TypedDict

import numpy as np

from quantas.models.reader import BasicReader

_CLAMPED_ELASTIC_MODULI = "SYMMETRIZED ELASTIC MODULI (kBar)"
_RELAXED_ELASTIC_MODULI = "TOTAL ELASTIC MODULI (kBar)"
_ATOMIC_MASS_UNIT_PER_ANGSTROM_CUBED = 1660.53906660

_POMASS_RE = re.compile(
    r"POMASS\s*=\s*([-+0-9.EeDd]+)\s*;\s*ZVAL",
    flags=re.IGNORECASE,
)
_IONS_PER_TYPE_RE = re.compile(
    r"ions\s+per\s+type\s*=\s*([0-9 \t]+)",
    flags=re.IGNORECASE,
)
_VOLUME_RE = re.compile(
    r"volume\s+of\s+cell\s*:\s*([-+0-9.EeDd]+)",
    flags=re.IGNORECASE,
)


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

        text = path.read_text(encoding="utf-8")
        lines = text.splitlines()
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
        self._data["density"] = self._read_density(text)
        self.completed = True

    @staticmethod
    def _read_density(text: str) -> float:
        """Return the final VASP cell density in kg m^-3 when available.

        VASP reports species masses through ``POMASS``, species populations
        through ``ions per type``, and one or more cell volumes during a run.
        The final reported volume is paired with the species-resolved mass.
        Missing or incomplete metadata leave density unavailable rather than
        invalidating an otherwise usable Elasticity input.
        """
        mass_matches = _POMASS_RE.findall(text)
        count_matches = _IONS_PER_TYPE_RE.findall(text)
        volume_matches = _VOLUME_RE.findall(text)
        if not mass_matches or not count_matches or not volume_matches:
            return 0.0

        try:
            masses = np.asarray(
                [
                    float(value.replace("D", "E").replace("d", "e"))
                    for value in mass_matches
                ],
                dtype=float,
            )
            counts = np.asarray(
                [int(value) for value in count_matches[0].split()],
                dtype=int,
            )
            volume = float(volume_matches[-1].replace("D", "E").replace("d", "e"))
        except ValueError:
            return 0.0

        if (
            masses.shape != counts.shape
            or masses.size == 0
            or np.any(~np.isfinite(masses))
            or np.any(masses <= 0.0)
            or np.any(counts <= 0)
            or not np.isfinite(volume)
            or volume <= 0.0
        ):
            return 0.0

        total_mass = float(np.dot(masses, counts))
        density = total_mass / volume * _ATOMIC_MASS_UNIT_PER_ANGSTROM_CUBED
        return float(density) if np.isfinite(density) and density > 0.0 else 0.0

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
