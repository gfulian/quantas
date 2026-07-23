# -*- coding: utf-8 -*-

"""CRYSTAL output reader for second-order elastic constants.

The reader extracts the final structure, volume, density, static pressure,
energy, and symmetrized stiffness matrix required by elasticity and
quasi-static thermoelastic workflows. CRYSTAL's ``PRESSURE`` keyword is
tracked explicitly because it determines whether hydrostatic pre-stress terms
are included in the reported elastic coefficients.
"""

from __future__ import annotations

from pathlib import Path
import re
from typing import TypeAlias, TypedDict

import numpy as np
from numpy.typing import NDArray

from quantas.core.geometry import analyze_symmetry
from quantas.interfaces.crystal.geometry import CrystalGeometryParser
from quantas.models.reader import BasicReader
from quantas.models.structures import CrystalStructure, SymmetryMetadata


FloatArray: TypeAlias = NDArray[np.float64]

_GEOMETRY_CONSISTENT = "GEOMETRY NOW FULLY CONSISTENT WITH THE GROUP"
_PRIMITIVE_CELL = "PRIMITIVE CELL - "
_ELASTICITY_OPTION = "ELASTCON OPTION"
_ELASTICITY_COMPLETED = "FINAL RESULTS START"
_ELASTICITY_CONSTANTS = "SYMMETRIZED ELASTIC CONSTANTS"
_ELASTIC_PRESSURE = "ELASTIC PROPERTIES AT PRESSURE"
_STRESS_PRESSURE = "PRESSURE IN GIGAPASCAL:"
_CELL_VOLUME = "VOLUME OF THE CELL:"
_CRYSTAL_DENSITY = "DENSITY OF THE CRYSTAL"

_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")
_ENERGY_RE = re.compile(
    r"TOTAL ENERGY\(DFT\)\(AU\)\(\s*\d+\)\s*"
    r"([-+]?\d+(?:\.\d*)?[EeDd][-+]?\d+)",
    flags=re.IGNORECASE,
)


class _ElasticityData(TypedDict):
    """Typed payload stored by :class:`CrystalElasticityReader`."""

    stiffness: FloatArray
    density: float
    volume: float
    energy: float
    pressure: float
    stress_pressure: float
    pressure_keyword_value: float
    prestress_applied: bool
    structure: CrystalStructure | None
    symmetry: SymmetryMetadata | None


def _float_values(line: str) -> list[float]:
    """Return all floating-point values occurring in one line.

    Parameters
    ----------
    line : str
        Text line to inspect.

    Returns
    -------
    list of float
        Floating-point values in source order.
    """
    return [
        float(value.replace("D", "E").replace("d", "e"))
        for value in _FLOAT_RE.findall(line)
    ]


class CrystalElasticityReader(BasicReader[None]):
    """Read a CRYSTAL second-order elastic-constants calculation.

    Parameters
    ----------
    filename : str, Path, or None, optional
        CRYSTAL text output file. When provided, the file is loaded during
        construction.
    symprec : float, optional
        Cartesian tolerance in angstrom used by spglib for final-structure
        symmetry analysis.
    angle_tolerance : float, optional
        Angular tolerance in degrees used by spglib. A negative value requests
        spglib's internal default.
    """

    def __init__(
        self,
        filename: str | Path | None = None,
        *,
        symprec: float = 1.0e-5,
        angle_tolerance: float = -1.0,
    ) -> None:
        super().__init__()
        self.symprec = float(symprec)
        self.angle_tolerance = float(angle_tolerance)
        self._data: _ElasticityData = self._empty_data()
        if filename is not None:
            self.load(filename)

    @staticmethod
    def _empty_data() -> _ElasticityData:
        """Return the initialized reader payload.

        Returns
        -------
        _ElasticityData
            Empty, typed reader data.
        """
        return {
            "stiffness": np.zeros((6, 6), dtype=np.float64),
            "density": 0.0,
            "volume": 0.0,
            "energy": np.nan,
            "pressure": np.nan,
            "stress_pressure": np.nan,
            "pressure_keyword_value": np.nan,
            "prestress_applied": False,
            "structure": None,
            "symmetry": None,
        }

    @property
    def stiffness(self) -> FloatArray:
        """Return the elastic stiffness matrix in Voigt notation, in GPa."""
        return self._data["stiffness"]

    @property
    def density(self) -> float:
        """Return the final crystal density in kg m^-3."""
        return self._data["density"]

    @property
    def volume(self) -> float:
        """Return the final primitive-cell volume in angstrom cubed."""
        return self._data["volume"]

    @property
    def energy(self) -> float:
        """Return the final static DFT energy in hartree."""
        return self._data["energy"]

    @property
    def pressure(self) -> float:
        """Return the pressure used in CRYSTAL's elastic correction, in GPa."""
        return self._data["pressure"]

    @property
    def stress_pressure(self) -> float:
        """Return the final pressure from the unstrained stress tensor, in GPa."""
        return self._data["stress_pressure"]

    @property
    def pressure_keyword_value(self) -> float:
        """Return the value supplied to CRYSTAL's ``PRESSURE`` keyword."""
        return self._data["pressure_keyword_value"]

    @property
    def prestress_applied(self) -> bool:
        """Return whether an explicit CRYSTAL ``PRESSURE`` keyword was found."""
        return self._data["prestress_applied"]

    @property
    def structure(self) -> CrystalStructure | None:
        """Return the final compact CRYSTAL structure, when available."""
        return self._data["structure"]

    @property
    def symmetry(self) -> SymmetryMetadata | None:
        """Return spglib symmetry metadata for the final structure."""
        return self._data["symmetry"]

    def load(self, filename: str | Path) -> None:
        """Read elastic and structural data from a CRYSTAL output file.

        Parameters
        ----------
        filename : str or Path
            CRYSTAL text output file.
        """
        path = Path(filename)
        self.completed = False
        self.error = None
        self._data = self._empty_data()

        if not self.is_elasticity_output(path):
            self.error = (
                "The file is not recognized as a CRYSTAL elastic-constants output."
            )
            return
        if not self.is_output_completed(path):
            self.error = "The CRYSTAL elastic-constants calculation is incomplete."
            return

        lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
        try:
            self._data["stiffness"] = self._read_stiffness(lines)
            self._data["pressure"] = self._read_elastic_pressure(lines)
            keyword_value = self._read_pressure_keyword(lines)
            if keyword_value is not None:
                self._data["pressure_keyword_value"] = keyword_value
                self._data["prestress_applied"] = True
            self._data["stress_pressure"] = self._read_last_scalar_after_marker(
                lines,
                _STRESS_PRESSURE,
                default=np.nan,
            )
            self._data["volume"] = self._read_last_scalar_after_marker(
                lines,
                _CELL_VOLUME,
                default=np.nan,
            )
            self._data["density"] = 1000.0 * self._read_final_density(lines)
            self._data["energy"] = self._read_final_energy(lines)
            structure = self._read_final_structure(path)
            self._data["structure"] = structure
            if structure is not None:
                self._data["symmetry"] = analyze_symmetry(
                    structure,
                    symprec=self.symprec,
                    angle_tolerance=self.angle_tolerance,
                )
                if not np.isfinite(self._data["volume"]):
                    self._data["volume"] = structure.volume
            if not np.isfinite(self._data["density"]):
                self._data["density"] = self.initial_crystal_density(path)
        except (OSError, ValueError, IndexError) as exc:
            self.error = str(exc)
            return

        if not np.isfinite(self.density) or self.density <= 0.0:
            self.error = "Unable to determine a positive final crystal density."
            return
        self.completed = True

    def is_elasticity_output(self, filename: str | Path) -> bool:
        """Return whether a CRYSTAL output contains an ELASTCON calculation.

        Parameters
        ----------
        filename : str or Path
            CRYSTAL text output file.

        Returns
        -------
        bool
            ``True`` when the ELASTCON section is present.
        """
        with Path(filename).open("r", encoding="utf-8") as stream:
            return any(_ELASTICITY_OPTION in line for line in stream)

    def is_output_completed(self, filename: str | Path) -> bool:
        """Return whether the CRYSTAL elastic calculation reached final results.

        Parameters
        ----------
        filename : str or Path
            CRYSTAL text output file.

        Returns
        -------
        bool
            ``True`` when the final-results marker is present.
        """
        with Path(filename).open("r", encoding="utf-8") as stream:
            return any(_ELASTICITY_COMPLETED in line for line in stream)

    def elasticity_start_line(self, filename: str | Path) -> int:
        """Return the zero-based line containing the first final stiffness row.

        Parameters
        ----------
        filename : str or Path
            CRYSTAL output file.

        Returns
        -------
        int
            Index of the first upper-triangular matrix row.

        Raises
        ------
        ValueError
            If the elastic-constants header is absent.
        """
        lines = Path(filename).read_text(encoding="utf-8").splitlines()
        indexes = [
            index for index, line in enumerate(lines) if _ELASTICITY_CONSTANTS in line
        ]
        if not indexes:
            raise ValueError("Elastic stiffness header not found in CRYSTAL output.")
        return indexes[-1] + 2

    def initial_crystal_density(self, filename: str | Path) -> float:
        """Return the initial primitive-cell density in kg m^-3.

        Parameters
        ----------
        filename : str or Path
            CRYSTAL output file.

        Returns
        -------
        float
            Initial density, or zero when unavailable.
        """
        consistent = False
        with Path(filename).open("r", encoding="utf-8") as stream:
            for line in stream:
                if _GEOMETRY_CONSISTENT in line:
                    consistent = True
                elif consistent and _PRIMITIVE_CELL in line:
                    values = _float_values(line)
                    if "DENSITY" in line.upper():
                        tail = line.upper().split("DENSITY", maxsplit=1)[-1]
                        density_values = _float_values(tail)
                        if density_values:
                            return float(density_values[0]) * 1000.0
                    if values:
                        index = -2 if len(values) >= 2 and values[-1] == 3.0 else -1
                        return float(values[index]) * 1000.0
        return 0.0

    @staticmethod
    def _read_stiffness(lines: list[str]) -> FloatArray:
        """Read the final symmetrized upper-triangular stiffness matrix.

        Parameters
        ----------
        lines : list of str
            CRYSTAL output lines.

        Returns
        -------
        ndarray
            Symmetric ``(6, 6)`` stiffness matrix in GPa.

        Raises
        ------
        ValueError
            If the stiffness section is absent or malformed.
        """
        indexes = [
            index for index, line in enumerate(lines) if _ELASTICITY_CONSTANTS in line
        ]
        if not indexes:
            raise ValueError("Elastic stiffness header not found in CRYSTAL output.")
        start = indexes[-1] + 2
        matrix: FloatArray = np.zeros((6, 6), dtype=np.float64)
        for row in range(6):
            values = _float_values(lines[start + row])
            expected = 6 - row
            if len(values) != expected:
                raise ValueError(
                    "Problem collecting elastic stiffness row "
                    f"{row + 1}: expected {expected} values, found {len(values)}."
                )
            matrix[row, row:] = values
        matrix += np.triu(matrix, 1).T
        return matrix

    @staticmethod
    def _read_elastic_pressure(lines: list[str]) -> float:
        """Read the pressure associated with the final elastic coefficients."""
        for line in reversed(lines):
            if _ELASTIC_PRESSURE in line:
                values = _float_values(line)
                if values:
                    return float(values[-1])
        return np.nan

    @staticmethod
    def _read_pressure_keyword(lines: list[str]) -> float | None:
        """Read an explicit ``PRESSURE`` keyword and its following value."""
        for index, line in enumerate(lines):
            if line.strip().upper() != "PRESSURE":
                continue
            for following in lines[index + 1 : index + 5]:
                values = _float_values(following)
                if values:
                    return float(values[0])
        return None

    @staticmethod
    def _read_last_scalar_after_marker(
        lines: list[str],
        marker: str,
        *,
        default: float,
    ) -> float:
        """Return the final scalar on a line containing ``marker``."""
        for line in reversed(lines):
            if marker in line:
                values = _float_values(line)
                if values:
                    return float(values[-1])
        return float(default)

    @staticmethod
    def _read_final_density(lines: list[str]) -> float:
        """Return the final crystal density in g cm^-3."""
        for line in reversed(lines):
            if _CRYSTAL_DENSITY in line:
                tail = line.split("=", maxsplit=1)[-1]
                values = _float_values(tail)
                if values:
                    return float(values[0])
        return np.nan

    @staticmethod
    def _read_final_energy(lines: list[str]) -> float:
        """Return the last CRYSTAL total DFT energy in hartree."""
        matches: list[float] = []
        for line in lines:
            match = _ENERGY_RE.search(line)
            if match is not None:
                matches.append(
                    float(match.group(1).replace("D", "E").replace("d", "e"))
                )
        if not matches:
            return np.nan
        return float(matches[-1])

    @staticmethod
    def _read_final_structure(path: Path) -> CrystalStructure | None:
        """Return the final optimized CRYSTAL structure when available."""
        parser = CrystalGeometryParser(path)
        optimized = parser.optimized_cells()
        if optimized:
            return optimized[-1]
        try:
            return parser.wavefunction_cell()
        except ValueError:
            pass
        try:
            return parser.initial_primitive_cell()
        except ValueError:
            return None
