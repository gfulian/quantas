# -*- coding: utf-8 -*-

"""Section-based CRYSTAL geometry parser shared by interface readers."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Iterable

import numpy as np
from numpy.typing import NDArray

from quantas.core.geometry.cells import lattice_from_parameters
from quantas.models.structures import CrystalStructure


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]

GEOMETRY_CONSISTENT = "GEOMETRY NOW FULLY CONSISTENT WITH THE GROUP"
GEOMETRY_WAVEFUNCTION = "GEOMETRY FOR WAVE FUNCTION -"
FINAL_OPTIMIZED_GEOMETRY = "FINAL OPTIMIZED GEOMETRY"
SUPERCELL_OPTION = "* SUPERCELL OPTION"
SUPERCELL_EXPANSION = "EXPANSION MATRIX OF PRIMITIVE CELL"
DIRECT_LATTICE_VECTORS = "DIRECT LATTICE VECTORS"
PRIMITIVE_TO_CRYSTALLOGRAPHIC = "TRANSFORMATION MATRIX PRIMITIVE-CRYSTALLOGRAPHIC CELL"

_ATOM_COUNT_RE = re.compile(r"ATOMS IN THE UNIT CELL:\s*(\d+)")
_SPACE_GROUP_RE = re.compile(r"SPACE\s+GROUP\s+N\.\s*:\s*(\d+)")
_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")


def _float_values(line: str) -> list[float]:
    """Extract floating-point values from one CRYSTAL output line."""
    return [
        float(value.replace("D", "E").replace("d", "e"))
        for value in _FLOAT_RE.findall(line)
    ]


def _is_int_token(token: str) -> bool:
    """Return whether a token is an integer literal."""
    try:
        int(token)
    except ValueError:
        return False
    return True


def _chemical_atomic_number(conventional: int) -> int:
    """Return the chemical atomic number encoded by CRYSTAL convention.

    CRYSTAL may add one or more hundreds to an atomic number to distinguish
    basis-set definitions or pseudopotential atoms. The chemical identity is
    therefore encoded by the last two decimal digits, while the original
    conventional number is retained separately in structure metadata.
    """
    return int(conventional) % 100


class CrystalGeometryParser:
    """Parse primitive, supercell, and optimized geometries from CRYSTAL output.

    The parser locates semantic section markers and tables instead of relying
    on fixed line offsets. This supports outputs with and without ``COORPRT``
    and provides one geometry implementation for phonon, QHA, elasticity, and
    future CRYSTAL interfaces.

    Parameters
    ----------
    source : str or Path or iterable of str
        CRYSTAL output path or pre-read output lines.
    """

    def __init__(self, source: str | Path | Iterable[str]) -> None:
        """Read CRYSTAL output lines and initialize the parser."""
        self.source: Path | None
        if isinstance(source, (str, Path)):
            self.source = Path(source)
            self.lines = self.source.read_text(
                encoding="utf-8", errors="strict"
            ).splitlines()
        else:
            self.source = None
            self.lines = [str(line).rstrip("\n") for line in source]

    @property
    def has_coorprt(self) -> bool:
        """Return whether the CRYSTAL input contains ``COORPRT``.

        Returns
        -------
        bool
            ``True`` when a line containing only the keyword is present.
        """
        return any(line.strip().upper() == "COORPRT" for line in self.lines)

    def expansion_matrix(self) -> IntArray:
        """Return the primitive-to-supercell expansion matrix.

        Returns
        -------
        ndarray
            Integer matrix with shape ``(3, 3)``. The identity is returned
            when no supercell expansion is present.

        Raises
        ------
        ValueError
            If an expansion marker is present but its matrix cannot be parsed.
        """
        index = self._find_first(SUPERCELL_EXPANSION)
        if index is None:
            return np.eye(3, dtype=np.int64)
        rows: list[list[int]] = []
        for line in self.lines[index + 1 : index + 8]:
            values = _float_values(line)
            if len(values) >= 3:
                rows.append([int(round(value)) for value in values[-3:]])
            if len(rows) == 3:
                break
        if len(rows) != 3:
            raise ValueError("Unable to parse CRYSTAL supercell expansion matrix")
        return np.asarray(rows, dtype=np.int64)

    def primitive_to_crystallographic(self) -> FloatArray | None:
        """Return CRYSTAL's primitive-to-crystallographic transformation.

        Returns
        -------
        ndarray or None
            Matrix with shape ``(3, 3)``, or ``None`` when CRYSTAL did not
            print the transformation.

        Raises
        ------
        ValueError
            If a transformation marker is present but fewer than nine values
            follow it.
        """
        index = self._find_first(PRIMITIVE_TO_CRYSTALLOGRAPHIC)
        if index is None:
            return None
        values: list[float] = []
        for line in self.lines[index + 1 : index + 6]:
            values.extend(_float_values(line))
            if len(values) >= 9:
                break
        if len(values) < 9:
            raise ValueError(
                "Unable to parse CRYSTAL primitive-to-crystallographic matrix"
            )
        return np.asarray(values[:9], dtype=np.float64).reshape(3, 3)

    def space_group_number(self) -> int | None:
        """Return the first explicit CRYSTAL space-group number.

        Returns
        -------
        int or None
            International space-group number, or ``None`` when absent.
        """
        for line in self.lines:
            match = _SPACE_GROUP_RE.search(line)
            if match is not None:
                return int(match.group(1))
        return None

    def initial_primitive_cell(self) -> CrystalStructure:
        """Return the primitive structure defined before supercell expansion.

        The method first prefers full ``PRIMITIVE CELL - CENTRING CODE``
        tables, including those printed by ``COORPRT`` or when reading
        ``fort.34``. If unavailable, it reconstructs the initial primitive
        structure from the earlier ``COORDINATES OF THE EQUIVALENT ATOMS``
        section.

        Returns
        -------
        CrystalStructure
            Primitive structure in the source crystallographic basis.

        Raises
        ------
        ValueError
            If no complete primitive geometry can be located.
        """
        supercell_index = self._find_first(SUPERCELL_OPTION)
        limit = len(self.lines) if supercell_index is None else supercell_index
        candidates: list[CrystalStructure] = []
        for index, line in enumerate(self.lines[:limit]):
            if "PRIMITIVE CELL - CENTRING CODE" in line:
                try:
                    candidates.append(
                        self._parse_centering_cell(index, label="initial primitive")
                    )
                except ValueError:
                    continue
        if candidates:
            return candidates[-1]
        equivalent = self._parse_equivalent_primitive(limit=limit)
        if equivalent is not None:
            return equivalent
        raise ValueError(
            "Unable to locate a primitive CRYSTAL geometry. Add COORPRT to "
            "the CRYSTAL input or provide an output containing the initial "
            "primitive-cell table."
        )

    def wavefunction_cell(self) -> CrystalStructure:
        """Return the geometry used to construct the wave function.

        Returns
        -------
        CrystalStructure
            Wave-function cell, which may be primitive or a phonon supercell.

        Raises
        ------
        ValueError
            If the wave-function geometry is absent or incomplete.
        """
        index = self._find_first(GEOMETRY_WAVEFUNCTION)
        if index is None:
            raise ValueError("No CRYSTAL wave-function geometry was found")
        return self._parse_geometry_after_marker(index, label="wavefunction cell")

    def optimized_cells(self) -> list[CrystalStructure]:
        """Return all final optimized geometries printed by a CRYSTAL QHA run.

        Returns
        -------
        list of CrystalStructure
            Optimized source cells in output order.
        """
        cells: list[CrystalStructure] = []
        for index, line in enumerate(self.lines):
            if FINAL_OPTIMIZED_GEOMETRY in line:
                cells.append(
                    self._parse_geometry_after_marker(
                        index,
                        label=f"optimized cell {len(cells)}",
                    )
                )
        return cells

    def _parse_geometry_after_marker(
        self,
        marker_index: int,
        *,
        label: str,
    ) -> CrystalStructure:
        """Parse a centring-code cell following a geometry marker."""
        for index in range(marker_index, min(marker_index + 30, len(self.lines))):
            if "CELL - CENTRING CODE" in self.lines[index]:
                return self._parse_centering_cell(index, label=label)
        raise ValueError(f"No complete cell table follows CRYSTAL section '{label}'")

    def _parse_centering_cell(self, index: int, *, label: str) -> CrystalStructure:
        """Parse one CRYSTAL centring-code geometry table."""
        parameters = self._cell_parameters_after(index)
        atom_line_index: int | None = None
        natoms: int | None = None
        for cursor in range(index, min(index + 40, len(self.lines))):
            match = _ATOM_COUNT_RE.search(self.lines[cursor])
            if match is not None:
                atom_line_index = cursor
                natoms = int(match.group(1))
                break
        if atom_line_index is None or natoms is None:
            raise ValueError("CRYSTAL cell table does not contain an atom count")

        numbers: list[int] = []
        conventional_numbers: list[int] = []
        positions: list[list[float]] = []
        last_atom_index = atom_line_index
        for cursor in range(atom_line_index + 1, len(self.lines)):
            tokens = self.lines[cursor].split()
            if (
                len(tokens) >= 7
                and _is_int_token(tokens[0])
                and tokens[1] in {"T", "F"}
                and _is_int_token(tokens[2])
            ):
                values = _float_values(self.lines[cursor])
                if len(values) >= 5:
                    conventional = int(tokens[2])
                    conventional_numbers.append(conventional)
                    numbers.append(_chemical_atomic_number(conventional))
                    positions.append([float(value) for value in values[-3:]])
                    last_atom_index = cursor
                    if len(numbers) == natoms:
                        break
            if cursor - atom_line_index > natoms + 100:
                break
        if len(numbers) != natoms:
            raise ValueError(
                f"CRYSTAL geometry declared {natoms} atoms but {len(numbers)} "
                "coordinate rows were parsed"
            )

        lattice = self._direct_lattice_after(last_atom_index)
        if lattice is None:
            lattice = lattice_from_parameters(*parameters)
        metadata = {
            "geometry_source": label,
            "cell_parameters": np.asarray(parameters, dtype=np.float64),
            "coorprt_present": self.has_coorprt,
            "crystal_conventional_atomic_numbers": np.asarray(
                conventional_numbers, dtype=np.int64
            ),
        }
        return CrystalStructure(
            lattice=lattice,
            fractional_positions=np.asarray(positions, dtype=np.float64),
            atomic_numbers=np.asarray(numbers, dtype=np.int64),
            label=label,
            metadata=metadata,
        )

    def _parse_equivalent_primitive(self, *, limit: int) -> CrystalStructure | None:
        """Parse the initial equivalent-atom primitive-cell section."""
        marker: int | None = None
        for index, line in enumerate(self.lines[:limit]):
            if "LATTICE PARAMETERS" in line and "- PRIMITIVE CELL" in line:
                marker = index
        if marker is None:
            return None
        parameters = self._cell_parameters_after(marker)
        coordinate_marker: int | None = None
        for index in range(marker, min(limit, marker + 80)):
            if "COORDINATES OF THE EQUIVALENT ATOMS" in self.lines[index]:
                coordinate_marker = index
                break
        if coordinate_marker is None:
            return None

        numbers: list[int] = []
        conventional_numbers: list[int] = []
        positions: list[list[float]] = []
        for line in self.lines[coordinate_marker + 1 : limit]:
            if "NUMBER OF SYMMETRY OPERATORS" in line:
                break
            tokens = line.split()
            if len(tokens) >= 8 and all(_is_int_token(token) for token in tokens[:4]):
                values = _float_values(line)
                if len(values) >= 7:
                    conventional = int(tokens[3])
                    conventional_numbers.append(conventional)
                    numbers.append(_chemical_atomic_number(conventional))
                    positions.append([float(value) for value in values[-3:]])
        if not numbers:
            return None
        return CrystalStructure(
            lattice=lattice_from_parameters(*parameters),
            fractional_positions=np.asarray(positions, dtype=np.float64),
            atomic_numbers=np.asarray(numbers, dtype=np.int64),
            label="initial primitive equivalent atoms",
            metadata={
                "geometry_source": "equivalent atoms",
                "cell_parameters": np.asarray(parameters, dtype=np.float64),
                "coorprt_present": self.has_coorprt,
                "crystal_conventional_atomic_numbers": np.asarray(
                    conventional_numbers, dtype=np.int64
                ),
            },
        )

    def _cell_parameters_after(
        self, index: int
    ) -> tuple[float, float, float, float, float, float]:
        """Return the first six-value cell-parameter row after an index."""
        for line in self.lines[index + 1 : index + 8]:
            values = _float_values(line)
            if len(values) >= 6:
                return tuple(float(value) for value in values[:6])  # type: ignore[return-value]
        raise ValueError("Unable to parse CRYSTAL lattice parameters")

    def _direct_lattice_after(self, index: int) -> FloatArray | None:
        """Return the first direct-lattice vector table after an atom table."""
        stop_markers = (
            FINAL_OPTIMIZED_GEOMETRY,
            GEOMETRY_WAVEFUNCTION,
            SUPERCELL_OPTION,
            "CRYSTALLOGRAPHIC CELL",
        )
        upper = min(len(self.lines), index + 500)
        marker: int | None = None
        for cursor in range(index + 1, upper):
            if any(value in self.lines[cursor] for value in stop_markers):
                break
            if DIRECT_LATTICE_VECTORS in self.lines[cursor]:
                marker = cursor
                break
        if marker is None:
            return None
        rows: list[list[float]] = []
        for line in self.lines[marker + 1 : marker + 12]:
            values = _float_values(line)
            if len(values) >= 3:
                rows.append([float(value) for value in values[-3:]])
            if len(rows) == 3:
                break
        if len(rows) != 3:
            return None
        return np.asarray(rows, dtype=np.float64)

    def _find_first(self, marker: str) -> int | None:
        """Return the first line index containing a marker."""
        for index, line in enumerate(self.lines):
            if marker in line:
                return index
        return None
