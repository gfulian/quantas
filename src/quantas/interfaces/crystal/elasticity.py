# -*- coding: utf-8 -*-

"""CRYSTAL output reader for second-order elastic constants.

The reader extracts the unstrained elastic-reference structure, volume,
density, static pressure, energy, and symmetrized stiffness matrix required by
elasticity and quasi-static thermoelastic workflows. Geometry and stress
records printed after the first elastic strain are deliberately excluded from
the reference state. CRYSTAL's ``PRESSURE`` and ``PRESSEOS`` keywords are
tracked explicitly because they determine whether hydrostatic pre-stress terms
are included in the reported elastic coefficients.
"""

from __future__ import annotations

from pathlib import Path
from typing import TypeAlias, TypedDict

import numpy as np
from numpy.typing import NDArray

from quantas.core.geometry import analyze_symmetry
from quantas.interfaces.crystal import markers, patterns
from quantas.interfaces.crystal.geometry import CrystalGeometryParser
from quantas.interfaces.crystal.output import CrystalOutputParser
from quantas.models.reader import BasicReader
from quantas.models.structures import CrystalStructure, SymmetryMetadata


FloatArray: TypeAlias = NDArray[np.float64]


class _ElasticityData(TypedDict):
    """Typed payload stored by :class:`CrystalElasticityReader`."""

    stiffness: FloatArray
    density: float
    volume: float
    energy: float
    scf_energy: float
    energy_provenance: dict[str, object]
    pressure: float
    stress_pressure: float
    pressure_keyword_value: float
    prestress_keyword: str | None
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
        for value in patterns.FLOAT_RE.findall(line)
    ]


class CrystalElasticityReader(BasicReader[None]):
    """Read a CRYSTAL second-order elastic-constants calculation.

    Parameters
    ----------
    filename : str, Path, or None, optional
        CRYSTAL text output file. When provided, the file is loaded during
        construction.
    symprec : float, optional
        Cartesian tolerance in angstrom used by spglib for elastic-reference
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
            "scf_energy": np.nan,
            "energy_provenance": {},
            "pressure": np.nan,
            "stress_pressure": np.nan,
            "pressure_keyword_value": np.nan,
            "prestress_keyword": None,
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
        """Return the elastic-reference crystal density in kg m^-3."""
        return self._data["density"]

    @property
    def volume(self) -> float:
        """Return the elastic-reference primitive-cell volume in angstrom cubed."""
        return self._data["volume"]

    @property
    def energy(self) -> float:
        """Return the elastic-reference total static energy in hartree.

        The value includes any a-posteriori correction explicitly included by
        CRYSTAL in a ``TOTAL ENERGY + ...`` record.  It therefore corresponds
        to the energy that should be used by static-E(V), HA/QHA, Kieffer, and
        quasi-static thermoelastic workflows.

        Returns
        -------
        float
            Resolved total static energy in hartree.
        """
        return self._data["energy"]

    @property
    def total_energy(self) -> float:
        """Return the resolved elastic-reference total energy.

        Returns
        -------
        float
            Total energy in hartree, including printed a-posteriori
            corrections when present.
        """
        return self._data["energy"]

    @property
    def scf_energy(self) -> float:
        """Return the uncorrected elastic-reference SCF energy.

        Returns
        -------
        float
            Electronic SCF energy in hartree.
        """
        return self._data["scf_energy"]

    @property
    def energy_provenance(self) -> dict[str, object]:
        """Return provenance describing resolution of the total energy.

        Returns
        -------
        dict
            Source markers, correction labels, SCF energy, and total-energy
            correction metadata.
        """
        return dict(self._data["energy_provenance"])

    @property
    def pressure(self) -> float:
        """Return the pressure used in CRYSTAL's elastic correction, in GPa."""
        return self._data["pressure"]

    @property
    def stress_pressure(self) -> float:
        """Return pressure from the elastic-reference stress tensor, in GPa."""
        return self._data["stress_pressure"]

    @property
    def pressure_keyword_value(self) -> float:
        """Return the value supplied to CRYSTAL ``PRESSURE`` or ``PRESSEOS``."""
        return self._data["pressure_keyword_value"]

    @property
    def prestress_keyword(self) -> str | None:
        """Return the CRYSTAL pre-stress keyword, when explicitly present."""
        return self._data["prestress_keyword"]

    @property
    def prestress_applied(self) -> bool:
        """Return whether CRYSTAL applied an explicit hydrostatic correction."""
        return self._data["prestress_applied"]

    @property
    def structure(self) -> CrystalStructure | None:
        """Return the unstrained elastic-reference structure, when available."""
        return self._data["structure"]

    @property
    def symmetry(self) -> SymmetryMetadata | None:
        """Return spglib symmetry metadata for the elastic-reference structure."""
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
            keyword = self._read_prestress_keyword(lines)
            if keyword is not None:
                keyword_name, keyword_value = keyword
                self._data["pressure_keyword_value"] = keyword_value
                self._data["prestress_keyword"] = keyword_name
                self._data["prestress_applied"] = True
            reference_lines = self._elastic_reference_lines(lines)
            self._data["volume"] = self._read_last_scalar_after_marker(
                reference_lines,
                markers.CELL_VOLUME,
                default=np.nan,
            )
            structure, geometry_start = self._read_reference_structure(
                reference_lines,
                volume=self._data["volume"],
            )
            pressure_lines = reference_lines[geometry_start:]
            self._data["stress_pressure"] = self._read_last_scalar_after_marker(
                pressure_lines,
                markers.STRESS_PRESSURE,
                default=np.nan,
            )
            self._data["density"] = 1000.0 * self._read_final_density(
                reference_lines
            )
            (
                self._data["scf_energy"],
                self._data["energy"],
                self._data["energy_provenance"],
            ) = self._read_reference_energies(reference_lines)
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
            self.error = (
                "Unable to determine a positive elastic-reference crystal density."
            )
            return
        self.completed = True

    def is_elasticity_output(self, filename: str | Path) -> bool:
        """Return whether a CRYSTAL output contains elastic constants.

        Parameters
        ----------
        filename : str or Path
            CRYSTAL text output file.

        Returns
        -------
        bool
            ``True`` when an elastic-producing operation or the standard
            symmetrized stiffness block is present.
        """
        with Path(filename).open("r", encoding="utf-8") as stream:
            for line in stream:
                if markers.ELASTICITY_CONSTANTS in line:
                    return True
                if any(
                    marker in line
                    for marker in markers.ELASTICITY_OPTION_MARKERS
                ):
                    return True
        return False

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
            return any(markers.ELASTICITY_RESULTS in line for line in stream)

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
            index for index, line in enumerate(lines) if markers.ELASTICITY_CONSTANTS in line
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
                if markers.GEOMETRY_CONSISTENT in line:
                    consistent = True
                elif consistent and markers.PRIMITIVE_CELL in line:
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
            index for index, line in enumerate(lines) if markers.ELASTICITY_CONSTANTS in line
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
            if markers.ELASTIC_PRESSURE in line:
                values = _float_values(line)
                if values:
                    return float(values[-1])
        return np.nan

    @staticmethod
    def _read_prestress_keyword(lines: list[str]) -> tuple[str, float] | None:
        """Read an explicit ``PRESSURE`` or ``PRESSEOS`` keyword and value."""
        for index, line in enumerate(lines):
            keyword = line.strip().upper()
            if keyword not in {"PRESSURE", "PRESSEOS"}:
                continue
            for following in lines[index + 1 : index + 5]:
                values = _float_values(following)
                if values:
                    return keyword, float(values[0])
        return None

    @staticmethod
    def _read_pressure_keyword(lines: list[str]) -> float | None:
        """Return the explicit CRYSTAL pre-stress value, when present."""
        keyword = CrystalElasticityReader._read_prestress_keyword(lines)
        return None if keyword is None else keyword[1]

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
            if markers.CRYSTAL_DENSITY in line:
                tail = line.split("=", maxsplit=1)[-1]
                values = _float_values(tail)
                if values:
                    return float(values[0])
        return np.nan

    @staticmethod
    def _read_final_energy(lines: list[str]) -> float:
        """Return the last resolved CRYSTAL total energy in hartree.

        This compatibility helper now follows the central CRYSTAL energy
        resolver: empirical DFT-D/gCP corrections are included when CRYSTAL
        prints them, and the SCF/DFT energy is used only when no correction is
        present.
        """
        _, total, _ = CrystalElasticityReader._read_reference_energies(lines)
        return total

    @staticmethod
    def _read_reference_energies(
        lines: list[str],
    ) -> tuple[float, float, dict[str, object]]:
        """Return SCF and total energies for the unstrained elastic state.

        Parameters
        ----------
        lines : list of str
            CRYSTAL output restricted to the unstrained reference state.

        Returns
        -------
        scf_energy, total_energy, provenance : tuple
            Electronic SCF energy, physical total energy, and parser
            provenance.  ``NaN`` values and an empty mapping are returned when
            no energy state can be resolved.
        """
        totals = CrystalOutputParser(lines).total_energies()
        if not totals:
            return np.nan, np.nan, {}
        total = totals[-1]
        metadata = dict(total.metadata)
        raw_scf_energy = metadata.get("scf_energy", total.value)
        scf_energy = (
            float(raw_scf_energy)
            if isinstance(raw_scf_energy, (int, float))
            else float(total.value)
        )
        metadata["selected_quantity"] = "total_energy"
        return scf_energy, float(total.value), metadata

    @staticmethod
    def _elastic_reference_lines(lines: list[str]) -> list[str]:
        """Return output records belonging to the unstrained elastic state.

        CRYSTAL may print optimized geometries and stresses while processing
        strained configurations, especially when ``COORPRT`` is active.  The
        first ``STRAIN MATRIX`` therefore marks the end of the reference state
        from which the elastic distortions are generated.

        Parameters
        ----------
        lines : list of str
            Complete CRYSTAL elastic output.

        Returns
        -------
        list of str
            Prefix ending immediately before the first strained configuration.
        """
        for index, line in enumerate(lines):
            if markers.ELASTIC_STRAIN_MATRIX in line:
                return lines[:index]
        return lines

    @staticmethod
    def _read_reference_structure(
        lines: list[str],
        *,
        volume: float,
    ) -> tuple[CrystalStructure | None, int]:
        """Return the geometry used as the unstrained elastic reference.

        The latest complete geometry before the first strain is preferred, but
        a reported elastic volume is used to reject unrelated geometry records.
        This keeps ``COORPRT`` output from selecting a subsequently optimized
        strained configuration.

        Parameters
        ----------
        lines : list of str
            CRYSTAL output restricted to the unstrained elastic-reference part.
        volume : float
            Primitive-cell volume reported by the elastic module.

        Returns
        -------
        structure, start_index : tuple
            Selected structure and the line at which its geometry block starts.
            ``(None, 0)`` is returned when no complete geometry is available.
        """
        parser = CrystalGeometryParser(lines)
        candidates: list[tuple[int, CrystalStructure]] = []

        optimized_indices = [
            index
            for index, line in enumerate(lines)
            if markers.FINAL_OPTIMIZED_GEOMETRY in line
        ]
        optimized = parser.optimized_cells()
        candidates.extend(zip(optimized_indices, optimized, strict=True))

        wave_index = next(
            (
                index
                for index, line in enumerate(lines)
                if markers.GEOMETRY_WAVEFUNCTION in line
            ),
            None,
        )
        if wave_index is not None:
            try:
                candidates.append((wave_index, parser.wavefunction_cell()))
            except ValueError:
                pass

        if not candidates:
            try:
                candidates.append((0, parser.initial_primitive_cell()))
            except ValueError:
                return None, 0

        candidates.sort(key=lambda item: item[0])
        if np.isfinite(volume):
            matching = [
                item
                for item in candidates
                if np.isclose(
                    item[1].volume,
                    volume,
                    rtol=2.0e-6,
                    atol=1.0e-6,
                )
            ]
            if matching:
                return matching[-1][1], matching[-1][0]
        return candidates[-1][1], candidates[-1][0]
