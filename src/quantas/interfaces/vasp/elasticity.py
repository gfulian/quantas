# -*- coding: utf-8 -*-

"""VASP OUTCAR reader for second-order elastic constants."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, TypedDict

import numpy as np
from numpy.typing import NDArray

from quantas.models.elastic_states import (
    ElasticTensorKind,
    PressureSource,
    PrestressProvenance,
)
from quantas.models.reader import BasicReader

from .output import VaspOutputParser

_CLAMPED_ELASTIC_MODULI = "SYMMETRIZED ELASTIC MODULI (kBar)"
_IONIC_ELASTIC_MODULI = "ELASTIC MODULI CONTR FROM IONIC RELAXATION (kBar)"
_RELAXED_ELASTIC_MODULI = "TOTAL ELASTIC MODULI (kBar)"
_ATOMIC_MASS_UNIT_PER_ANGSTROM_CUBED = 1660.53906660
_VASP_VOIGT_LABELS = ("XX", "YY", "ZZ", "XY", "YZ", "ZX")
_QUANTAS_VOIGT_LABELS = ("XX", "YY", "ZZ", "YZ", "ZX", "XY")
_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"

_POMASS_RE = re.compile(
    rf"POMASS\s*=\s*({_FLOAT})\s*;\s*ZVAL",
    flags=re.IGNORECASE,
)
_IONS_PER_TYPE_RE = re.compile(
    r"ions\s+per\s+type\s*=\s*([0-9 \t]+)",
    flags=re.IGNORECASE,
)
_VOLUME_RE = re.compile(
    rf"volume\s+of\s+cell\s*:\s*({_FLOAT})",
    flags=re.IGNORECASE,
)
_EXTERNAL_PRESSURE_RE = re.compile(
    rf"external\s+pressure\s*=\s*({_FLOAT})\s*kB",
    flags=re.IGNORECASE,
)
_IN_KB_RE = re.compile(
    rf"^\s*in\s+kB\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s+"
    rf"({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s*$",
    flags=re.IGNORECASE,
)


class _ElasticityData(TypedDict):
    """Typed payload stored by the VASP elasticity interface reader."""

    stiffness: NDArray[np.float64]
    clamped_stiffness: NDArray[np.float64]
    ionic_relaxation_stiffness: NDArray[np.float64] | None
    relaxed_stiffness: NDArray[np.float64] | None
    density: float
    reference_volume_angstrom3: float | None
    reference_stress_gpa: NDArray[np.float64] | None
    reference_pressure_gpa: float | None
    ibrion: int | None
    isif: int | None
    selected_level: str
    source: Path | None


class VASPElasticityReader(BasicReader[None]):
    """Read second-order elastic constants from a VASP calculation.

    The reader accepts a VASP calculation directory, its ``OUTCAR`` path, or
    ``vasprun.xml`` when a sibling ``OUTCAR`` is available.  VASP's clamped-ion,
    ionic-relaxation, and total elastic-moduli blocks are preserved separately.
    The relaxed-ion ``TOTAL ELASTIC MODULI`` block is selected when present;
    otherwise the clamped-ion block is used.

    VASP reports elastic components in the order ``XX YY ZZ XY YZ ZX``.  Quantas
    uses ``11 22 33 23 13 12`` (``XX YY ZZ YZ ZX XY``), so both tensor axes are
    reordered explicitly by label before values are exposed.  Values are converted
    from kbar to GPa.

    The reference stress/pressure is collected from the first unstrained stress
    record in the OUTCAR using VASP's sign convention (positive compression).  No
    CRYSTAL-specific finite-prestress transformation is applied.  Current VASP
    documentation identifies the reported moduli as finite-difference strain--stress
    derivatives / energy second derivatives, but does not establish equivalence to
    Quantas' named finite-prestress conventions.  The reader therefore records the
    tensor kind as :class:`~quantas.models.elastic_states.ElasticTensorKind.UNKNOWN`
    until that convention is validated independently.

    Parameters
    ----------
    filename : str or pathlib.Path or None, optional
        VASP run directory, ``OUTCAR``, or sibling ``vasprun.xml`` to load.
    """

    def __init__(self, filename: str | Path | None = None) -> None:
        super().__init__()
        self._data = self._empty_data()
        if filename is not None:
            self.load(filename)

    @staticmethod
    def _empty_data() -> _ElasticityData:
        """Return an empty reader payload."""
        zeros = np.zeros((6, 6), dtype=np.float64)
        return {
            "stiffness": zeros.copy(),
            "clamped_stiffness": zeros.copy(),
            "ionic_relaxation_stiffness": None,
            "relaxed_stiffness": None,
            "density": 0.0,
            "reference_volume_angstrom3": None,
            "reference_stress_gpa": None,
            "reference_pressure_gpa": None,
            "ibrion": None,
            "isif": None,
            "selected_level": "unavailable",
            "source": None,
        }

    @property
    def stiffness(self) -> NDArray[np.float64]:
        """Return the selected elastic stiffness matrix in GPa."""
        return self._data["stiffness"].copy()

    @property
    def clamped_stiffness(self) -> NDArray[np.float64]:
        """Return the clamped-ion VASP stiffness matrix in GPa."""
        return self._data["clamped_stiffness"].copy()

    @property
    def ionic_relaxation_stiffness(self) -> NDArray[np.float64] | None:
        """Return the ionic-relaxation contribution in GPa when available."""
        matrix = self._data["ionic_relaxation_stiffness"]
        return None if matrix is None else matrix.copy()

    @property
    def relaxed_stiffness(self) -> NDArray[np.float64] | None:
        """Return the total relaxed-ion VASP stiffness matrix in GPa when available."""
        matrix = self._data["relaxed_stiffness"]
        return None if matrix is None else matrix.copy()

    @property
    def selected_level(self) -> str:
        """Return ``"relaxed_ion"`` or ``"clamped_ion"`` for the selected tensor."""
        return self._data["selected_level"]

    @property
    def density(self) -> float:
        """Return the reference-cell density in kg m^-3 when available."""
        return float(self._data["density"])

    @property
    def reference_volume_angstrom3(self) -> float | None:
        """Return the unstrained reference-cell volume in angstrom cubed."""
        volume = self._data["reference_volume_angstrom3"]
        return None if volume is None else float(volume)

    @property
    def reference_stress_gpa(self) -> NDArray[np.float64] | None:
        """Return the unstrained VASP stress tensor in GPa, positive in compression."""
        stress = self._data["reference_stress_gpa"]
        return None if stress is None else stress.copy()

    @property
    def reference_pressure_gpa(self) -> float | None:
        """Return the unstrained hydrostatic pressure in GPa, positive in compression."""
        pressure = self._data["reference_pressure_gpa"]
        return None if pressure is None else float(pressure)

    @property
    def ibrion(self) -> int | None:
        """Return the effective ``IBRION`` value parsed from OUTCAR when available."""
        return self._data["ibrion"]

    @property
    def isif(self) -> int | None:
        """Return the effective ``ISIF`` value parsed from OUTCAR when available."""
        return self._data["isif"]

    @property
    def tensor_kind(self) -> ElasticTensorKind:
        """Return the current conservative Quantas tensor-convention classification."""
        return ElasticTensorKind.RAW_STRESS_STRAIN

    @property
    def prestress(self) -> PrestressProvenance:
        """Return pressure and finite-prestress provenance for the selected tensor."""
        pressure = self.reference_pressure_gpa
        if pressure is None:
            return PrestressProvenance(
                tensor_kind=ElasticTensorKind.RAW_STRESS_STRAIN,
                pressure_source=PressureSource.UNAVAILABLE,
            )
        return PrestressProvenance(
            tensor_kind=ElasticTensorKind.RAW_STRESS_STRAIN,
            pressure_gpa=pressure,
            pressure_source=PressureSource.OUTPUT_STRESS,
        )

    @property
    def metadata(self) -> dict[str, Any]:
        """Return VASP elasticity parsing and convention provenance."""
        return {
            "interface": "vasp",
            "source": None if self._data["source"] is None else str(self._data["source"]),
            "selected_level": self.selected_level,
            "source_voigt_order": "XX YY ZZ XY YZ ZX",
            "quantas_voigt_order": "11 22 33 23 13 12",
            "tensor_definition": "vasp-finite-difference-strain-stress",
            "tensor_kind": self.tensor_kind.value,
            "quantas_prestress_correction_applied": False,
            "finite_prestress_semantics": "raw-stress-strain-requires-pressure-adjustment",
            "stress_sign_convention": "positive-compression",
            "reference_stress_gpa": (
                None
                if self.reference_stress_gpa is None
                else self.reference_stress_gpa.tolist()
            ),
            "ibrion": self.ibrion,
            "isif": self.isif,
            "reference_pressure_gpa": self.reference_pressure_gpa,
            "reference_volume_angstrom3": self.reference_volume_angstrom3,
        }

    def load(self, filename: str | Path) -> None:
        """Read elastic constants from a VASP run directory or OUTCAR.

        Parameters
        ----------
        filename : str or pathlib.Path
            Calculation directory, ``OUTCAR``, or sibling ``vasprun.xml``.

        Notes
        -----
        Filesystem errors propagate.  Scientific format errors set ``completed`` to
        ``False`` and populate ``error``.
        """
        self.completed = False
        self.error = None
        self._data = self._empty_data()
        path = self._resolve_outcar(filename)
        self._data["source"] = path
        text = path.read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()

        try:
            clamped = self._read_elastic_table(lines, _CLAMPED_ELASTIC_MODULI)
        except ValueError as exc:
            self.error = str(exc)
            return
        if clamped is None:
            self.error = (
                f"'{path}' does not appear to be a VASP output containing "
                "elastic moduli."
            )
            return

        try:
            ionic = self._read_elastic_table(lines, _IONIC_ELASTIC_MODULI)
            relaxed = self._read_elastic_table(lines, _RELAXED_ELASTIC_MODULI)
        except ValueError as exc:
            self.error = str(exc)
            return

        self._data["clamped_stiffness"] = clamped
        self._data["ionic_relaxation_stiffness"] = ionic
        self._data["relaxed_stiffness"] = relaxed
        if relaxed is not None:
            self._data["stiffness"] = relaxed.copy()
            self._data["selected_level"] = "relaxed_ion"
        else:
            self._data["stiffness"] = clamped.copy()
            self._data["selected_level"] = "clamped_ion"

        stress, pressure = self._read_reference_stress(lines)
        self._data["reference_stress_gpa"] = stress
        self._data["reference_pressure_gpa"] = pressure
        volume = self._read_reference_volume(filename, text)
        self._data["reference_volume_angstrom3"] = volume
        self._data["density"] = self._read_density(text, volume=volume)
        self._data["ibrion"] = self._read_integer_parameter(lines, "IBRION")
        self._data["isif"] = self._read_integer_parameter(lines, "ISIF")
        self.completed = True

    @staticmethod
    def _resolve_outcar(source: str | Path) -> Path:
        """Resolve a VASP run directory or sibling file to ``OUTCAR``.

        Parameters
        ----------
        source : str or pathlib.Path
            Directory, ``OUTCAR``, or ``vasprun.xml`` path.

        Returns
        -------
        pathlib.Path
            Existing OUTCAR path.

        Raises
        ------
        FileNotFoundError
            If the source or required OUTCAR does not exist.
        ValueError
            If a non-VASP filename is supplied.
        """
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(f"VASP elasticity source does not exist: {path}")
        if path.is_dir():
            outcar = path / "OUTCAR"
        elif path.name == "OUTCAR":
            outcar = path
        elif path.name.casefold() == "vasprun.xml":
            outcar = path.parent / "OUTCAR"
        else:
            # Preserve compatibility with historical callers that used renamed
            # OUTCAR excerpts in tests or local workflows.
            outcar = path
        if not outcar.is_file():
            raise FileNotFoundError(f"VASP elasticity source requires OUTCAR: {outcar}")
        return outcar

    @staticmethod
    def _read_elastic_table(
        lines: list[str], header: str
    ) -> NDArray[np.float64] | None:
        """Return one VASP elastic-moduli table in Quantas Voigt order.

        Parameters
        ----------
        lines : list of str
            OUTCAR lines.
        header : str
            Exact VASP table heading.

        Returns
        -------
        ndarray or None
            ``6 x 6`` matrix in GPa and Quantas order, or ``None`` when absent.

        Raises
        ------
        ValueError
            If a present table has malformed labels, shape, or values.
        """
        indices = [index for index, line in enumerate(lines) if header in line]
        if not indices:
            return None
        start = indices[-1]
        direction_index = next(
            (
                index
                for index in range(start + 1, min(start + 6, len(lines)))
                if lines[index].strip().startswith("Direction")
            ),
            None,
        )
        if direction_index is None:
            raise ValueError(f"Malformed VASP elastic-moduli header: {header}")
        column_labels = tuple(lines[direction_index].split()[1:])
        if set(column_labels) != set(_VASP_VOIGT_LABELS) or len(column_labels) != 6:
            raise ValueError(
                f"Unexpected VASP elastic-moduli column labels for {header}: "
                f"{' '.join(column_labels)}"
            )

        rows: dict[str, list[float]] = {}
        for line in lines[direction_index + 1 : direction_index + 12]:
            parts = line.split()
            if not parts or parts[0] not in _VASP_VOIGT_LABELS:
                continue
            if len(parts) != 7:
                raise ValueError(f"Malformed VASP elastic-moduli row: {line.strip()}")
            try:
                rows[parts[0]] = [_float(value) for value in parts[1:]]
            except ValueError as exc:
                raise ValueError(
                    f"Invalid VASP elastic-moduli value in row {parts[0]}"
                ) from exc
        if set(rows) != set(_VASP_VOIGT_LABELS):
            raise ValueError(f"Incomplete VASP elastic-moduli table: {header}")

        source = np.asarray(
            [
                [
                    rows[row][column_labels.index(column)]
                    for column in _VASP_VOIGT_LABELS
                ]
                for row in _VASP_VOIGT_LABELS
            ],
            dtype=np.float64,
        )
        source /= 10.0  # kbar -> GPa
        order = [_VASP_VOIGT_LABELS.index(label) for label in _QUANTAS_VOIGT_LABELS]
        matrix = source[np.ix_(order, order)]
        if not np.all(np.isfinite(matrix)) or not np.allclose(
            matrix, matrix.T, rtol=0.0, atol=1.0e-8
        ):
            raise ValueError(f"VASP elastic-moduli table is not finite and symmetric: {header}")
        return matrix

    @staticmethod
    def _read_reference_stress(
        lines: list[str],
    ) -> tuple[NDArray[np.float64] | None, float | None]:
        """Return the first unstrained VASP stress and hydrostatic pressure.

        VASP prints stress components as ``XX YY ZZ XY YZ ZX`` in kbar with
        positive values denoting compression.  The first stress block precedes the
        finite lattice distortions used for elastic constants and therefore
        represents the reference state.
        """
        for line in lines:
            match = _IN_KB_RE.match(line)
            if match is None:
                continue
            values = np.asarray([_float(value) for value in match.groups()]) / 10.0
            xx, yy, zz, xy, yz, zx = values
            stress = np.asarray(
                [[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]], dtype=np.float64
            )
            return stress, float(np.trace(stress) / 3.0)

        for line in lines:
            match = _EXTERNAL_PRESSURE_RE.search(line)
            if match is not None:
                return None, _float(match.group(1)) / 10.0
        return None, None

    @staticmethod
    def _read_reference_volume(source: str | Path, text: str) -> float | None:
        """Return the unstrained VASP cell volume in angstrom cubed.

        A sibling ``vasprun.xml`` is preferred because it retains the lattice at
        full precision.  Standalone OUTCAR inputs fall back to the first reported
        cell volume.
        """
        path = Path(source)
        xml = (path / "vasprun.xml") if path.is_dir() else path.parent / "vasprun.xml"
        if path.name.casefold() == "vasprun.xml":
            xml = path
        if xml.is_file():
            try:
                return float(VaspOutputParser(xml).initial_structure().volume)
            except (FileNotFoundError, ValueError):
                pass
        matches = _VOLUME_RE.findall(text)
        if not matches:
            return None
        try:
            return _float(matches[0])
        except ValueError:
            return None

    @staticmethod
    def _read_density(text: str, *, volume: float | None = None) -> float:
        """Return the VASP reference-cell density in kg m^-3 when available.

        Parameters
        ----------
        text : str
            OUTCAR text containing masses and species populations.
        volume : float or None, optional
            Preferred reference volume in angstrom cubed.  When unavailable, the
            first OUTCAR cell volume is used.

        Returns
        -------
        float
            Density in kg m^-3, or zero when required metadata are unavailable.

        Notes
        -----
        Later ``IBRION=6`` volume records include finite trial distortions and must
        not redefine the reference state.
        """
        mass_matches = _POMASS_RE.findall(text)
        count_matches = _IONS_PER_TYPE_RE.findall(text)
        volume_matches = _VOLUME_RE.findall(text)
        if not mass_matches or not count_matches or (volume is None and not volume_matches):
            return 0.0

        try:
            masses = np.asarray([_float(value) for value in mass_matches], dtype=float)
            counts = np.asarray(
                [int(value) for value in count_matches[0].split()], dtype=int
            )
            reference_volume = _float(volume_matches[0]) if volume is None else float(volume)
        except ValueError:
            return 0.0

        if (
            masses.shape != counts.shape
            or masses.size == 0
            or np.any(~np.isfinite(masses))
            or np.any(masses <= 0.0)
            or np.any(counts <= 0)
            or not np.isfinite(reference_volume)
            or reference_volume <= 0.0
        ):
            return 0.0

        total_mass = float(np.dot(masses, counts))
        density = total_mass / reference_volume * _ATOMIC_MASS_UNIT_PER_ANGSTROM_CUBED
        return float(density) if np.isfinite(density) and density > 0.0 else 0.0

    @staticmethod
    def _read_integer_parameter(lines: list[str], name: str) -> int | None:
        """Return one integer VASP run parameter from OUTCAR when available."""
        pattern = re.compile(rf"^\s*{re.escape(name)}\s*=\s*(-?\d+)\b")
        for line in lines:
            match = pattern.match(line)
            if match is not None:
                return int(match.group(1))
        return None

    def is_elasticity_output(self, filename: str | Path) -> bool:
        """Return whether a VASP source contains the clamped-ion elastic table.

        Parameters
        ----------
        filename : str or pathlib.Path
            VASP calculation directory, OUTCAR, or sibling ``vasprun.xml``.

        Returns
        -------
        bool
            ``True`` when the clamped-ion elastic-moduli table is present.
        """
        path = self._resolve_outcar(filename)
        with path.open("r", encoding="utf-8", errors="replace") as stream:
            return any(_CLAMPED_ELASTIC_MODULI in line for line in stream)

    def elasticity_start_line(self, filename: str | Path) -> int:
        """Return the first row of the preferred VASP stiffness table.

        Parameters
        ----------
        filename : str or pathlib.Path
            VASP calculation directory, OUTCAR, or sibling ``vasprun.xml``.

        Returns
        -------
        int
            Zero-based line index of the first row in the preferred elastic table.

        Raises
        ------
        ValueError
            If the source contains no recognized elastic-moduli table.

        Notes
        -----
        This compatibility helper retains the historical reader API.  New code
        should use :attr:`stiffness` and the explicit clamped/relaxed properties.
        """
        path = self._resolve_outcar(filename)
        clamped_line: int | None = None
        relaxed_line: int | None = None
        with path.open("r", encoding="utf-8", errors="replace") as stream:
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


def _float(value: str) -> float:
    """Parse a VASP floating-point token including Fortran ``D`` exponents."""
    return float(value.replace("D", "E").replace("d", "e"))


__all__ = ["VASPElasticityReader"]
