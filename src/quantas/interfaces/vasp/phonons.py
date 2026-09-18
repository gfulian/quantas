# -*- coding: utf-8 -*-

"""Read primitive-cell Gamma phonons from VASP run outputs."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Any

import numpy as np
from numpy.typing import NDArray

from quantas.core.chemistry import atomic_mass, number2symbol
from quantas.core.geometry import analyze_symmetry, reduce_to_primitive_cell
from quantas.interfaces.vasp.output import VaspOutputParser
from quantas.models.phonons import PhononModeData
from quantas.models.reader import BasicReader
from quantas.models.structures import (
    CellNormalization,
    StructureReconstructionDiagnostics,
    StructureVolumeSeries,
)


FloatArray = NDArray[np.float64]

_MODE_RE = re.compile(
    r"^\s*(?P<index>\d+)\s+(?P<kind>f(?:/i)?)\s*=\s*"
    r"(?P<thz>[-+0-9.EeDd]+)\s+THz\s+"
    r"[-+0-9.EeDd]+\s+2PiTHz\s+"
    r"(?P<wavenumber>[-+0-9.EeDd]+)\s+cm-1\s+"
    r"(?P<mev>[-+0-9.EeDd]+)\s+meV\s*$"
)
_MODE_HEADER = "Eigenvectors and eigenvalues of the dynamical matrix"
_TRANSLATION_SCORE_MIN = 0.90
_TRANSLATION_SEPARATION_MIN = 0.20


class VaspPhononReader(BasicReader[None]):
    """Read one primitive-cell VASP Gamma-point phonon calculation.

    The current interface deliberately supports only the Gamma point of the
    primitive cell used by VASP.  It does not reconstruct primitive-cell
    dispersions from force constants computed in a supercell and it does not
    parse VASP 6 ``LPHON_DISPERSION``/``QPOINTS`` output.  A source cell that is
    reducible to more than one primitive-cell repetition is therefore rejected
    rather than silently interpreting folded supercell modes as a primitive
    Gamma-only calculation.

    VASP reports normalized eigenvectors of the mass-weighted dynamical matrix.
    Quantas stores the same unit-norm, mass-weighted representation used by the
    backend-neutral :class:`~quantas.models.phonons.PhononModeData` contract.
    The three rigid translations are identified from their projection onto the
    mass-weighted translational subspace, not from an arbitrary frequency
    cutoff.  Their small numerical Gamma frequencies are retained in metadata
    but set to exactly zero in the thermodynamic frequency array.

    Parameters
    ----------
    source : str or pathlib.Path
        VASP calculation directory, ``vasprun.xml``, or sibling ``OUTCAR``.
    symprec : float, optional
        Cartesian spglib tolerance in angstrom used to verify that the source
        cell is primitive.

    Raises
    ------
    ValueError
        If the run is incomplete, contains no readable dynamical matrix, does
        not describe a primitive-cell Gamma calculation, or has inconsistent
        frequencies/eigenvectors.
    ImportError
        If spglib is unavailable for primitive-cell validation.
    """

    def __init__(self, source: str | Path, *, symprec: float = 1.0e-5) -> None:
        super().__init__()
        self.symprec = float(symprec)
        self._data: dict[str, Any] = {}
        self.load(source)

    def load(self, filename: str | Path) -> None:
        """Load and validate one VASP Gamma-point phonon run.

        Parameters
        ----------
        filename : str or pathlib.Path
            VASP calculation directory, ``vasprun.xml``, or sibling ``OUTCAR``.
        """
        self.completed = False
        self.error = None
        self._data = {}
        try:
            self._load(Path(filename))
        except (OSError, ValueError, ImportError) as exc:
            self.error = str(exc)
            return
        self.completed = True

    def _load(self, source: Path) -> None:
        parser = VaspOutputParser(source)
        if parser.termination().status.value != "normal":
            raise ValueError("VASP phonon source did not terminate normally")

        structure = parser.initial_structure()
        reduction = reduce_to_primitive_cell(
            structure,
            symprec=self.symprec,
            no_idealize=True,
        )
        if reduction.repetitions != 1:
            raise ValueError(
                "VASP phonon source is not a primitive-cell Gamma calculation; "
                "direct dispersion/folded-supercell reconstruction from VASP "
                "outputs is not implemented yet"
            )

        document = parser.document
        dynmats = document.root.findall(".//dynmat")
        if len(dynmats) != 1:
            raise ValueError(
                "VASP Gamma phonon reader requires exactly one dynamical-matrix "
                f"record; found {len(dynmats)}"
            )
        dynmat = dynmats[0]
        eigenvectors = _xml_eigenvectors(dynmat, structure.natoms)
        raw_frequencies = _outcar_gamma_frequencies(
            document.outcar_text,
            expected_modes=3 * structure.natoms,
        )
        if eigenvectors.shape[0] != raw_frequencies.size:
            raise ValueError(
                "VASP phonon frequencies and dynamical-matrix eigenvectors "
                "contain different mode counts"
            )

        atom_symbols = tuple(number2symbol(int(z)) for z in structure.atomic_numbers)
        eigenvectors = _normalized_eigenvectors(eigenvectors)
        translation_scores = _translation_projection_scores(
            eigenvectors,
            atom_symbols,
        )
        translation_indices = _translation_indices(translation_scores)
        frequencies = raw_frequencies.copy()
        frequencies[list(translation_indices)] = 0.0

        steps = parser.ionic_steps()
        if not steps:
            raise ValueError("VASP phonon source contains no reference electronic state")
        reference_step = steps[0]
        reference_energy = reference_step.energies.sigma_zero_energy
        symmetry = analyze_symmetry(structure, symprec=self.symprec)

        mode_data = PhononModeData(
            frequencies=frequencies[np.newaxis, :],
            eigenvectors=eigenvectors[np.newaxis, :, :, :],
            atom_symbols=atom_symbols,
            frequency_unit="cm^-1",
            eigenvector_normalization="mass-weighted-unit",
            metadata={
                "interface": "vasp",
                "vasp_version": parser.version,
                "ibrion": document.parameter("IBRION"),
                "source": str(parser.source.directory),
                "q_point": [0.0, 0.0, 0.0],
                "source_eigenvectors": "vasprun.xml/dynmat/eigenvectors",
                "source_frequencies": "OUTCAR Gamma dynamical-matrix block",
                "source_eigenvector_coordinates": "Cartesian",
                "raw_frequencies_cm^-1": raw_frequencies.tolist(),
                "translation_indices": list(translation_indices),
                "translation_projection_scores": translation_scores.tolist(),
                "translation_policy": "mass-weighted-rigid-translation-subspace",
                "dispersion_support": "not_implemented",
            },
        )

        normalization = CellNormalization(
            basis="primitive",
            source_basis="vasp-primitive",
            expansion_matrix=np.eye(3, dtype=np.int64),
            repetitions=1,
            source_atoms=structure.natoms,
            normalized_atoms=structure.natoms,
        )
        diagnostics = StructureReconstructionDiagnostics(
            status="exact",
            source_atoms=structure.natoms,
            reconstructed_atoms=structure.natoms,
            expected_repetitions=1,
            minimum_replica_count=1,
            maximum_replica_count=1,
            maximum_translation_residual=0.0,
            rms_translation_residual=0.0,
            message="VASP source cell verified as primitive; source basis preserved.",
        )
        structure_series = StructureVolumeSeries(
            reference=structure,
            lattices=structure.lattice[np.newaxis, :, :],
            fractional_positions=structure.fractional_positions[np.newaxis, :, :],
            volumes=np.asarray([structure.volume], dtype=np.float64),
            normalization=normalization,
            symmetry=symmetry,
            diagnostics=(diagnostics,),
            orientation="vasp",
            reference_index=0,
            metadata={
                "interface": "vasp",
                "gamma_only": True,
                "dispersion_support": "not_implemented",
            },
        )

        self._data = {
            "source": parser.source,
            "version": parser.version,
            "ibrion": document.parameter("IBRION"),
            "structure": structure,
            "structure_series": structure_series,
            "energy": float(reference_energy.value),
            "energy_provenance": {
                "selected_quantity": "sigma_zero_energy",
                "source_marker": "vasprun.xml:first_reference_state:e_0_energy",
                "source_ionic_step": int(reference_step.index),
                "corrections": [],
                "unit": reference_energy.unit,
            },
            "mode_data": mode_data,
            "phonons": {0: frequencies.copy()},
            "raw_frequencies": raw_frequencies,
            "translation_indices": translation_indices,
            "translation_scores": translation_scores,
        }

    @property
    def natom(self) -> int:
        """Return the primitive-cell atom count."""
        return int(self._data["structure"].natoms)

    @property
    def lattice(self) -> FloatArray:
        """Return the primitive VASP lattice in angstrom."""
        return np.asarray(self._data["structure"].lattice, dtype=np.float64).copy()

    @property
    def volume(self) -> float:
        """Return the primitive-cell volume in cubic angstrom."""
        return float(self._data["structure"].volume)

    @property
    def energy(self) -> float:
        """Return the reference ``energy(sigma->0)`` in eV."""
        return float(self._data["energy"])

    @property
    def scf_energy(self) -> float:
        """Return ``NaN`` because no distinct HA SCF energy is selected."""
        return float("nan")

    @property
    def energy_provenance(self) -> dict[str, Any]:
        """Return VASP reference-energy provenance."""
        return dict(self._data["energy_provenance"])

    @property
    def nphonon(self) -> int:
        """Return the number of Gamma phonon branches."""
        return int(self.natom * 3)

    @property
    def qpoints(self) -> int:
        """Return one because this interface is currently Gamma-only."""
        return 1

    @property
    def qcoords(self) -> FloatArray:
        """Return the single Gamma q-point in fractional reciprocal coordinates."""
        return np.zeros((1, 3), dtype=np.float64)

    @property
    def qcoords_fractional(self) -> FloatArray:
        """Return the single Gamma q-point in fractional reciprocal coordinates."""
        return self.qcoords

    @property
    def q_position_source(self) -> str:
        """Return the explicit source convention for the Gamma coordinate."""
        return "vasp-gamma-only"

    @property
    def weights(self) -> FloatArray:
        """Return the unit integration weight for Gamma."""
        return np.ones(1, dtype=np.float64)

    @property
    def shrinkf(self) -> NDArray[np.int64]:
        """Return identity reciprocal shrinking factors for Gamma-only data."""
        return np.ones(3, dtype=np.int64)

    @property
    def dim(self) -> NDArray[np.int64]:
        """Return the identity phonon-supercell matrix."""
        return np.eye(3, dtype=np.int64)

    @property
    def phonons(self) -> dict[int, FloatArray]:
        """Return Gamma phonon frequencies in wavenumbers."""
        return {0: np.asarray(self._data["phonons"][0], dtype=np.float64).copy()}

    def phonons_array(self) -> FloatArray:
        """Return Gamma phonon frequencies as one q-point array.

        Returns
        -------
        ndarray
            Frequencies with shape ``(1, 3*natom)`` in ``cm^-1``.
        """
        return np.asarray([self._data["phonons"][0]], dtype=np.float64)

    @property
    def mode_data(self) -> PhononModeData:
        """Return Gamma frequencies and unit-norm mass-weighted eigenvectors."""
        return self._data["mode_data"]

    @property
    def structure_series(self) -> StructureVolumeSeries:
        """Return the one-volume primitive structural path."""
        return self._data["structure_series"]

    @property
    def units(self) -> dict[str, str]:
        """Return explicit VASP phonon input units."""
        return {
            "energy": "eV",
            "volume": "angstrom^3",
            "frequency": "cm^-1",
            "length": "angstrom",
        }

    @property
    def translation_indices(self) -> tuple[int, int, int]:
        """Return zero-based indices classified as rigid translations."""
        values = self._data["translation_indices"]
        return int(values[0]), int(values[1]), int(values[2])

    @property
    def translation_projection_scores(self) -> FloatArray:
        """Return rigid-translation subspace projection scores for all modes."""
        return np.asarray(self._data["translation_scores"], dtype=np.float64).copy()


def _outcar_gamma_frequencies(text: str | None, *, expected_modes: int) -> FloatArray:
    """Return signed Gamma frequencies in wavenumbers from one OUTCAR block."""
    if text is None:
        raise ValueError(
            "VASP Gamma phonon parsing currently requires OUTCAR for signed "
            "frequency labels"
        )
    start = text.rfind(_MODE_HEADER)
    if start < 0:
        raise ValueError("VASP OUTCAR contains no dynamical-matrix mode block")
    frequencies: list[float] = []
    indices: list[int] = []
    for line in text[start:].splitlines()[2:]:
        match = _MODE_RE.match(line)
        if match is None:
            if frequencies and len(frequencies) >= expected_modes:
                break
            continue
        indices.append(int(match.group("index")))
        value = float(match.group("wavenumber").replace("D", "E").replace("d", "e"))
        if match.group("kind") == "f/i":
            value = -abs(value)
        frequencies.append(value)
        if len(frequencies) == expected_modes:
            break
    if indices != list(range(1, expected_modes + 1)):
        raise ValueError(
            "VASP OUTCAR Gamma phonon modes are missing, duplicated, or out of order"
        )
    return np.asarray(frequencies, dtype=np.float64)


def _xml_eigenvectors(node: Any, natoms: int) -> FloatArray:
    """Return VASP dynamical-matrix eigenvectors with shape ``(3N, N, 3)``."""
    varray = node.find("varray[@name='eigenvectors']")
    if varray is None:
        raise ValueError("VASP vasprun.xml dynmat contains no eigenvectors")
    rows: list[list[float]] = []
    for row in varray.findall("v"):
        if row.text is None:
            raise ValueError("VASP dynamical-matrix eigenvector row is empty")
        rows.append([float(value) for value in row.text.split()])
    expected_modes = 3 * natoms
    array = np.asarray(rows, dtype=np.float64)
    if array.shape != (expected_modes, expected_modes):
        raise ValueError(
            "VASP dynamical-matrix eigenvectors must have shape "
            f"({expected_modes}, {expected_modes}); observed {array.shape}"
        )
    return array.reshape(expected_modes, natoms, 3)


def _normalized_eigenvectors(eigenvectors: FloatArray) -> FloatArray:
    """Return unit-norm VASP mass-weighted dynamical-matrix eigenvectors."""
    array = np.asarray(eigenvectors, dtype=np.float64)
    flat = array.reshape(array.shape[0], -1)
    norms = np.linalg.norm(flat, axis=1)
    if np.any(~np.isfinite(norms)) or np.any(norms <= np.finfo(np.float64).eps):
        raise ValueError("VASP phonon eigenvector has zero or invalid norm")
    return array / norms[:, np.newaxis, np.newaxis]


def _translation_projection_scores(
    eigenvectors: FloatArray,
    atom_symbols: tuple[str, ...],
) -> FloatArray:
    """Return projection of each mode onto the mass-weighted translation space."""
    masses = np.asarray([atomic_mass(symbol) for symbol in atom_symbols], dtype=np.float64)
    natoms = len(atom_symbols)
    basis = np.zeros((3, natoms, 3), dtype=np.float64)
    root_mass = np.sqrt(masses)
    for axis in range(3):
        basis[axis, :, axis] = root_mass
        basis[axis] /= np.linalg.norm(basis[axis])
    flat_modes = np.asarray(eigenvectors, dtype=np.float64).reshape(eigenvectors.shape[0], -1)
    flat_basis = basis.reshape(3, -1)
    projections = flat_modes @ flat_basis.T
    return np.sqrt(np.sum(np.square(projections), axis=1))


def _translation_indices(scores: FloatArray) -> tuple[int, int, int]:
    """Return three modes robustly identified as rigid translations."""
    values = np.asarray(scores, dtype=np.float64)
    if values.ndim != 1 or values.size < 3:
        raise ValueError("VASP phonon translation classification requires at least 3 modes")
    order = np.argsort(values)[::-1]
    selected = np.sort(order[:3])
    minimum = float(np.min(values[selected]))
    if minimum < _TRANSLATION_SCORE_MIN:
        raise ValueError(
            "VASP Gamma phonons do not contain three well-resolved rigid "
            f"translations (minimum projection score {minimum:.6f})"
        )
    if values.size > 3:
        fourth = float(values[order[3]])
        if minimum - fourth < _TRANSLATION_SEPARATION_MIN:
            raise ValueError(
                "VASP Gamma rigid translations are not sufficiently separated "
                "from optical modes"
            )
    return int(selected[0]), int(selected[1]), int(selected[2])


__all__ = ["VaspPhononReader"]
