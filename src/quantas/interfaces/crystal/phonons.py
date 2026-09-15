# -*- coding: utf-8 -*-

"""Read CRYSTAL phonon outputs for Quantas HA and QHA workflows."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Sequence

import numpy as np

from quantas.core.geometry import (
    analyze_symmetry,
    primitive_lattice_from_supercell,
    reconstruct_primitive_structure,
    supercell_repetitions,
)
from quantas.interfaces.crystal import markers, patterns
from quantas.interfaces.crystal.geometry import CrystalGeometryParser
from quantas.interfaces.crystal.output import CrystalOutputParser
from quantas.interfaces.crystal.phonon_modes import CrystalPhononModeParser
from quantas.models.reader import BasicReader
from quantas.models.structures import (
    CellNormalization,
    StructureReconstructionDiagnostics,
    StructureVolumeSeries,
)


class CrystalPhononReader(BasicReader):
    """Read a CRYSTAL phonon calculation for Quantas thermodynamics.

    The reader supports Gamma-only calculations, explicit phonon supercells, and
    ``SCELPHONO`` dispersion calculations. It preserves CRYSTAL-specific parsing at
    the interface boundary while exposing normalized energies, structures,
    q-point metadata, and frequencies to HA/QHA consumers. When a supercell is
    present, :attr:`energy`, :attr:`lattice`, and :attr:`natom` follow the physical
    normalization implied by the CRYSTAL calculation; :attr:`structure_series`
    provides the corresponding primitive structural representation and provenance.

    Parameters
    ----------
    crystal_output : str or pathlib.Path or None, optional
        CRYSTAL output file to load immediately. If ``None``, create an empty
        reader and call :meth:`load` later.

    Notes
    -----
    Ordinary recognition or completeness failures are reported through
    :attr:`~quantas.models.reader.BasicReader.error` and leave ``completed=False``.
    Malformed scientific records that cannot be interpreted unambiguously may
    raise ``ValueError`` from the specialized parser helpers."""

    _is_supercell = False
    _is_scelphono = False

    def __init__(self, crystal_output=None):
        """
        Initialize the CRYSTAL phonon reader.

        Parameters
        ----------
        crystal_output : str or pathlib.Path or None, optional
            CRYSTAL output file to load immediately. If ``None``, an empty
            reader is created and :meth:`load` can be called later.
        """
        super().__init__()
        self._is_supercell = False
        self._is_scelphono = False
        self._data = self._empty_data()

        if crystal_output is not None:
            self.load(crystal_output)
        return

    @staticmethod
    def _empty_data():
        """
        Return an empty data dictionary for a CRYSTAL phonon calculation.

        Returns
        -------
        dict
            Mutable per-instance storage used by the reader.
        """
        return {
            "unitcell": {},
            "supercell": {},
            "expansion": np.identity(3, dtype=int),
            "energy": 0.0,
            "scf_energy": np.nan,
            "energy_provenance": {},
            "kpoints": 1,
            "qpoints": 1,
            "qcoords": {},
            "nphonon": 0,
            "phonons": {},
            "weights": {},
            "shrinkf": np.ones(3, dtype=int),
            "q_position_source": "unavailable",
            "structure_series": None,
            "mode_data": None,
            "source_path": None,
        }

    def load(self, file):
        """Read one CRYSTAL phonon output into the reader state.

        The existing state is reset before parsing. Geometry, cell normalization,
        q-point metadata, authoritative total energy, phonon frequencies, and compact
        structural provenance are collected in one pass.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL phonon output file.

        Notes
        -----
        An unrecognized or incomplete phonon output is represented by ``error`` and
        ``completed=False``. Parser inconsistencies that would require guessing, such
        as contradictory q-point tables, are allowed to propagate as explicit errors."""
        self._data = self._empty_data()
        self._data["source_path"] = Path(file)
        geometry = CrystalGeometryParser(file)

        if not self.is_frequency_calculation(file):
            self.error = (
                "The file is not recognized as a CRYSTAL phonon output."
            )
            return

        self.supercell_on = self.is_supercell(file)
        self.scelphono_on = self.is_phonon_dispersion(file)
        #
        # Collect system information
        if self.supercell_on:
            self.supercell = self.set_wf_cell(file)
            self.unitcell = self.set_init_cell(file)
            self.dim = self.set_expansion(file)
            self._data["unitcell"]["lattice"] = primitive_lattice_from_supercell(
                self._data["supercell"]["lattice"],
                self.dim,
            )
            self.qpoints, self.qcoords, self.weights, self.shrinkf = self.set_q_mesh(
                file
            )

        else:
            self.unitcell = self.set_wf_cell(file)
            self.supercell = self.set_wf_cell(file)
            self.qcoords = [[0.0, 0.0, 0.0]]
            self.weights = [1]
            self._data["q_position_source"] = "gamma"
        # Collect energy values and phonons
        self._data["energy"] = self.set_energy(file)
        (
            self._data["scf_energy"],
            self._data["energy_provenance"],
        ) = self._resolve_energy_provenance(file, self._data["energy"])
        self.phonons = self.set_phonons(file)
        self._data["structure_series"] = self._build_structure_series(geometry)
        #
        self._check(file)
        return

    def _check(self, file):
        """ """
        if self.energy == 0.0:
            self.error = "No unit cell energy in {0}".format(file)
            return
        if not np.any(self.lattice):
            self.error = "No unit cell lattice in {0}".format(file)
            return
        if not bool(self.phonons):
            self.error = "No phonon data in {0}".format(file)
            return
        self.completed = True
        return

    def is_frequency_calculation(self, file):
        """Return whether an output contains a CRYSTAL frequency calculation.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL output file.

        Returns
        -------
        bool
            ``True`` when the CRYSTAL frequency-calculation marker is present."""
        with open(file, "r") as f:
            for line in f:
                if markers.FREQUENCY_CALCULATION in line:
                    return True
            return False

    def is_supercell(self, file):
        """Return whether CRYSTAL used an explicit phonon supercell.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL phonon output file.

        Returns
        -------
        bool
            ``True`` when the ``SUPERCELL`` option is present in the output."""
        with open(file, "r") as f:
            for line in f:
                if markers.SUPERCELL_OPTION in line:
                    return True
            return False

    @property
    def supercell_on(self):
        """
        Return the flag that tells if the input file is related to a supercell.
        """
        return self._is_supercell

    @supercell_on.setter
    def supercell_on(self, bool_value):
        """Store whether an explicit phonon supercell is active.

        Parameters
        ----------
        bool_value : bool
            Supercell-state flag."""
        self._is_supercell = bool_value
        return

    def is_phonon_dispersion(self, file):
        """Return whether CRYSTAL used ``SCELPHONO`` phonon dispersion.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL phonon output file.

        Returns
        -------
        bool
            ``True`` when the ``SCELPHONO`` option is present."""
        with open(file, "r") as f:
            for line in f:
                if markers.SCELPHONO_OPTION in line:
                    return True
            return False

    @property
    def scelphono_on(self):
        """
        Return the flag that tells if the input file is related to phonon
        dispersion relations calculation.
        """
        return self._is_scelphono

    @scelphono_on.setter
    def scelphono_on(self, bool_value):
        """Store whether the calculation uses ``SCELPHONO`` dispersion.

        Parameters
        ----------
        bool_value : bool
            Dispersion-state flag."""
        self._is_scelphono = bool_value
        return

    def is_hessian_interpolated(self, file):
        """Return whether CRYSTAL used Hessian interpolation.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL phonon output file.

        Returns
        -------
        bool
            ``True`` when the Hessian-interpolation marker is present."""
        with open(file, "r") as f:
            for line in f:
                if markers.HESSIAN_INTERPOLATION in line:
                    return True
            return False

    @property
    def natom(self):
        """Return the atom count represented by each stored phonon spectrum.

        For ``SCELPHONO`` dispersion and primitive Gamma calculations this is the
        primitive-cell atom count. For an explicit supercell Gamma calculation it is
        the supercell atom count, matching the number of frequencies actually stored."""
        if self.supercell_on:
            if self.scelphono_on:
                return self._data["unitcell"]["natom"]
            else:
                return self._data["supercell"]["natom"]
        else:
            return self._data["unitcell"]["natom"]

    @property
    def lattice(self):
        """Return the lattice associated with the stored phonon spectrum.

        Returns the primitive lattice for ``SCELPHONO`` dispersion and primitive
        Gamma calculations, and the explicit supercell lattice for a supercell Gamma
        calculation. Lattice vectors are expressed in angstrom."""
        if self.supercell_on:
            if self.scelphono_on:
                return self._data["unitcell"]["lattice"]
            else:
                return self._data["supercell"]["lattice"]
        else:
            return self._data["unitcell"]["lattice"]

    @property
    def volume(self):
        """Return the volume associated with :attr:`lattice` in ``angstrom^3``."""
        if self.supercell_on:
            if self.scelphono_on:
                return np.linalg.det(self._data["unitcell"]["lattice"])
            else:
                return np.linalg.det(self._data["supercell"]["lattice"])
        else:
            return np.linalg.det(self._data["unitcell"]["lattice"])

    @property
    def energy(self):
        """
        Return the total unit-cell or supercell energy used by the phonon state.

        CRYSTAL ``CENTRAL POINT`` energies are treated as the physical total
        energy and therefore include any empirical corrections already
        propagated by CRYSTAL.

        Returns
        -------
        float
            Total energy in hartree after the reader's cell normalization.
        """
        if self.supercell_on:
            if self.scelphono_on:
                return self._data["energy"] / self.kpoints
            else:
                return self._data["energy"]
        else:
            return self._data["energy"] / self.kpoints

    @property
    def total_energy(self):
        """Return the total phonon-reference energy.

        Returns
        -------
        float
            Total energy in hartree after the reader's cell normalization.
        """
        return self.energy

    @property
    def scf_energy(self):
        """Return the uncorrected SCF energy associated with the reference.

        Returns
        -------
        float
            Electronic SCF energy in hartree after the reader's cell
            normalization, or ``NaN`` when the central point could not be
            associated unambiguously with an SCF state.
        """
        value = float(self._data["scf_energy"])
        if not np.isfinite(value):
            return np.nan
        if self.supercell_on:
            if self.scelphono_on:
                return value / self.kpoints
            return value
        return value / self.kpoints

    @property
    def energy_provenance(self) -> dict[str, object]:
        """Return provenance for SCF and total phonon-reference energies.

        Returns
        -------
        dict
            Source markers, correction labels, and SCF/total matching
            diagnostics.
        """
        return dict(self._data["energy_provenance"])

    @property
    def nphonon(self):
        """Return the number of stored phonon branches per q-point.

        The value is ``3 * natom`` using the same primitive/supercell normalization as
        :attr:`natom`."""
        if self.supercell_on:
            if self.scelphono_on:
                return self._data["unitcell"]["natom"] * 3
            else:
                return self._data["supercell"]["natom"] * 3
        else:
            return self._data["unitcell"]["natom"] * 3

    @property
    def unitcell(self):
        """
        Return the crystal unit cell in tuple format.
        """
        return (
            self._data["unitcell"]["natom"],
            self._data["unitcell"]["numbers"],
            self._data["unitcell"]["positions"],
            self._data["unitcell"]["lattice"],
        )

    @unitcell.setter
    def unitcell(self, cell_data):
        """Store primitive/unit-cell data in the historical tuple layout.

        Parameters
        ----------
        cell_data : tuple
            ``(natom, numbers, positions, lattice)`` with fractional positions and
            lattice vectors in angstrom."""
        self._data["unitcell"]["natom"] = cell_data[0]
        self._data["unitcell"]["numbers"] = cell_data[1]
        self._data["unitcell"]["positions"] = cell_data[2]
        self._data["unitcell"]["lattice"] = cell_data[3]
        return

    @property
    def supercell(self):
        """
        Return the crystal supercell in tuple format.
        """
        return (
            self._data["supercell"]["natom"],
            self._data["supercell"]["numbers"],
            self._data["supercell"]["positions"],
            self._data["supercell"]["lattice"],
        )

    @supercell.setter
    def supercell(self, cell_data):
        """Store supercell data in the historical tuple layout.

        Parameters
        ----------
        cell_data : tuple
            ``(natom, numbers, positions, lattice)`` with fractional positions and
            lattice vectors in angstrom."""
        self._data["supercell"]["natom"] = cell_data[0]
        self._data["supercell"]["numbers"] = cell_data[1]
        self._data["supercell"]["positions"] = cell_data[2]
        self._data["supercell"]["lattice"] = cell_data[3]
        return

    @property
    def dim(self):
        """
        Return the expansion matrix employed to build the supercell.
        """
        return self._data["expansion"]

    @dim.setter
    def dim(self, expansion):
        """Store the primitive-to-supercell expansion matrix.

        Parameters
        ----------
        expansion : array-like
            ``(3, 3)`` expansion matrix. A copy is retained by the reader."""
        self._data["expansion"] = expansion.copy()
        return

    @property
    def units(self) -> dict[str, str]:
        """Return the physical units exposed by this CRYSTAL reader.

        Returns
        -------
        dict
            Energy, volume, frequency, and structural length units.
        """
        return {
            "energy": "Ha",
            "volume": "angstrom^3",
            "frequency": "cm^-1",
            "length": "angstrom",
        }

    @property
    def kpoints(self):
        """Return the number of primitive-cell repetitions in the expansion.

        The value is the rounded determinant of the ``(3, 3)`` expansion matrix and is
        used to normalize supercell quantities where appropriate."""
        return int(np.around(np.linalg.det(self._data["expansion"]), 0))

    @property
    def qpoints(self):
        """
        Return the number of sampled **q**-points.
        """
        return self._data["qpoints"]

    @qpoints.setter
    def qpoints(self, value: int):
        """Store the number of sampled phonon q-points.

        Parameters
        ----------
        value : int
            Number of q-points."""
        self._data["qpoints"] = value
        return

    @property
    def qcoords(self):
        """Return CRYSTAL q-point coordinate numerators keyed by q-point index.

        Use :attr:`qcoords_fractional` for primitive reciprocal fractional coordinates.
        For Gamma-only calculations the sole coordinate is ``(0, 0, 0)``."""
        return self._data["qcoords"]

    @qcoords.setter
    def qcoords(self, array):
        """Store CRYSTAL q-point coordinate numerators.

        Parameters
        ----------
        array : array-like
            Coordinate rows ordered by q-point index, with shape ``(qpoints, 3)``."""
        for i in range(len(array)):
            self._data["qcoords"][i] = array[i]
        return

    @property
    def weights(self):
        """Return CRYSTAL q-point integration weights keyed by q-point index."""
        return self._data["weights"]

    @weights.setter
    def weights(self, array):
        """Store q-point integration weights.

        Parameters
        ----------
        array : array-like
            Weight values ordered by q-point index."""
        for i in range(len(array)):
            self._data["weights"][i] = array[i]
        return

    @property
    def shrinkf(self):
        """Return the three CRYSTAL reciprocal-space shrinking factors."""
        return self._data["shrinkf"]

    @shrinkf.setter
    def shrinkf(self, array):
        """Store the CRYSTAL reciprocal shrinking factors.

        Parameters
        ----------
        array : array-like
            Three shrinking factors. A copy is retained by the reader."""
        self._data["shrinkf"] = array.copy()
        return

    @property
    def q_position_source(self) -> str:
        """Return the origin of the q-point coordinates.

        Returns
        -------
        str
            Parser provenance for the q-point coordinate metadata.
        """
        return str(self._data.get("q_position_source", "unavailable"))

    @property
    def qcoords_fractional(self) -> np.ndarray:
        """Return fractional primitive reciprocal q-point coordinates.

        Returns
        -------
        numpy.ndarray
            Array with shape ``(qpoints, 3)``.  CRYSTAL prints integer
            coordinate numerators and three shrinking factors; this property
            performs the component-wise division once, at the interface
            boundary.

        Raises
        ------
        ValueError
            If the shrinking factors are missing or non-positive.
        """
        shrinkf = np.asarray(self.shrinkf, dtype=np.float64)
        if shrinkf.shape != (3,) or np.any(shrinkf <= 0.0):
            raise ValueError("phonon shrinking factors must contain three positives")
        return self.qcoords_array() / shrinkf[np.newaxis, :]

    @property
    def phonons(self):
        """Return phonon frequencies keyed by q-point index.

        Each value contains ``nphonon`` frequencies in ``cm^-1``."""
        return self._data["phonons"]

    @phonons.setter
    def phonons(self, dictionary):
        """Store phonon-frequency blocks keyed by q-point index.

        Parameters
        ----------
        dictionary : dict
            Mapping whose values contain frequencies in ``cm^-1``."""
        self._data["phonons"] = dictionary
        return

    def qcoords_array(self):
        """
        Return q-point coordinates as an ordered NumPy array.

        Returns
        -------
        numpy.ndarray
            Array with shape ``(qpoints, 3)``.
        """
        return np.asarray(
            [self._data["qcoords"][i] for i in range(self.qpoints)],
            dtype=np.float64,
        )

    def weights_array(self):
        """
        Return q-point weights as an ordered NumPy array.

        Returns
        -------
        numpy.ndarray
            Array with shape ``(qpoints,)``.
        """
        return np.asarray(
            [self._data["weights"][i] for i in range(self.qpoints)],
            dtype=np.float64,
        )

    def phonons_array(self):
        """
        Return phonon frequencies as an ordered NumPy array.

        Returns
        -------
        numpy.ndarray
            Array with shape ``(qpoints, nphonon)``.
        """
        return np.asarray(
            [self._data["phonons"][i] for i in range(self.qpoints)],
            dtype=np.float64,
        )

    @property
    def mode_data(self):
        """Return parsed mass-weighted phonon eigenvectors when available.

        Returns
        -------
        PhononModeData or None
            Frequencies and unit-norm eigenvectors, or ``None`` when CRYSTAL
            did not print eigenvectors in the source output. The data are
            parsed lazily and cached because HA thermodynamics do not require
            eigenvectors unless mode continuity is inspected.

        Raises
        ------
        ValueError
            If CRYSTAL prints an incomplete or inconsistent eigenvector block.
        """
        cached = self._data.get("mode_data")
        if cached is not None:
            return cached
        source = self._data.get("source_path")
        if source is None:
            return None
        parsed = CrystalPhononModeParser(source).parse(self.nphonon)
        self._data["mode_data"] = parsed
        return parsed

    @property
    def structure_series(self):
        """Return the compact primitive structural series when available.

        Returns
        -------
        StructureVolumeSeries or None
            One-volume primitive structural path and its normalization data.
        """
        return self._data.get("structure_series")

    def _build_structure_series(self, geometry: CrystalGeometryParser):
        """Build compact primitive structural metadata for this calculation."""
        reference = geometry.initial_primitive_cell()
        source = geometry.wavefunction_cell()
        expansion = np.asarray(self.dim, dtype=np.int64)
        repetitions = supercell_repetitions(expansion)
        if repetitions > 1 and source.natoms == reference.natoms * repetitions:
            primitive, diagnostics = reconstruct_primitive_structure(
                source,
                expansion,
                reference,
            )
        elif source.natoms == reference.natoms:
            primitive = source
            diagnostics = StructureReconstructionDiagnostics(
                status="exact",
                source_atoms=source.natoms,
                reconstructed_atoms=source.natoms,
                expected_repetitions=1,
                minimum_replica_count=1,
                maximum_replica_count=1,
                maximum_translation_residual=0.0,
                rms_translation_residual=0.0,
                message="Source and compact primitive cells coincide.",
            )
        else:
            raise ValueError(
                "CRYSTAL source and primitive atom counts are inconsistent "
                f"({source.natoms} versus {reference.natoms} x {repetitions})"
            )
        symmetry = analyze_symmetry(primitive, symprec=1.0e-5)
        basis = (
            "phonon_supercell"
            if self.supercell_on and not self.scelphono_on
            else "primitive"
        )
        normalization = CellNormalization(
            basis=basis,
            source_basis=("wavefunction_supercell" if repetitions > 1 else "primitive"),
            expansion_matrix=expansion,
            repetitions=repetitions,
            source_atoms=source.natoms,
            normalized_atoms=self.natom,
        )
        return StructureVolumeSeries(
            reference=primitive,
            lattices=np.asarray([primitive.lattice], dtype=np.float64),
            fractional_positions=np.asarray(
                [primitive.fractional_positions],
                dtype=np.float64,
            ),
            volumes=np.asarray([primitive.volume], dtype=np.float64),
            normalization=normalization,
            symmetry=symmetry,
            primitive_to_crystallographic=(geometry.primitive_to_crystallographic()),
            diagnostics=(diagnostics,),
            orientation="crystal",
            reference_index=0,
            metadata={
                "interface": "crystal",
                "coorprt_present": geometry.has_coorprt,
                "space_group_number_from_output": geometry.space_group_number(),
            },
            source_lattices=np.asarray([source.lattice], dtype=np.float64),
            source_fractional_positions=(source.fractional_positions.copy(),),
        )

    def set_init_cell(self, file):
        """Return the primitive input geometry parsed from CRYSTAL output.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL output file.

        Returns
        -------
        tuple
            Atom count, atomic numbers, fractional positions, and lattice.
        """
        structure = CrystalGeometryParser(file).initial_primitive_cell()
        return (
            structure.natoms,
            structure.atomic_numbers.copy(),
            structure.fractional_positions.copy(),
            structure.lattice.copy(),
        )

    def set_wf_cell(self, file):
        """Return the CRYSTAL wave-function geometry.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL output file.

        Returns
        -------
        tuple
            Atom count, atomic numbers, fractional positions, and lattice.
        """
        structure = CrystalGeometryParser(file).wavefunction_cell()
        return (
            structure.natoms,
            structure.atomic_numbers.copy(),
            structure.fractional_positions.copy(),
            structure.lattice.copy(),
        )

    def set_expansion(self, file):
        """Return the CRYSTAL primitive-to-supercell expansion matrix.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL phonon output file.

        Returns
        -------
        numpy.ndarray
            ``(3, 3)`` expansion matrix printed by CRYSTAL. The determinant gives the
            number of primitive-cell repetitions represented by the supercell."""
        sline = self._get_start_line(file, markers.SUPERCELL_EXPANSION) + 1

        with open(file, "r") as f:
            data = f.readlines()

        expansion = np.zeros((3, 3), dtype=float)

        for i in range(3):
            line = data[sline + i].split()
            del line[0]
            expansion[i] = np.asarray(line, dtype=float)
        return expansion

    def set_q_mesh(
        self, file: str | Path
    ) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
        """Read the phonon q-point mesh used by CRYSTAL.

        The preferred source is the explicit table headed by ``K WEIGHT
        COORD``.  Coordinates in that table are integer numerators and are
        converted to fractional reciprocal coordinates only by
        :attr:`qcoords_fractional`, using the three printed shrinking factors.
        The per-q-point ``DISPERSION K POINT NUMBER`` records are parsed as an
        independent consistency check and as a fallback for older outputs.

        Parameters
        ----------
        file : str or pathlib.Path
            Path of the CRYSTAL phonon output.

        Returns
        -------
        tuple
            Number of q-points, integer coordinate numerators, q-point
            weights, and shrinking factors.

        Raises
        ------
        ValueError
            If the printed q-point tables are incomplete or mutually
            inconsistent.
        """
        with open(file, "r", encoding="utf-8", errors="replace") as stream:
            lines = stream.readlines()

        table = self._parse_dispersion_qpoint_table(lines)
        markers = self._parse_dispersion_qpoint_markers(lines)

        if table is not None:
            qcoords, qweights, shrinkf = table
            if markers is not None:
                marker_coords, marker_weights = markers
                if marker_coords.shape != qcoords.shape or not np.array_equal(
                    marker_coords, qcoords
                ):
                    raise ValueError(
                        "CRYSTAL q-point list and dispersion sections use "
                        "different coordinates or ordering"
                    )
                if not np.allclose(marker_weights, qweights, rtol=0.0, atol=1.0e-12):
                    raise ValueError(
                        "CRYSTAL q-point list and dispersion sections use "
                        "different weights"
                    )
            self._data["q_position_source"] = "crystal-dispersion-table"
            return len(qweights), qcoords, qweights, shrinkf

        if markers is not None:
            qcoords, qweights = markers
            marker_shrinkf = self._parse_shrinking_factors(lines)
            if marker_shrinkf is None:
                raise ValueError(
                    "CRYSTAL dispersion q-points were found, but their "
                    "shrinking factors are missing"
                )
            self._data["q_position_source"] = "crystal-dispersion-sections"
            return len(qweights), qcoords, qweights, marker_shrinkf

        qpoints, qcoords, qweights, shrinkf = self._set_q_mesh_legacy(lines)
        self._data["q_position_source"] = "crystal-legacy-qmesh"
        return qpoints, qcoords, qweights, shrinkf

    @staticmethod
    def _parse_dispersion_qpoint_table(
        lines: Sequence[str],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
        """Parse CRYSTAL's explicit thermodynamic q-point table.

        Parameters
        ----------
        lines : sequence of str
            Complete output file split into lines.

        Returns
        -------
        tuple or None
            Integer coordinates, weights, and shrinking factors, or ``None``
            when the table is absent.

        Raises
        ------
        ValueError
            If a table header is present but its records are malformed.
        """
        header_index = next(
            (
                index
                for index, line in enumerate(lines)
                if "K       WEIGHT       COORD" in line
            ),
            None,
        )
        if header_index is None:
            return None

        record_pattern = re.compile(
            r"^\s*\*?\s*(\d+)\s+"
            r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+"
            r"([+-]?\d+)\s+([+-]?\d+)\s+([+-]?\d+)"
        )
        records = []
        shrinkf = None
        for line in lines[header_index + 1 :]:
            if "WITH SHRINKING FACTORS" in line:
                shrinkf = CrystalPhononReader._shrinking_factors_from_line(line)
                break
            match = record_pattern.match(line)
            if match is None:
                continue
            records.append(
                (
                    int(match.group(1)),
                    float(match.group(2)),
                    int(match.group(3)),
                    int(match.group(4)),
                    int(match.group(5)),
                )
            )

        if not records:
            raise ValueError("CRYSTAL q-point table contains no readable records")
        if shrinkf is None:
            raise ValueError("CRYSTAL q-point table has no shrinking factors")

        expected = list(range(1, len(records) + 1))
        indices = [record[0] for record in records]
        if indices != expected:
            raise ValueError("CRYSTAL q-point table indices are not consecutive")

        weights = np.asarray([record[1] for record in records], dtype=np.float64)
        coordinates = np.asarray(
            [[record[2], record[3], record[4]] for record in records],
            dtype=np.float64,
        )
        return coordinates, weights, shrinkf

    @staticmethod
    def _parse_dispersion_qpoint_markers(
        lines: Sequence[str],
    ) -> tuple[np.ndarray, np.ndarray] | None:
        """Parse q-point coordinates repeated before each dispersion block.

        Parameters
        ----------
        lines : sequence of str
            Complete output file split into lines.

        Returns
        -------
        tuple or None
            Integer coordinates and weights ordered by q-point index, or
            ``None`` when no dispersion records are present.

        Raises
        ------
        ValueError
            If records are duplicated, incomplete, or out of order.
        """
        pattern = re.compile(
            r"DISPERSION K POINT NUMBER\s+(\d+)\s+"
            r"COORD:\s+\w\(\s*([+-]?\d+)\s+([+-]?\d+)\s+"
            r"([+-]?\d+)\s*\)\s+WEIGHT:\s*"
            r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))"
        )
        records: dict[int, tuple[int, int, int, float]] = {}
        for line in lines:
            match = pattern.search(line)
            if match is None:
                continue
            index = int(match.group(1))
            value = (
                int(match.group(2)),
                int(match.group(3)),
                int(match.group(4)),
                float(match.group(5)),
            )
            if index in records and records[index] != value:
                raise ValueError(
                    f"CRYSTAL q-point {index} is printed with inconsistent data"
                )
            records[index] = value

        if not records:
            return None
        expected = list(range(1, max(records) + 1))
        if sorted(records) != expected:
            raise ValueError("CRYSTAL dispersion q-point indices are incomplete")
        coordinates = np.asarray(
            [
                [records[index][0], records[index][1], records[index][2]]
                for index in expected
            ],
            dtype=np.float64,
        )
        weights = np.asarray(
            [records[index][3] for index in expected],
            dtype=np.float64,
        )
        return coordinates, weights

    @staticmethod
    def _shrinking_factors_from_line(line: str) -> np.ndarray:
        """Return the three integer shrinking factors from one CRYSTAL line."""
        match = re.search(
            r"IS1\s*=\s*(\d+)\s+IS2\s*=\s*(\d+)\s+"
            r"IS3\s*=\s*(\d+)",
            line,
        )
        if match is None:
            raise ValueError("Unable to parse CRYSTAL shrinking factors")
        values = np.asarray(
            [int(match.group(1)), int(match.group(2)), int(match.group(3))],
            dtype=np.float64,
        )
        if np.any(values <= 0.0):
            raise ValueError("CRYSTAL shrinking factors must be positive")
        return values

    @staticmethod
    def _parse_shrinking_factors(lines: Sequence[str]) -> np.ndarray | None:
        """Find shrinking factors anywhere in a CRYSTAL phonon output."""
        for line in lines:
            if "WITH SHRINKING FACTORS" in line:
                return CrystalPhononReader._shrinking_factors_from_line(line)
        return None

    def _set_q_mesh_legacy(
        self, data: Sequence[str]
    ) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
        """Retain support for historical INTERPHESS output layouts."""
        hess = True
        sline = next(
            (index for index, line in enumerate(data) if markers.HESSIAN_INTERPOLATION in line),
            None,
        )
        if sline is None:
            hess = False
            sline = next(
                (index for index, line in enumerate(data) if markers.SCELPHONO_QPOINTS in line),
                None,
            )
        if sline is None:
            raise ValueError("Unable to locate CRYSTAL phonon q-point metadata")

        qpoints = int(data[sline].split()[-4])
        qcoords = np.zeros((qpoints, 3), dtype=np.float64)
        qmesh = np.zeros(3, dtype=np.float64)
        qweight = np.zeros(qpoints, dtype=np.float64)
        if hess:
            sline += 9
        else:
            table_line = next(
                (
                    index
                    for index, line in enumerate(data)
                    if "K       WEIGHT       COORD" in line
                ),
                None,
            )
            if table_line is None:
                raise ValueError("Unable to locate CRYSTAL q-point table")
            sline = table_line + 1

        for index in range(qpoints):
            fields = data[sline + index].split()
            qweight[index] = float(fields[2])
            qcoords[index] = np.asarray(fields[3:6], dtype=np.float64)

        shrink_line = data[sline + qpoints].split()
        qmesh[0] = float(shrink_line[6])
        qmesh[1] = float(shrink_line[9])
        qmesh[2] = float(shrink_line[12])
        return qpoints, qcoords, qweight, qmesh

    def set_energy(self, file):
        """Return the authoritative CRYSTAL central-point phonon energy.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL phonon output file.

        Returns
        -------
        float
            Central-point total energy in hartree, in the cell normalization printed
            by CRYSTAL. The public :attr:`energy` property applies the reader's
            primitive/supercell normalization.

        Raises
        ------
        ValueError
            If the central-point energy cannot be identified unambiguously."""
        with open(file, "r") as f:
            for line in f:
                match = patterns.CENTRAL_POINT_RE.search(line)
                if match is not None:
                    return float(
                        match.group("energy")
                        .replace("D", "E")
                        .replace("d", "e")
                    )

        raise ValueError("CRYSTAL central-point energy not found")

    @staticmethod
    def _resolve_energy_provenance(
        file: str | Path,
        central_point_energy: float,
    ) -> tuple[float, dict[str, object]]:
        """Associate a phonon central point with its SCF/total energy state.

        The ``CENTRAL POINT`` value remains authoritative for phonon input
        generation because it is the energy attached by CRYSTAL to the
        undisplaced reference configuration.  When the generic CRYSTAL parser
        can match that value to an earlier SCF state, Quantas additionally
        records the uncorrected SCF energy and the empirical-correction
        provenance without changing the numerical value used by HA/QHA.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL phonon output.
        central_point_energy : float
            Reference energy parsed from ``CENTRAL POINT`` in hartree.

        Returns
        -------
        scf_energy, provenance : tuple
            Matched uncorrected SCF energy and provenance mapping.  The SCF
            value is ``NaN`` when no total-energy state matches the central
            point within printing precision.
        """
        total_records = CrystalOutputParser(file).total_energies()
        provenance: dict[str, object] = {
            "selected_quantity": "total_energy",
            "source_marker": "CENTRAL POINT",
            "central_point_energy_hartree": float(central_point_energy),
            "matched_scf_state": False,
            "corrections": [],
        }
        if not total_records:
            return np.nan, provenance

        # In a standalone CRYSTAL FREQCALC output the first SCF state is the
        # undisplaced input configuration from which the central point is
        # constructed.  Do not search later displaced states merely because
        # one happens to have a numerically similar energy.
        matched = total_records[0]
        if not np.isclose(
            float(matched.value),
            float(central_point_energy),
            rtol=1.0e-11,
            atol=1.0e-8,
        ):
            return np.nan, provenance

        metadata = dict(matched.metadata)
        raw_scf_energy = metadata.get("scf_energy", matched.value)
        scf_energy = (
            float(raw_scf_energy)
            if isinstance(raw_scf_energy, (int, float))
            else float(matched.value)
        )
        correction_labels = metadata.get("corrections", ())
        if not isinstance(correction_labels, (list, tuple)):
            correction_labels = ()
        provenance.update(
            {
                "matched_scf_state": True,
                "scf_energy_hartree": scf_energy,
                "resolved_total_energy_hartree": float(matched.value),
                "resolved_total_source_marker": metadata.get("source_marker"),
                "corrections": list(correction_labels),
                "total_correction_energy_hartree": metadata.get(
                    "total_correction_energy"
                ),
            }
        )
        return scf_energy, provenance

    def set_phonons(self, file):
        """Parse phonon-frequency blocks from a CRYSTAL output file.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL frequency-calculation output file.

        Returns
        -------
        dict[int, numpy.ndarray]
            Mapping from sequential q-point block index to a one-dimensional
            array of ``self.nphonon`` frequencies in cm^-1.

        Raises
        ------
        OSError
            If the output file cannot be opened.
        ValueError
            If a parsed mode index or frequency is not numeric.
        """
        phonons = {}
        band_counter = 0

        with open(file, "r") as f:
            data = f.readlines()

        for i in range(len(data)):
            if markers.FREQUENCY_HEADER in data[i]:
                finished = False
                band = np.zeros(self.nphonon, dtype=float)
                line_counter = 0
                freq_counter = 0

                while not finished or freq_counter != self.nphonon:
                    line = data[i + line_counter + 2].split()

                    if len(line) == 0:
                        finished = True
                        continue

                    mode_n1 = int(line[0][:-1])
                    mode_n2 = int(line[1])

                    for j in range(mode_n2 - mode_n1 + 1):
                        band[freq_counter] = float(line[3])
                        freq_counter += 1

                    line_counter += 1

                phonons[band_counter] = band
                band_counter += 1

                if band_counter == self.qpoints:
                    break

        return phonons

    def _get_start_line(self, file, search_string):
        found = False
        sline = 0
        with open(file, "r") as f:
            for line in f:
                if search_string in line:
                    found = True
                    break
                else:
                    sline += 1
        if found:
            return sline
        else:
            return None
