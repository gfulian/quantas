# -*- coding: utf-8 -*-

"""Read CRYSTAL phonon outputs for Quantas HA and QHA workflows."""

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
from quantas.models.reader import BasicReader
from quantas.models.structures import (
    CellNormalization,
    StructureReconstructionDiagnostics,
    StructureVolumeSeries,
)


class CrystalPhononReader(BasicReader):
    """
    This class reads and stores phonon properties of crystals obtained
    through the CRYSTAL14/17 code.

    Parameters
    ----------

    file: str
        Path to the CRYSTAL output file.

    """

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
            "kpoints": 1,
            "qpoints": 1,
            "qcoords": {},
            "nphonon": 0,
            "phonons": {},
            "weights": {},
            "shrinkf": np.ones(3, dtype=int),
            "q_position_source": "unavailable",
            "structure_series": None,
        }

    def load(self, file):
        """
        Read and store the information collected in the CRYSTAL output file.

        Parameters
        ----------

        file: str
            Path to the CRYSTAL output file.

        """
        self._data = self._empty_data()
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
        """
        This method checks if the CRYSTAL14/17 output file is related to
        a FREQCALC calculation.

        Parameters
        ----------

        file: str
            Path to the CRYSTAL14/17 output file.


        Returns
        -------

        bool
            Returns True if the output is correct, otherwise False.

        """
        with open(file, "r") as f:
            for line in f:
                if markers.FREQUENCY_CALCULATION in line:
                    return True
            return False

    def is_supercell(self, file):
        """
        This method checks if the CRYSTAL14/17 output file is related to
        a FREQCALC calculation using the a supercell approach.

        Parameters
        ----------

        file: str
            Path to the CRYSTAL14/17 output file.


        Returns
        -------

        bool
            Returns True if the output is correct, otherwise False.

        """
        with open(file, "r") as f:
            for line in f:
                if markers.SUPERCELL_OPTION in line:
                    return True
            return False

    @property
    def supercell_on(self):
        """
        Get the flag that tells if the input file is related to a supercell.
        """
        return self._is_supercell

    @supercell_on.setter
    def supercell_on(self, bool_value):
        """
        Set the flag that tells if the input file is related to a supercell.
        """
        self._is_supercell = bool_value
        return

    def is_phonon_dispersion(self, file):
        """
        This method checks if the CRYSTAL14/17 output file is related to
        a FREQCALC calculation using the SCELPHONO keyword (phonon
        dispersion relations).

        Parameters
        ----------

        file: str
            Path to the CRYSTAL14/17 output file.


        Returns
        -------

        bool
            Returns True if the output is correct, otherwise False.

        """
        with open(file, "r") as f:
            for line in f:
                if markers.SCELPHONO_OPTION in line:
                    return True
            return False

    @property
    def scelphono_on(self):
        """
        Get the flag that tells if the input file is related to phonon
        dispersion relations calculation.
        """
        return self._is_scelphono

    @scelphono_on.setter
    def scelphono_on(self, bool_value):
        """
        Set the flag that tells if the input file is related to phonon
        dispersion relations calculation.
        """
        self._is_scelphono = bool_value
        return

    def is_hessian_interpolated(self, file):
        """
        This method checks if FREQCALC calculation used Hessian interpolation
        scheme.

        Parameters
        ----------

        file: str
            Path to the CRYSTAL14/17 output file.


        Returns
        -------

        bool
            Returns True if the output is correct, otherwise False.

        """
        with open(file, "r") as f:
            for line in f:
                if markers.HESSIAN_INTERPOLATION in line:
                    return True
            return False

    @property
    def natom(self):
        """
        Get the number of atoms in the unit cell (if phonon dispersion
        relations or if :math:`\\Gamma`-point frequencies) or in the
        supercell.
        """
        if self.supercell_on:
            if self.scelphono_on:
                return self._data["unitcell"]["natom"]
            else:
                return self._data["supercell"]["natom"]
        else:
            return self._data["unitcell"]["natom"]

    @property
    def lattice(self):
        """
        Get the unit cell (if phonon dispersion relations or if
        :math:`\\Gamma`-point frequencies) or the supercell lattice vectors.
        """
        if self.supercell_on:
            if self.scelphono_on:
                return self._data["unitcell"]["lattice"]
            else:
                return self._data["supercell"]["lattice"]
        else:
            return self._data["unitcell"]["lattice"]

    @property
    def volume(self):
        """
        Get the unit cell (if phonon dispersion relations or if
        :math:`\\Gamma`-point frequencies) or the supercell volume.
        """
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
        Get the unit cell (if phonon dispersion relations or if
        :math:`\\Gamma`-point frequencies) or the supercell energy.
        """
        if self.supercell_on:
            if self.scelphono_on:
                return self._data["energy"] / self.kpoints
            else:
                return self._data["energy"]
        else:
            return self._data["energy"] / self.kpoints

    @property
    def nphonon(self):
        """
        Get the number of frequencies per band in unit cell (if phonon
        dispersion relations or if :math:`\\Gamma`-point frequencies)
        or in the supercell.
        """
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
        Get the crystal unit cell in tuple format.
        """
        return (
            self._data["unitcell"]["natom"],
            self._data["unitcell"]["numbers"],
            self._data["unitcell"]["positions"],
            self._data["unitcell"]["lattice"],
        )

    @unitcell.setter
    def unitcell(self, cell_data):
        """
        Set the crystal unit cell in tuple format.
        """
        self._data["unitcell"]["natom"] = cell_data[0]
        self._data["unitcell"]["numbers"] = cell_data[1]
        self._data["unitcell"]["positions"] = cell_data[2]
        self._data["unitcell"]["lattice"] = cell_data[3]
        return

    @property
    def supercell(self):
        """
        Get the crystal supercell in tuple format.
        """
        return (
            self._data["supercell"]["natom"],
            self._data["supercell"]["numbers"],
            self._data["supercell"]["positions"],
            self._data["supercell"]["lattice"],
        )

    @supercell.setter
    def supercell(self, cell_data):
        """
        Set the crystal supercell in tuple format.
        """
        self._data["supercell"]["natom"] = cell_data[0]
        self._data["supercell"]["numbers"] = cell_data[1]
        self._data["supercell"]["positions"] = cell_data[2]
        self._data["supercell"]["lattice"] = cell_data[3]
        return

    @property
    def dim(self):
        """
        Get the expansion matrix employed to build the supercell.
        """
        return self._data["expansion"]

    @dim.setter
    def dim(self, expansion):
        """
        Set the expansion matrix employed to build the supercell.
        """
        self._data["expansion"] = expansion.copy()
        return

    @property
    def kpoints(self):
        """
        Get the number of sampled *k*-points, determined from the expansion
        matrix.
        """
        return int(np.around(np.linalg.det(self._data["expansion"]), 0))

    @property
    def qpoints(self):
        """
        Get the number of sampled **q**-points.
        """
        return self._data["qpoints"]

    @qpoints.setter
    def qpoints(self, value: int):
        """
        Set the number of sampled **q**-points.
        """
        self._data["qpoints"] = value
        return

    @property
    def qcoords(self):
        """
        Get the coordinates of sampled **q**-points, in dict format.
        """
        return self._data["qcoords"]

    @qcoords.setter
    def qcoords(self, array):
        """
        Set the coordinates of sampled **q**-points, in dict format.
        """
        for i in range(len(array)):
            self._data["qcoords"][i] = array[i]
        return

    @property
    def weights(self):
        """
        Get the weights of each phonon band, in dict format.
        """
        return self._data["weights"]

    @weights.setter
    def weights(self, array):
        """
        Set the weights of each phonon band.
        """
        for i in range(len(array)):
            self._data["weights"][i] = array[i]
        return

    @property
    def shrinkf(self):
        """
        Get the Hessian interpolation mesh used in the INTERPHESS keyword.
        """
        return self._data["shrinkf"]

    @shrinkf.setter
    def shrinkf(self, array):
        """
        Set the Hessian interpolation mesh used in the INTERPHESS keyword.
        """
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
        """
        Get the phonon bands, in dict format.
        """
        return self._data["phonons"]

    @phonons.setter
    def phonons(self, dictionary):
        """
        Set the phonon bands, in dict format.
        """
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
        """
        This method sets the expasion matrix used to build the supercell.

        Thus, it could be:

          - a unit cell, if the CRYSTAL output is related to a
            :math:`\\Gamma`-point frequency calculation, or

          - a supercell, if either SUPERCELL of SCELPHONO were employed.

        Parameters
        ----------

        file: str
            Path of the CRYSTAL14/17 output file.

        Returns
        -------

        expansion: ndarray
            :math:`3 \\times 3` array of the expansion matrix.

        """
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
        """
        This method sets the energy of the cell. It reads the value from the
        central point (equilibrium) of displacement.

        Parameters
        ----------

        file: str
            Path of the CRYSTAL14/17 output file.

        Returns
        -------

        float
            Energy of the crystal cell.

        """
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

    def set_phonons(self, file):
        """ """
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
