# -*- coding: utf-8 -*-

"""Read CRYSTAL quasi-harmonic outputs and structural volume paths."""

import numpy as np

from quantas.core.geometry import (
    analyze_symmetry,
    reconstruct_primitive_structure,
    supercell_repetitions,
)
from quantas.interfaces.crystal.geometry import CrystalGeometryParser
from quantas.models.reader import BasicReader
from quantas.models.structures import (
    CellNormalization,
    StructureVolumeSeries,
)

# Geometry strings
geometry_consistent = "GEOMETRY NOW FULLY CONSISTENT WITH THE GROUP"
geometry_sym_changed = "SYMMETRY CHANGED DURING OLD RUN :"
geometry_wf = "GEOMETRY FOR WAVE FUNCTION - "
primitive_cell = "PRIMITIVE CELL - "
supercell_option = " * SUPERCELL OPTION"
supercell_expansion = "EXPANSION MATRIX OF PRIMITIVE CELL"
lattice_vectors = "DIRECT LATTICE VECTORS CARTESIAN COMPONENTS"
# FREQCALC-specific strings
frequency_calculation = "EIGENVALUES (EIGV) OF THE MASS WEIGHTED HESSIAN"
frequency_calculation += " MATRIX AND HARMONIC"
scelphono_option = "PHONON FREQUENCIES AT A SET OF K POINTS BY USING "
scelphono_option += "A SUPERCELL"
hess_interpolation = "ACTIVATED INTERPOLATION OF THE HESSIAN UP TO"
scelphono_qpoints = "THAT PERMITS THE CALCULATION OF MODES AT"
central_energy = "    CENTRAL POINT"
frequency_header = "MODES         EIGV          FREQUENCIES     IRREP"
# QHA-specific strings
qha_header = "QUASI-HARMONIC APPROXIMATION"
qha_freq = "FREQUENCY #"
qha_ev = "SORTING VOLUMES/ENERGIES"
qha_opt = "FINAL OPTIMIZED GEOMETRY"
qha_restart = "READING DATA FROM RESTART UNIT"


class CrystalQHAReader(BasicReader):
    """
    Reader for CRYSTAL17 output file obtained via the QHA keyword.

    .. seealso:: CRYSTAL17 tutorial on QHA_ calculation.

    .. _QHA: http://tutorials.crystalsolutions.eu/tutorial.html?td=Tutorial_QHA&tf=QHA

    Parameters
    ----------

    crystal_output: str
        Path to the CRYSTAL output file.

    """

    _is_supercell = False
    _is_restarted = False

    def __init__(self, crystal_output=None):
        """
        Initialize the CRYSTAL QHA reader.

        Parameters
        ----------
        crystal_output : str or pathlib.Path or None, optional
            CRYSTAL QHA output file to load immediately. If ``None``, an empty
            reader is created and :meth:`load` can be called later.
        """
        super().__init__()
        self._is_supercell = False
        self._is_restarted = False
        self._data = self._empty_data()

        if crystal_output is not None:
            self.load(crystal_output)
        return

    @staticmethod
    def _empty_data():
        """
        Return an empty data dictionary for a CRYSTAL QHA calculation.

        Returns
        -------
        dict
            Mutable per-instance storage used by the reader.
        """
        return {
            "points": 0.0,
            "unitcell": [],
            "supercell": [],
            "expansion": np.identity(3, dtype=int),
            "energy": 0.0,
            "kpoints": 1,
            "qpoints": 1,
            "qcoords": {},
            "nphonon": 0,
            "phonons": {},
            "weights": {},
            "shrinkf": np.ones(3, dtype=int),
            "q_position_source": "unavailable-crystal-qha-supercell-modes",
            "structure_series": None,
        }

    def load(self, file):
        """
        Read and store information from a CRYSTAL QHA output file.

        Parameters
        ----------
        file : str or pathlib.Path
            Path to the CRYSTAL QHA output file.
        """
        self.completed = False
        self.error = None
        self._data = self._empty_data()
        geometry = CrystalGeometryParser(file)

        if not self.is_qha(file):
            self.error = (
                "The file is not recognized as a CRYSTAL QHA output."
            )
            return

        self.points = self.set_qha_points(file)
        if self.points < 1:
            self.error = "The CRYSTAL QHA output contains no volume series."
            return

        elif self.points < 4:
            self.error = "Insufficient number of unit cell volumes explored"
            return

        self.restarted_on = self.is_restarted(file)
        if self.restarted_on:
            self.error = "Restarted CRYSTAL QHA calculations are not supported."
            return

        self.supercell_on = self.is_supercell(file)

        if self.supercell_on:
            self.dim = self.set_expansion(file)

        self.qpoints, self.qcoords, self.weights = self.set_dummy_qmesh()
        ordered_volumes = self.set_volume(file)
        cells = self.set_unit_cells(file)
        self._data["unitcell"] = self.set_ordered_cells(cells, ordered_volumes)
        self._data["structure_series"] = self._build_structure_series(
            geometry,
            ordered_volumes,
        )

        self._data["energy"] = self.set_energy(file)

        self.phonons = self.set_phonons(file)

        self.completed = True
        return

    def is_qha(self, file):
        """
        This method checks if the CRYSTAL14/17 output file is related to
        a QHA calculation.

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
                if qha_header in line:
                    return True
            return False
        return

    def is_restarted(self, file):
        """
        This method checks if the CRYSTAL14/17 output file is related to
        a restarted QHA calculation.

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
                if qha_restart in line:
                    return True
            return False
        return

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
                if supercell_option in line:
                    return True
            return False

    @property
    def results(self):
        """
        Get the data collected from the CRYSTAL17 output with `dict` type.
        """
        return self._data

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

    @property
    def restarted_on(self):
        """
        Get the flag that tells if the input file is related to a restarted
        calculation.
        """
        return self._is_restarted

    @restarted_on.setter
    def restarted_on(self, bool_value):
        """
        Set the flag that tells if the input file is related to a restarted
        calculation.
        """
        self._is_restarted = bool_value
        return

    @property
    def points(self):
        """
        Get the number of unit cell volumes explored in QHA analysis.
        """
        return self._data["points"]

    @points.setter
    def points(self, value: int):
        """
        Set the number of unit cell volumes explored in QHA analysis.
        """
        self._data["points"] = value
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
    def natom(self):
        """
        Get the number of atoms in the unit cell.
        """
        return int(self._data["unitcell"][0]["natom"] / self.kpoints)

    @property
    def kpoints(self):
        """
        Get the number of sampled *k*-points, determined from the expansion
        matrix.
        """
        return int(np.around(np.linalg.det(self._data["expansion"]), 0))

    @property
    def energy(self):
        """
        Get the unit cell (if phonon dispersion relations or if
        :math:`\\Gamma`-point frequencies) or the supercell energy.
        """
        return self._data["energy"] / self.kpoints

    @property
    def volume(self):
        """
        Get the unit cell volumes
        """
        volumes = np.zeros(self.points, dtype=float)
        for i in range(self.points):
            volumes[i] = np.linalg.det(self._data["unitcell"][i]["lattice"])
        return volumes / self.kpoints

    @property
    def nphonon(self):
        """
        Get the number of frequencies per band in unit cell (if phonon
        dispersion relations or if :math:`\\Gamma`-point frequencies)
        or in the supercell.
        """
        return self.natom * 3

    @property
    def qpoints(self):
        """
        Get the number of sampled **q**-points. For QHA, this is the same as
        the number of *k*-points.
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
        """Return q-point provenance for native CRYSTAL QHA modes.

        Returns
        -------
        str
            Marker explaining that the output does not label supercell
            eigenmodes by primitive-cell q point.
        """
        return str(self._data["q_position_source"])

    @property
    def qcoords_fractional(self) -> None:
        """Return unavailable primitive q-point labels as ``None``.

        Returns
        -------
        None
            Native CRYSTAL QHA output follows supercell modes with volume but
            does not print a reliable mapping from each mode block to a
            primitive-cell q point.
        """
        return None

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
            Array with shape ``(qpoints, nphonon, points)``.
        """
        return np.asarray(
            [self._data["phonons"][i] for i in range(self.qpoints)],
            dtype=np.float64,
        )

    @property
    def structure_series(self):
        """Return the compact primitive structural path.

        Returns
        -------
        StructureVolumeSeries or None
            Primitive lattices, coordinates, symmetry, and reconstruction
            diagnostics for all QHA volumes.
        """
        return self._data.get("structure_series")

    def _build_structure_series(
        self,
        geometry: CrystalGeometryParser,
        ordered_volumes,
    ):
        """Reduce optimized QHA source cells to a primitive structural path."""
        reference = geometry.initial_primitive_cell()
        expansion = np.asarray(self.dim, dtype=np.int64)
        repetitions = supercell_repetitions(expansion)
        source_cells = [
            self._structure_from_mapping(item, label=f"QHA source cell {index}")
            for index, item in enumerate(self._data["unitcell"])
        ]
        primitive_cells = []
        diagnostics = []
        for source in source_cells:
            primitive, diagnostic = reconstruct_primitive_structure(
                source,
                expansion,
                reference,
            )
            primitive_cells.append(primitive)
            diagnostics.append(diagnostic)
        primitive_volumes = np.asarray(
            [cell.volume for cell in primitive_cells],
            dtype=np.float64,
        )
        normalized_volumes = np.asarray(ordered_volumes, dtype=np.float64) / repetitions
        if not np.allclose(
            primitive_volumes,
            normalized_volumes,
            rtol=0.0,
            atol=2.0e-5,
        ):
            raise ValueError(
                "reconstructed primitive volumes do not match the CRYSTAL QHA "
                "normalization"
            )
        reference_index = int(np.argmin(np.abs(primitive_volumes - reference.volume)))
        compact_reference = primitive_cells[reference_index]
        symmetry = analyze_symmetry(compact_reference, symprec=1.0e-5)
        normalization = CellNormalization(
            basis="primitive",
            source_basis=("qha_supercell" if repetitions > 1 else "primitive"),
            expansion_matrix=expansion,
            repetitions=repetitions,
            source_atoms=source_cells[0].natoms,
            normalized_atoms=compact_reference.natoms,
        )
        return StructureVolumeSeries(
            reference=compact_reference,
            lattices=np.asarray(
                [cell.lattice for cell in primitive_cells],
                dtype=np.float64,
            ),
            fractional_positions=np.asarray(
                [cell.fractional_positions for cell in primitive_cells],
                dtype=np.float64,
            ),
            volumes=primitive_volumes,
            normalization=normalization,
            symmetry=symmetry,
            primitive_to_crystallographic=(geometry.primitive_to_crystallographic()),
            diagnostics=tuple(diagnostics),
            orientation="crystal",
            reference_index=reference_index,
            metadata={
                "interface": "crystal-qha",
                "coorprt_present": geometry.has_coorprt,
                "space_group_number_from_output": geometry.space_group_number(),
            },
            source_lattices=np.asarray(
                [cell.lattice for cell in source_cells],
                dtype=np.float64,
            ),
            source_fractional_positions=tuple(
                cell.fractional_positions.copy() for cell in source_cells
            ),
        )

    @staticmethod
    def _structure_from_mapping(mapping, *, label):
        """Convert the historical QHA cell mapping to a structure object."""
        from quantas.models.structures import CrystalStructure

        return CrystalStructure(
            lattice=np.asarray(mapping["lattice"], dtype=np.float64),
            fractional_positions=np.asarray(mapping["positions"], dtype=np.float64),
            atomic_numbers=np.asarray(mapping["numbers"], dtype=np.int64),
            label=label,
        )

    def set_qha_points(self, file):
        """
        This method sets the number of points (number of unit cell volumes)
        explored during the QHA analysis.

        They are taken from the number of lines in the EOS output section.

        Parameters
        ----------

        file: str
            Path of the CRYSTAL14/17 output file.

        Returns
        -------

        points: int
            Number of unit cell volumes considered.

        """
        sline = self._get_start_line(file, qha_ev) + 4

        with open(file, "r") as f:
            data = f.readlines()

        points = 0
        for i in range(sline, sline + 100):
            if len(data[i].split()) != 2:
                break
            else:
                points += 1
        return points

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
        sline = self._get_start_line(file, supercell_expansion) + 1

        with open(file, "r") as f:
            data = f.readlines()

        expansion = np.zeros((3, 3), dtype=float)

        for i in range(3):
            line = data[sline + i].split()
            del line[0]
            expansion[i] = np.asarray(line, dtype=float)
        return expansion

    def set_energy(self, file):
        """
        This method sets the energy values for each unit cell.

        They are taken from the number of lines in the EOS output section.

        Parameters
        ----------

        file: str
            Path of the CRYSTAL14/17 output file.

        Returns
        -------

        energy: ndarray
            Unit cell energy values with `float` type.

        """
        sline = self._get_start_line(file, qha_ev) + 4

        with open(file, "r") as f:
            data = f.readlines()

        energy = np.zeros(self.points, dtype=float)
        for i in range(self.points):
            energy[i] = float(data[sline + i].split()[1])

        return energy

    def set_volume(self, file):
        """
        This method sets the unit cell volume values, retrieved from
        the number of lines in the EOS output section. They are used to
        re-order the unit cell data.

        Parameters
        ----------

        file: str
            Path of the CRYSTAL14/17 output file.

        Returns
        -------

        volume: ndarray
            Unit cell volumes with `float` type.

        """
        sline = self._get_start_line(file, qha_ev) + 4

        with open(file, "r") as f:
            data = f.readlines()

        volume = np.zeros(self.points, dtype=float)
        for i in range(self.points):
            volume[i] = float(data[sline + i].split()[0])

        return volume

    def set_unit_cells(self, file):
        """Return optimized source cells printed by the CRYSTAL QHA run.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL QHA output file.

        Returns
        -------
        list of dict
            Historical cell mappings used by the QHA reader.
        """
        cells = []
        for structure in CrystalGeometryParser(file).optimized_cells():
            cells.append(
                {
                    "natom": structure.natoms,
                    "numbers": structure.atomic_numbers.copy(),
                    "positions": structure.fractional_positions.copy(),
                    "lattice": structure.lattice.copy(),
                }
            )
        return cells

    def set_ordered_cells(self, cells, ordered_volumes):
        """
        This method reorders the collected unit cell data according to the
        sorted volumes reported in the CRYSTAL output.

        Parameters
        ----------

        cells: list
            List containing the unit cell data, each element with `dict`
            type.

        ordered_volumes: ndarray
            Array containing the sorted unit cell volumes.

        Returns
        -------

        ordered_cells: list
            List of ordered unit cell data by increasing unit cell volume.

        """
        indexes = []
        for i in range(len(cells)):
            volume = np.linalg.det(cells[i]["lattice"])
            for j in range(len(cells)):
                if np.isclose(volume, ordered_volumes[j]):
                    indexes.append(j)

        ordered_cells = []
        for i in range(len(cells)):
            ordered_cells.append(cells[indexes.index(i)])

        return ordered_cells

    def set_optimized_cell(self, file, idx: int):
        """Return an optimized geometry beginning at a known marker index.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL QHA output file.
        idx : int
            Line index of ``FINAL OPTIMIZED GEOMETRY``.

        Returns
        -------
        tuple
            Atom count, atomic numbers, fractional positions, and lattice.
        """
        parser = CrystalGeometryParser(file)
        structure = parser._parse_geometry_after_marker(  # noqa: SLF001
            idx,
            label="optimized cell",
        )
        return (
            structure.natoms,
            structure.atomic_numbers.copy(),
            structure.fractional_positions.copy(),
            structure.lattice.copy(),
        )

    def set_phonons(self, file):
        """ """
        phonons = {}
        nfreq = self.natom * self.qpoints * 3
        phonon_matrix = np.zeros((nfreq, self.points), dtype=float)

        with open(file, "r") as f:
            data = f.readlines()

        freq_idx = []
        for i in range(len(data)):
            if qha_freq in data[i]:
                if data[i].split()[0] == "FREQUENCY":
                    freq_idx.append(i + 3)

        for i in range(nfreq):
            for j in range(self.points):
                frequency_line = data[freq_idx[i] + j].split()
                phonon_matrix[i, j] = float(frequency_line[1])

        bands = phonon_matrix.reshape(self.qpoints, self.natom * 3, self.points)

        for i in range(self.qpoints):
            phonons[i] = bands[i].copy()

        return phonons

    def set_dummy_qmesh(self):
        """Create equal-weight storage blocks for native CRYSTAL QHA modes.

        CRYSTAL follows the supercell Gamma eigenmodes with volume, but the
        final QHA table does not label those modes by primitive-cell q point.
        Quantas retains ``det(expansion)`` equal-weight blocks to preserve the
        historical frequency-array shape and thermodynamic normalization. The
        placeholder zero coordinates are never exported as physical q-point
        labels; the generated YAML marks their positions as unavailable.

        Returns
        -------
        tuple
            Number of equal-weight blocks, placeholder coordinates, and unit
            weights.
        """
        qpoints = self.kpoints
        qcoords = {}
        qweights = {}
        for i in range(qpoints):
            qcoords[i] = np.zeros(3, dtype=float)
            qweights[i] = 1.0
        return qpoints, qcoords, qweights

    def _get_start_line(self, file, search_string):
        sline = 0
        with open(file, "r") as f:
            for line in f:
                if search_string in line:
                    break
                else:
                    sline += 1
        return sline
