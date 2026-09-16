# -*- coding: utf-8 -*-

"""Read CRYSTAL quasi-harmonic outputs and structural volume paths."""

from __future__ import annotations

import numpy as np

from quantas.core.geometry import (
    analyze_symmetry,
    reconstruct_primitive_structure,
    supercell_repetitions,
)
from quantas.interfaces.crystal import markers
from quantas.interfaces.crystal.geometry import CrystalGeometryParser
from quantas.models.reader import BasicReader
from quantas.models.structures import (
    CellNormalization,
    StructureVolumeSeries,
)


class CrystalQHAReader(BasicReader):
    """Read a native CRYSTAL quasi-harmonic calculation.

    The reader collects the volume-energy series, optimized structures, and
    volume-resolved phonon frequencies printed by CRYSTAL's native QHA workflow.
    Source supercells are reduced to a compact primitive structural path for
    Quantas, while the original expansion matrix and reconstruction provenance are
    retained.

    CRYSTAL follows supercell Gamma eigenmodes with volume but does not provide a
    reliable primitive-cell q-point label for every stored mode block. The reader
    therefore preserves equal-weight historical storage blocks while exposing
    ``qcoords_fractional=None`` and explicit q-position provenance rather than
    inventing physical q-point coordinates.

    Parameters
    ----------
    crystal_output : str or pathlib.Path or None, optional
        CRYSTAL QHA output to load immediately. If ``None``, create an empty
        reader and call :meth:`load` later."""

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
            "mode_continuity": "unknown",
            "mode_continuity_metadata": {},
        }

    def load(self, file):
        """Read one native CRYSTAL QHA output into the reader state.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL QHA output file.

        Notes
        -----
        The reader rejects restarted native-QHA outputs and requires at least four
        volume states. Recognition and supported-workflow failures are stored in
        ``error`` with ``completed=False``; structural inconsistencies that would make
        primitive reconstruction ambiguous raise explicitly."""
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
        if self.has_verified_mode_continuity(file):
            self._data["mode_continuity"] = "verified"
            self._data["mode_continuity_metadata"] = {
                "method": "crystal-qha",
                "source": "crystal",
            }

        self.completed = True
        return

    def is_qha(self, file):
        """Return whether an output contains a native CRYSTAL QHA calculation.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL output file.

        Returns
        -------
        bool
            ``True`` when the native QHA header is present."""
        with open(file, "r") as f:
            for line in f:
                if markers.QHA_HEADER in line:
                    return True
            return False
        return

    def has_verified_mode_continuity(self, file) -> bool:
        """Return whether CRYSTAL reports completed QHA mode continuity.

        Parameters
        ----------
        file : str or pathlib.Path
            Native CRYSTAL QHA output.

        Returns
        -------
        bool
            ``True`` when the final CRYSTAL continuity table is present.
        """
        with open(file, "r") as stream:
            return any(markers.QHA_CONTINUITY_FOUND in line for line in stream)

    @property
    def mode_continuity(self) -> str:
        """Return the continuity status established by native CRYSTAL QHA."""
        return str(self._data.get("mode_continuity", "unknown"))

    @property
    def mode_continuity_metadata(self) -> dict[str, object]:
        """Return provenance for native CRYSTAL QHA mode continuity."""
        return dict(self._data.get("mode_continuity_metadata", {}))

    def is_restarted(self, file):
        """Return whether a native CRYSTAL QHA calculation was restarted.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL QHA output file.

        Returns
        -------
        bool
            ``True`` when the QHA restart marker is present."""
        with open(file, "r") as f:
            for line in f:
                if markers.QHA_RESTART in line:
                    return True
            return False
        return

    def is_supercell(self, file):
        """Return whether the QHA phonons use an explicit CRYSTAL supercell.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL QHA output file.

        Returns
        -------
        bool
            ``True`` when the ``SUPERCELL`` option is present."""
        with open(file, "r") as f:
            for line in f:
                if markers.SUPERCELL_OPTION in line:
                    return True
            return False

    @property
    def results(self):
        """Return the internal parsed CRYSTAL QHA payload.

        Notes
        -----
        This property exposes the historical reader mapping for compatibility. New
        workflow code should prefer the typed/public properties such as :attr:`energy`,
        :attr:`volume`, :attr:`phonons_array`, and :attr:`structure_series`."""
        return self._data

    @property
    def supercell_on(self):
        """
        Return the flag that tells if the input file is related to a supercell.
        """
        return self._is_supercell

    @supercell_on.setter
    def supercell_on(self, bool_value):
        """Store whether the native QHA run uses a phonon supercell.

        Parameters
        ----------
        bool_value : bool
            Supercell-state flag."""
        self._is_supercell = bool_value
        return

    @property
    def restarted_on(self):
        """
        Return the flag that tells if the input file is related to a restarted
        calculation.
        """
        return self._is_restarted

    @restarted_on.setter
    def restarted_on(self, bool_value):
        """Store whether the native QHA output is a restarted run.

        Parameters
        ----------
        bool_value : bool
            Restart-state flag."""
        self._is_restarted = bool_value
        return

    @property
    def points(self):
        """
        Return the number of unit cell volumes explored in QHA analysis.
        """
        return self._data["points"]

    @points.setter
    def points(self, value: int):
        """Store the number of QHA volume states.

        Parameters
        ----------
        value : int
            Number of sampled volume states."""
        self._data["points"] = value
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
    def natom(self):
        """Return the primitive-cell atom count represented by the QHA modes.

        The source optimized cells may be supercells; their atom count is divided by
        the expansion determinant to obtain the primitive normalization."""
        return int(self._data["unitcell"][0]["natom"] / self.kpoints)

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
        """Return the number of primitive-cell repetitions in the QHA supercell.

        The historical property name is retained for compatibility; numerically this
        is the rounded determinant of the expansion matrix, not an electronic
        Brillouin-zone k-point count."""
        return int(np.around(np.linalg.det(self._data["expansion"]), 0))

    @property
    def energy(self):
        """Return primitive-normalized QHA static energies in hartree.

        The native QHA table stores source-cell energies. Quantas divides the full
        ``(points,)`` series by the expansion determinant so energies and primitive
        volumes share the same normalization."""
        return self._data["energy"] / self.kpoints

    @property
    def volume(self):
        """Return primitive-normalized QHA volumes in ``angstrom^3``.

        Returns
        -------
        numpy.ndarray
            Array with shape ``(points,)`` ordered consistently with the native QHA
            energy and phonon series."""
        volumes = np.zeros(self.points, dtype=float)
        for i in range(self.points):
            volumes[i] = np.linalg.det(self._data["unitcell"][i]["lattice"])
        return volumes / self.kpoints

    @property
    def nphonon(self):
        """Return the number of primitive phonon branches per storage block.

        The value is ``3 * natom`` after primitive-cell normalization."""
        return self.natom * 3

    @property
    def qpoints(self):
        """Return the number of equal-weight native-QHA storage blocks.

        For a supercell calculation this equals the expansion determinant for
        historical thermodynamic normalization. These blocks are not reliable
        primitive-cell q-point labels; :attr:`qcoords_fractional` therefore returns
        ``None``."""
        return self._data["qpoints"]

    @qpoints.setter
    def qpoints(self, value: int):
        """Store the number of native-QHA frequency blocks.

        Parameters
        ----------
        value : int
            Number of equal-weight storage blocks."""
        self._data["qpoints"] = value
        return

    @property
    def qcoords(self):
        """Return placeholder coordinates for native-QHA storage blocks.

        The zero vectors preserve the historical array layout only and must not be
        interpreted as physical primitive-cell q-point labels. See
        :attr:`q_position_source`."""
        return self._data["qcoords"]

    @qcoords.setter
    def qcoords(self, array):
        """Store placeholder coordinates for native-QHA blocks.

        Parameters
        ----------
        array : array-like
            Coordinate rows ordered by storage-block index. These placeholders are not
            physical primitive q-point labels."""
        for i in range(len(array)):
            self._data["qcoords"][i] = array[i]
        return

    @property
    def weights(self):
        """Return equal integration weights for native-QHA storage blocks."""
        return self._data["weights"]

    @weights.setter
    def weights(self, array):
        """Store integration weights for native-QHA blocks.

        Parameters
        ----------
        array : array-like
            Weight values ordered by storage-block index."""
        for i in range(len(array)):
            self._data["weights"][i] = array[i]
        return

    @property
    def shrinkf(self):
        """Return the historical reciprocal-space shrinking-factor storage."""
        return self._data["shrinkf"]

    @shrinkf.setter
    def shrinkf(self, array):
        """Store the historical reciprocal shrinking factors.

        Parameters
        ----------
        array : array-like
            Three shrinking factors. A copy is retained by the reader."""
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
        """Return volume-resolved frequencies keyed by storage-block index.

        Each value has shape ``(nphonon, points)`` and contains frequencies in
        ``cm^-1`` ordered along the QHA volume path."""
        return self._data["phonons"]

    @phonons.setter
    def phonons(self, dictionary):
        """Store volume-resolved phonon-frequency blocks.

        Parameters
        ----------
        dictionary : dict
            Mapping from storage-block index to arrays with shape
            ``(nphonon, points)`` in ``cm^-1``."""
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
        """Return the number of volume states in the native QHA table.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL QHA output file.

        Returns
        -------
        int
            Number of contiguous volume-energy rows printed by CRYSTAL."""
        sline = self._get_start_line(file, markers.QHA_VOLUME_ENERGY_TABLE) + 4

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
        """Return the CRYSTAL primitive-to-supercell expansion matrix.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL QHA output file.

        Returns
        -------
        numpy.ndarray
            ``(3, 3)`` expansion matrix. Its determinant defines the number of
            primitive-cell repetitions used to normalize source energies and volumes."""
        sline = self._get_start_line(file, markers.SUPERCELL_EXPANSION) + 1

        with open(file, "r") as f:
            data = f.readlines()

        expansion = np.zeros((3, 3), dtype=float)

        for i in range(3):
            line = data[sline + i].split()
            del line[0]
            expansion[i] = np.asarray(line, dtype=float)
        return expansion

    def set_energy(self, file):
        """Return source-cell energies from the native QHA table.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL QHA output file.

        Returns
        -------
        numpy.ndarray
            Array with shape ``(points,)`` containing energies in hartree as printed
            by CRYSTAL. The public :attr:`energy` property divides these values by the
            supercell repetition count to expose primitive-normalized energies."""
        sline = self._get_start_line(file, markers.QHA_VOLUME_ENERGY_TABLE) + 4

        with open(file, "r") as f:
            data = f.readlines()

        energy = np.zeros(self.points, dtype=float)
        for i in range(self.points):
            energy[i] = float(data[sline + i].split()[1])

        return energy

    def set_volume(self, file):
        """Return source-cell volumes from the native QHA table.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL QHA output file.

        Returns
        -------
        numpy.ndarray
            Array with shape ``(points,)`` containing source volumes in
            ``angstrom^3``. The public :attr:`volume` property exposes the corresponding
            primitive-normalized volumes."""
        sline = self._get_start_line(file, markers.QHA_VOLUME_ENERGY_TABLE) + 4

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
        """Order optimized CRYSTAL cells by the native QHA volume sequence.

        Parameters
        ----------
        cells : list of dict
            Parsed optimized source-cell mappings.
        ordered_volumes : array-like
            Source-cell volumes in the order printed by the native QHA summary.

        Returns
        -------
        list of dict
            Cell mappings reordered to match ``ordered_volumes``."""
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
        """Parse volume-resolved phonon frequencies from CRYSTAL QHA output.

        Parameters
        ----------
        file : str or pathlib.Path
            CRYSTAL output produced by a native QHA calculation.

        Returns
        -------
        dict[int, numpy.ndarray]
            Mapping from equal-weight phonon block index to an array with shape
            ``(self.natom * 3, self.points)`` containing frequencies in cm^-1
            along the sampled volume path.

        Raises
        ------
        OSError
            If the output file cannot be opened.
        ValueError
            If a parsed frequency value is not numeric.
        """
        phonons = {}
        nfreq = self.natom * self.qpoints * 3
        phonon_matrix = np.zeros((nfreq, self.points), dtype=float)

        with open(file, "r") as f:
            data = f.readlines()

        freq_idx = []
        for i in range(len(data)):
            if markers.QHA_FREQUENCY in data[i]:
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
