# -*- coding: utf-8 -*-
"""Read Phonopy YAML data and optional VASP energies for Quantas."""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from scipy import constants as cs

from quantas.core.chemistry.symbols import symbol2number as s2n
from quantas.core.physics.units import convert_frequency as cf
from quantas.models.reader import BasicReader

h = cs.Planck
eV = cs.value("electron volt")


class PhonopyReader(BasicReader[None]):
    """Read a Phonopy mesh and its associated displacement structure.

    The main YAML file must be accompanied by a displacement YAML file with
    the same stem and the suffix ``_disp.yaml``. When a VASP XML file with the
    same stem is available, its electronic energy is included in the parsed
    data. Other Phonopy backends can still provide phonon data, but their energy
    must be supplied separately before an HA or QHA calculation.

    Parameters
    ----------
    phonopy_output : str or Path
        Phonopy YAML file to read.
    """

    def __init__(self, phonopy_output: str | Path) -> None:
        """Initialize the reader and load the requested Phonopy dataset."""
        super().__init__()
        self._data: dict[str, Any] = {}
        self.load(phonopy_output)

    def load(self, file: str | Path) -> None:
        """Read one Phonopy YAML file and its companion data.

        Parameters
        ----------
        file : str or Path
            Phonopy YAML file.
        """
        self._data = {
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
        }

        file = str(file)
        dispfile = os.path.splitext(file)[0] + "_disp.yaml"
        if not os.path.exists(dispfile):
            self.error = (
                f"Required Phonopy displacement file not found: {dispfile}"
            )
            return

        vasprun = os.path.splitext(file)[0] + ".xml"
        if os.path.exists(vasprun):
            root = ET.parse(vasprun).getroot()
            energy_results = root.findall("calculation/energy/i")
            for energy_elem in energy_results:
                if energy_elem.attrib["name"] == "e_wo_entrp":
                    if energy_elem.text is not None:
                        self._data["energy"] = float(energy_elem.text)

        if os.path.exists(file):
            yaml_error, exception, idata = self.load_yaml(file)
            if yaml_error or idata is None:
                self.error = exception or "Unable to parse phonopy YAML input"
                return
        else:
            self.error = " File not found!"
            return

        self.unitcell = self.set_unit_cell(idata)

        error, exception, dispdata = self.load_yaml(dispfile)
        if error or dispdata is None:
            self.error = exception or "Unable to parse phonopy displacement YAML"
            return
        self.supercell = self.set_displacements(dispdata)

        self.set_phonons(idata)
        self.completed = True

    @property
    def unitcell(self) -> tuple[Any, Any, Any, Any]:
        """Return unit-cell atom count, lattice, positions, and atomic numbers."""
        return (
            self._data["unitcell"]["natom"],
            self._data["unitcell"]["lattice"],
            self._data["unitcell"]["positions"],
            self._data["unitcell"]["numbers"],
        )

    @unitcell.setter
    def unitcell(self, cell_data: tuple[Any, Any, Any, Any]) -> None:
        """Store unit-cell data in the reader's historical tuple layout."""
        self._data["unitcell"]["natom"] = cell_data[0]
        self._data["unitcell"]["numbers"] = cell_data[1]
        self._data["unitcell"]["positions"] = cell_data[2]
        self._data["unitcell"]["lattice"] = cell_data[3]

    @property
    def supercell(self) -> tuple[Any, Any, Any]:
        """Return supercell lattice, positions, and atomic numbers."""
        return (
            self._data["supercell"]["lattice"],
            self._data["supercell"]["positions"],
            self._data["supercell"]["numbers"],
        )

    @supercell.setter
    def supercell(self, cell_data: tuple[Any, Any, Any, Any]) -> None:
        """Store supercell data in the reader's historical tuple layout."""
        self._data["supercell"]["natom"] = cell_data[0]
        self._data["supercell"]["numbers"] = cell_data[1]
        self._data["supercell"]["positions"] = cell_data[2]
        self._data["supercell"]["lattice"] = cell_data[3]

    @property
    def dim(self) -> Any:
        """Return the Phonopy supercell expansion matrix."""
        return self._data["expansion"]

    @dim.setter
    def dim(self, matrix: Any) -> None:
        """Store a copy of the Phonopy supercell expansion matrix."""
        self._data["expansion"] = matrix.copy()

    @property
    def natom(self) -> int:
        """Return the unit-cell atom count."""
        return self._data["unitcell"]["natom"]

    @property
    def lattice(self) -> Any:
        """Return the unit-cell lattice matrix."""
        return self._data["unitcell"]["lattice"]

    @property
    def volume(self) -> float:
        """Return the signed determinant of the unit-cell lattice."""
        return np.linalg.det(self._data["unitcell"]["lattice"])

    @property
    def energy(self) -> float:
        """Return the electronic energy read from the companion VASP XML file."""
        return self._data["energy"]

    @property
    def units(self) -> dict[str, str]:
        """Return the physical units exposed by this Phonopy reader.

        The current interface reads static energies from a companion VASP XML
        file and converts Phonopy frequencies from THz to wavenumbers.

        Returns
        -------
        dict
            Energy, volume, frequency, and structural length units.
        """
        return {
            "energy": "eV",
            "volume": "angstrom^3",
            "frequency": "cm^-1",
            "length": "angstrom",
        }

    @property
    def kpoints(self) -> int:
        """Return the number of primitive-cell repetitions in the supercell."""
        return int(np.around(np.linalg.det(self._data["expansion"]), 0))

    @property
    def qpoints(self) -> int:
        """Return the number of sampled phonon q-points."""
        return self._data["qpoints"]

    @property
    def nphonon(self) -> int:
        """Return the number of phonon branches per q-point."""
        return self._data["nphonon"]

    @property
    def phonons(self) -> Any:
        """Return phonon frequencies indexed by q-point."""
        return self._data["phonons"]

    @property
    def qcoords(self) -> Any:
        """Return fractional q-point coordinates."""
        return self._data["qcoords"]

    @property
    def weights(self) -> Any:
        """Return q-point integration weights."""
        return self._data["weights"]

    @property
    def shrinkf(self) -> Any:
        """Return the historical reciprocal-space shrinking factors."""
        return self._data["shrinkf"]

    def load_yaml(
        self,
        file: str | Path,
    ) -> tuple[bool, str | None, dict[str, Any] | None]:
        """Read a YAML mapping without raising parser errors.

        Parameters
        ----------
        file : str or Path
            YAML file to read.

        Returns
        -------
        tuple
            Error flag, optional error description, and parsed mapping.
        """
        idata = None
        error = False
        exception = None
        with open(file, "r", encoding="utf-8") as stream:
            try:
                idata = yaml.safe_load(stream)
            except yaml.YAMLError:
                error = True
        return error, exception, idata

    def set_displacements(
        self,
        idata: dict[str, Any],
    ) -> tuple[int, Any, Any, Any]:
        """Extract supercell geometry from a Phonopy displacement mapping."""
        self.dim = np.asarray(idata["supercell_matrix"])

        natom = len(idata["supercell"]["points"])
        numbers = np.zeros(natom, dtype=int)
        positions = np.zeros((natom, 3), dtype=int)
        lattice = idata["supercell"]["lattice"]

        for index in range(natom):
            symbol = idata["supercell"]["points"][index]["symbol"]
            numbers[index] = s2n(symbol)
            positions[index] = idata["supercell"]["points"][index]["coordinates"]

        return natom, numbers, positions, lattice

    def set_unit_cell(self, idata: dict[str, Any]) -> tuple[int, Any, Any, Any]:
        """Extract unit-cell geometry from a Phonopy mesh mapping."""
        natom = idata["natom"]
        numbers = np.zeros(natom, dtype=int)
        positions = np.zeros((natom, 3), dtype=int)
        lattice = idata["lattice"]

        for index in range(natom):
            symbol = idata["points"][index]["symbol"]
            numbers[index] = s2n(symbol)
            positions[index] = idata["points"][index]["coordinates"]

        return natom, numbers, positions, lattice

    def set_phonons(self, idata: dict[str, Any]) -> None:
        """Store q-points, weights, and frequencies from a Phonopy mapping."""
        nf = 3 * idata["natom"]
        nq = idata["nqpoint"]

        self._data["nphonon"] = nf
        self._data["qpoints"] = nq

        for index in range(nq):
            qdata = idata["phonon"][index]
            self._data["qcoords"][index] = qdata["q-position"]
            self._data["weights"][index] = qdata["weight"]
            band = np.zeros(nf, dtype=float)
            for branch in range(nf):
                band[branch] = cf(
                    qdata["band"][branch]["frequency"],
                    "THz",
                    "cm^-1",
                )
            self._data["phonons"][index] = band
