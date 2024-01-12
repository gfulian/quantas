# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

import numpy as np

from quantas.core.reader import BasicReader

# Geometry strings
geometry_consistent = 'GEOMETRY NOW FULLY CONSISTENT WITH THE GROUP'
geometry_sym_changed = 'SYMMETRY CHANGED DURING OLD RUN :'
geometry_wf = 'GEOMETRY FOR WAVE FUNCTION - '
primitive_cell = 'PRIMITIVE CELL - '
supercell_option = ' * SUPERCELL OPTION'
supercell_expansion = 'EXPANSION MATRIX OF PRIMITIVE CELL'
lattice_vectors = 'DIRECT LATTICE VECTORS CARTESIAN COMPONENTS'
# FREQCALC-specific strings
frequency_calculation = 'EIGENVALUES (EIGV) OF THE MASS WEIGHTED HESSIAN'
frequency_calculation += ' MATRIX AND HARMONIC'
scelphono_option = 'PHONON FREQUENCIES AT A SET OF K POINTS BY USING '
scelphono_option += 'A SUPERCELL'
hess_interpolation = 'ACTIVATED INTERPOLATION OF THE HESSIAN UP TO'
scelphono_qpoints = 'THAT PERMITS THE CALCULATION OF MODES AT'
central_energy = '    CENTRAL POINT'
frequency_header = 'MODES         EIGV          FREQUENCIES     IRREP'
# QHA-specific strings
qha_header = 'QUASI-HARMONIC APPROXIMATION'
qha_freq = 'FREQUENCY #'
qha_ev = 'SORTING VOLUMES/ENERGIES'
qha_opt = 'FINAL OPTIMIZED GEOMETRY'
qha_restart = 'READING DATA FROM RESTART UNIT'
# ELASTCON-specific strings
soec_option = 'ELASTCON OPTION'
soec_completed = 'FINAL RESULTS START'
soec_ec = 'SYMMETRIZED ELASTIC CONSTANTS'


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

    def __init__(self, crystal_output):
        BasicReader.__init__(self, crystal_output)
        return

    def load(self, file):
        """
        Read and store the information collected in the CRYSTAL output file.

        Parameters
        ----------

        file: str
            Path to the CRYSTAL output file.

        """
        self._data = {
            'unitcell': {},
            'supercell': {},
            'expansion': np.identity(3, dtype=int),
            'energy': 0.,
            'kpoints': 1,
            'qpoints': 1,
            'qcoords': {},
            'nphonon': 0,
            'phonons': {},
            'weights': {},
            'shrinkf': np.ones(3, dtype=int)
            }

        if not self.is_frequency_calculation(file):
            self.error = 'This does not appear to be a CRYSTAL14/17 '
            self.error += 'output \nrelated to phonon calculations'
            return

        self.supercell_on = self.is_supercell(file)
        self.scelphono_on = self.is_phonon_dispersion(file)
        #
        # Collect system information
        if self.supercell_on:
            self.supercell = self.set_wf_cell(file)
            self.unitcell = self.set_init_cell(file)
            self.dim = self.set_expansion(file)
            self._data['unitcell']['lattice'] = np.dot(
                self._data['supercell']['lattice'], np.linalg.inv(self.dim))
            self.qpoints, self.qcoords, self.weights, self.shrinkf = \
                          self.set_q_mesh(file)

        else:
            self.unitcell = self.set_wf_cell(file)
            self.supercell = self.set_wf_cell(file)
            self.qcoords = [[0., 0., 0.]]
            self.weights = [1]
        # Collect energy values and phonons
        self._data['energy'] = self.set_energy(file)
        self.phonons = self.set_phonons(file)
        #
        self._check(file)
        return

    def _check(self, file):
        """
        """
        if self.energy == 0.:
            self.error = 'No unit cell energy in {0}'.format(file)
            return
        if not np.any(self.lattice):
            self.error = 'No unit cell lattice in {0}'.format(file)
            return
        if not bool(self.phonons):
            self.error = 'No phonon data in {0}'.format(file)
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
        with open(file, 'r') as f:
            for line in f:
                if frequency_calculation in line:
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
        with open(file, 'r') as f:
            for line in f:
                if supercell_option in line:
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
        with open(file, 'r') as f:
            for line in f:
                if scelphono_option in line:
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
        with open(file, 'r') as f:
            for line in f:
                if hess_interpolation in line:
                    return True
            return False

    @property
    def natom(self):
        """
        Get the number of atoms in the unit cell (if phonon dispersion
        relations or if :math:`\Gamma`-point frequencies) or in the
        supercell.
        """
        if self.supercell_on:
            if self.scelphono_on:
                return self._data['unitcell']['natom']
            else:
                return self._data['supercell']['natom']
        else:
            return self._data['unitcell']['natom']

    @property
    def lattice(self):
        """
        Get the unit cell (if phonon dispersion relations or if
        :math:`\Gamma`-point frequencies) or the supercell lattice vectors.
        """
        if self.supercell_on:
            if self.scelphono_on:
                return self._data['unitcell']['lattice']
            else:
                return self._data['supercell']['lattice']
        else:
            return self._data['unitcell']['lattice']

    @property
    def volume(self):
        """
        Get the unit cell (if phonon dispersion relations or if
        :math:`\Gamma`-point frequencies) or the supercell volume.
        """
        if self.supercell_on:
            if self.scelphono_on:
                return np.linalg.det(self._data['unitcell']['lattice'])
            else:
                return np.linalg.det(self._data['supercell']['lattice'])
        else:
            return np.linalg.det(self._data['unitcell']['lattice'])

    @property
    def energy(self):
        """
        Get the unit cell (if phonon dispersion relations or if
        :math:`\Gamma`-point frequencies) or the supercell energy.
        """
        if self.supercell_on:
            if self.scelphono_on:
                return self._data['energy'] / self.kpoints
            else:
                return self._data['energy']
        else:
            return self._data['energy'] / self.kpoints

    @property
    def nphonon(self):
        """
        Get the number of frequencies per band in unit cell (if phonon
        dispersion relations or if :math:`\Gamma`-point frequencies)
        or in the supercell.
        """
        if self.supercell_on:
            if self.scelphono_on:
                return self._data['unitcell']['natom'] * 3
            else:
                return self._data['supercell']['natom'] * 3
        else:
            return self._data['unitcell']['natom'] * 3

    @property
    def unitcell(self):
        """
        Get the crystal unit cell in tuple format.
        """
        return (
            self._data['unitcell']['natom'],
            self._data['unitcell']['numbers'],
            self._data['unitcell']['positions'],
            self._data['unitcell']['lattice']
            )

    @unitcell.setter
    def unitcell(self,  cell_data):
        """
        Set the crystal unit cell in tuple format.
        """
        self._data['unitcell']['natom'] = cell_data[0]
        self._data['unitcell']['numbers'] = cell_data[1]
        self._data['unitcell']['positions'] = cell_data[2]
        self._data['unitcell']['lattice'] = cell_data[3]
        return

    @property
    def supercell(self):
        """
        Get the crystal supercell in tuple format.
        """
        return (
            self._data['supercell']['natom'],
            self._data['supercell']['numbers'],
            self._data['supercell']['positions'],
            self._data['supercell']['lattice']
            )

    @supercell.setter
    def supercell(self, cell_data):
        """
        Set the crystal supercell in tuple format.
        """
        self._data['supercell']['natom'] = cell_data[0]
        self._data['supercell']['numbers'] = cell_data[1]
        self._data['supercell']['positions'] = cell_data[2]
        self._data['supercell']['lattice'] = cell_data[3]
        return

    @property
    def dim(self):
        """
        Get the expansion matrix employed to build the supercell.
        """
        return self._data['expansion']

    @dim.setter
    def dim(self, expansion):
        """
        Set the expansion matrix employed to build the supercell.
        """
        self._data['expansion'] = expansion.copy()
        return

    @property
    def kpoints(self):
        """
        Get the number of sampled *k*-points, determined from the expansion
        matrix.
        """
        return int(np.around(np.linalg.det(self._data['expansion']), 0))

    @property
    def qpoints(self):
        """
        Get the number of sampled **q**-points.
        """
        return self._data['qpoints']

    @qpoints.setter
    def qpoints(self, value: int):
        """
        Set the number of sampled **q**-points.
        """
        self._data['qpoints'] = value
        return

    @property
    def qcoords(self):
        """
        Get the coordinates of sampled **q**-points, in dict format.
        """
        return self._data['qcoords']

    @qcoords.setter
    def qcoords(self, array):
        """
        Set the coordinates of sampled **q**-points, in dict format.
        """
        for i in range(len(array)):
            self._data['qcoords'][i] = array[i]
        return

    @property
    def weights(self):
        """
        Get the weights of each phonon band, in dict format.
        """
        return self._data['weights']

    @weights.setter
    def weights(self, array):
        """
        Set the weights of each phonon band.
        """
        for i in range(len(array)):
            self._data['weights'][i] = array[i]
        return

    @property
    def shrinkf(self):
        """
        Get the Hessian interpolation mesh used in the INTERPHESS keyword.
        """
        return self._data['shrinkf']

    @shrinkf.setter
    def shrinkf(self, array):
        """
        Set the Hessian interpolation mesh used in the INTERPHESS keyword.
        """
        self._data['shrinkf'] = array.copy()
        return

    @property
    def phonons(self):
        """
        Get the phonon bands, in dict format.
        """
        return self._data['phonons']

    @phonons.setter
    def phonons(self, dictionary):
        """
        Set the phonon bands, in dict format.
        """
        self._data['phonons'] = dictionary
        return

    def set_init_cell(self, file):
        """
        This method sets the geometry of the crystal to that provided as
        input.

        Parameters
        ----------

        file: str
            Path of the CRYSTAL14/17 output file.

        Returns
        -------

        natom: int
            Number of atoms in the crystal cell

        numbers: ndarray
            Array containing atomic numbers with `int` type.

        positions: ndarray
            2D array containing atomic fractional coordinates with `float`
            type.

        """
        sline = self._get_start_line(file, geometry_consistent)

        with open(file, 'r') as f:
            data = f.readlines()

        natoms = int(data[sline+7].split()[-1])
        numbers = np.zeros(natoms, dtype=int)
        positions = np.zeros((natoms, 3), dtype=np.float64)
        lattice = np.zeros((3, 3), dtype=np.float64)

        for i in range(natoms):
            atomline = data[sline+10+i].split()
            numbers[i] = int(atomline[2])
            for j in range(3):
                positions[i][j] = np.float64(atomline[4+j])

        return natoms, numbers, positions, lattice

    def set_wf_cell(self, file):
        """
        This method sets the geometry of the crystal to that employed for the
        construction of the wave function.

        Thus, it could be:

          - a unit cell, if the CRYSTAL output is related to a
            :math:`\Gamma`-point frequency calculation, or

          - a supercell, if either SUPERCELL of SCELPHONO were employed.

        Parameters
        ----------

        file: str
            Path of the CRYSTAL14/17 output file.

        Returns
        -------

        natom: int
            Number of atoms in the crystal cell

        numbers: ndarray
            Array containing atomic numbers with `int` type.

        positions: ndarray
            2D array containing atomic fractional coordinates with `float`
            type.

        lattice: ndarray
            2D array containing cell lattice vectors.

        """
        sline = self._get_start_line(file, geometry_wf)

        with open(file, 'r') as f:
            data = f.readlines()

        natoms = int(data[sline+8].split()[-1])
        numbers = np.zeros(natoms, dtype=int)
        positions = np.zeros((natoms, 3), dtype=np.float64)
        lattice = np.zeros((3, 3), dtype=np.float64)

        for i in range(natoms):
            atomline = data[sline+11+i].split()
            numbers[i] = int(atomline[2])
            for j in range(3):
                positions[i][j] = np.float64(atomline[4+j])

        for i in range(sline + 11 + natoms, len(data)):
            if lattice_vectors in data[i]:
                for j in range(3):
                    lattice[j] = np.asarray(data[i+2+j].split(),
                                            dtype=np.float64)
                break

        return natoms, numbers, positions, lattice

    def set_expansion(self, file):
        """
        This method sets the expasion matrix used to build the supercell.

        Thus, it could be:

          - a unit cell, if the CRYSTAL output is related to a
            :math:`\Gamma`-point frequency calculation, or

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

        with open(file, 'r') as f:
            data = f.readlines()

        expansion = np.zeros((3, 3), dtype=float)

        for i in range(3):
            line = data[sline+i].split()
            del(line[0])
            expansion[i] = np.asarray(line, dtype=float)
        return expansion

    def set_q_mesh(self, file):
        """
        This method sets the **q**-points mesh used to sample phonon
        properties.

        Parameters
        ----------

        file: str
            Path of the CRYSTAL14/17 output file.

        Returns
        -------

        qpoints: int
            Number of sampled **q**-points.

        qcoords: ndarray
            2D array containing the coordinates of the sampled **q**-points
            with `float` type.

        qweight: ndarray
            Array containing the weights of each phonon band with `float`
            type.

        qmesh: ndarray
            Array containing the Hessian interpolation mesh that was set with
            the INTERPHESS keyword in CRYSTAL

        """
        # Try to search for the Hessian interpolation
        hess = True
        sline = self._get_start_line(file, hess_interpolation)

        if sline is None:
            hess = False
            sline = self._get_start_line(file, scelphono_qpoints)

        with open(file, 'r') as f:
            data = f.readlines()

        qpoints = int(data[sline].split()[-4])
        qcoords = np.zeros((qpoints, 3), dtype=float)
        qmesh = np.zeros(3, dtype=float)
        qweight = np.zeros(qpoints, dtype=float)

        if hess:
            sline += 9
        else:
            for i in range(len(data)):
                if 'K       WEIGHT       COORD' in data[i]:
                    sline = i+1

        for i in range(qpoints):
            line = data[sline+i].split()
            qweight[i] = float(line[2])
            for j in range(3):
                qcoords[i, j] = float(line[j+3])

        sline += qpoints
        shrink_line = data[sline].split()
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
        sline = self._get_start_line(file, central_energy)

        with open(file, 'r') as f:
            data = f.readlines()

        return float(data[sline].split()[2])

    def set_phonons(self, file):
        """
        """
        phonons = {}
        band_counter = 0

        with open(file, 'r') as f:
            data = f.readlines()

        for i in range(len(data)):

            if frequency_header in data[i]:

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
                band_counter +=1

                if band_counter == self.qpoints:
                    break

        return phonons

    def _get_start_line(self, file, search_string):
        found = False
        sline = 0
        with open(file, 'r') as f:
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

    _data = {
        'points': 0.,
        'unitcell': [],
        'supercell': [],
        'expansion': np.identity(3, dtype=int),
        'energy': 0.,
        'kpoints': 1,
        'qpoints': 1,
        'qcoords': {},
        'nphonon': 0,
        'phonons': {},
        'weights': {},
        'shrinkf': np.ones(3, dtype=int)
        }

    def __init__(self, crystal_output):
        BasicReader.__init__(self, crystal_output)
        return

    def load(self, file):
        """

        """
        if not self.is_qha(file):
            self.error = 'This does not appear to be a CRYSTAL14/17 '
            self.error += 'output \nrelated to QHA calculations'
            return

        self.points = self.set_qha_points(file)
        if self.points < 1:
            self.error = 'This QHA file seems to not contain any volume'
            self.error += 'information'
            return

        elif self.points < 4:
            self.error = 'Insufficient number of unit cell volumes explored'
            return

        self.restarted_on = self.is_restarted(file)
        if self.restarted_on:
            self.error = 'At present, restarted calculations are not '
            self.error += 'supported.\nSorry for the inconvenience.'
            return

        self.supercell_on = self.is_supercell(file)

        if self.supercell_on:
            self.dim = self.set_expansion(file)

        self.qpoints, self.qcoords, self.weights = self.set_dummy_qmesh()
        ordered_volumes = self.set_volume(file)
        cells = self.set_unit_cells(file)
        self._data['unitcell'] = self.set_ordered_cells(cells,
                                                        ordered_volumes)

        self._data['energy'] = self.set_energy(file)

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
        with open(file, 'r') as f:
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
        with open(file, 'r') as f:
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
        with open(file, 'r') as f:
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
        return self._data['points']

    @points.setter
    def points(self, value: int):
        """
        Set the number of unit cell volumes explored in QHA analysis.
        """
        self._data['points'] = value
        return

    @property
    def dim(self):
        """
        Get the expansion matrix employed to build the supercell.
        """
        return self._data['expansion']

    @dim.setter
    def dim(self, expansion):
        """
        Set the expansion matrix employed to build the supercell.
        """
        self._data['expansion'] = expansion.copy()
        return

    @property
    def natom(self):
        """
        Get the number of atoms in the unit cell.
        """
        return int(self._data['unitcell'][0]['natom'] / self.kpoints)

    @property
    def kpoints(self):
        """
        Get the number of sampled *k*-points, determined from the expansion
        matrix.
        """
        return int(np.around(np.linalg.det(self._data['expansion']), 0))

    @property
    def energy(self):
        """
        Get the unit cell (if phonon dispersion relations or if
        :math:`\Gamma`-point frequencies) or the supercell energy.
        """
        return self._data['energy'] / self.kpoints

    @property
    def volume(self):
        """
        Get the unit cell volumes
        """
        volumes = np.zeros(self.points, dtype=float)
        for i in range(self.points):
            volumes[i] = np.linalg.det(self._data['unitcell'][i]['lattice'])
        return volumes / self.kpoints

    @property
    def nphonon(self):
        """
        Get the number of frequencies per band in unit cell (if phonon
        dispersion relations or if :math:`\Gamma`-point frequencies)
        or in the supercell.
        """
        return self.natom * 3

    @property
    def qpoints(self):
        """
        Get the number of sampled **q**-points. For QHA, this is the same as
        the number of *k*-points.
        """
        return self._data['qpoints']

    @qpoints.setter
    def qpoints(self, value: int):
        """
        Set the number of sampled **q**-points.
        """
        self._data['qpoints'] = value
        return

    @property
    def qcoords(self):
        """
        Get the coordinates of sampled **q**-points, in dict format.
        """
        return self._data['qcoords']

    @qcoords.setter
    def qcoords(self, array):
        """
        Set the coordinates of sampled **q**-points, in dict format.
        """
        for i in range(len(array)):
            self._data['qcoords'][i] = array[i]
        return

    @property
    def weights(self):
        """
        Get the weights of each phonon band, in dict format.
        """
        return self._data['weights']

    @weights.setter
    def weights(self, array):
        """
        Set the weights of each phonon band.
        """
        for i in range(len(array)):
            self._data['weights'][i] = array[i]
        return

    @property
    def shrinkf(self):
        """
        Get the Hessian interpolation mesh used in the INTERPHESS keyword.
        """
        return self._data['shrinkf']

    @shrinkf.setter
    def shrinkf(self, array):
        """
        Set the Hessian interpolation mesh used in the INTERPHESS keyword.
        """
        self._data['shrinkf'] = array.copy()
        return

    @property
    def phonons(self):
        """
        Get the phonon bands, in dict format.
        """
        return self._data['phonons']

    @phonons.setter
    def phonons(self, dictionary):
        """
        Set the phonon bands, in dict format.
        """
        self._data['phonons'] = dictionary
        return

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

        with open(file, 'r') as f:
            data = f.readlines()

        points = 0
        for i in range(sline, sline+100):
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
            :math:`\Gamma`-point frequency calculation, or

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

        with open(file, 'r') as f:
            data = f.readlines()

        expansion = np.zeros((3, 3), dtype=float)

        for i in range(3):
            line = data[sline+i].split()
            del(line[0])
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

        with open(file, 'r') as f:
            data = f.readlines()

        energy = np.zeros(self.points, dtype=float)
        for i in range(self.points):
            energy[i] = float(data[sline+i].split()[1])

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

        with open(file, 'r') as f:
            data = f.readlines()

        volume = np.zeros(self.points, dtype=float)
        for i in range(self.points):
            volume[i] = float(data[sline+i].split()[0])

        return volume

    def set_unit_cells(self, file):
        """
        This method sets the optimized crystal cells.

        Parameters
        ----------

        file: str
            Path of the CRYSTAL14/17 output file.

        Returns
        -------

        cells: dict
            Dictionary containing the unit cells.

        """
        with open(file, 'r') as f:
            data = f.readlines()

        cells = []
        for i in range(len(data)):
            if qha_opt in data[i]:
                cell = {}
                cell_data = self.set_optimized_cell(file, i)
                cell['natom'] = cell_data[0]
                cell['numbers'] = cell_data[1].copy()
                cell['positions']= cell_data[2].copy()
                cell['lattice']= cell_data[3].copy()
                cells.append(cell)
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
            volume = np.linalg.det(cells[i]['lattice'])
            for j in range(len(cells)):
                if np.isclose(volume, ordered_volumes[j]):
                    indexes.append(j)

        ordered_cells = []
        for i in range(len(cells)):
            ordered_cells.append(cells[indexes.index(i)])

        return ordered_cells

    def set_optimized_cell(self, file, idx: int):
        """
        This method sets the geometry of the crystal with optimized geometry.

        Parameters
        ----------

        file: str
            Path of the CRYSTAL14/17 output file.

        Returns
        -------

        natom: int
            Number of atoms in the crystal cell

        numbers: ndarray
            Array containing atomic numbers with `int` type.

        positions: ndarray
            2D array containing atomic fractional coordinates with `float`
            type.

        lattice: ndarray
            2D array containing cell lattice vectors.

        """
        with open(file, 'r') as f:
            data = f.readlines()

        natoms = int(data[idx+8].split()[-1])
        numbers = np.zeros(natoms, dtype=int)
        positions = np.zeros((natoms, 3), dtype=np.float64)
        lattice = np.zeros((3, 3), dtype=np.float64)

        for i in range(natoms):
            atomline = data[idx+11+i].split()
            numbers[i] = int(atomline[2])
            for j in range(3):
                positions[i][j] = np.float64(atomline[4+j])

        for i in range(idx + 11 + natoms, len(data)):
            if lattice_vectors in data[i]:
                for j in range(3):
                    lattice[j] = np.asarray(data[i+2+j].split(),
                                            dtype=np.float64)
                break

        return natoms, numbers, positions, lattice

    def set_phonons(self, file):
        """
        """
        phonons = {}
        nfreq = self.natom * self.qpoints * 3
        phonon_matrix = np.zeros((nfreq, self.points), dtype=float)

        with open(file, 'r') as f:
            data = f.readlines()

        freq_idx = []
        for i in range(len(data)):
            if qha_freq in data[i]:
                if data[i].split()[0] == 'FREQUENCY':
                    freq_idx.append(i+3)

        for i in range(nfreq):
            for j in range(self.points):
                frequency_line = data[freq_idx[i] + j].split()
                phonon_matrix[i, j] = float(frequency_line[1])

        bands = phonon_matrix.reshape(self.qpoints, self.natom*3, self.points)

        for i in range(self.qpoints):
            phonons[i] = bands[i].copy()

        return phonons

    def set_dummy_qmesh(self):
        """
        """
        qpoints = self.kpoints
        qcoords = {}
        qweights = {}
        for i in range(qpoints):
            qcoords[i] = np.zeros(3, dtype=float)
            qweights[i] = 1.
        return qpoints, qcoords, qweights

    def _get_start_line(self, file, search_string):
        sline = 0
        with open(file, 'r') as f:
            for line in f:
                if search_string in line:
                    break
                else:
                    sline += 1
        return sline


class CrystalSOECReader(BasicReader):
    """
    """

    _data = {
        'stiffness': np.zeros((6,6), dtype=float),
        'density': 0.
        }

    def __init__(self, crystal_soec_output):
        """
        """
        BasicReader.__init__(self, crystal_soec_output)
        return

    def load(self, file):
        """
        """
        if not self.is_soec_output(file):
            self.error = 'This does not appear to be a CRYSTAL14/17 '
            self.error += 'output \nrelated to elastic constants calculations'
            return

        if not self.is_output_completed(file):
            self.error = 'This CRYSTAL14/17 SOEC calculation has not been '
            self.error += 'completed'
            return

        self._data['density'] = self.crystal_density_init(file)

        sline = self.soec_start_line(file)

        with open(file, 'r') as f:
            data = f.readlines()

        for i in range(6):
            line = data[i+sline].replace('|','').split()
            for j in range(6-i):
                try:
                    self._data['stiffness'][i, j+i] = float(line[j])
                except ValueError:
                    self.error = 'Problem in collecting the SOEC matrix '
                    self.error += 'element {} {}\n'.format(i+1, j+i+1)
                    self.error += 'Please, check the CRYSTAL output file'
                    return

        self._data['stiffness'] += np.triu(self._data['stiffness'], 1).T

        self.completed = True
        return

    @property
    def stiffness(self):
        """
        The elastic stiffness tensor in Voigt notation.
        """
        return self._data['stiffness']

    @property
    def density(self):
        """
        The crystal density, expressed in kg m^-3.
        """
        return self._data['density']

    def is_soec_output(self, file):
        """
        This method checks if the CRYSTAL14/17 output file is related to
        an ELASTCON calculation.

        Parameters
        ----------

        file: str
            Path to the CRYSTAL14/17 output file.


        Returns
        -------

        bool
            Returns True if the output is correct, otherwise False.

        """
        with open(file, 'r') as f:
            for line in f:
                if soec_option in line:
                    return True
            return False

    def is_output_completed(self, file):
        """
        This method checks if the CRYSTAL14/17 output file is related to a
        completed ELASTCON calculation.

        Parameters
        ----------

        file: str
            Path to the CRYSTAL14/17 output file.


        Returns
        -------

        bool
            Returns True if the simulation was succesfully concluded,
            otherwise False.

        """
        with open(file, 'r') as f:
            for line in f:
                if soec_completed in line:
                    return True
            return False

    def soec_start_line(self, file):
        """
        This method returns the line of the ouput where the elastic
        stiffness matrix should be present.

        Parameters
        ----------

        file: str
            Path to the CRYSTAL14/17 output file.


        Returns
        -------

        int
            Line number in the CRYSTAL14/17 output file corresponding to the
            first row of the elastic constants tensor in Voigt notation.

        """
        counter = 0
        with open(file, 'r') as f:
            for line in f:
                if soec_ec in line:
                    counter += 2
                    break
                else:
                    counter += 1
        return counter

    def crystal_density_init(self, file):
        """
        This method returns the density of the crystal from the geometry
        provided as input in the CRYSTAL14/17 output file.

        Parameters
        ----------

        file: str
            Path to the CRYSTAL14/17 output file.


        Returns
        -------

        float
            Crystal density, expressed in kg m^-3.

        """
        density = 0.
        consistent = False
        with open(file, 'r') as f:
            for line in f:
                if geometry_consistent in line:
                    consistent = True
                if primitive_cell in line:
                    if consistent:
                        density = float(line.split()[-2]) * 1000.
                        break
                    else:
                        continue
        return density
