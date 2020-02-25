# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

import os
import yaml
import numpy as np
from scipy import constants as cs
import xml.etree.ElementTree as ET

from quantas.core.reader import BasicReader

from quantas.utils.chemistry.symbols import symbol2number as s2n
from quantas.utils.physics.units import convert_frequency as cf

h = cs.Planck
eV = cs.value(u'electron volt')

class PhonopyReader(BasicReader):
    """
    This class reads and stores phonon properties of crystals obtained
    through the phonopy code.

    .. note::

        There must be other two files in the same folder of the input
        file, namely the phonopy displacement output
        (``phonopy_disp.yaml``, or the deprecated ``disp.yaml``)
        and the VASP xml output. The ``phonopy_disp.yaml`` must be renamed
        with same base name of the input file, terminating ``_disp.yaml``.
        The ``vasprun.xml`` must be renamed with the same
        root name of the phonopy output file.

    .. note::

        At the moment, only calculation performed using VASP are fully
        supported by Quantas. You can still employ this software to
        create an input for (Q)HA calculations, but energy information
        will be missing in the generated input file. In this case,
        you have to provide the missing data by hand.
        Sorry for the inconvenience.

    """

    def __init__(self, phonopy_output):
        BasicReader.__init__(self, phonopy_output)
        # Data dictionary with default values
        return

    def load(self, file):
        """
        This method reads the phonopy output file (.yaml) provided as input.

        Parameters
        ----------

        file: str
            phonopy output file (.yaml)

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
        #
        # Check if the 'phonopy_disp.yaml' exists
        #
        dispfile = os.path.splitext(file)[0] + '_disp.yaml'
        if not os.path.exists(dispfile):
            self.error = 'Phonopy displacement file (phonopy_disp.yaml) not '
            self.error += 'available!\n'
            self.error += 'Will now exit\n'
            return
        #
        # Check if the VASP xml ouput exists
        #
        vasprun = os.path.splitext(file)[0] + '.xml'
        if os.path.exists(vasprun):
            root = ET.parse(vasprun).getroot()
            energy_results = root.findall('calculation/energy/i')
            for energy_elem in energy_results:
                if energy_elem.attrib['name'] == 'e_wo_entrp':
                    self._data['energy'] = float(energy_elem.text)

        if os.path.exists(file):
            yaml_error, exception, idata = self.load_yaml(file)
            if yaml_error:
                self.error = exception
                return
        else:
            self.error = ' File not found!'.format(file)
            return

        self.unitcell = self.set_unit_cell(idata)

        error, exception, dispdata = self.load_yaml(dispfile)
        if error:
            self.error = exception
            return
        self.supercell = self.set_displacements(dispdata)

        self.set_phonons(idata)

        self.completed = True
        return

    @property
    def unitcell(self):
        return (self._data['unitcell']['natom'],
                self._data['unitcell']['lattice'],
                self._data['unitcell']['positions'],
                self._data['unitcell']['numbers'])

    @unitcell.setter
    def unitcell(self, cell_data):
        self._data['unitcell']['natom'] = cell_data[0]
        self._data['unitcell']['numbers'] = cell_data[1]
        self._data['unitcell']['positions'] = cell_data[2]
        self._data['unitcell']['lattice'] = cell_data[3]
        return

    @property
    def supercell(self):
        return (self._data['supercell']['lattice'],
                self._data['supercell']['positions'],
                self._data['supercell']['numbers'])

    @supercell.setter
    def supercell(self, cell_data):
        self._data['supercell']['natom'] = cell_data[0]
        self._data['supercell']['numbers'] = cell_data[1]
        self._data['supercell']['positions'] = cell_data[2]
        self._data['supercell']['lattice'] = cell_data[3]
        return

    @property
    def dim(self):
        return self._data['expansion']

    @dim.setter
    def dim(self, matrix):
        self._data['expansion'] = matrix.copy()
        return

    @property
    def natom(self):
        return self._data['unitcell']['natom']

    @property
    def lattice(self):
        return self._data['unitcell']['lattice']

    @property
    def volume(self):
        return np.linalg.det(self._data['unitcell']['lattice'])

    @property
    def energy(self):
        return self._data['energy']

    @property
    def kpoints(self):
        return int(np.around(np.linalg.det(self._data['expansion']), 0))

    @property
    def qpoints(self):
        return self._data['qpoints']

    @property
    def nphonon(self):
        return self._data['nphonon']

    @property
    def phonons(self):
        return self._data['phonons']

    @property
    def qcoords(self):
        return self._data['qcoords']

    @property
    def weights(self):
        return self._data['weights']

    @property
    def shrinkf(self):
        return self._data['shrinkf']

    def load_yaml(self, file):
        idata = None
        error = False
        exception = None
        with open(file, 'r') as stream:
            try:
                idata = yaml.safe_load(stream)
            except yaml.YAMLError as exception:
                error = True
            return error, exception, idata

    def set_displacements(self, idata):
        """
        """
        self.dim = np.asarray(idata['supercell_matrix'])

        natom = len(idata['supercell']['points'])
        numbers = np.zeros(natom, dtype=int)
        positions = np.zeros((natom, 3), dtype=int)
        lattice = idata['supercell']['lattice']

        for i in range(natom):
            symbol = idata['supercell']['points'][i]['symbol']
            numbers[i] = s2n(symbol)
            positions[i] = idata['supercell']['points'][i]['coordinates']

        return natom, numbers, positions, lattice

    def set_unit_cell(self, idata):
        """
        """
        natom = idata['natom']
        numbers = np.zeros(natom, dtype=int)
        positions = np.zeros((natom, 3), dtype=int)
        lattice = idata['lattice']

        for i in range(natom):
            symbol = idata['points'][i]['symbol']
            numbers[i] = s2n(symbol)
            positions[i] = idata['points'][i]['coordinates']

        return natom, numbers, positions, lattice

    def set_phonons(self, idata):
        """
        """
        nf = 3 * idata['natom']
        nq = idata['nqpoint']

        self._data['nphonon'] = nf
        self._data['qpoints'] = nq

        for i in range(nq):
            qdata = idata['phonon'][i]
            self._data['qcoords'][i] = qdata['q-position']
            self._data['weights'][i] = qdata['weight']
            band = np.zeros(nf, dtype=float)
            for j in range(nf):
                band[j] = cf(qdata['band'][j]['frequency'], 'THz', 'cm^-1')
            self._data['phonons'][i] = band
        return
