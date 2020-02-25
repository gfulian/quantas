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
import yaml

from quantas.core.reader import BasicReader


class QHAInputFileReader(BasicReader):
    """
    This class is responsible of parsing an input file in YAML format. It
    checks the basic file format and collect data.

    After a successfull data collection, it checks if the required data is
    fine for (Q)HA calculations.

    Attributes
    ----------

    completed: bool
        Flag that is True when the input file was completely read, otherwise
        it should be False.

    error: None or str
        Error encountered during data handling and storing.

    data: dict
        Dictionary containing the data collected from the input.
    """

    def __init__(self, qha_input):
        BasicReader.__init__(self, qha_input)
        return

    def load(self, filename):
        """
        Read the input filename and collect the information stored in it.

        Parameters
        ----------

        filename: str
            Path to the input file
        
        """
        with open(filename, 'r') as stream:
            try:
                self._data = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                if hasattr(exc, 'problem_mark'):
                    self.error = 'Error in parsing YAML file: \n'
                    self.error += '  - '+ str(exc.problem_mark) + '\n'
                    self.error += '  - '+ str(exc.problem) + '\n'
                    self.error += '\nPlease, correct data and retry.'
                else:
                    self.error = 'Error in parsing YAML file: unknown problem.'
                    self.error += '\n\nPlease, correct data and retry.'
            except UnicodeDecodeError:
                self.error = "The provided input is not a YAML file."
        self._check()
        return

    def _check(self):
        if isinstance(self._data, type(None)):
            return
        if not 'volume' in self._data:
            self.error = 'No volume values in input file'
            return
        if not 'energy' in self._data:
            self.error = 'No energy values in input file'
            return
        if not 'phonon' in self._data:
            self.error = 'No phonon data in input file'
            return
        self.completed = True
        return

    @staticmethod
    def _get_float(item):
        try:
            return np.float64(item)
        except ValueError:
            return np.array(item[0].split(), dtype=np.float64)

    @staticmethod
    def _get_int(item):
        try:
            return np.int64(item)
        except ValueError:
            return np.array(item[0].split(), dtype=np.int64)

    @staticmethod
    def _get_string(item):
        return item.rstrip()

    @property
    def data(self):
        return self._data

    @property
    def jobname(self):
        if 'job' in self._data:
            return self._get_string(self._data['job'])
        else:
            return 'Unknown'

    @property
    def natoms(self):
        return self._get_int(self._data['natom'])

    @property
    def kpoints(self):
        return int(np.around(np.linalg.det(self._data['supercell']),0))

    @property
    def qpoints(self):
        return self._get_int(self._data['qpoints'])

    @property
    def nvol(self):
        """ Number of volumes in the input """
        if isinstance(self.volume, np.ndarray):
            return self.volume.shape[0]
        else:
            return 1

    @property
    def volume(self):
        return self._get_float(self._data['volume'])

    @property
    def energy(self):
        return self._get_float(self._data['energy'])

    @property
    def frequencies(self):
        matrix = np.zeros((self.qpoints, self.natoms*3, self.nvol),
                          dtype=np.float64)
        for i, q in enumerate(self._data['phonon']):
            for j, frequency in enumerate(q['band']):
                matrix[i][j] = np.array(
                    self._get_float(frequency['frequency']),
                    dtype=np.float64
                    )
        return matrix

    @property
    def weights(self):
        w = np.zeros(self.qpoints,dtype=np.float64)
        for i, q in enumerate( self._data['phonon'] ):
            w[i] = self._get_float(q['weight'])
        return w

    @property
    def total_q_points(self):
        return self.weights.sum()
