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
from numpy import linalg

from quantas.core.reader import BasicReader


class SOECInputFileReader(BasicReader):
    """
    This class is responsible of parsing an input file (extension .dat)
    in a raw way, meaning that it checks the basic file format.

    Attributes
    ----------
    data:
        Dictionary containing the SOEC data.

    completed: bool
        Flag used to assert if the input file was succesfully read.

    error: str
        Text explaining the error encountered in the input file.
    """

    completed = False
    error = ''

    _data ={
        'job': 'Unknown',
        'stiffness': None,
        'density': 0.
        }

    def __init__(self, soec_input):
        """ Class constructor.

        Parameters
        ----------
        soec_input: str
            Path of the SOEC input file.
        """
        BasicReader.__init__(self, soec_input)

        self._load(soec_input)
        return

    @property
    def jobname(self):
        """ This method returns the name of the current job. """
        return self._data['job']

    @property
    def stiffness(self):
        """ This method returns the stiffness matrix. """
        return self._data['stiffness']

    @property
    def density(self):
        """ This method returns the density of the crystal. """
        return self._data['density']

    def _load(self, filename):
        """
        This method parses the input file, looking for the SOEC data
        """

        with open(filename, 'r') as f:
            data = f.readlines()

        if len(data) < 6:
            self.error = 'File length is too short to contain '
            self.error += 'the SOEC matrix'
            return

        # Does the first line contain the job name?
        line = data[0]
        try:
            np.float64(line.strip()[0])
            start = 0
        except ValueError:
            self._data['job'] = line.rstrip()
            start = 1

        # Collect the SOEC matrix...
        matstr = [ data[i+start] for i in range(6) ]
        # ... and convert it to float numbers
        matrix = np.zeros((6,6), dtype=np.float64)
        try:
            for i in range(6):
                line = matstr[i].split()
                n = len(line)
                for j in range(n):
                    matrix[i][j+(6-n)] = line[j]
        except ValueError:
            self.error = 'SOEC components should be float numbers'
            return
        else:
            # Take its from from the string matrix
            form = []
            for i in range(6):
                form.append(len(matstr[i].split()))
            # Check matrix symmetry and make it 6x6
            if form == [6,5,4,3,2,1]:  # upper triangular
                matrix = matrix + np.triu(matrix, 1).transpose()
            if form == [1,2,3,4,5,6]:  # lower triangular
                matrix = matrix + np.tril(matrix, 1).transpose()
            if linalg.norm(matrix - matrix.transpose()) > 1e-3:
                self.error = 'The SOEC matrix should be symmetric '
                self.error += 'or triangular'
                return
            else:
                # Store the matrix
                self._data['stiffness'] = matrix.copy()

        # Is the density present?
        dline = start+6
        try:
            self._data['density'] = np.float64(data[dline].split()[0])
        except IndexError:
            pass
        except ValueError:
            self.error = 'Density value should be a float number'
            self.completed = False
            return
        self.completed = True
        return
