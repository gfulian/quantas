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

soec_total = 'ELASTIC MODULI  (kBar)'
soec_clamped = 'SYMMETRIZED ELASTIC MODULI (kBar)'
soec_relaxed = 'TOTAL ELASTIC MODULI (kBar)'


class VASPSOECReader(BasicReader):
    """
    """

    _data = {
        'stiffness': np.zeros((6,6), dtype=float),
        'density': 0.
        }

    def __init__(self, vasp_outcar):
        """
        """
        BasicReader.__init__(self, vasp_outcar)
        return

    def load(self, file):
        """
        """
        if not self.is_soec_output(file):
            self.error = "'{}' does not appear to be a VASP output\n".format(
                file)
            self.error += 'related to elastic moduli calculation'
            return

        sline = self.soec_start_line(file)

        with open(file, 'r') as f:
            data = f.readlines()

        for i in range(6):
            line = data[i+sline].split()
            for j in range(6):
                try:
                    self._data['stiffness'][i, j] = float(line[j+1]) / 10.
                except ValueError:
                    self.error = 'Problem in collecting the SOEC matrix '
                    self.error += 'element {} {}\n'.format(i+1, j+1)
                    self.error += 'Please, check the OUTCAR file'
                    return

        self._data['stiffness'][[3, 5]] = self._data['stiffness'][[5, 3]]
        self._data['stiffness'][:, [3, 5]] = self._data['stiffness'][:, [5, 3]]

        self.completed = True
        return

    @property
    def stiffness(self):
        """
        The elastic stiffness tensor in Voigt notation, expressed in GPa.
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
        This method checks if the VASP OUTCAR contains elastic moduli.

        Parameters
        ----------

        file: str
            Path to the VASP OUTCAR file.


        Returns
        -------

        bool
            Returns True if the output is correct, otherwise False.

        """
        with open(file, 'r') as f:
            for line in f:
                if soec_clamped in line:
                    return True
            return False

    def soec_start_line(self, file):
        """
        This method returns the line of the ouput where the elastic
        stiffness matrix should be present.

        It returns the line corresponding to the total (clamped ion + ionic
        relaxation) elastic moduli if it is present in the OUTCAR, otherwise
        it returns that related to the clamped ion.

        Parameters
        ----------

        file: str
            Path to the VASP OUTCAR file.


        Returns
        -------

        int
            Line number in the VASP OUTCAR file corresponding to the
            first row of the elastic constants tensor in Voigt notation.

        """
        relaxed = False
        counter = 0
        with open(file, 'r') as f:
            for line in f:
                if soec_clamped in line:
                    clamped_line = counter + 3
                if soec_relaxed in line:
                    relaxed = True
                    relaxed_line = counter + 3
                else:
                    counter += 1
        if relaxed:
            return relaxed_line
        else:
            return clamped_line
