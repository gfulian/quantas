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


class EoSInputFileReader(BasicReader):
    """
    This class is responsible of parsing an input file (extension .dat)
    in a raw way, meaning that it checks the basic file format.

    Attributes
    ----------
    data:
        Dictionary containing the input data.

    completed: bool
        Flag used to assert if the input file was succesfully read.

    error: str
        Text explaining the error encountered in the input file.

    """

    _data ={
        'job': 'Unknown',   # Job name
        'v': None,          # Volume values
        'sigmav': None,     # Error in volume measures
        'p': None,          # Pressure values
        'sigmap': None,     # Error in pressure measures
        'a': None,          # a lattice parameter values
        'sigmaa': None,     # Error in a-axis measures
        'b': None,          # b lattice parameter values
        'sigmab': None,     # Error in b-axis measures
        'c': None,          # c lattice parameter values
        'sigmac': None,     # Error in a-axis measures
        }

    _formats = [
        'v', 'sigmav',
        'p', 'sigmap',
        'a', 'sigmaa',
        'b', 'sigmab',
        'c', 'sigmac'
        ]

    _formats = {
        'p': ['p'],
        'v': ['v'],
        'a': ['a'],
        'b': ['b'],
        'c': ['c'],
        'sigmap': ['sigmap', 'sigp'],
        'sigmav': ['sigmav', 'sigv'],
        'sigmaa': ['sigmaa', 'siga'],
        'sigmab': ['sigmab', 'sigb'],
        'sigmac': ['sigmac', 'sigc'],
        }

    _format = {}

    def __init__(self, eos_input):
        """
        Constructor method of the class.

        Parameters
        ----------

        eos_input: str
            Path of the input file.

        """
        BasicReader.__init__(self, eos_input)
        return

    @property
    def jobname(self):
        return self._data['job']

    @property
    def format(self):
        return self._format

    @property
    def data(self):
        return self._data

    @property
    def volume(self):
        return self._data['v']

    @property
    def sigma_v(self):
        return self._data['sigmav']

    @property
    def pressure(self):
        return self._data['p']

    @property
    def sigma_p(self):
        return self._data['sigmap']

    @property
    def a(self):
        return self._data['a']

    @property
    def sigma_a(self):
        return self._data['sigmaa']

    @property
    def b(self):
        return self._data['b']

    @property
    def sigma_b(self):
        return self._data['sigmab']

    @property
    def c(self):
        return self._data['c']

    @property
    def sigma_c(self):
        return self._data['sigmac']

    def is_comment(self, text):
        """
        Check if a parsed line is a comment or an empty line.

        Parameters
        ----------

        text: str
            Line of text to be checked.

        Returns
        -------

        bool
            True if the line is empty or a comment, otherwise False.

        """
        if (len(text.strip()) == 0 or text.strip()[0] == '#'):
            return True
        else:
            return False

    def load(self, eos_input):
        """
        Read the input file and store the data that are contained.

        Parameters
        ----------

        eos_input: str
            Path to the EoS input file.

        """
        with open(eos_input, 'r') as f:
            data = f.readlines()

        for i in range(len(data)):
            line = data[i]
            # If empty or commented (#) lines are found, skip them
            if self.is_comment(line):
                continue

            keyword = line.split()[0]
            if keyword.lower() == 'job':
                self.set_job_name(data, i+1)
            if keyword.lower() == 'format':
                self.set_data_format(data, i+1)
            if keyword.lower() == 'data':
                self.set_data(data, i)

        if not isinstance(self.error, type(None)):
            self.completed = False
        else:
            self.completed = True
        return

    def set_job_name(self, inpdata, idx):
        """
        Set the name of the job related to the provided input file.

        Parameters
        ----------

        inpdata: list
            List containing the input file data.

        idx: int
            Index of the line where the name of the job is supposed to be.

        """
        while self._data['job'] == 'Unknown':
            line = inpdata[idx].rstrip('\n')
            if (len(line) == 0 or line.strip()[0] == '#'):
                idx += 1
            else:
                self._data['job'] = line
        return

    def set_data_format(self, inpdata, idx):
        """
        Set the format used in the input file to provide pressure and
        structural data.

        Parameters
        ----------

        inpdata: list
            List containing the input file data.

        idx: int
            Index of the line where the data format is supposed to be.

        """
        found = False
        while not found:
            line = inpdata[idx]
            if self.is_comment(line):
                idx += 1
                continue
            else:
                found = True
                dformat = inpdata[idx].split()
                for i in range(len(dformat)):
                    for key, items in self._formats.items():
                        if dformat[i].lower() in items:
                            self._format[dformat[i].lower()] = i
        return

    def set_data(self, inpdata, idx):
        """
        Set the data with those reported in the input file .

        Parameters
        ----------

        inpdata: list
            List containing the input file data.

        idx: int
            Index of the line where the data format is supposed to be.

        """
        if self._format['p'] == None:
            self.error = 'No pressure data provided in input file'
            return

        data = {}
        for quantity in self._format:
            data[quantity] = []

        stop = False
        while not stop:
            idx += 1
            if idx >= len(inpdata):
                stop = True
                continue

            line = inpdata[idx]
            if self.is_comment(line):
                continue

            values = line.split()
            if len(values) == 0:
                stop = True
                continue

            for quantity, n in self._format.items():
                data[quantity].append(float(values[n]))

        for quantity in self._format:
            self._data[quantity] = np.asarray(data[quantity])
        return
