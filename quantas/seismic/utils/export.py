# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

import h5py
import os
import click
import numpy as np
from quantas.cmdline.utils.messages import echo
from quantas.cmdline.utils.messages import confirm

# ---------------------------------------------------------------------------#
# Formatting constants
# ---------------------------------------------------------------------------#
dl = '{: >8.6f} {: >8.6f} {: >9.5f} {: >8.3f} {: >8.4f} {: >8.4f} {: >8.4f} '
dl += '{: >9.5f} {: >8.3f} {: >9.5f} {: >9.5f} {: >9.5f} {: >9.4f} {: >9.4f}\n'


class QuantasSeismicExport(object):
    """
    """
    data = {}
    error = None

    def __init__(self):
        """ Constructor method """
        return

    def load(self, filename):
        """

        """
        self.basename = os.path.splitext(filename)[0]
        with h5py.File(filename, 'r') as f:
            try:
                self.info = f.attrs['info']
            except KeyError:
                self.error = 'This does not appear to be a Quantas output file'
                return

            try:
                f['primary'][:]
            except KeyError:
                self.error = 'This does not appear to be a Quantas '
                self.error += 'seismic output file'
                return

            self.keys = ['ssecondary', 'fsecondary', 'primary']
            for key in f.keys():
                self.data[key] = {}
                self.data[key]['values'] = f[key][:]

        return

    def run(self):
        echo(self.info)
        is_running = True

        while is_running:
            echo('')
            selection = click.prompt(
                'Select data to export', type=click.Choice(self.keys),
                show_choices=True)
            outfile = self.basename + '_' + selection + '.dat'

            first = True
            for key in self.data.keys():
                if selection in key:
                    data = self.data[key]['values']

            if os.path.exists(outfile):
                if confirm("File '{}' exists. Overwrite it?".format(outfile)):
                    pass
                else:
                    outfile = click.prompt('Please, enter a file name')

            self.new_info = self.info.splitlines()
            self.new_info[9] = 'This file contains the following dataset:'

            header = ''
            for line in self.new_info:
                if ' data ' in line:
                    if selection not in line:
                        continue
                header += line
                header += '\n'

            n, _ = data.shape

            with open(outfile, 'w') as f:
                f.write(header)
                f.write('\n\n')
                for i in range(n):
                    line = data[i]
                    f.write(dl.format(*line))

            echo("Data successfully exported to '{}'".format(outfile))

            is_running = confirm('Would you like to export other data?')
        return

