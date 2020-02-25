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

from .table import SOECTableFile


class QuantasSOECExport(object):
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
                f['stiffness'][:]
            except KeyError:
                self.error = 'This does not appear to be a Quantas SOEC output'
                self.error += ' file'
                return

            self.keys = ['E', 'LC', 'G', 'Nu']
            for key in f.keys():
                if not 'polar_' in key:
                    continue
                else:
                    if 'waves' in key and not 'waves' in self.keys:
                        self.keys.append('waves')
                    self.data[key] = {}
                    self.data[key]['values'] = f[key][:]
                    self.data[key]['unit'] = f[key].attrs['unit']
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
                    if first:
                        degree = self.data[key]['values'][:,0]
                        data = self.data[key]['values'][:,1:]
                        unit = self.data[key]['unit']
                        first = False
                    else:
                        data = np.hstack((data, self.data[key]['values'][:,1:]))

            if os.path.exists(outfile):
                if confirm("File '{}' exists. Overwrite it?".format(outfile)):
                    pass
                else:
                    outfile = click.prompt('Please, enter a file name')

            tabled_data = SOECTableFile()
            tabled_data.write(outfile, degree, data, selection, unit)
                
            echo("Data successfully exported to '{}'".format(outfile))
            
            is_running = confirm('Would you like to export other data?')
        return
