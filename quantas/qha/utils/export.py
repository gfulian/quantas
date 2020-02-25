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
from quantas.cmdline.utils.messages import echo
from quantas.cmdline.utils.messages import confirm

from .table import QHATableFile


class QuantasQHAExport(object):
    """
    """
    x = {}
    y = {}
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
                self.x['values'] = f['T'][:]
                self.x['unit'] = f['T'].attrs['unit']
                self.y['values'] = f['P'][:]
                self.y['unit'] = f['P'].attrs['unit']
            except KeyError:
                self.error = 'This does not appear to be a Quantas HA output '
                self.error += 'file'
                return

            self.keys = []
            for key in f.keys():
                if key == 'T' or key == 'P':
                    continue
                else:
                    self.keys.append(key)
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

            if os.path.exists(outfile):
                if confirm("File '{}' exists. Overwrite it?".format(outfile)):
                    pass
                else:
                    outfile = click.prompt('Please, enter a file name')

            tabled_data = QHATableFile()
            tabled_data.write(
                outfile,
                self.x['values'], 'T', self.x['unit'],
                self.y['values'], 'P', self.y['unit'],
                self.data[selection]['values'],
                selection,
                self.data[selection]['unit'],
                )
            echo("Data successfully exported to '{}'".format(outfile))
            
            is_running = confirm('Would you like to export other data?')
        return
    
