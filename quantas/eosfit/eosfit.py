# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

"""
This module provides the Equation of State calculator object employed by
Quantas, alongside specific functions used by interactive runs.
"""
import numpy as np

# Basic claculator
from quantas.core.calculator import BasicCalculator

# File reader
from quantas.IO.eos_reader import EoSInputFileReader

from .utils.interactive import is_available
from .utils.interactive import confirm
from .utils.interactive import select_eos
from .utils.interactive import select_data
from .utils.interactive import get_label
from .utils.interactive import set_fixed_parameters
from .utils.interactive import set_parameters
from .utils.interactive import set_weights

# Equation of State formulations for P-V fits
from .utils.pv_eos import Murnaghan
from .utils.pv_eos import BirchMurnaghan
from .utils.pv_eos import NaturalStrain
from .utils.pv_eos import Vinet
from .utils.pv_eos import Tait


class EoSFitCalculator(BasicCalculator):
    """
    This class performs Equation of State (EoS) fitting on provided
    experimental or theoretical *P-V* data.

    """
    def __init__(self, settings):
        """ Constructor method of the class. """
        BasicCalculator.__init__(self)
        # units
        self._vunit = settings['lenght_unit']
        self._punit = settings['pressure_unit']
        self.log = settings['logfile']
        return

    def read_input(self, filename):
        """
        This method read the provided input file and stores a dictionary
        containing the necessary data.

        Parameters
        ----------
        filename: str
            Complete path of the inputfile.

        Returns
        -------
        dict
            Dictionary containing the input data.

        """
        
        self.echo('Reading input file: {0}'.format(filename))
        self.echo_debug('Instantiate file reader')
        data = EoSInputFileReader(filename)
        #
        # Error check
        if not data.completed:
            return data.error

        quantities = 4
        for q in ['v', 'a', 'b', 'c']:
            if not isinstance(data._data[q], np.ndarray):
                quantities -= 1
        if quantities == 0:
            error = 'No structural data provided in input file'
            return error

        self.input = data
        return

    def report_input_data(self):
        """
        This method write to the output stream the main information about
        the input.
        """
        data = self.input
        self.echo('')
        self.echo('Job: ' + data.jobname)
        self.echo('')
        self.echo('{0: ^10}  {1: ^10}  {2: ^10}  {3: ^10}  {4: ^10}'.format(
            'P (GPa)', 'V (A^3)', 'a (A)', 'b (A)', 'c (A)'))
        self.echo('{0:-^10}  {0:-^10}  {0:-^10}  {0:-^10}  {0:-^10}'.format(''))
        for i in range(data.pressure.shape[0]):
            msg = '{0: ^10.3f}'.format(data.pressure[i])
            if is_available(data.volume):
                msg += '  {0: ^10.5f}'.format(data.volume[i])
            else:
                msg += '  {0: ^10}'.format(str('--'))
            if is_available(data.a):
                msg += '  {0: ^10.5f}'.format(data.a[i])
            else:
                msg += '  {0: ^10}'.format(str('--'))
            if is_available(data.b):
                msg += '  {0: ^10.5f}'.format(data.b[i])
            else:
                msg += '  {0: ^10}'.format(str('--'))
            if is_available(data.c):
                msg += '  {0: ^10.5f}'.format(data.c[i])
            else:
                msg += '  {0: ^10}'.format(str('--'))
            self.echo(msg)
        self.echo('')
        # Report weights
        if is_available(data.sigma_p):
            self.echo(' - Found weights for pressure')
        if is_available(data.sigma_v):
            self.echo(' - Found weights for volume')
        if is_available(data.sigma_a):
            self.echo(' - Found weights for a-axis')
        if is_available(data.sigma_b):
            self.echo(' - Found weights for b-axis')
        if is_available(data.sigma_c):
            self.echo(' - Found weights for c-axis')
        self.echo('')
        return

    def export_hdf5(self, filename=None):
        """ Dummy method, as this calculator does not export data. """
        return

    def run(self):
        """
        This method starts the (interactive) calculation.
        """
        data = self.input.data

        is_running = True
        eos_defaults = None
        
        self.echo('{0}{1:-^78}{0}'.format(
            '#', 'Interactive equation of state fitting started'))
        while is_running:
            self.echo('')
            sel_data, data_name, values, linearized = select_data(data)
            self.echo('')
            labels = get_label(sel_data)
            change_data = False
            do_init_guess = True

            while not change_data:
                eos, order, eos_defaults = select_eos(eos_defaults)
                change_eos = False

                if do_init_guess:
                    guess = eos.guess(values, self.input.pressure)
                    do_init_guess = False

                while not change_eos:
                    parameters = set_parameters(guess, labels)
                    fix = set_fixed_parameters(labels)
                    self.echo('')
                    weights = set_weights(sel_data, data)

                    results = eos.fit(
                        values, self.input.pressure, guess=guess,
                        sigmaV=weights[1], sigmaP=weights[0], order=order,
                        linear=linearized, fixed=fix
                        )

                    self.report_fit(
                        results, values, data_name, eos.__repr__(), order,
                        weights, fix, labels
                        )

                    if confirm('Update the parameters'):
                        guess = results.beta[:]
                        self.echo('')

                    change_data = not confirm('Would you like to continue '
                                              'fitting on the same data?')
                    self.echo('')
                    if change_data:
                        change_eos = True
                    else:
                        change_eos = confirm('Would you like to change '
                                             'the EoS?')
                        self.echo('')

            is_running = confirm("Would you like to fit other data?")

        self.echo('{0}{1:-^78}{0}'.format(
            '#', 'Interactive equation of state fitting finished'))
        self.completed = True    
        return

    def report_fit(self, results, values, value_name, eos_name, eos_order,
                   weights, fix, labels):
        """
        """
        implieds = {
            2: [1, 2],
            3: [2],
            4: []
            }

        if value_name == 'Volume':
            value_label = value_name + ' (' + self._vunit + '^3)'
        else:
            value_label = value_name + ' (' + self._vunit + ')'

        self.echo('')
        self.echo('{}{:-^78}{}'.format('#','','#'))
        self.echo('{}{: ^78}{}'.format(
            '#', 'EOS results for ' + value_name + ' fitting', '#'))
        self.echo('{}{:-^78}{}'.format('#','','#'))
        self.echo('')
        self.echo('{: >15}: {}'.format('EOS formulation', eos_name))
        self.echo('{: >15}: {}'.format('EOS order', eos_order))
        self.echo('')

        if isinstance(weights[0], np.ndarray):
            self.echo('Results weighted for pressure')
        else:
            self.echo('Results not weighted for pressure')
        if isinstance(weights[1], np.ndarray):
            self.echo('Results weighted for {}'.format(value_name))
        else:
            self.echo('Results not weighted for for {}'.format(value_name))
        self.echo('')

        self.echo('Fitting parameters:')
        for i in range(len(results.beta)):
            res_row = '  {0:6} = {1: 10.5f} +/- {2: 10.5f}'.format(
                labels[i], results.beta[i], results.sd_beta[i])
            if fix[i] == 0:
                res_row += ' (FIXED)'
            if i in implieds[eos_order]:
                res_row +=  ' [IMPLIED]'
            self.echo(res_row)
        self.echo('')

        self.echo('Covariance matrix:')
        for i in range(len(results.beta)):
            self.echo('  {0: 12.5e} {1: 12.5e} {2: 12.5e} {3: 12.5e}'.format(
                *results.cov_beta[i]))
        self.echo('')
        self.echo('Residual variance: {0:6.4f}'.format(results.res_var))
        self.echo('')
        self.echo('Exited ODR loop because:')
        for i in range(len(results.stopreason)):
            self.echo('  {}. {}'.format(i+1, results.stopreason[i]))
        self.echo('')
        self.echo('Calculated data:\n')
        self.echo('{: ^15} {: ^15} {: ^15} {: ^15} {: ^15}'.format(
            value_label, 'eps_'+value_name[0], 'Pressure (GPa)',
            'P_obs (GPa)', 'eps_P'))
        self.echo('{0:-^15} {0:-^15} {0:-^15} {0:-^15} {0:-^15}'.format(''))
        for i in range(len(values)):
            line = ''
            line += '{: ^15.5f} '.format(values[i])
            line += '{: ^ 15.5f} '.format(results.delta[i])
            line += '{: ^15.5f} '.format(self.input.pressure[i])
            line += '{: ^15.5f} '.format(results.y[i])
            line += '{: ^ 15.5f}'.format(results.eps[i])
            self.echo(line)
        self.echo('')
        self.echo('{}{:-^78}{}'.format('#','','#'))
        self.echo('')
        return
