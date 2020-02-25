# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" The HA command. """

import os
import sys
import click
import numpy as np

from quantas.cmdline.utils.messages import quantas_title
from quantas.cmdline.utils.messages import quantas_warning
from quantas.cmdline.utils.messages import quantas_error
from quantas.cmdline.utils.messages import quantas_finish
from quantas.cmdline.utils.messages import init_logfile
from quantas.cmdline.utils.messages import echo
from quantas.cmdline.utils.messages import echo_error
from quantas.cmdline.utils.messages import echo_highlight
from quantas.cmdline.utils.messages import confirm

from quantas.citations import biblio_header
from quantas.citations import biblio_footer
from quantas.citations import quantas_citation
from quantas.citations import qha_citation

from ..harmonic import HACalculator

@click.command('ha')
@click.argument('filename', type=click.Path(exists=True))
@click.option('-o', '--outfile', help='Output file where data will be stored, '
              'without extension.', default=None, metavar='out_file',
              show_default='Input file base name')
@click.option('-T', '--temperature', nargs=3, type=float,
              metavar='min max step',
              help='Temperature range provided as a tuple.',
              default=(298.15, 298.15, 1.), show_default='(298.15, 298.15, 1.)')
@click.option('--eunit', help='Measurement unit for energy values.',
              type=click.Choice(['Ha', 'eV', 'Ry'], case_sensitive=True),
              default='Ha', show_default='Ha')
@click.option('--vunit', help='Measurement unit for volume values.',
              type=click.Choice(['A', 'bohr'], case_sensitive=True),
              default='A', show_default='A')
@click.option('--funit', help='Measurement unit for phonon frequency values.',
              type=click.Choice(['cm-1', 'THz', 'Hz'], case_sensitive=True),
              default='cm-1', show_default='cm^-1')
@click.option('--tunit', help='Measurement unit for temperature values.',
              type=click.Choice(['K', 'C'], case_sensitive=True),
              default='K', show_default='K')
@click.option('-q', '--quiet', 'silent', is_flag=True, default=False,
              help='Output will not be printed on screen.',
              )
@click.option('-d', '--debug', is_flag=True, help='Activate debug option.',
              default=False)
def ha_calculation(filename, outfile, temperature, eunit, vunit, funit, tunit, silent,
                   debug):
    """ Harmonic Approximation calculation.

    This command requires a file (FILENAME) that will be read to provide the
    input data for the calculations.
    """
    
    outbase = os.path.splitext(filename)[0]
    if isinstance(outfile, type(None)):
        outfile = outbase + '_HA.hdf5'
    else:
        outfile += '.hdf5'
    logfile = os.path.splitext(outfile)[0] + '.txt'
    init_logfile(logfile)

    echo_highlight(quantas_title(), logfile, silent=silent)

    if not os.path.isfile(filename):
        echo_error(quantas_error(), logfile, bold=True)
        echo_error('{} is not a file'.format(filename), logfile)
        return

    runtime_settings = {
        'trange': temperature,
        'energy_unit': eunit,
        'lenght_unit': vunit,
        'frequency_unit': 'cm^-1' if funit == 'cm-1' else funit,
        'temperature_unit': tunit,
        'debug': debug,
        'silent': silent,
        'logfile': logfile,
        }

    echo(print_settings(runtime_settings), logfile, silent=silent)

    calculator = HACalculator(runtime_settings)

    error = calculator.read_input(filename)

    if not isinstance(error, type(None)):
        echo_error(quantas_error(), bold=True)
        echo_error(error, logfile)
        return

    calculator.report_input_data()
    calculator.run()

    if not calculator.completed:
        echo_error(quantas_error(), logfile, bold=True)
        echo_error(calculator.error, logfile)
        return

    if os.path.exists(outfile):
        ans = confirm('Output file {} exists. '
                      'Would you like to overwrite it?'.format(outfile))
        if ans:
            calculator.export_hdf5(outfile)
        else:
            echo('Results not saved', logfile)
    else:
        calculator.export_hdf5(outfile)

    echo_highlight(biblio_header(), logfile, silent=silent)
    echo_highlight(quantas_citation(), logfile, silent=silent)
    echo_highlight(qha_citation(), logfile, silent=silent)
    echo_highlight(biblio_footer(), logfile, silent=silent)

    echo_highlight(quantas_finish(), logfile, silent=silent)
    return


def print_settings(settings):
    """
    This function returns the harmonic approximation settings .

    Returns
    -------

    text: str
        Pretty-printed settings for the current Quantas run.

    """
    text = '\nCalculator: Harmonic Approximation\n'

    text += '\nTemperature settings\n'
    text += '-------------------------------------\n'
    text += ' - {:5} {: 8.2f} {}\n'.format(
        'min:', settings['trange'][0], settings['temperature_unit'])
    text += ' - {:5} {: 8.2f} {}\n'.format(
        'max:', settings['trange'][1], settings['temperature_unit'])
    text += ' - {:5} {: 8.2f} {}\n'.format(
        'step:', settings['trange'][2], settings['temperature_unit'])

    text += '\nMeasurement units\n'
    text += '-------------------------------------\n'
    text += ' - {:12} {}\n'.format('energy:', settings['energy_unit'])
    text += ' - {:12} {}\n'.format('lenght:', settings['lenght_unit'])
    text += ' - {:12} {}\n'.format('frequency:', settings['frequency_unit'])
    text += ' - {:12} {}\n'.format('temperature:', settings['temperature_unit'])
    return text
    
