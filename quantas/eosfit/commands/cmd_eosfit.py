# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" The Equation of state fitting command. """

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

from quantas.citations import biblio_header
from quantas.citations import biblio_footer
from quantas.citations import quantas_citation
from quantas.citations import eos_citation

from ..eosfit import EoSFitCalculator

@click.command('eosfit')
@click.argument('filename', type=click.Path(exists=True))
@click.option('-o', '--outfile', help='Output file for EoS log, '
              'without extension.', default=None, metavar='out_file',
              show_default='Input file base name')
@click.option('--punit', help='Measurement unit for pressure values.',
              type=click.Choice(['GPa'], case_sensitive=True),
              default='GPa', show_default='GPa')
@click.option('--vunit', help='Measurement unit for volume values.',
              type=click.Choice(['A'], case_sensitive=True),
              default='A', show_default='A')
def eos_calculation(filename, outfile, punit, vunit):
    """ Equation of state (EoS) fitting.

    This command requires a file (FILENAME) that will be read to provide the
    input data for the calculations.
    """
    
    outbase = os.path.splitext(filename)[0]
    if isinstance(outfile, type(None)):
        logfile = outbase + '.log'
    else:
        logfile += '.log'
    init_logfile(logfile)

    echo_highlight(quantas_title(), logfile)

    if not os.path.isfile(filename):
        echo_error(quantas_error(), logfile, bold=True)
        echo_error('{} is not a file'.format(filename), logfile)
        return

    runtime_settings = {
        'pressure_unit': punit,
        'lenght_unit': vunit,
        'logfile': logfile,
        }

    echo(print_settings(runtime_settings), logfile)

    calculator = EoSFitCalculator(runtime_settings)

    error = calculator.read_input(filename)

    if not isinstance(error, type(None)):
        echo_error(quantas_error(), bold=True)
        echo_error(error, logfile)
        return

    calculator.report_input_data()
    calculator.run()

    echo_highlight(biblio_header(), logfile)
    echo_highlight(quantas_citation(), logfile)
    echo_highlight(eos_citation(), logfile)
    echo_highlight(biblio_footer(), logfile)

    echo_highlight(quantas_finish(), logfile)
    return


def print_settings(settings):
    """
    This function returns the harmonic approximation settings .

    Returns
    -------

    text: str
        Pretty-printed settings for the current Quantas run.

    """
    text = '\nCalculator: Equation of state (EoS) fitting\n'

    text += '\nMeasurement units\n'
    text += '-------------------------------------\n'
    text += ' - {:12} {}\n'.format('pressure:', settings['pressure_unit'])
    text += ' - {:12} {}\n'.format('lenght:', settings['lenght_unit'])
    return text
