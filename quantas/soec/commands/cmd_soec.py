# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" The SOEC command. """

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
from quantas.citations import soec_citation

from ..soec import SOECCalculator

@click.command('soec')
@click.argument('filename', type=click.Path(exists=True))
@click.option('-o', '--outfile', help='Output file where data will be stored, '
              'without extension.', default=None, metavar='out_file',
              show_default='Input file base name')
@click.option('--punit', help='Measurement unit for pressure values.',
              type=click.Choice(['GPa'], case_sensitive=True),
              default='GPa', show_default='GPa')
@click.option('--polar', help='Calculate elastic properties on (xy), (xz) and'
              ' (yz) planes.', is_flag=True, default=False)
@click.option('-p', '--plot', help='Create polar plots of the results.',
              is_flag=True, default=False)
@click.option('--dpi', help='Resolution (DPI) of the output figures.',
              default=80, show_default=80)
@click.option('-q', '--quiet', 'silent', is_flag=True, default=False,
              help='Output will not be printed on screen.',
              )
@click.option('-d', '--debug', is_flag=True, help='Activate debug option.',
              default=False)
def soec_calculation(filename, outfile, punit, polar, plot, dpi, silent, debug):
    """ Second-order elastic moduli analisys.

    This command requires a file (FILENAME) that will be read to provide the
    input data for the calculations.
    """
    
    outbase = os.path.splitext(filename)[0]
    if isinstance(outfile, type(None)):
        outfile = outbase + '_SOEC.hdf5'
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
        'pressure_unit': punit,
        'debug': debug,
        'silent': silent,
        'polar': polar,
        'plotting': plot,
        'dpi': dpi,
        'logfile': logfile,
        'outfig': os.path.splitext(outfile)[0],
        }

    echo(print_settings(runtime_settings), logfile, silent=silent)

    calculator = SOECCalculator(runtime_settings)

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
    echo_highlight(soec_citation(), logfile, silent=silent)
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
    text = '\nCalculator: Second-order elastic moduli analysis\n'

    text += '\nMeasurement units\n'
    text += '-------------------------------------\n'
    text += ' - {:12} {}\n'.format('pressure:', settings['pressure_unit'])
    return text
