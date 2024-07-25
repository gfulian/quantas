# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" The Seismic command. """

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
from quantas.citations import christoffel_citation

from ..seismic import SeismicCalculator

@click.command('seismic')
@click.argument('filename', type=click.Path(exists=True))
@click.option('-o', '--outfile', help='Output file where data will be stored, '
              'without extension.', default=None, metavar='out_file',
              show_default='Input file base name')
@click.option('--ntheta', help='Number of polar angular points.',
              default=180, show_default=180)
@click.option('--nphi', help='Number of azimuthal angular points.',
              default=4*180, show_default='4*ntheta')
@click.option('--punit', help='Measurement unit for pressure values.',
              type=click.Choice(['GPa'], case_sensitive=True),
              default='GPa', show_default='GPa')
@click.option('-p', '--plot', help='Create 3D and 2D plots of the results.',
              is_flag=True, default=False)
@click.option('--no-2D', help='Do not create 2D plots of the results.',
              is_flag=True, default=False)
@click.option('--no-3D', help='Do not create 3D plots of the results.',
              is_flag=True, default=False)
@click.option('--proj', help='Type of 2D projection', default='eqar',
              type=click.Choice(['eqar', 'stereo'], case_sensitive=False),
              show_default='eqar')
@click.option('--dpi', help='Resolution (DPI) of the output figures.',
              default=80, show_default=80)
@click.option('--only-plot', help='Plot previous results.',
              is_flag=True, default=False)
@click.option('-q', '--quiet', 'silent', is_flag=True, default=False,
              help='Output will not be printed on screen.',
              )
@click.option('-d', '--debug', is_flag=True, help='Activate debug option.',
              default=False)
def seismic_calculation(filename, outfile, ntheta, nphi, punit, plot, no_2d,
                        no_3d, proj, dpi, only_plot, silent, debug):
    """ Calculation of 3D seismic wave velocities from elastic moduli, solving
    the Christoffel's equation.

    This command requires a file (FILENAME) that will be read to provide the
    input data for the calculations.
    """
    
    outbase = os.path.splitext(filename)[0]
    if isinstance(outfile, type(None)):
        outfile = outbase + '_SEISMIC.hdf5'
    else:
        outfile += '.hdf5'
    logfile = os.path.splitext(outfile)[0] + '.txt'
    init_logfile(logfile)

    echo_highlight(quantas_title(), logfile, silent=silent)

    if not os.path.isfile(filename):
        echo_error(quantas_error(), logfile, bold=True)
        echo_error('{} is not a file'.format(filename), logfile)
        return

    if ntheta != 180:
        if nphi == 720:
            nphi = 4 * ntheta

    runtime_settings = {
        'ntheta': int(ntheta),
        'nphi': int(nphi),
        'pressure_unit': punit,
        'debug': debug,
        'silent': silent,
        'plotting': plot,
        '2D': not no_2d,
        '3D': not no_3d,
        'projection': proj,
        'dpi': dpi,
        'logfile': logfile,
        'outfig': os.path.splitext(outfile)[0],
        }

    echo(print_settings(runtime_settings), logfile, silent=silent)

    calculator = SeismicCalculator(runtime_settings)

    if only_plot:
        error = calculator.plot(filename)
        if not isinstance(error, type(None)):
            echo_error(quantas_error(), bold=True)
            echo_error(error, logfile)
            return

    else:
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

        if calculator.plotting:
            calculator.plot()
            
    echo_highlight(biblio_header(), logfile, silent=silent)
    echo_highlight(quantas_citation(), logfile, silent=silent)
    echo_highlight(christoffel_citation(), logfile, silent=silent)
    echo_highlight(biblio_footer(), logfile, silent=silent)

    echo_highlight(quantas_finish(), logfile, silent=silent)
    return


def print_settings(settings):
    """
    This function returns the settings used during the current seismic wave
    calculation.

    Returns
    -------

    text: str
        Pretty-printed settings for the current Quantas run.

    """
    text = "\nCalculator: wave velocities from Christoffel's equation\n"

    text += '\nNumber of angular points\n'
    text += '-------------------------------------\n'
    text += ' - {:12} {}\n'.format('ntheta:', settings['ntheta'])
    text += ' - {:12} {}\n'.format('nphi:', settings['nphi'])

    text += '\nMeasurement units\n'
    text += '-------------------------------------\n'
    text += ' - {:12} {}\n'.format('pressure:', settings['pressure_unit'])

    text += '\nPlotting\n'
    text += '-------------------------------------\n'
    text += ' - {:12} {}\n'.format('requested:', settings['plotting'])
    if settings['plotting']:
        text += ' - {:12} {}\n'.format('dpi:', settings['dpi'])
        text += ' - {:12} {}\n'.format('3D plots:', settings['3D'])
        text += ' - {:12} {}\n'.format('2D plots:', settings['2D'])
        if settings['2D']:
            text += ' - {:12} {}\n'.format(
                'projection:',
                'Lambert equal area' if settings['projection']=='eqar'
                else 'stereographic')
    return text
