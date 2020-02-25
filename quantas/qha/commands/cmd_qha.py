# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" The QHA command. """

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

from ..qha import QHACalculator

eos_formulations = {
    'M': 'Murnaghan',
    'BM': '3rd-order Birch-Murnaghan',
    'PT': '3rd-order Poirier-Tarantola',
    'V': '3rd-order Vinet'
    }

eos_functions = {
    'M': 'murnaghan',
    'BM': 'birchmurnaghan',
    'PT': 'pouriertarantola',
    'V': 'vinet'
    }


@click.command('qha')
@click.argument('filename', type=click.Path(exists=True))
@click.option('-o', '--outfile', help='Output file where data will be stored, '
              'without extension.', default=None, metavar='out_file',
              show_default='Input file base name')
@click.option('-s', '--scheme', 'qha_scheme',
              help='QHA scheme, select between frequency (freq) or '
              'thermodynamic (td) interpolation.',
              type=click.Choice(['freq', 'td'], case_sensitive=False),
              default='td', show_default='td')
@click.option('-m', '--minimization', 'minimization',
              help='Volume minimization scheme, select between polynomial '
              'roots (poly) or equation of state fitting (eos).',
              type=click.Choice(['poly', 'eos'], case_sensitive=False),
              default='poly', show_default='poly')
@click.option('-T', '--temperature', nargs=3, type=float,
              metavar='min max step',
              help='Temperature range provided as a tuple.',
              default=(298.15, 298.15, 1.), show_default='(298.15, 298.15, 1.)')
@click.option('-P', '--pressure', nargs=3, type=float,
              metavar='min max step',
              help='Pressure range provided as a tuple.',
              default=(0., 0., 1.), show_default='(0., 0., 1.)')
@click.option('--fdeg', help='Set the degree of the polynomials used to fit'
              ' phonon frequency values.', type=click.IntRange(2, 5),
              default=3, show_default=3)
@click.option('--edeg', help='Set the degree of the polynomials used to fit'
              ' energy values.', type=click.IntRange(2, 5),
              default=3, show_default=3)
@click.option('--eos', help='Set the equation of state formulation.',
              type=click.Choice(['M', 'BM', 'PT', 'V'], case_sensitive=True),
              default='BM', show_default='BM')
@click.option('--eunit', help='Measurement unit for energy values.',
              type=click.Choice(['Ha', 'eV', 'Ry'], case_sensitive=True),
              default='Ha', show_default='Ha')
@click.option('--vunit', help='Measurement unit for volume values.',
              type=click.Choice(['A', 'bohr'], case_sensitive=True),
              default='A', show_default='A')
@click.option('--funit', help='Measurement unit for phonon frequency values.',
              type=click.Choice(['cm-1', 'THz', 'Hz'], case_sensitive=True),
              default='cm-1', show_default='cm-1')
@click.option('--tunit', help='Measurement unit for temperature values.',
              type=click.Choice(['K', 'C'], case_sensitive=True),
              default='K', show_default='K')
@click.option('--punit', help='Measurement unit for pressure values.',
              type=click.Choice(['GPa', 'kbar'], case_sensitive=True),
              default='GPa', show_default='GPa')
@click.option('-q', '--quiet', 'silent', is_flag=True, default=False,
              help='Output will not be printed on screen.',
              )
@click.option('-d', '--debug', is_flag=True, help='Activate debug option.',
              default=False)
def qha_calculation(filename, outfile, qha_scheme, minimization, temperature,
                    pressure, fdeg, edeg, eos, eunit, vunit, funit, tunit,
                    punit, silent, debug):
    """Quasi-Harmonic Approximation calculation.

    This command requires a file (FILENAME) that will be read to provide the
    input data for the calculations.
    """
    outbase = os.path.splitext(filename)[0]
    if isinstance(outfile, type(None)):
        outfile = outbase + '_QHA.hdf5'
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
        'qha_scheme': qha_scheme,
        'minimization': minimization,
        'edeg': edeg,
        'fdeg': fdeg,
        'eos': eos,
        'eos_function': eos_functions[eos],
        'trange': temperature,
        'prange': pressure,
        'energy_unit': eunit,
        'lenght_unit': vunit,
        'frequency_unit': 'cm^-1' if funit == 'cm-1' else funit,
        'temperature_unit': tunit,
        'pressure_unit': punit,
        'debug': debug,
        'silent': silent,
        'logfile': logfile,
        }

    echo(print_settings(runtime_settings), logfile, silent=silent)

    calculator = QHACalculator(runtime_settings)

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

    calculator.report_results()

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
    text = '\nCalculator: Quasi-Harmonic Approximation\n'

    text += '\nQuasi-Harmonic Approximation approach\n'
    text += '-------------------------------------\n'
    text += ' - scheme: '
    if settings['qha_scheme'] == 'freq':
        text += 'frequency interpolation\n'
    elif settings['qha_scheme'] == 'td':
        text += 'thermodynamics interpolation\n'
    text += ' - volume minimization: '
    if settings['minimization'] == 'eos':
        text += 'Equation of State (EoS)\n'
        text += ' - EoS formulation: '
        text += '{}\n'.format(eos_formulations[settings['eos']])
    if settings['minimization'] == 'poly':
        text += 'polynomial functions\n'

    text += '\nPolynomial fitting settings (degree)\n'
    text += '-------------------------------------\n'
    text += ' - {:10} {}\n'.format('energy:', settings['edeg'])
    text += ' - {:10} {}\n'.format('frequency:', settings['fdeg'])

    text += '\nTemperature settings\n'
    text += '-------------------------------------\n'
    text += ' - {:5} {: 8.2f} {}\n'.format(
        'min:', settings['trange'][0], settings['temperature_unit'])
    text += ' - {:5} {: 8.2f} {}\n'.format(
        'max:', settings['trange'][1], settings['temperature_unit'])
    text += ' - {:5} {: 8.2f} {}\n'.format(
        'step:', settings['trange'][2], settings['temperature_unit'])

    text += '\nPressure settings\n'
    text += '-------------------------------------\n'
    text += ' - {:5} {: 8.2f} {}\n'.format(
        'min:', settings['prange'][0], settings['pressure_unit'])
    text += ' - {:5} {: 8.2f} {}\n'.format(
        'max:', settings['prange'][1], settings['pressure_unit'])
    text += ' - {:5} {: 8.2f} {}\n'.format(
        'step:', settings['prange'][2], settings['pressure_unit'])

    text += '\nMeasurement units\n'
    text += '-------------------------------------\n'
    text += ' - {:12} {}\n'.format('energy:', settings['energy_unit'])
    text += ' - {:12} {}\n'.format('lenght:', settings['lenght_unit'])
    text += ' - {:12} {}\n'.format('frequency:', settings['frequency_unit'])
    text += ' - {:12} {}\n'.format('temperature:', settings['temperature_unit'])
    text += ' - {:12} {}\n'.format('pressure:', settings['pressure_unit'])
    return text
