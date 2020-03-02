# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" The SOEC input file generator command. """

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
from quantas.cmdline.utils.messages import prompt

from quantas.citations import biblio_header
from quantas.citations import biblio_footer
from quantas.citations import quantas_citation

from quantas.IO.soec_writer import SOECInputCreator

interfaces = [
    'crystal',
    'vasp'
    ]

@click.command('soec')
@click.argument('filename', type=click.Path(exists=True))
@click.option('-o', '--outfile', help='Output file where data will be stored, '
              'without extension.', default=None, metavar='out_file',
              show_default='Input file base name')
@click.option('-i', '--interface', help='Interface for ab initio codes.',
              type=click.Choice(interfaces, case_sensitive=True),
              default='crystal', show_default='crystal')
def soec_inpgen(filename, outfile, interface):
    """ Input generator for second-order elastic moduli analisys.

    This command requires a file (FILENAME) that will be read to provide the
    input data for the calculations.
    """
    outbase = os.path.splitext(filename)[0]
    if isinstance(outfile, type(None)):
        outfile = outbase + '_soec_input.dat'
    else:
        outfile += '.dat'

    echo_highlight(quantas_title())

    if not os.path.isfile(filename):
        echo_error(quantas_error(), bold=True)
        echo_error('{} is not a file'.format(filename))
        return

    generator = SOECInputCreator(interface)

    try:
        completed, error = generator.read(filename)
    except KeyError:
        echo_error(quantas_error(), bold=True)
        echo_error(
            "File '{}' does not appear as a valid input file".format(filename))
        return
    except UnicodeDecodeError:
        echo_error(quantas_error(), bold=True)
        echo_error(
            "File '{}' is in binary format".format(filename))
        return

    if not completed:
        echo_error(quantas_error(), bold=True)
        echo_error(error)
        return

    echo("Preparing the input file for Quantas: '{}'".format(outfile))
    msg = '\nPlease, enter a short description for the input file: '
    jobname = prompt(msg, default='Unknown', show_default=False)
    generator.write(outfile, jobname)

    echo_highlight(biblio_header())
    echo_highlight(quantas_citation())
    echo_highlight(biblio_footer())
    echo_highlight(quantas_finish())
    return
