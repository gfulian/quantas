# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" The HA input file generator command. """

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

from quantas.IO.qha_writer import QHAInputCreator

from quantas.citations import biblio_header
from quantas.citations import biblio_footer
from quantas.citations import quantas_citation

interfaces = [
    'crystal',
    'crystal-qha',
    'phonopy'
    ]

@click.command('ha')
@click.argument('filename', type=click.Path(exists=True))
@click.option('-o', '--outfile', help='Output file where data will be stored, '
              'without extension.', default='quantas_ha', metavar='out_file',
              show_default='quantas_ha')
@click.option('-l', '--list', 'is_list', is_flag=True, default=False,
              show_default=False, help='Input files provided as a file list')
@click.option('-r', '--ref', type=click.INT, default=0, show_default=0,
              help='Reference file for (Q)HA input')
@click.option('-i', '--interface', help='Interface for ab initio codes.',
              type=click.Choice(interfaces, case_sensitive=True),
              default='crystal', show_default='crystal')
def ha_inpgen(filename, outfile, is_list, ref, interface):
    """ Input generator for (Quasi-)Harmonic Approximation calculations.

    This command requires a file (FILENAME) that will be read to provide the
    input data for the input generation.
    """
    
    outbase = os.path.splitext(filename)[0]
    if isinstance(outfile, type(None)):
        outfile = outbase + '_ha_input.yaml'
    else:
        outfile += '.yaml'

    echo_highlight(quantas_title())

    if not os.path.isfile(filename):
        echo_error(quantas_error(), bold=True)
        echo_error("'{}' is not a file".format(filename))
        return

    generator = QHAInputCreator(interface)

    try:
        completed, error = generator.read(filename, is_list, ref, sym)
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
    if interface != 'crystal-qha':
        if ref > len(generator.phondata):
            echo_error(quantas_error(), bold=True)
            echo_error('Invalid reference provided')
            return
        generator.write(outfile, ref)
    else:
        generator.write(outfile)

    echo_highlight(biblio_header())
    echo_highlight(quantas_citation())
    echo_highlight(biblio_footer())
    echo_highlight(quantas_finish())
    return
