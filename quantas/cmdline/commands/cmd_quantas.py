# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

import os
import sys
import click

from .cmd_inpgen import inpgen
from .cmd_export import export
from quantas.harmonic.commands.cmd_ha import ha_calculation
from quantas.qha.commands.cmd_qha import qha_calculation
from quantas.soec.commands.cmd_soec import soec_calculation
from quantas.eosfit.commands.cmd_eosfit import eos_calculation

from ..utils.general import MostSimilarCommandGroup
from ..utils.general import print_version
from quantas.cmdline.utils.messages import quantas_title
from quantas.cmdline.utils.messages import quantas_error
from quantas.cmdline.utils.messages import echo_highlight
from quantas.cmdline.utils.messages import echo_error

@click.group(cls=MostSimilarCommandGroup,
             context_settings={'help_option_names': ['-h', '--help']})
@click.option('-v', '--version', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True,
              help='Show the software version and exit.')
@click.pass_context
def cli(ctx):
    """\b
________                       __
\_____  \  __ _______    _____/  |______    ______
 /  / \  \|  |  \__  \  /    \   __\__  \  /  ___/
/   \_/.  \  |  // __ \|   |  \  |  / __ \_\___ \\
\_____\ \_/____/(____  /___|  /__| (____  /____  >
       \__>          \/     \/          \/     \/ 
    """
    return


cli.add_command(inpgen)
cli.add_command(export)
cli.add_command(ha_calculation)
cli.add_command(qha_calculation)
cli.add_command(soec_calculation)
cli.add_command(eos_calculation)
