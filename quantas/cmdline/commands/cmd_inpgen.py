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

from ..utils.general import MostSimilarCommandGroup

from quantas.harmonic.commands.cmd_inpgen import ha_inpgen
from quantas.soec.commands.cmd_inpgen import soec_inpgen


@click.group(cls=MostSimilarCommandGroup,
             context_settings={'help_option_names': ['-h', '--help']})
def inpgen():
    """ Generate inputs file for Quantas from input files."""
    return

inpgen.add_command(ha_inpgen)
inpgen.add_command(soec_inpgen)
