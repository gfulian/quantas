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
from ..utils.general import print_version

from quantas.harmonic.commands.cmd_export import ha_export
from quantas.qha.commands.cmd_export import qha_export
from quantas.soec.commands.cmd_export import soec_export

@click.group(cls=MostSimilarCommandGroup,
             context_settings={'help_option_names': ['-h', '--help']})
def export():
    """ Export results from binary format (HDF5) to text format. """
    return

export.add_command(ha_export)
export.add_command(qha_export)
export.add_command(soec_export)
