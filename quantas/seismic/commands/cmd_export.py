# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" The Seismic export command. """

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

from ..utils.export import QuantasSeismicExport


@click.command('seismic')
@click.argument('filename', type=click.Path(exists=True))
def seismic_export(filename):
    """ Export results obtained from the Seismic analysis.

    This command requires an HDF5 file (FILENAME) that contains the results
    from second-order elastic moduli analysis, which will be exported.
    """
    echo_highlight(quantas_title())

    if not os.path.isfile(filename):
        echo_error(quantas_error(), bold=True)
        echo_error("'{}' is not a file".format(filename))
        return

    export = QuantasSeismicExport()
    export.load(filename)

    if not isinstance(export.error, type(None)):
        echo_error(quantas_error(), bold=True)
        echo_error(export.error)
        return

    export.run()

    echo_highlight(biblio_header())
    echo_highlight(quantas_citation())
    echo_highlight(biblio_footer())

    echo_highlight(quantas_finish())
    return
