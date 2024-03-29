# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

import numpy as np
from quantas.interfaces.crystal import CrystalSOECReader
from quantas.interfaces.vasp import VASPSOECReader


class SOECInputCreator(object):
    """
    This class creates input files for the second-order elastic constants
    anslysis implemented in Quantas.
    """
    
    interface_filter = {
        'crystal': CrystalSOECReader,
        'vasp': VASPSOECReader,
        }

    soec_data = None
    
    interface = None

    def __init__(self, interface=None):
        """ Constructor method. """
        self.interface_flag = interface
        return

    def read(self, filename):
        """
        Read the ab initio output file to collect the elastic moduli data.

        Parameters
        ----------

        filename: str
            Path to the output file.
        """
        interface = self.interface_filter[self.interface_flag]
        data = interface(filename)

        if not data.completed:
            return False, data.error
        else:
            self.soec_data = data

        return True, None

    def write(self, outfile, jobname='Unknown'):
        """
        Crate the SOEC input for Quantas using the collected results from
        ab initio output file.

        Parameters
        ----------

        outfile: str
            Path of the Quantas SOEC input file.

        jobname: str, optional
            Job name.

        """
        data = []
        data.append(jobname)
        for i in range(6):
            soec_line = ''
            for j in range(6):
                    soec_line += '{: > 11.4f} '.format(
                        self.soec_data.stiffness[i, j])
            data.append(soec_line)
        data.append('{:9.1f}'.format(self.soec_data.density))

        with open(outfile, 'w') as f:
            f.write('\n'.join(data))
        return
