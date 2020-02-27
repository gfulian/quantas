# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" Test cases for the harmonic approximation calculation run """

import os
import unittest
import tempfile
import numpy as np

from quantas.IO.soec_reader import SOECInputFileReader
from quantas.soec.soec import SOECCalculator

class SOECCalculationTest(unittest.TestCase):

    valid_soec_input = '''Periclase (MgO) SOEC
   300.7440    101.7180    101.7180      0.0000      0.0000      0.0000 
   101.7180    300.7440    101.7180      0.0000      0.0000      0.0000 
   101.7180    101.7180    300.7440      0.0000      0.0000      0.0000 
     0.0000      0.0000      0.0000    169.0550      0.0000      0.0000 
     0.0000      0.0000      0.0000      0.0000    169.0550      0.0000 
     0.0000      0.0000      0.0000      0.0000      0.0000    169.0550 
   3513.0
'''
    
    def test_run(self):
        
        runtime_settings = {
            'pressure_unit': 'GPa',
            'debug': False,
            'silent': True,
            'plotting': False,
            'dpi': 80,
            'logfile': None,
            'outfig': None,
            }
        calculator = SOECCalculator(runtime_settings)

        reference = {
            'avgs': np.array(
                [[168.06, 330.99225746655424, 141.2382, 0.17175189667722418],
                 [168.06, 314.06592076074, 132.1228040490129, 0.18853790237540158],
                 [168.06, 322.5892778210904, 136.68050202450644, 0.1800852098247745]]
                ),
            'eigv': np.array([169.055,
                              169.055,
                              169.055,
                              199.026,
                              199.026,
                              504.18000000000006]
                             )
            }
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
            filename = tmp.name
            basename = os.path.split(filename)[1]
            tmp.write(self.valid_soec_input)
            tmp.flush()
            calculator.read_input(filename)
            tmp.close()
            os.unlink(tmp.name)

        calculator.run()
        res = calculator.results

        self.assertTrue(calculator.completed, 'Run should be completed')
        for key, values in res.items():
            if key in reference:
                self.assertTrue(np.isclose(values, reference[key]).all(),
                                "Mismatch for the '{}' results")
        return
        
        
        
        
