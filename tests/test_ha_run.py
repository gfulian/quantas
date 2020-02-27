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

from quantas.IO.qha_reader import QHAInputFileReader
from quantas.harmonic.harmonic import HACalculator

class HarmonicCalculationTest(unittest.TestCase):

    valid_qha_input = '''job: MgO Gamma-point frequencies
natom:   2
supercell:
- [      1,      0,      0 ]
- [      0,      1,      0 ]
- [      0,      0,      1 ]
qpoints: 1
volume: [      18.89686185 ]
energy: [  -2.754625427009E+02 ]
phonon:
- q-position: [    0.0000000,    0.0000000,    0.0000000 ]
  weight: 1    
  band:
  - # 1
    frequency: [     0.0000000000 ]
  - # 2
    frequency: [     0.0000000000 ]
  - # 3
    frequency: [     0.0000000000 ]
  - # 4
    frequency: [   381.0964000000 ]
  - # 5
    frequency: [   381.0964000000 ]
  - # 6
    frequency: [   381.0964000000 ]
'''
    
    def test_run(self):
        
        runtime_settings = {
            'trange': (298.15, 298.15, 1.),
            'energy_unit': 'Ha',
            'lenght_unit': 'A',
            'frequency_unit': 'cm^-1',
            'temperature_unit': 'K',
            'debug': False,
            'silent': True,
            'logfile': None,
            }
        calculator = HACalculator(runtime_settings)

        reference = {
            'U0': -275.462542700920,
            'Uzp': 0.002604604211,
            'Uth': 0.000984626351,
            'S': 0.004947231852,
            'Cv': 0.007221334436,
            'TS': 0.001475017177
            }
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
            filename = tmp.name
            basename = os.path.split(filename)[1]
            tmp.write(self.valid_qha_input)
            tmp.flush()
            calculator.read_input(filename)
            tmp.close()
            os.unlink(tmp.name)

        calculator.run()
        res = calculator.results

        self.assertTrue(calculator.completed, 'Run should be completed')
        for key, value in res.items():
            if key in reference:
                self.assertTrue(np.isclose(value, reference[key]),
                                '{} should be equal to {: 12.8f}'.format(
                                    key, reference[key]))

        return
        
        
        
        
