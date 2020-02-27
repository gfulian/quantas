# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" Test cases for the input file readers. """

import os
import unittest
import tempfile
import numpy as np

from quantas.IO.qha_reader import QHAInputFileReader

class QuasiHarmonicReaderTest(unittest.TestCase):

    valid_qha_input = '''
job: MgO Gamma-point frequencies
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
    invalid_qha_input = '''
job: MgO Gamma-point frequencies
natom:   2      
supercell:
- [      1,      0,      0 ]
- [      0,      1,      0 ]
- [      0,      0,      1 ]
qpoints: 1
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
    def test_reading(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
            filename = tmp.name
            basename = os.path.split(filename)[1]
            tmp.write(self.valid_qha_input)
            tmp.flush()
            reader = QHAInputFileReader(filename)

        self.assertTrue(reader.completed, 'Reading process should be completed')
        self.assertEqual(reader.jobname, 'MgO Gamma-point frequencies',
                         "Jobname should be 'MgO Gamma-point frequencies'")
        self.assertEqual(reader.kpoints, 1, 'Should be 1')
        self.assertEqual(reader.qpoints, 1, 'Should be 1')
        self.assertEqual(reader.volume, np.array([18.89686185]),
                         'Volume should be [18.89686185]')
        self.assertEqual(type(reader.volume), type(np.array([0])),
                         'Volume should be numpy.ndarray')
        self.assertEqual(reader.energy, np.array([-2.754625427009E+02]),
                         'Energy should be [-2.754625427009E+02]')
        self.assertEqual(type(reader.energy), type(np.array([0])),
                         'Energy should be numpy.ndarray')
        self.assertEqual(reader.frequencies.shape, (1, 6, 1),
                         'Frequency should have shape (1, 6, 1)')
            
        tmp.close()
        os.unlink(tmp.name)

    def test_reading_with_error(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
            filename = tmp.name
            basename = os.path.split(filename)[1]
            tmp.write(self.invalid_qha_input)
            tmp.flush()
            reader = QHAInputFileReader(filename)

        self.assertFalse(reader.completed,
                         'Reading process should be not completed')
        self.assertEqual(reader.error, 'No volume values in input file',
                         "Error should be 'No volume values in input file'")
            
        tmp.close()
        os.unlink(tmp.name)


if __name__ == '__main__':
    unittest.main()
