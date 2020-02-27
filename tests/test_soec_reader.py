# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" Test cases for the SOEC input file reader. """

import os
import unittest
import tempfile
import numpy as np

from quantas.IO.soec_reader import SOECInputFileReader


class SOECReaderTest(unittest.TestCase):

    valid_soec_input = '''Periclase (MgO) SOEC
   300.7440    101.7180    101.7180      0.0000      0.0000      0.0000 
   101.7180    300.7440    101.7180      0.0000      0.0000      0.0000 
   101.7180    101.7180    300.7440      0.0000      0.0000      0.0000 
     0.0000      0.0000      0.0000    169.0550      0.0000      0.0000 
     0.0000      0.0000      0.0000      0.0000    169.0550      0.0000 
     0.0000      0.0000      0.0000      0.0000      0.0000    169.0550 
   3513.0
'''
    invalid_soec_input = '''
Periclase (MgO) SOEC
   300.7440    101.7180    101.7180      0.0000      0.0000      0.0000  
   101.7180    101.7180    300.7440      0.0000      0.0000      0.0000 
     0.0000      0.0000      0.0000    169.0550      0.0000      0.0000 
     0.0000      0.0000      0.0000      0.0000    169.0550      0.0000 
     0.0000      0.0000      0.0000      0.0000      0.0000    169.0550 
   3513.0
'''
    def test_reading(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
            filename = tmp.name
            basename = os.path.split(filename)[1]
            tmp.write(self.valid_soec_input)
            tmp.flush()
            reader = SOECInputFileReader(filename)

        self.assertTrue(reader.completed, 'Reading process should be completed')
        self.assertEqual(reader.jobname, 'Periclase (MgO) SOEC',
                         "Jobname should be 'Periclase (MgO) SOEC'")
        self.assertEqual(reader.density, 3513.0, 'Should be 3513.0')
        self.assertEqual(type(reader.stiffness), type(np.array([0])),
                         'Stiffness matrix should be numpy.ndarray')
        self.assertEqual(reader.stiffness.shape, (6, 6),
                         'Stiffness matrix should have shape (6, 6)')
            
        tmp.close()
        os.unlink(tmp.name)

    def test_reading_with_error(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
            filename = tmp.name
            basename = os.path.split(filename)[1]
            tmp.write(self.invalid_soec_input)
            tmp.flush()
            with self.assertRaises(IndexError):
                reader = SOECInputFileReader(filename)
                self.assertFalse(
                    reader.completed,
                    'Reading process should be not completed')
                self.assertEqual(
                    reader.error,
                    'The SOEC matrix does not have the expected (6, 6) '
                    'shape.\nPlease, check the input file and retry.',
                    "Error should be related to the shape of the matrix")
            
        tmp.close()
        os.unlink(tmp.name)


if __name__ == '__main__':
    unittest.main()
