# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

"""
Test case for the generation of input for (quasi-)harmonic approximation
calculations.
"""

import os
import unittest
import tempfile
import numpy as np

from quantas.IO.soec_reader import SOECInputFileReader
from quantas.IO.soec_writer import SOECInputCreator

class QuasiHarmonicInputGeneration(unittest.TestCase):

    crystal_output = '''
[ ... omissis ... ]

 GEOMETRY NOW FULLY CONSISTENT WITH THE GROUP

[ ... omissis ... ]

 GEOMETRY FOR WAVE FUNCTION - DIMENSIONALITY OF THE SYSTEM    3
 (NON PERIODIC DIRECTION: LATTICE PARAMETER FORMALLY SET TO 500)
 *******************************************************************************
 LATTICE PARAMETERS (ANGSTROMS AND DEGREES) - BOHR = 0.5291772083 ANGSTROM
 PRIMITIVE CELL - CENTRING CODE 5/0 VOLUME=    18.896862 - DENSITY  3.513 g/cm^3
         A              B              C           ALPHA      BETA       GAMMA 
     2.98975014     2.98975014     2.98975014    60.000000  60.000000  60.000000
 *******************************************************************************
 ATOMS IN THE ASYMMETRIC UNIT    2 - ATOMS IN THE UNIT CELL:    2
     ATOM                 X/A                 Y/B                 Z/C    
 *******************************************************************************
      1 T  12 MG    0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      2 T   8 O    -5.000000000000E-01 -5.000000000000E-01 -5.000000000000E-01

 TRANSFORMATION MATRIX PRIMITIVE-CRYSTALLOGRAPHIC CELL
 -1.0000  1.0000  1.0000  1.0000 -1.0000  1.0000  1.0000  1.0000 -1.0000

 *******************************************************************************
 CRYSTALLOGRAPHIC CELL (VOLUME=         75.58744740)
         A              B              C           ALPHA      BETA       GAMMA 
     4.22814520     4.22814520     4.22814520    90.000000  90.000000  90.000000

 COORDINATES IN THE CRYSTALLOGRAPHIC CELL
     ATOM                 X/A                 Y/B                 Z/C    
 *******************************************************************************
      1 T  12 MG    0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      2 T   8 O    -5.000000000000E-01 -5.000000000000E-01 -5.000000000000E-01

 T = ATOM BELONGING TO THE ASYMMETRIC UNIT

[ ... omissis ... ]

 ************************************************************************
 ************************************************************************
 * ELASTCON OPTION, W.F. PERGER, B. CIVALLERI, R. DOVESI                *
 * REFERENCES TO BE QUOTED WHEN USING THIS MODULE:                      *
 *                                                                      *
 * W.F. Perger, C. Criswell, B. Civalleri, R. Dovesi                    *
 * Ab-initio calculation of elastic constants of crystalline            *
 * systems with the CRYSTAL code                                        *
 * Comp. Phys. Comm. 180 (2009) 1753-1759                               *
 *                                                                      *
 * IT IS ASSUMED THAT THE INPUT IS FROM AN OPTIMIZED RUN                *

 CONVERGENCE CRITERION FOR ENERGY, TOLDEE =          1.0000E-08
 CONVERGENCE CRITERION FOR RMS GRADIENT, TOLDEG =    3.0000E-04
 CONVERGENCE CRITERION FOR RMS DISPLACEMENT, TOLDEX  6.2000E-04
 MAGNITUDE OF STEPSIZE FOR ADIMENSIONAL STRAINS =    1.0000E-02
 NUMBER OF POINTS (TOTAL) IN NUMERICAL 2ND DERIVATIVE=        3

[ ... omissis ... ]

 ************************************************************************
 ************************************************************************
                           FINAL RESULTS START
 ************************************************************************
 ************************************************************************

 THE CALCULATION HAS BEEN PERFORMED WITH   3 POINTS AND
 A STEP OF     0.01000
 THIS PERMITS TO PERFORM A FITTING UP TO SECOND   ORDER

 DATA FOR MAXIMUM NUMBER OF POINTS AND ORDER OF FIT

 SYMMETRIZED ELASTIC CONSTANTS FOR CUBIC        CASE, IN GPa

 |   300.744  101.718  101.718    0.000    0.000    0.000 |
 |            300.744  101.718    0.000    0.000    0.000 |
 |                     300.744    0.000    0.000    0.000 |
 |                              169.055    0.000    0.000 |
 |                                       169.055    0.000 |
 |                                                169.055 |


 ELASTIC MODULI (COMPLIANCE TENSOR), IN TPa^-1

 |    4.0108  -1.0137  -1.0137   0.0000   0.0000   0.0000 |
 |             4.0108  -1.0137   0.0000   0.0000   0.0000 |
 |                      4.0108   0.0000   0.0000   0.0000 |
 |                               5.9152   0.0000   0.0000 |
 |                                        5.9152   0.0000 |
 |                                                 5.9152 |


 BULK MODULUS K, SHEAR MODULUS G, YOUNG MODULUS E AND
 POISSON RATIO v (ALL IN GPa) ACCORDING TO VOIGT-REUSS-HILL

    K_V     G_V     K_R     G_R     K_H     G_H    E_H      v_H

  168.06  141.24  168.06  132.12  168.06  136.68  322.59   0.180

 SEISMIC VELOCITIES BY CHRISTOFFEL EQUATION (km/s)

       WAVE VECTOR                Vp        Vs1       Vs2

 [ 0.000  0.000  1.000]          9.252     6.937     6.937
 [ 0.000  1.000  0.000]          9.252     6.937     6.937
 [ 1.000  0.000  0.000]          9.252     6.937     6.937
 [ 1.000  1.000  0.000]         10.266     6.937     5.322
 [ 1.000  0.000  1.000]         10.266     6.937     5.322
 [ 0.000  1.000  1.000]         10.266     6.937     5.322
 [ 1.000  1.000  1.000]         10.583     5.910     5.910

[ ... omissis ... ]
'''

    vasp_outcar = '''
[ ... omissis ... ]

  ELASTIC MODULI  (kBar)
 Direction    XX          YY          ZZ          XY          YZ          ZX
 --------------------------------------------------------------------------------
 XX        2784.1164    778.9562    682.9594     -0.5893      0.0000      0.0000
 YY         782.2267   2789.2637    687.5498     -0.0871      0.0000      0.0000
 ZZ         706.2948    707.4161   3371.6084     -0.0002      0.0000      0.0000
 XY          -3.7580     -4.2092     -4.3732   1001.9607      0.0000      0.0000
 YZ          -0.5399     -0.4014     -0.3660     -0.7450    938.7658      0.4406
 ZX           0.0426     -0.1471     -0.1427     -0.0347     -0.9174    938.8021
 --------------------------------------------------------------------------------


  SYMMETRIZED ELASTIC MODULI (kBar)
 Direction    XX          YY          ZZ          XY          YZ          ZX
 --------------------------------------------------------------------------------
 XX        2784.1164    780.5915    694.6271     -2.1737      0.0000      0.0000
 YY         780.5915   2789.2637    697.4829     -2.1481      0.0000      0.0000
 ZZ         694.6271    697.4829   3371.6084     -2.1867      0.0000      0.0000
 XY          -2.1737     -2.1481     -2.1867   1001.9607      0.0000      0.0000
 YZ           0.0000      0.0000      0.0000      0.0000    938.7658     -0.2384
 ZX           0.0000      0.0000      0.0000      0.0000     -0.2384    938.8021
 --------------------------------------------------------------------------------

[ ... omissis ... ]

 ELASTIC MODULI CONTR FROM IONIC RELAXATION (kBar)
 Direction    XX          YY          ZZ          XY          YZ          ZX
 --------------------------------------------------------------------------------
 XX       -1306.2075   -283.4437    -53.5551     -2.6423      0.0000     -0.0000
 YY        -283.4437  -1300.9683    -38.5893     -3.6141      0.0000      0.0000
 ZZ         -53.5551    -38.5893  -1509.3462      4.6626     -0.0000     -0.0000
 XY          -2.6423     -3.6141      4.6626   -513.2031     -0.0000     -0.0000
 YZ           0.0000      0.0000      0.0000     -0.0000   -485.8520      0.7341
 ZX          -0.0000     -0.0000     -0.0000      0.0000      0.7341   -483.0727
 --------------------------------------------------------------------------------


 TOTAL ELASTIC MODULI (kBar)
 Direction    XX          YY          ZZ          XY          YZ          ZX
 --------------------------------------------------------------------------------
 XX        1477.9089    497.1477    641.0720     -4.8160      0.0000     -0.0000
 YY         497.1477   1488.2954    658.8937     -5.7622      0.0000      0.0000
 ZZ         641.0720    658.8937   1862.2622      2.4759     -0.0000     -0.0000
 XY          -4.8160     -5.7622      2.4759    488.7576     -0.0000     -0.0000
 YZ           0.0000      0.0000      0.0000     -0.0000    452.9138      0.4957
 ZX          -0.0000     -0.0000     -0.0000      0.0000      0.4957    455.7293
 --------------------------------------------------------------------------------

[ ... omissis ... ]
'''
    def test_input_from_CRYSTAL(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as cryout:
            filename = cryout.name
            basename = os.path.split(filename)[1]
            cryout.write(self.crystal_output)
            cryout.flush()
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
                outfile = tmp.name
                outbase = os.path.split(outfile)[1]
                generator = SOECInputCreator(interface='crystal')
                generator.read(filename)
                generator.write(outfile, 'MgO Elastic moduli analysis')

                reader = SOECInputFileReader(outfile)

                tmp.close()
                os.unlink(tmp.name)

            cryout.close()
            os.unlink(cryout.name)

        self.assertTrue(reader.completed, 'Reading process should be completed')
        self.assertEqual(reader.jobname, 'MgO Elastic moduli analysis',
                         "Jobname should be 'MgO Elastic moduli analysis'")
        self.assertEqual(reader.density, 3513.0, 'Should be 3513.0')
        self.assertEqual(type(reader.stiffness), type(np.array([0])),
                         'Stiffness matrix should be numpy.ndarray')
        self.assertEqual(reader.stiffness.shape, (6, 6),
                         'Stiffness matrix should have shape (6, 6)')
        return

    def test_input_from_VASP(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as outcar:
            filename = outcar.name
            basename = os.path.split(filename)[1]
            outcar.write(self.vasp_outcar)
            outcar.flush()
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
                outfile = tmp.name
                outbase = os.path.split(outfile)[1]
                generator = SOECInputCreator(interface='vasp')
                generator.read(filename)
                generator.write(outfile, 'OHAp Elastic moduli analysis')

                reader = SOECInputFileReader(outfile)

                tmp.close()
                os.unlink(tmp.name)

            outcar.close()
            os.unlink(outcar.name)

        self.assertTrue(reader.completed, 'Reading process should be completed')
        self.assertEqual(reader.jobname, 'OHAp Elastic moduli analysis',
                         "Jobname should be 'OHAp Elastic moduli analysis'")
        self.assertEqual(reader.density, 0.0, 'Should be 0.0')
        self.assertEqual(type(reader.stiffness), type(np.array([0])),
                         'Stiffness matrix should be numpy.ndarray')
        self.assertEqual(reader.stiffness.shape, (6, 6),
                         'Stiffness matrix should have shape (6, 6)')
        return
