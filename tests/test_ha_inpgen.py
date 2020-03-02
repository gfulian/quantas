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

from quantas.IO.qha_reader import QHAInputFileReader
from quantas.IO.qha_writer import QHAInputCreator

class QuasiHarmonicInputGeneration(unittest.TestCase):

    crystal_output = '''
[ ... omissis ... ]
 EEEEEEEEEE STARTING  DATE 26 02 2020 TIME 10:02:56.6
 MGO Bulk - B3LYP/POB-TZVP - Geometry optimization 3rd run                       

 CRYSTAL CALCULATION
 (INPUT ACCORDING TO THE INTERNATIONAL TABLES FOR X-RAY CRYSTALLOGRAPHY)
 CRYSTAL FAMILY                       :  CUBIC       
 CRYSTAL CLASS  (GROTH - 1921)        :  CUBIC HEXAKISOCTAHEDRAL              

 SPACE GROUP (CENTROSYMMETRIC)        :  F M 3 M         

 LATTICE PARAMETERS  (ANGSTROMS AND DEGREES) - CONVENTIONAL CELL
        A           B           C        ALPHA        BETA       GAMMA
     4.22815     4.22815     4.22815    90.00000    90.00000    90.00000


 NUMBER OF IRREDUCIBLE ATOMS IN THE CONVENTIONAL CELL:    2

 INPUT COORDINATES

 ATOM AT. N.              COORDINATES
   1  12     0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
   2   8     5.000000000000E-01  5.000000000000E-01  5.000000000000E-01

 *******************************************************************************

 << INFORMATION >>: FROM NOW ON, ALL COORDINATES REFER TO THE PRIMITIVE CELL

 *******************************************************************************

 LATTICE PARAMETERS  (ANGSTROMS AND DEGREES) - PRIMITIVE CELL
       A          B          C        ALPHA      BETA     GAMMA        VOLUME
    2.98975    2.98975    2.98975    60.00000  60.00000  60.00000     18.896862

 COORDINATES OF THE EQUIVALENT ATOMS (FRACTIONAL UNITS)

 N. ATOM EQUIV AT. N.          X                  Y                  Z

   1   1   1   12 MG    0.00000000000E+00  0.00000000000E+00  0.00000000000E+00

   2   2   1    8 O    -5.00000000000E-01 -5.00000000000E-01 -5.00000000000E-01

 NUMBER OF SYMMETRY OPERATORS         :   48
 *******************************************************************************
 * GEOMETRY EDITING - INPUT COORDINATES ARE GIVEN IN ANGSTROM         
 *******************************************************************************

 GEOMETRY NOW FULLY CONSISTENT WITH THE GROUP

 FRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQFRQ

                       FREQUENCY CALCULATION

 INFORMATION **** NOOPTGEOM **** NO OPTIMIZATION IS PERFORMED BEFORE FREQUENCIES CALCULATION
 INFORMATION **** INPFREQ ****  PRINTING OF EIGENVECTORS
 INFORMATION **** INPFREQ ****  ANALYSIS OF THE VIBRATIONAL MODES
 INFORMATION **** INPFREQ **** NEW DEFAULT - IMPOSING ECKART CONDITIONS TO THE HESSIAN

 GCALCO - MAX INDICES DIRECT LATTICE VECTOR    14    14    14
 NO.OF VECTORS CREATED 6999 STARS  105 RMAX    59.79196 BOHR

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

 ****   48 SYMMOPS - TRANSLATORS IN FRACTIONAL UNITS
 **** MATRICES AND TRANSLATORS IN THE CRYSTALLOGRAPHIC REFERENCE FRAME
   V INV                    ROTATION MATRICES                   TRANSLATORS
   1   1  1.00  0.00  0.00  0.00  1.00  0.00  0.00  0.00  1.00  0.00  0.00  0.00
   2   2  0.00  1.00  0.00  1.00  0.00  0.00 -1.00 -1.00 -1.00  0.00  0.00  0.00
   3   3 -1.00 -1.00 -1.00  0.00  0.00  1.00  0.00  1.00  0.00  0.00  0.00  0.00
   4   4  0.00  0.00  1.00 -1.00 -1.00 -1.00  1.00  0.00  0.00  0.00  0.00  0.00
   5   6  0.00  0.00  1.00  1.00  0.00  0.00  0.00  1.00  0.00  0.00  0.00  0.00
   6   5  0.00  1.00  0.00  0.00  0.00  1.00  1.00  0.00  0.00  0.00  0.00  0.00
   7   8  1.00  0.00  0.00  0.00  0.00  1.00 -1.00 -1.00 -1.00  0.00  0.00  0.00
   8   7  1.00  0.00  0.00 -1.00 -1.00 -1.00  0.00  1.00  0.00  0.00  0.00  0.00
   9  10 -1.00 -1.00 -1.00  0.00  1.00  0.00  1.00  0.00  0.00  0.00  0.00  0.00
  10   9  0.00  0.00  1.00  0.00  1.00  0.00 -1.00 -1.00 -1.00  0.00  0.00  0.00
  11  12  0.00  1.00  0.00 -1.00 -1.00 -1.00  0.00  0.00  1.00  0.00  0.00  0.00
  12  11 -1.00 -1.00 -1.00  1.00  0.00  0.00  0.00  0.00  1.00  0.00  0.00  0.00
  13  13  0.00 -1.00  0.00 -1.00  0.00  0.00  0.00  0.00 -1.00  0.00  0.00  0.00
  14  14 -1.00  0.00  0.00  0.00 -1.00  0.00  1.00  1.00  1.00  0.00  0.00  0.00
  15  16  0.00  0.00 -1.00  1.00  1.00  1.00  0.00 -1.00  0.00  0.00  0.00  0.00
  16  15  1.00  1.00  1.00  0.00  0.00 -1.00 -1.00  0.00  0.00  0.00  0.00  0.00
  17  17 -1.00  0.00  0.00  0.00  0.00 -1.00  0.00 -1.00  0.00  0.00  0.00  0.00
  18  18  0.00  0.00 -1.00  0.00 -1.00  0.00 -1.00  0.00  0.00  0.00  0.00  0.00
  19  21  0.00  0.00 -1.00 -1.00  0.00  0.00  1.00  1.00  1.00  0.00  0.00  0.00
  20  22  1.00  1.00  1.00 -1.00  0.00  0.00  0.00 -1.00  0.00  0.00  0.00  0.00
  21  19  0.00 -1.00  0.00  1.00  1.00  1.00 -1.00  0.00  0.00  0.00  0.00  0.00
  22  20  0.00 -1.00  0.00  0.00  0.00 -1.00  1.00  1.00  1.00  0.00  0.00  0.00
  23  23  1.00  1.00  1.00  0.00 -1.00  0.00  0.00  0.00 -1.00  0.00  0.00  0.00
  24  24 -1.00  0.00  0.00  1.00  1.00  1.00  0.00  0.00 -1.00  0.00  0.00  0.00
  25  25 -1.00  0.00  0.00  0.00 -1.00  0.00  0.00  0.00 -1.00  0.00  0.00  0.00
  26  26  0.00 -1.00  0.00 -1.00  0.00  0.00  1.00  1.00  1.00  0.00  0.00  0.00
  27  27  1.00  1.00  1.00  0.00  0.00 -1.00  0.00 -1.00  0.00  0.00  0.00  0.00
  28  28  0.00  0.00 -1.00  1.00  1.00  1.00 -1.00  0.00  0.00  0.00  0.00  0.00
  29  30  0.00  0.00 -1.00 -1.00  0.00  0.00  0.00 -1.00  0.00  0.00  0.00  0.00
  30  29  0.00 -1.00  0.00  0.00  0.00 -1.00 -1.00  0.00  0.00  0.00  0.00  0.00
  31  32 -1.00  0.00  0.00  0.00  0.00 -1.00  1.00  1.00  1.00  0.00  0.00  0.00
  32  31 -1.00  0.00  0.00  1.00  1.00  1.00  0.00 -1.00  0.00  0.00  0.00  0.00
  33  34  1.00  1.00  1.00  0.00 -1.00  0.00 -1.00  0.00  0.00  0.00  0.00  0.00
  34  33  0.00  0.00 -1.00  0.00 -1.00  0.00  1.00  1.00  1.00  0.00  0.00  0.00
  35  36  0.00 -1.00  0.00  1.00  1.00  1.00  0.00  0.00 -1.00  0.00  0.00  0.00
  36  35  1.00  1.00  1.00 -1.00  0.00  0.00  0.00  0.00 -1.00  0.00  0.00  0.00
  37  37  0.00  1.00  0.00  1.00  0.00  0.00  0.00  0.00  1.00  0.00  0.00  0.00
  38  38  1.00  0.00  0.00  0.00  1.00  0.00 -1.00 -1.00 -1.00  0.00  0.00  0.00
  39  40  0.00  0.00  1.00 -1.00 -1.00 -1.00  0.00  1.00  0.00  0.00  0.00  0.00
  40  39 -1.00 -1.00 -1.00  0.00  0.00  1.00  1.00  0.00  0.00  0.00  0.00  0.00
  41  41  1.00  0.00  0.00  0.00  0.00  1.00  0.00  1.00  0.00  0.00  0.00  0.00
  42  42  0.00  0.00  1.00  0.00  1.00  0.00  1.00  0.00  0.00  0.00  0.00  0.00
  43  45  0.00  0.00  1.00  1.00  0.00  0.00 -1.00 -1.00 -1.00  0.00  0.00  0.00
  44  46 -1.00 -1.00 -1.00  1.00  0.00  0.00  0.00  1.00  0.00  0.00  0.00  0.00
  45  43  0.00  1.00  0.00 -1.00 -1.00 -1.00  1.00  0.00  0.00  0.00  0.00  0.00
  46  44  0.00  1.00  0.00  0.00  0.00  1.00 -1.00 -1.00 -1.00  0.00  0.00  0.00
  47  47 -1.00 -1.00 -1.00  0.00  1.00  0.00  0.00  0.00  1.00  0.00  0.00  0.00
  48  48  1.00  0.00  0.00 -1.00 -1.00 -1.00  0.00  0.00  1.00  0.00  0.00  0.00

 DIRECT LATTICE VECTORS CARTESIAN COMPONENTS (ANGSTROM)
          X                    Y                    Z
   0.000000000000E+00   0.211407260000E+01   0.211407260000E+01
   0.211407260000E+01   0.000000000000E+00   0.211407260000E+01
   0.211407260000E+01   0.211407260000E+01   0.000000000000E+00


 CARTESIAN COORDINATES - PRIMITIVE CELL
 *******************************************************************************
 *      ATOM          X(ANGSTROM)         Y(ANGSTROM)         Z(ANGSTROM)
 *******************************************************************************
    1    12 MG    0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
    2     8 O     2.114072600000E+00  2.114072600000E+00  2.114072600000E+00

[ ... omissis ... ]

 THERE ARE NO SYMMETRY ALLOWED DIRECTIONS
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SYMM        TELAPSE        0.22 TCPU        0.10
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT INT_SCREEN  TELAPSE        0.23 TCPU        0.10

 *******************************************************************************
 *                                                                             *
 *                                                                             *
 *         FFFFF  RRRR   EEEE   QQQ   U   U  EEEE  N   N   CCC  Y   Y          *
 *         F      R   R  E     Q   Q  U   U  E     NN  N  C      Y Y           *
 *         FFF    RRRR   EEEE  Q   Q  U   U  EEEE  N N N  C       Y            *
 *         F      R   R  E     Q  QQ  U   U  E     N  NN  C       Y            *
 *         F      R   R  EEEE   QQ Q   UUU   EEEE  N   N   CCC    Y            *
 *                                                                             *
 *                                                                             *
 * CALCULATION OF PHONON FREQUENCIES AT THE GAMMA POINT.                       *
 *                                                                             *
 * SYMMETRY IS EXPLOITED TO BUILD THE TOTAL HESSIAN MATRIX.                    *
 * (F. PASCALE PHD THESIS TURIN-PARIS 2002)                                    *
 *                                                                             *
 *******************************************************************************


 * INTENSITIES COMPUTED VIA THE BERRY PHASE APPROACH                           *
 *                                                                             *
 * REFERENCES TO BE QUOTED WHEN USING THIS MODULE:                             *
 *                                                                             *
 * F. Pascale, C.M. Zicovich-Wilson, F. Lopez, B. Civalleri                    *
 * R. Orlando, R. Dovesi                                                       *
 * The calculation of the vibration frequencies of crystalline                 *
 * compounds and its implementation in the CRYSTAL code                        *
 * J. Comput. Chem. 25 (2004) 888-897                                          *
 *                                                                             *
 * C.M. Zicovich-Wilson, F. Pascale, C. Roetti, V.R. Saunders,                 *
 * R. Orlando, R. Dovesi                                                       *
 * The calculation of the vibration frequencies of alpha-quartz:               *
 * the effect of hamiltonian and basis set                                     *
 * J. Comput. Chem. 25 (2004) 1873-1881                                        *
 *******************************************************************************



 ATOMS ISOTOPIC MASS (AMU) FOR FREQUENCY CALCULATION 

    1 MG     23.9850    2 O      15.9949

 INFORMATION CONCERNING THE SCF+GRADIENT CALCULATIONS REQUIRED FOR
 GENERATING FREQUENCIES. IN PRINCIPLE 3N+1 SCF + GRADIENT
 CALCULATIONS  ARE REQUIRED;
 FOR EACH OF THEM THE REMAINING POINT SYMMETRY IS INDICATED.
 POINT SYMMETRY PERMITS TO GENERATE GRADIENTS FOR DISPLACEMENT B 
 STARTING FROM THE GRADIENT GENERATED BY DISPLACEMENT A.

   N   LABEL SYMBOL DISPLACEMENT     SYM. 

   1     EQUILIBRIUM GEOMETRY           48

   2      1   MG           DX        GENERATED BY TRANSLATIONAL INVARIANCE
   3      1   MG           DY        GENERATED BY TRANSLATIONAL INVARIANCE
   4      1   MG           DZ        GENERATED BY TRANSLATIONAL INVARIANCE
   5      2   O            DX           8
   6      2   O            DY        GENERATED FROM LINE X WITH OP   6
   7      2   O            DZ        GENERATED FROM LINE X WITH OP   5
 USE OF RESIDUAL SYMMETRY AFTER DISPLACEMENT

 NUMERICAL GRADIENT COMPUTED WITH A SINGLE DISPLACEMENT (+DX) FOR EACH
 CARTESIAN COORDINATE WITH RESPECT TO THE EQUILIBRIUM CONFIGURATION
 DX=  0.003 ANGSTROM
 NUMBER OF IRREDUCIBLE ATOMS                     2
 NUMBER OF SCF+GRADIENT CALCULATIONS             2

  ATOM  SYMOP  ORDER
    1     48      6
    2     48      6

 ATOM  : IRREDUCIBLE ATOM
 SYMOP : NUMBER OF SYMMETRY OPERATORS THAT DOESN'T MOVE THE IRREDUCIBLE ATOM
 ORDER : MAXIMUM ORDER AMONG THE OPERATORS OF THE IRREDUCIBLE ATOM

 *******************************************************************************


 GCALCO - MAX INDICES DIRECT LATTICE VECTOR    14    14    14
 NO.OF VECTORS CREATED 6999 STARS  105 RMAX    59.79196 BOHR

 CAPPA:IS1 12;IS2 12;IS3 12; K PTS MONK NET  72; SYMMOPS:K SPACE  48;G SPACE  48

 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT gordsh1     TELAPSE        0.29 TCPU        0.10

 MATRIX SIZE: P(G)  128375, F(G)   22615, P(G) IRR    4600, F(G) IRR    2292
 MAX G-VECTOR INDEX FOR 1- AND 2-ELECTRON INTEGRALS 679

 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT INPUT       TELAPSE        0.31 TCPU        0.12

 THERE ARE NO SYMMETRY ALLOWED DIRECTIONS
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SYMM        TELAPSE        0.31 TCPU        0.12
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT INT_SCREEN  TELAPSE        0.32 TCPU        0.12
 INFORMATION **** EXCBUF **** EXCH. BIPO BUFFER: WORDS USED =      561465

 DFT PARAMETERS

     ATOM       ELECTRONS   NET CHARGE   R(ANGSTROM)
   1  12  MG     12.0000      0.0000     1.60000000
   2   8  O       8.0000      0.0000     0.74000000

 SIZE OF GRID=       1296
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT MAKE_GRID2  TELAPSE        0.75 TCPU        0.56
 BECKE WEIGHT FUNCTION
 RADSAFE =     2.00
 TOLERANCES - DENSITY:10**- 6; POTENTIAL:10**- 9; GRID WGT:10**-14

 RADIAL INTEGRATION  - INTERVALS (POINTS,UPPER LIMIT):            1( 75,  4.0*R)

 ANGULAR INTEGRATION - INTERVALS (ACCURACY LEVEL [N. POINTS] UPPER LIMIT):
  1(  4[  86]   0.2)  2(  8[ 194]   0.5)  3( 12[ 350]   0.9)  4( 16[ 974]   3.5)
  5( 12[ 350]9999.0)
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT INT_CALC    TELAPSE        0.80 TCPU        0.60

 *******************************************************************************
 MGO Bulk - B3LYP/POB-TZVP - Geometry optimization 3rd run                       
 CRYSTAL - SCF - TYPE OF CALCULATION :  RESTRICTED CLOSED SHELL
 *******************************************************************************

 CAPPA:IS1 12;IS2 12;IS3 12; K PTS MONK NET  72; SYMMOPS:K SPACE  48;G SPACE  48

 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SDIK        TELAPSE        0.80 TCPU        0.61

[ ... omissis ... ]

 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFG      TELAPSE      197.20 TCPU      196.43
 INFORMATION **** EXCPOG **** EXCH. BIPO BUFFER LENGTH (WORDS) =     2245860
 INFORMATION **** GENPOG **** BIPO BUFFER LENGTH (WORDS) =      348500
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SHELXG      TELAPSE      206.52 TCPU      205.73

 CARTESIAN FORCES IN HARTREE/BOHR (ANALYTICAL)
   ATOM                     X                   Y                   Z
   1  12             3.384609451742E-16  3.384614387620E-16  5.677628646160E-16
   2   8            -3.338350159050E-16 -3.338355094927E-16 -5.677628646160E-16

 RESULTANT FORCE     4.625929269271E-18  4.625929269271E-18  0.000000000000E+00

 THERE ARE NO SYMMETRY ALLOWED DIRECTIONS
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT TOTGRA      TELAPSE      220.49 TCPU      219.66

 HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
                   FORCE CONSTANT MATRIX - NUMERICAL ESTIMATE
 HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
 MAX ABS(DGRAD): MAXIMUM ABSOLUTE GRADIENT DIFFERENCE WITH RESPECT TO
                 THE CENTRAL POINT
 DE:             ENERGY DIFFERENCE WITH RESPECT TO THE CENTRAL POINT
                 (DE IS EXPECTED TO BE POSITIVE FOR ALL DISPLACEMENTS)

   ATOM      MAX ABS(DGRAD)       TOTAL ENERGY (AU)  N.CYC      DE       SYM
    CENTRAL POINT             -2.754625427009E+02    37     0.0000E+00    48
    2 O  DX   2.9899E-04      -2.754625418373E+02     6     8.6363E-07     8
    2 O  DY                   GENERATED FROM A PREVIOUS LINE
    2 O  DZ                   GENERATED FROM A PREVIOUS LINE

 GCALCO - MAX INDICES DIRECT LATTICE VECTOR    14    14    14
 NO.OF VECTORS CREATED 6999 STARS  105 RMAX    59.79196 BOHR

 CAPPA:IS1 12;IS2 12;IS3 12; K PTS MONK NET  72; SYMMOPS:K SPACE  48;G SPACE  48

 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT gordsh1     TELAPSE      380.69 TCPU      379.41

 MATRIX SIZE: P(G)  128375, F(G)   22615, P(G) IRR    4600, F(G) IRR    2292
 MAX G-VECTOR INDEX FOR 1- AND 2-ELECTRON INTEGRALS 679

 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT INPUT       TELAPSE      380.70 TCPU      379.41

 THERE ARE NO SYMMETRY ALLOWED DIRECTIONS
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT SYMM        TELAPSE      380.70 TCPU      379.42
 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT INT_SCREEN  TELAPSE      380.71 TCPU      379.42

 ATOMIC BORN CHARGE TENSOR (UNITS OF e, ELECTRON CHARGE).
 DYNAMIC CHARGE = 1/3 * TRACE .

 ATOM   1 MG DYNAMIC CHARGE     2.012923

              1           2           3
   1     2.0129E+00  8.4145E-30  8.4145E-30
   2     6.3109E-30  2.0129E+00  1.2622E-29
   3     8.4145E-30  1.2622E-29  2.0129E+00

 ATOM   2 O  DYNAMIC CHARGE    -2.012923

              1           2           3
   1    -2.0129E+00 -8.4145E-30 -8.4145E-30
   2    -6.3109E-30 -2.0129E+00 -1.2622E-29
   3    -8.4145E-30 -1.2622E-29 -2.0129E+00

 +++ SYMMETRY ADAPTION OF VIBRATIONAL MODES +++

 SYMMETRY INFORMATION:
 K-LITTLE GROUP: CLASS TABLE, CHARACTER TABLE.
 IRREP-(DIMENSION, NO. IRREDUCIBLE SETS)
 [WARNINGS: (1) ONLY ACTIVE IRREPS ARE GENERATED AND LISTED.
            (2) ONLY RELEVANT CLASSES ARE CONSIDERED IN THE CHARACTER TABLE
            (3) SYMBOLS MAY NOT FULLY COINCIDE WITH THOSE FROM TEXT BOOKS.]

 (P, D, RP, RD, STAND FOR PAIRING, DOUBLING, REAL PAIRING AND REAL DOUBLING 
 OF THE IRREPS (SEE MANUAL))

 CLASS  | GROUP OPERATORS (SEE SYMMOPS KEYWORD)
 --------------------------------------------------------------------
 C2     |   2;   3;   4;
 C3     |   5;  11;   7;   9;   6;  12;  10;   8;
 C2'    |  13;  14;  17;  18;  23;  24;
 C4     |  15;  16;  19;  20;  21;  22;
 I      |  25;
 SGH    |  26;  27;  28;
 S6     |  29;  35;  31;  33;  30;  36;  34;  32;
 SGD    |  37;  38;  41;  42;  47;  48;
 S4     |  39;  40;  43;  44;  45;  46;

 IRREP/CLA      E     C2     C3    C2'     C4      I    SGH     S6    SGD     S4
 -------------------------------------------------------------------------------
  MULTIP |      1      3      8      6      6      1      3      8      6      6
 -------------------------------------------------------------------------------
    Fu   |   3.00  -1.00   0.00  -1.00   1.00  -3.00   1.00   0.00   1.00  -1.00
 

 Fu -(3,   2);

 BORN CHARGE VECTOR IN THE BASIS OF NORMAL MODES (UNITS OF e*M_E**(-1/2) ).
 e AND M_E ARE UNITS OF ELECTRON CHARGE AND MASS, RESPECTIVELY.

 MODE         X              Y              Z
   1    0.13002E-17    0.41026E-47    0.41026E-47
   2    0.00000E+00    0.00000E+00    0.00000E+00
   3    0.00000E+00    0.00000E+00    0.00000E+00
   4   -0.63622E-31   -0.15220E-01   -0.95433E-31
   5   -0.15220E-01   -0.47717E-31   -0.63622E-31
   6   -0.63622E-31   -0.95433E-31   -0.15220E-01

 HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
 VIBRATIONAL CONTRIBUTIONS TO THE STATIC DIELECTRIC TENSOR (OSCILLATOR
 STRENGTHS) ARE PURE NUMBERS. THEY ARE COMPUTED FOR EACH nth MODE AS:

   f_(n,ij) = 1 / (4 * pi * eps0) * 4 * pi / V * Z_(n,i) * Z_(n,j) / nu_n**2

 WHERE:
  1/(4*pi*eps0)  1 A.U. [M*L**3*T**(-4)*C**(-2)]
  V              CELL VOLUME (BOHR**3) [L**3]
  Z_(n,i)        ith COMPONENT OF BORN CHARGE VECTOR IN THE BASIS
                 OF NORMAL MODES ( e*M_E**(-1/2) ) [C*T*M**(-1/2)]
  nu_n           FREQUENCY (HARTREE) [T**(-1)]
  e, M_E         UNITS OF ELECTRON CHARGE AND MASS, RESPECTIVELY
  M, L, T, C     MASS, LENGTH, TIME, CURRENT, RESPECTIVELY

 HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

 MODE      CARTESIAN AXES SYSTEM

   4       0.00000    0.00000    0.00000
           0.00000    7.57072    0.00000
           0.00000    0.00000    0.00000

   5       7.57072    0.00000    0.00000
           0.00000    0.00000    0.00000
           0.00000    0.00000    0.00000

   6       0.00000    0.00000    0.00000
           0.00000    0.00000    0.00000
           0.00000    0.00000    7.57072

 SUM TENSOR OF THE VIBRATIONAL CONTRIBUTIONS TO THE' STATIC DIELECTRIC     TENSOR

           7.57072    0.00000    0.00000
           0.00000    7.57072    0.00000
           0.00000    0.00000    7.57072

 INTEGRATED IR INTENSITIES, IN UNITS OF KM/MOL, ARE COMPUTED UNDER THE
 HYPOTHESIS OF ISOTROPIC RESPONSE (I.E. POWDER SAMPLE):

   INTENS_n = 1 / (4 * pi * eps0) * pi * N_AV / 3 / c**2 *
                             * d_n * ( Z_(n,x)**2 + Z_(n,y)**2 + Z_(n,z)**2 )
            = 0.17770712E+07 * d_n * ( Z_(n,x)**2 + Z_(n,y)**2 + Z_(n,z)**2 )

 WHERE:
  1/(4*pi*eps0)  1 A.U. [M*L**3*T**(-4)*C**(-2)]
  N_AV           AVOGADRO'S NUMBER [QM**(-1)]
  c              SPEED OF LIGHT [L*T**(-1)]
  d_n            DEGENERACY OF THE MODE
  Z_(n,x)        xth COMPONENT OF BORN CHARGE VECTOR IN THE BASIS
                 OF NORMAL MODES ( e*M_E**(-1/2) ) [C*T*M**(-1/2)]
  M,L,T,C,QM     MASS, LENGTH, TIME, CURRENT, QUANTITY OF MATTER, RESPECTIVELY

 HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
 EIGENVALUES (EIGV) OF THE MASS WEIGHTED HESSIAN MATRIX AND HARMONIC TRANSVERSE
 OPTICAL (TO) FREQUENCIES. IRREP LABELS REFER TO SYMMETRY REPRESENTATION
 ANALYSIS; A AND I INDICATE WHETHER THE MODE IS ACTIVE OR INACTIVE,
 RESPECTIVELY, FOR IR AND RAMAN; INTEGRATED IR INTENSITIES IN BRACKETS.

 CONVERSION FACTORS FOR FREQUENCIES:
      1 CM**(-1) =   0.4556335E-05 HARTREE
      1 THZ      =   0.3335641E+02 CM**(-1)

 HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    MODES         EIGV          FREQUENCIES     IRREP  IR   INTENS    RAMAN
             (HARTREE**2)   (CM**-1)     (THZ)             (KM/MOL)
    1-   3    0.0000E+00      0.0000    0.0000  (Fu )   A (     0.00)   I
    4-   6    0.3015E-05    381.0964   11.4250  (Fu )   A (  1234.92)   I

 HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
 VIBRATIONAL ANALYSIS:
 1) EACH PAIR OF BONDED ATOMS (I.E. WITHIN THEIR VAN DER WAALS DISTANCE) A AND
 B IS EXAMINED TO SEE IF THERE IS A LARGE RELATIVE MOTION BETWEEN THEM.
 2) IF SO, THE AB MOTION IS DECOMPOSED IN THREE COMPONENTS: ALONG A-B(LONG),
 ON THE PLANE CONTAING A THIRD ATOM C (ANG) AND OUT OF THE PLANE (OUT).
 LONG+ANG+OUT=1.
 3) THE MODES ARE CLASSIFIED SO AS: LONG > 0.85 -> STRETCHING (S);
 ANG > 0.85 -> BENDING (B); OTHERWISE -> OTHER (O). TYPE (?) MEANS THE MODE IS
 LIKELY TO BE ROTATIONAL. THE POTENTIAL ENERGY CONTRIBUTION (PEC) IS IN KJ/MOL.
 HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  MODE FRQ(CM**-1) IRREP TYP     A      B  LONG      C  ANG     PEC
    1      0.0000  (Fu )
    2      0.0000  (Fu )
    3      0.0000  (Fu )
    4    381.0964  (Fu )
                         (B)    2 O    1 MG(0.0)   2 O (1.0)   -0.000
                         (S)    2 O    1 MG(1.0)               -0.000
                         (B)    2 O    1 MG(0.0)   2 O (1.0)   -0.000
                         (B)    2 O    1 MG(0.0)   2 O (1.0)   -0.000
    5    381.0964  (Fu )
                         (S)    2 O    1 MG(1.0)                0.000
                         (B)    2 O    1 MG(0.0)   2 O (1.0)    0.000
                         (B)    2 O    1 MG(0.0)   2 O (1.0)    0.000
                         (B)    2 O    1 MG(0.0)   2 O (1.0)    0.000
    6    381.0964  (Fu )
                         (B)    2 O    1 MG(0.0)   2 O (1.0)    0.000
                         (B)    2 O    1 MG(0.0)   2 O (1.0)    0.000
                         (S)    2 O    1 MG(1.0)                0.000
                         (S)    2 O    1 MG(1.0)                0.000

 NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES (IN BOHR)

 FREQ(CM**-1)      0.00      0.00      0.00    381.10    381.10    381.10

 AT.   1 MG X     0.1582    0.0000    0.0000    0.0000   -0.0726    0.0000
            Y     0.0000    0.1582    0.0000   -0.0726    0.0000    0.0000
            Z     0.0000    0.0000    0.1582    0.0000    0.0000   -0.0726
 AT.   2 O  X     0.1582    0.0000    0.0000    0.0000    0.1089    0.0000
            Y     0.0000    0.1582    0.0000    0.1089    0.0000    0.0000
            Z     0.0000    0.0000    0.1582    0.0000    0.0000    0.1089

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
                generator = QHAInputCreator(interface='crystal')
                generator.read(filename, False, reference=0, use_symm=False)
                generator.write_from_multiple_files(
                    outfile, 'MgO Gamma-point frequencies', ref=0)

                reader = QHAInputFileReader(outfile)

                tmp.close()
                os.unlink(tmp.name)

            cryout.close()
            os.unlink(cryout.name)

        self.assertTrue(reader.completed,
                        'Reading process should be completed')
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
        return
                
