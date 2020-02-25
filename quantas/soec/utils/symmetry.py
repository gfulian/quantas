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
import itertools as it

Cij_symmetry = {
   'cubic': np.array([[1, 7, 7, 0, 0, 0],
                      [7, 1, 7, 0, 0, 0],
                      [7, 7, 1, 0, 0, 0],
                      [0, 0, 0, 4, 0, 0],
                      [0, 0, 0, 0, 4, 0],
                      [0, 0, 0, 0, 0, 4]]),

   'hexagonal': np.array([[1, 7, 8, 9, 10, 0],
                          [7, 1, 8, 0,-9, 0],
                          [8, 8, 3, 0, 0, 0],
                          [9, -9, 0, 4, 0, 0],
                          [10, 0, 0, 0, 4, 0],
                          [0, 0, 0, 0, 0, 6]]),

   'trigonal_high': np.array([[1, 7, 8, 9, 10, 0],
                              [7, 1, 8, 0,-9, 0],
                              [8, 8, 3, 0, 0, 0],
                              [9, -9, 0, 4, 0, 0],
                              [10, 0, 0, 0, 4, 0],
                              [0, 0, 0, 0, 0, 6]]),

   'trigonal_low': np.array([[1,  7,  8,  9,  10,  0 ],
                             [7,  1,  8, -9, -10,  0 ],
                             [8,  8,  3,  0,   0,  0 ],
                             [9, -9,  0,  4,   0, -10],
                             [10,-10, 0,  0,   4,  9 ],
                             [0,  0,  0, -10 , 9,  6 ]]),

   'tetragonal_high': np.array([[1, 7, 8, 0, 0, 0],
                                [7, 1, 8, 0, 0, 0],
                                [8, 8, 3, 0, 0, 0],
                                [0, 0, 0, 4, 0, 0],
                                [0, 0, 0, 0, 4, 0],
                                [0, 0, 0, 0, 0, 6]]),

   'tetragonal_low': np.array([[1, 7, 8, 0, 0, 11],
                               [7, 1, 8, 0, 0, -11],
                               [8, 8, 3, 0, 0, 0],
                               [0, 0, 0, 4, 0, 0],
                               [0, 0, 0, 0, 4, 0],
                               [11, -11, 0, 0, 0, 6]]),

   'orthorhombic':    np.array([[ 1,  7,  8,  0,  0,  0],
                                [ 7,  2, 12,  0,  0,  0],
                                [ 8, 12,  3,  0,  0,  0],
                                [ 0,  0,  0,  4,  0,  0],
                                [ 0,  0,  0,  0,  5,  0],
                                [ 0,  0,  0,  0,  0,  6]]),

   'monoclinic': np.array([[ 1,  7,  8,  0,  10,  0],
                           [ 7,  2, 12,  0, 14,  0],
                           [ 8, 12,  3,  0, 17,  0],
                           [ 0,  0,  0,  4,  0,  20],
                           [10, 14, 17,  0,  5,  0],
                           [ 0,  0,  0, 20,  0,  6]]),

    'triclinic': np.array([[ 1,  7,  8,  9,  10, 11],
                           [ 7,  2, 12,  13, 14, 15],
                           [ 8, 12,  3,  16, 17, 18],
                           [ 9, 13, 16,  4,  19, 20],
                           [10, 14, 17, 19,  5,  21],
                           [11, 15, 18, 20,  21, 6 ]]),
   }

def check_symmetry(m):
    """
    This function return the crystal system related to the input SOEC matrix.

    Arguments
    ---------
    m: ndarray(dtype=float, dim=2)
        2D-array of stiffness matrix in Voigt's notation (6x6).

    Returns
    -------
    symm: string
        String description of the crystal system.
    """

    for s in Cij_symmetry.keys():
        if run_check(m, s) is True:
            symm = s
            break
        else:
            continue
    return symm


def run_check(m, symmetry):
    tol = 1.e-3

    # Make mapping
    Cij_map = {}
    for i in range(6):
        for j in range(i,6):
            Cij_map[(i,j)] = Cij_symmetry[symmetry][i,j]

    # Add the lower triangle to Cij_map, e.g. C21 = C12
    for (i1,i2) in Cij_map.copy().keys():
        Cij_map[(i2,i1)] = Cij_map[(i1,i2)]

    # Check the symmetry of the SOEC matrix according to the map
    for i in range(0,22):
        l = []
        for key in Cij_map.keys():
            # Collect the combinations
            if Cij_map[key] == i:
                l.append(key)
        if len(l) != 0:
            if i == 0:
                for index in l:
                    if abs(m[index]) > tol:
                        return False
                    else:
                        continue
            else:
                for k in it.combinations(l, 2):
                    if abs(m[k[0]] - m[k[1]]) > tol:
                        return False
                    else:
                        continue

    # Additional check for hexagonal system
    if symmetry == 'hexagonal':
        if abs(m[5][5]-((m[0][0]-m[0][1])/2.)) > tol:
            return False
    return True
