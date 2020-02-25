# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" Derivative module. """
import numpy as np

class LenError(Exception):
    pass

def derivative(X, Y):
    """ This function returns the derivativeof two 1-D arrays.

    Parameters
    ----------

    X: ndarray(float, ndim=1)
        Array of independent variable.
    Y: ndarray(float, ndim=1)
        Array of dependent variable

    Returns
    -------

    ndarray(float, ndim=1)
        Array of derivative dX/dY

    """
    if X.shape[0] <= 1 and Y.shape[0] <= 1:
        raise LenError
    else:
        return np.gradient(Y) / np.gradient(X)
