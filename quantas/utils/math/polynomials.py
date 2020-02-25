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
This module contains functions used to handle polynomial fitting.
"""
import numpy as np
from .fast_math import multi_interpolate, multi_R2

def polyfit(X, Y, degree):
    """ Polynomial regression based on numpy.polynomial.polynomial
    """
    if len(Y.shape) == 1:
        return np.polynomial.polynomial.polyfit(X, Y, degree)
    else:
        return np.polynomial.polynomial.polyfit(X, Y.T, degree).T

def single_R_squared(X, Y, pars):
    """ Calculates the R-squared values related to a polynomial fit """
    p = np.polynomial.polynomial.Polynomial(pars)
    yhat = p(X)
    ybar = np.sum(Y)/(len(X))
    ssreg = np.sum((yhat - ybar)**2)
    sstot = np.sum((Y - ybar)**2)
    return ssreg/sstot

def R_squared(X, Y, pars):
    """ Calculates the R-squared values related to a polynomial fit """
    if len(Y.shape) == 1:
        return single_R_squared(X, Y, pars)
    else:
        return multi_R2(X, Y, pars)

def interpolate(X, pars):
    """
    """
    if len(pars.shape) == 1 and type(X) != np.ndarray:
        values = np.polynomial.polynomial.polyval(X, pars)
    elif len(pars.shape) == 1 and type(X) == np.ndarray:
        values = np.polynomial.polynomial.polyval(X, pars)
    elif len(pars.shape) == 2:
        values = multi_interpolate(X, pars)
    elif len(pars.shape) == 3:
        values = np.zeros( (pars.shape[0], pars.shape[1], X.shape[0]),
                           dtype = np.float64)
        for i in range(pars.shape[0]):
            values[i] = multi_interpolate(X, pars[i])
    return values

def find_polynomial_minimum(parameters, addendum=None):
    """
    """
    f = np.polynomial.polynomial.Polynomial(parameters)
    df = f.deriv()
    if addendum is not None:
        df += addendum
    crit_points = df.roots()
    real_points = crit_points[crit_points.imag==0].real
    test = f.deriv(2)(real_points)
    minima = real_points[test > 0]
    return minima
