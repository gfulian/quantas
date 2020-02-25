# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

cimport cython
from cython.parallel import prange

import numpy as np

cdef extern from "math.h" nogil:
    double sum(double x)
    double pow(double x, double y)
    double sin(double x)
    double cos(double x)

cdef double sInterp(double x, double n, double coeff) nogil:
    return coeff * pow(x, n)

@cython.boundscheck(False)
@cython.wraparound(False)
def multi_interpolate_array(double[::1] X, double[:,::1] coeffs):
    cdef Py_ssize_t nv = coeffs.shape[0]
    cdef Py_ssize_t nc = coeffs.shape[1]
    cdef Py_ssize_t nx = X.shape[0]
    cdef Py_ssize_t i, j, k

    result = np.zeros( (nv,nx), dtype=np.float64 )
    cdef double[:,::1] result_view = result

    for i in prange(nv, nogil=True):
        for j in range(nc):
            for k in range(nx):
                result_view[i, k] += sInterp(
                    X[k], j, coeffs[i,j])
    
    return result

@cython.boundscheck(False)
@cython.wraparound(False)
def multi_interpolate_scalar(double x, double[:,::1] coeffs):
    cdef Py_ssize_t nv = coeffs.shape[0]
    cdef Py_ssize_t nc = coeffs.shape[1]
    cdef Py_ssize_t i, j

    result = np.zeros( nv, dtype=np.float64 )
    cdef double[::1] result_view = result

    for i in prange(nv, nogil=True):
        for j in range(nc):
            result_view[i] += sInterp(x, j, coeffs[i,j])
    
    return result


def multi_interpolate(X, coeffs):
    if type(X) != np.ndarray:
        return multi_interpolate_scalar(X, coeffs)
    else:
        return multi_interpolate_array(X, coeffs)


@cython.boundscheck(False)
@cython.wraparound(False)
def multi_R2(double[::1] X, double[:,::1] Y, double[:,::1] coeffs):
    cdef Py_ssize_t nx = X.shape[0]
    cdef Py_ssize_t ny = Y.shape[0]
    cdef Py_ssize_t nc = coeffs.shape[1]

    cdef Py_ssize_t i, j

    yhat = multi_interpolate_array(X, coeffs)
    ybar = np.zeros( ny, dtype=np.float64 )
    ssreg = np.zeros( ny, dtype=np.float64 )
    sstot = np.zeros( ny, dtype=np.float64 )
    R2 = np.zeros( ny, dtype=np.float64 )
    
    cdef double[::1] ybar_view = ybar
    cdef double[::1] ssreg_view = ssreg
    cdef double[::1] sstot_view = sstot
    cdef double[::1] R2_view = R2
    cdef double[:,::1] Y_view = Y
    cdef double[:,::1] yhat_view = yhat

    for i in prange(ny, nogil=True):
        for j in range(nx):
            ybar_view[i] += Y_view[i,j]/nx

    for i in prange(ny, nogil=True):
        for j in range(nx):
            ssreg_view[i] += pow(yhat_view[i,j]-ybar_view[i],2.)
            sstot_view[i] += pow(Y_view[i,j]-ybar_view[i],2.)

    for i in prange(ny, nogil=True):
        if sstot_view[i] == 0.:
            R2_view[i] = 0.
        else:
            R2_view[i] = ssreg_view[i]/sstot_view[i]
    return R2


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef vector(double theta, double phi):
    vec = np.zeros( 3, dtype=np.float64 )
    cdef double[::1] v = vec
    v[0] = sin(theta)*cos(phi)
    v[1] = sin(theta)*sin(phi)
    v[2] = cos(theta)
    return vec

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cofactor(double[:,::1] mat):
    cof = np.zeros( (3,3), dtype=np.float64 )
    cdef double[:,::1] c = cof
    cdef double[:,::1] m = mat
    c[0][0] = m[1][1]*m[2][2] - m[1][2]*m[2][1]
    c[0][1] = m[1][2]*m[2][0] - m[1][0]*m[2][2]
    c[0][2] = m[1][0]*m[2][1] - m[1][1]*m[2][0]
    c[1][0] = m[0][2]*m[2][1] - m[0][1]*m[2][2]
    c[1][1] = m[0][0]*m[2][2] - m[0][2]*m[2][0]
    c[1][2] = m[0][1]*m[2][0] - m[0][0]*m[2][1]
    c[2][0] = m[0][1]*m[1][2] - m[0][2]*m[1][1]
    c[2][1] = m[0][2]*m[1][0] - m[0][0]*m[1][2]
    c[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0]
    return

