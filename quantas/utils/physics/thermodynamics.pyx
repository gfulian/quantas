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

import scipy.constants as cs
import numpy as np

cdef double NA = cs.Avogadro

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef enthalpy(double[:,::1] V, double[::1] p, double[:,::1] U):
    """
    This functions calculates the enthalpy (H) of the system
    according to the formula:

    .. math:: H\\big(T,P \\big) = U\\big(T,P \\big) + pV\\big(T,P \\big)

    Attributes
    ----------
    V: ndarray
         2D array of volume values (in m^3).
    p: ndarray
         2D array of pressure values (in Pa).
    U: ndarray
         2D array of internal energy (in kJ/mol).
        
    Returns
    -------
    H: ndarray
         2D array of enthalpy values (in kJ/mol).

    """
    cdef Py_ssize_t n = V.shape[0]
    cdef Py_ssize_t m = V.shape[1]
    cdef int i, j
    cdef double[::1] p_v = p
    cdef double[:,::1] V_v = V
    cdef double[:,::1] U_v = U

    result = np.zeros( (n,m), dtype=np.float64 )
    cdef double[:,::1] result_v = result
    
    for i in prange(n, nogil=True):
        for j in range(m):
            result_v[i,j] = U[i,j] + p[j]*V[i,j]*NA/1000
    
    return result

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef gibbs(double[:,::1] V, double[::1] p, double[:,::1] F):
    """
    This functions calculates the Gibbs free energy (G) of the system
    according to the formula:

    .. math:: G \\big(T,P \\big) = F\\big(T,P \\big) + pV\\big(T,P \\big)

    Attributes
    ----------
    V: ndarray
         2D array of volume values (in m^3).
    p: ndarray
         2D array of pressure values (in Pa).
    F: ndarray
         2D array of Helmholtz free energy (in kJ/mol).
        
    Returns
    -------
    H: ndarray
         2D array of Gibbs free energy values (in kJ/mol).

    """
    cdef Py_ssize_t n = V.shape[0]
    cdef Py_ssize_t m = V.shape[1]
    cdef int i, j
    cdef double[::1] p_v = p
    cdef double[:,::1] V_v = V
    cdef double[:,::1] F_v = F

    result = np.zeros( (n,m), dtype=np.float64 )
    cdef double[:,::1] result_v = result
    
    for i in prange(n, nogil=True):
        for j in range(m):
            result_v[i,j] = F[i,j] + p[j]*V[i,j]*NA/1000.
    
    return result


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef adiabatic_bulk_modulus(double[::1] T, double[:,::1] V, double[:,::1] KT,
                             double[:,::1] alpha, double[:,::1] Cv, double pf):
    """
    This method calculates the adiabatic bulk modulus (KS) according to:

    .. math::

        K_S = K_T + \\frac{\\alpha_V^2 V^2 K_T}{C_V}
        

    Attributes
    ----------
    T: ndarray
        1D array of temperature values (in Kelvin).
    V: ndarray
        2D array of volume values (in m^3).
    KT: ndarray
        2D array of isothermal bulk modulus at P(V)-T condition.
    alpha: ndarray
        2D array of thermal expansion coefficient at P(V)-T condition (in K^-1).
    Cv: ndarray
        2D array of isochoric heat capacity at P(V)-T condition (in J mol^-1).
    factor: double
        factor to convert the from and to Pa.

    Returns
    -------
    ndarray
        Adiabatic bulk modulus (KS).
        """
    cdef Py_ssize_t n = V.shape[0] 
    cdef Py_ssize_t m = V.shape[1]
    cdef int i, j
    cdef double[::1] T_v = T
    cdef double[:,::1] V_v = V
    cdef double[:,::1] KT_v = KT
    cdef double[:,::1] alpha_v = alpha
    cdef double[:,::1] Cv_v = Cv
    cdef double factor = pf

    KS =  np.zeros( (n,m), dtype=np.float64 )
    cdef double[:,::1] KS_v = KS

    for i in prange(n, nogil=True):
        for j in range(m):
            if T_v[i] == 0. or Cv_v[i,j] == 0.:
                KS_v[i,j] = KT_v[i,j]
            else:
                KS_v[i,j] = KT_v[i,j] + ( alpha_v[i,j] * alpha_v[i,j] * \
                    (factor * KT_v[i,j] * factor * KT_v[i,j]) * \
                    V_v[i,j] * T_v[i] / (Cv_v[i,j]/NA)
                    )/factor
    return KS


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef gruneisen_parameter(double[:,::1] V, double[:,::1] KT,
                          double[:,::1] alpha, double[:,::1] Cv):
    """
    This method calculates the the Grüneisen parameter by

    .. math::

       \\gamma = \\frac{ \\alpha K_T V }{ C_V }.
        

    Attributes
    ----------
    V: ndarray
        2D array of volume values (in m^3).
    KT: ndarray
        2D array of isothermal bulk modulus at P(V)-T condition(in Pa).
    alpha: ndarray
        2D array of thermal expansion coefficient at P(V)-T condition (in K^-1).
    Cv: ndarray
        2D array of isochoric heat capacity at P(V)-T condition (in J mol^-1).

    Returns
    -------
    ndarray
        Grüneisen parameters (gamma).
    """
    cdef Py_ssize_t n = V.shape[0] 
    cdef Py_ssize_t m = V.shape[1]
    
    cdef int i, j
    cdef double[:,::1] V_v = V
    cdef double[:,::1] KT_v = KT
    cdef double[:,::1] alpha_v = alpha
    cdef double[:,::1] Cv_v = Cv

    gamma =  np.zeros( (n,m), dtype=np.float64 )
    cdef double[:,::1] gamma_v = gamma

    for i in prange(n, nogil=True):
        for j in range(m):
            if Cv_v[i,j] == 0.:
                gamma_v[i,j] = 0.
            else:
                gamma_v[i,j] = V_v[i,j] * KT_v[i,j] * alpha_v[i,j] /\
                               (Cv_v[i,j]/NA)
    return gamma
