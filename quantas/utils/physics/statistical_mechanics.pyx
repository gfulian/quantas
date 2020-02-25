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

cdef double H = cs.Planck
cdef double KB = cs.Boltzmann
cdef double NA = cs.Avogadro

cdef extern from "math.h" nogil:
    double log(double x)
    double exp(double x)
    double expm1(double x)
    double pow(double x, double y)


cdef inline double ho_zpe(double temperature, double omega) nogil:
    """ Zero-point energy of single harmonic oscillator. """
    if omega <= 0.:
        return 0.
    else:
        return pow(10., -3.) * NA * H * (omega * 0.5)

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef zero_point_energy(double[::1] temperature, double[:, :, ::1] band, double[::1] weights):
    """
    Calculate the zero-point energy for a complete phonon band structure as:

    .. math::

        U_{zp}\\big(V \\big) = \\sum_{\\vec{k}} \\sum_{i=0}^{3N}
        \\frac{1}{2} h \\nu_i(\\vec{k})

    with :math:`\\vec{k}` the sampled *k*-points, :math:`N` the number of
    atoms in the considered unit cell and :math:`h` the Planck constant.

    .. note::

        Zero-point energy is independent on temperature.

    Parameters
    ----------

    temperature: ndarray
        Array of temperature values with `float` type.

    band: ndarray(ndim=3)
        Array of phonon band frequencies with `float` type.

    weights: ndarray
        Array of weights for each phonon band with `float` type.

    Returns
    -------

    result: ndarray
        Array of zero-point energy values with `float` type.

    """
    cdef Py_ssize_t nt = temperature.shape[0]
    cdef Py_ssize_t nb = band.shape[0]
    cdef Py_ssize_t nf = band.shape[1]
    cdef Py_ssize_t nv = band.shape[2]
    cdef int i, j, k, n

    result = np.zeros( (nt,nv), dtype=np.float64 )
    cdef double[:,::1] result_view = result


    for i in prange(nt, nogil=True):
        for j in range(nb):
            for k in range(nf):
                for n in range(nv):
                    result_view[i, n] += ho_zpe(temperature[i], band[j,k,n]) * \
                                         weights[j]

    return result


cdef inline double ho_eth(double temperature, double omega) nogil:
    """ Thermal internal energy of single harmonic oscillator. """
    cdef double x, n, l
    if omega <= 0. or temperature == 0.:
        return 0.
    else:
        n = H * omega
        x = n / (KB * temperature)
        return pow(10., -3.) * NA * (n/expm1(x))


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef thermal_energy(double[::1] temperature, double[:, :, ::1] band, double[::1] weights):
    """
    Calculate the thermal internal energy for a complete phonon band
    structure as:

    .. math::

        U_{th}\\big(T, V \\big) = \\sum_{\\vec{k}} \\sum_{i=0}^{3N}
        \\frac{\\nu_i(\\vec{k})}
        {e^{\\frac{h \\nu_i(\\vec{k})}{k_B T}} -1}

    with :math:`\\vec{k}` the sampled *k*-points, :math:`N` the number of
    atoms in the considered unit cell, :math:`h` the Planck constant and
    :math:`k_B` the Boltzmann constant.

    Parameters
    ----------

    temperature: ndarray
        Array of temperature values with `float` type.

    band: ndarray(ndim=3)
        Array of phonon band frequencies with `float` type.

    weights: ndarray
        Array of weights for each phonon band with `float` type.

    Returns
    -------

    result: ndarray
        Array of thermal internal energy values with `float` type.

    """
    cdef Py_ssize_t nt = temperature.shape[0]
    cdef Py_ssize_t nb = band.shape[0]
    cdef Py_ssize_t nf = band.shape[1]
    cdef Py_ssize_t nv = band.shape[2]
    cdef int i, j, k, n

    result = np.zeros( (nt,nv), dtype=np.float64 )
    cdef double[:,::1] result_view = result


    for i in prange(nt, nogil=True):
        for j in range(nb):
            for k in range(nf):
                for n in range(nv):
                    result_view[i, n] += ho_eth(temperature[i], band[j,k,n]) * \
                                         weights[j]

    return result

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef internal_energy(double[::1] U0, double[::1] Uzp, double[:,::1] Uth):
    """
    Calculate the total internal energy of the system as:

    .. math::

        U\\big(T, V \\big) = U_{zp}\\big(V \\big) + U_{th}\\big(T, V \\big) +
        U_0\\big(V \\big)

    Returns
    -------

    result: ndarray
        Array of total internal energy values with `float` type.

    """
    cdef Py_ssize_t nt = Uth.shape[0]
    cdef Py_ssize_t nv = U0.shape[0]
    cdef int i, j

    cdef double[::1] U0_v = U0
    cdef double[::1] Uzp_v = Uzp
    cdef double[:,::1] Uth_v = Uth

    result = np.zeros( (nt,nv), dtype=np.float64 )
    cdef double[:,::1] result_view = result


    for i in prange(nt, nogil=True):
        for j in range(nv):
            result_view[i, j] += U0_v[j] + Uzp_v[j] + Uth_v[i,j]

    return result

cdef double ho_S(double temperature, double omega) nogil:
    """ Entropy of single harmonic oscillator. """
    cdef double x, n, l
    if omega <= 0. or temperature == 0.:
        return 0.
    else:
        x = H * omega / (KB * temperature)
        n = 1./expm1(x)
        l = log(1. - exp(-x))
        return NA * KB * (n * x - l)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef entropy(double[::1] temperature, double[:, :, ::1] band, double[::1] weights):
    """
    Calculate the entropy for a complete phonon band structure according to:

    .. math::

        S\\big(T, V \\big) = N_A k_B \\sum_{\\vec{k}} \\sum_{i=0}^{3N}{
        \\Bigg[
        \\frac{1}{k_B T} \\frac{h \\nu_i(\\vec{k})}
        {e^{\\frac{h \\nu_i(\\vec{k})}{k_B T}} - 1} -
        ln \\big(1 - e^{- \\frac{h \\nu_i(\\vec{k})}{k_B T}} \\big)
        \\Bigg]}

    with :math:`\\vec{k}` the sampled *k*-points, :math:`N` the number of
    atoms in the considered unit cell, :math:`h` the Planck constant and
    :math:`k_B` the Boltzmann constant.

    Parameters
    ----------

    temperature: ndarray
        Array of temperature values with `float` type.

    band: ndarray(ndim=3)
        Array of phonon band frequencies with `float` type.

    weights: ndarray
        Array of weights for each phonon band with `float` type.

    Returns
    -------
    result: ndarray(ndim=2)
        2D matrix containing the entropy values (in J/mol).

    """
    cdef Py_ssize_t nt = temperature.shape[0]
    cdef Py_ssize_t nb = band.shape[0]
    cdef Py_ssize_t nf = band.shape[1]
    cdef Py_ssize_t nv = band.shape[2]
    cdef int i, j, k, n

    result = np.zeros( (nt,nv), dtype=np.float64 )
    cdef double[:,::1] result_view = result

    for i in prange(nt, nogil=True):
        for j in range(nb):
            for k in range(nf):
                for n in range(nv):
                    result_view[i, n] += ho_S(temperature[i], band[j,k,n]) * \
                                         weights[j]

    return result


cdef inline double ho_F(double temperature, double omega) nogil:
    """ Vibrational free energy of single harmonic oscillator. """
    cdef double x, hw, kt
    if omega <= 0.:
        return 0.

    hw = H * omega
    if temperature == 0.:
        return pow(10.,-3) * NA * (0.5 * hw)
    else:
        x = H * omega / (KB * temperature)
        kt = KB * temperature
        return pow(10.,-3) * NA * (0.5 * hw + kt * log(1.-exp(-x)))


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef vibrational_free_energy(double[::1] temperature, double[:, :, ::1] band, double[::1] weights):
    """
    Calculate the vibrational free energy for a complete phonon band
    structure as:

    .. math::

        F_{vib}^{QHA}\\big(T, V \\big) = k_B T \\sum_{\\vec{k}}
        \\sum_{i=0}^{3N}
        \\Bigg[ ln \\big(1 - e^{\\frac{h \\nu_i(\\vec{k})}{k_B T}} \\big) \\Bigg]

    with :math:`\\vec{k}` the sampled *k*-points, :math:`N` the number of
    atoms in the considered unit cell, :math:`h` the Planck constant and
    :math:`k_B` the Boltzmann constant.

    Parameters
    ----------

    temperature: ndarray
        Array of temperature values with `float` type.

    band: ndarray(ndim=3)
        Array of phonon band frequencies with `float` type.

    weights: ndarray
        Array of weights for each phonon band with `float` type.

    Returns
    -------

    result: ndarray
        Array of vibrational free energy values with `float` type.

    """
    cdef Py_ssize_t nt = temperature.shape[0]
    cdef Py_ssize_t nb = band.shape[0]
    cdef Py_ssize_t nf = band.shape[1]
    cdef Py_ssize_t nv = band.shape[2]
    cdef int i, j, k, n

    result = np.zeros( (nt,nv), dtype=np.float64 )
    cdef double[:,::1] result_view = result

    for i in prange(nt, nogil=True):
        for j in range(nb):
            for k in range(nf):
                for n in range(nv):
                    result_view[i, n] += ho_F(temperature[i], band[j,k,n]) * \
                                         weights[j]

    return result


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef free_energy(double[::1] U0, double[:,::1] Fvib):
    cdef Py_ssize_t nt = Fvib.shape[0]
    cdef Py_ssize_t nv = U0.shape[0]
    cdef int i, j

    cdef double[::1] U0_v = U0
    cdef double[:,::1] Fvib_v = Fvib

    result = np.zeros( (nt,nv), dtype=np.float64 )
    cdef double[:,::1] result_view = result


    for i in prange(nt, nogil=True):
        for j in range(nv):
            result_view[i, j] += U0_v[j] + Fvib[i,j]

    return result


cdef inline double ho_Cv(double temperature, double omega) nogil:
    """ Isochoric heat capacity of single harmonic oscillator. """
    cdef double x, e, n
    if omega <= 0. or temperature == 0.:
        return 0.
    else:
        x = H * omega / (KB * temperature)
        e = exp(x)
        n = pow((x/expm1(x)),2)
        if n == 0.:
            return 0.
        else:
            return NA * KB * e * n

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef isochoric_heat_capacity(double[::1] temperature, double[:, :, ::1] band, double[::1] weights):
    cdef Py_ssize_t nt = temperature.shape[0]
    cdef Py_ssize_t nb = band.shape[0]
    cdef Py_ssize_t nf = band.shape[1]
    cdef Py_ssize_t nv = band.shape[2]
    cdef int i, j, k, n

    result = np.zeros( (nt,nv), dtype=np.float64 )
    cdef double[:,::1] result_view = result

    for i in prange(nt, nogil=True):
        for j in range(nb):
            for k in range(nf):
                for n in range(nv):
                    result_view[i, n] += ho_Cv(temperature[i], band[j,k,n]) * \
                                         weights[j]

    return result
