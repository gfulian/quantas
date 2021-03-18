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
from scipy.integrate import quad
from scipy.constants import Avogadro as NA
from scipy.constants import Planck as h
from scipy.constants import Boltzmann as kb

class Kieffer(object):

    def __init__(self, frequencies, cutoff=1.e-10):
        """ Constructor method for the Kieffer's calculator.

        Parameters
        ----------

        frequencies: ndarray
            Array of acoustic frequencies (in Hz).
        """
        self.cutoff = cutoff
        self.acofreq = frequencies
        return

    @property
    def acofreq(self):
        """ Acoustic frequencies stored in the class.

        Returns
        -------
        ndarray(dtype=float, ndim=1)
            Array containing the acoustic frequency values (in Hz).
        """
        return self._acofreq

    @acofreq.setter
    def acofreq(self, frequencies):
        """ Acoustic frequencies stored in the class.

        Parameters
        ----------
        frequencies: ndarray(dtype=float, ndim=1)
            Array containing the acoustic frequency values (in Hz).
        """
        self._acofreq = np.asarray(frequencies)
        return

    @property
    def acofreq_exp(self):
        """ Read-only property.

        Returns
        -------
        ndarray(dtype=float, ndim=1)
            Array of the exponents without temperature.
        """
        return self._acofreq * h / kb

    def helmholtz(self, temperature):
        """ Calculate the acoustic contribution to the Helmholtz free energy
        according to the Kieffer's model.

        Parameters
        ----------

        temperature: float
            Temperature value at which the contribution is calculated.

        Return
        ------

        value: float
            Acoustic contribution to the Helmholtz free energy.
        """
        value = 0.
        for acofreq in self.acofreq_exp:
            value += self._helmholtz_integral(temperature, acofreq)
        return value

    def _helmholtz_integral(self, temperature, xnti):
        """
        """
        def helmholtz_function(x, temperature, xnti):
            xi = xnti / temperature
            num = np.power((np.arcsin(x/xi)),2)*np.log(1-np.exp(-x))
            den = np.power((np.power(xi, 2.) - np.power(x, 2.)), 0.5)
            value = num / den
            return value

        if temperature == 0.:
            return 0.
        wmin = 1.e-6
        wmax = xnti/temperature
        function  = lambda x: helmholtz_function(x, temperature, xnti)
        integral, err = quad(function, wmin, wmax, epsrel=self.cutoff)
        factor = 3 * temperature * NA * kb * np.power(2./np.pi, 3.)

        return integral * factor

    def heat_capacity(self, temperature):
        """ Calculate the acoustic contribution to the isochoric (constant
        volume) heat capacity according to the Kieffer's model.

        Parameters
        ----------

        temperature: float
            Temperature value at which the contribution is calculated.

        Return
        ------

        value: float
            Acoustic contribution to the isochoric heat capacity.
        """
        value = 0.
        for acofreq in self.acofreq_exp:
            value += self._heat_capacity_integral(temperature, acofreq)
        return value

    def _heat_capacity_integral(self, temperature, xnti):
        """
        """
        def heat_capacity_function(x, temperature, xnti):
            """ Function for the acoustic contribution to the heat
            capacity according to the Kieffer's model.

            Parameters
            ----------

            x: float
                Current value of the term hv_i/(k_B T).
            temperature: float
                Temperature value (in K).
            xnti: float
                Maximum value of the term hv_i/(k_B T).

            Return
            ------
            value: float
                Heat capacity value.
            """
            xi = xnti / temperature
            num = np.power(np.arcsin(x/xi), 2.) * np.power(x, 2.)
            num *= np.exp(x)
            den = np.power((np.power(xi, 2.) - np.power(x, 2.)), 0.5)
            den *= np.power(np.exp(x) - 1., 2.)
            value = num / den
            return value

        if temperature == 0.:
            return 0.
        wmin = 1.e-6
        wmax = xnti/temperature
        function  = lambda x: heat_capacity_function(x, temperature, xnti)
        integral, err = quad(function, wmin, wmax, epsrel=self.cutoff)

        return integral * 3 * NA * kb * np.power(2./np.pi, 3.)

    def entropy(self, temperature):
        """ Calculate the acoustic contribution to entropy according to
        the Kieffer's model.

        Parameters
        ----------

        temperature: float
            Temperature value at which the contribution is calculated.

        Return
        ------

        value: float
            Acoustic contribution to entropy.
        """
        value = 0.
        for acofreq in self.acofreq_exp:
            value += self._entropy_integral(temperature, acofreq)
        return value

    def _entropy_integral(self, temperature, xnti):
        """
        """
        def entropy_function(x, temperature, xnti):
            """ Function for the acoustic contribution to the entropy
            according to the Kieffer's model.

            Parameters
            ----------

            x: float
                Current value of the term hv_i/(k_B T).
            temperature: float
                Temperature value (in K).
            xnti: float
                Maximum value of the term hv_i/(k_B T).

            Return
            ------
            value: float
                Entropy value.
            """
            xi = xnti / temperature
            # Calculate the first addendum
            num = np.power(np.arcsin(x/xi), 2.) *x
            den = np.power((np.power(xi, 2.) - np.power(x, 2.)), 0.5)
            den *= np.power(np.exp(x) - 1., 2.)
            first_term = num / den

            # Calculate the second addendum
            num = np.power(np.arcsin(x/xi), 2.)
            den = np.power((np.power(xi, 2.) - np.power(x, 2.)), 0.5)
            second_term = num * np.log(1. - np.exp(-x)) / den

            value = first_term - second_term
            return value

        if temperature == 0.:
            return 0.
        wmin = 1.e-6
        wmax = xnti/temperature
        function  = lambda x: entropy_function(x, temperature, xnti)
        integral, err = quad(function, wmin, 0.999*wmax, epsrel=self.cutoff)

        return integral * 3 * NA * kb * np.power(2./np.pi, 3.)

