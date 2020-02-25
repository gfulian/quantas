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
This module contains some volume-integrated equation of state (EoS)
formulations that can be used to fit :math:`E(V)` data.

..note:

    :math:`E(V)` could be either electronic (static) energy or Helmholtz free
    energy

"""
import numpy as np
from scipy.optimize import curve_fit

from quantas.eosfit.utils.pv_eos import Murnaghan
from quantas.eosfit.utils.pv_eos import BirchMurnaghan
from quantas.eosfit.utils.pv_eos import NaturalStrain
from quantas.eosfit.utils.pv_eos import Vinet


class EnergyEOS(object):
    """
    This class holds all the volume-integrated equations of state used by
    Quantas.

    """

    def __init__(self):
        """ Constructor method."""
        return

    def _set_eos_type(self, eos):
        eos_values = {
            'murnaghan': 'self.murnaghan',
            'birchmurnaghan': 'self.birchmurnaghan',
            'vinet': 'self.vinet',
            'pouriertarantola': 'self.pouriertarantola'
            }
        chosen_eos = None
        for name, string in eos_values.items():
            if eos == name:
                chosen_eos = string
        return chosen_eos

    def _set_pv_eos_type(self, eos):
        eos_values = {
            'murnaghan': Murnaghan,
            'birchmurnaghan': BirchMurnaghan,
            'vinet': Vinet,
            'pouriertarantola': NaturalStrain
            }
        chosen_eos = None
        for name, string in eos_values.items():
            if eos == name:
                chosen_eos = string
        return chosen_eos

    def fit(self, eostag, V, E):
        """ This method is called to perform the equation of state fit, using
        one of the implemented formulations.

        Parameters
        ----------
        eostag: str
            Acronym of the EOS formulation.
        V: ndarray(dtype=float, ndim=1)
            Array of volume values.
        E: ndarray(dtype=float, ndim=1)
            Array of energy values.

        Returns
        -------
        eos_parameters: ndarray
            Array of the fitted EOS parameters .


        .. note::

            Available EOS formulations (:code:`eostag`) are: Murnaghan (M),
            third-order Birch-Murnaghan (BM3), Vinet (V) and Poirier-Tarantola
            (PT).
        """
        # Set the EOS function
        eos = self._set_eos_type(eostag)
        # Perform a guess on the parameters
        p0 = self.guess(V, E)
        # Perform the EOS fit using the guess as initial parameters
        eos_parameters, eos_cov = curve_fit(
            eval(eos), V, E, p0, ftol=1.e-15, xtol=1.e-15)
        # EoS errors
        eos_errors = np.sqrt(np.diag(eos_cov))
        return eos_parameters, eos_errors

    def guess(self, V, E):
        """ Parabolic fit of the :math:`E(V)` data in order to provide an
        initial guess for the EOS fit.

        .. math::

            E = aV^2 +bV + c

        With the fitting parameters, it is possible to obtain:

        .. math::

            V_0 = -b/\\big(2a\\big)

        .. math::

            E_0 = aV_0^2 +bV_0 + c

        .. math::

            K_0 = 2aV_0

        .. math::

            K^{\\prime} = 4

        Parameters
        ----------

        V: ndarray
            Array of volume values with `float` type.

        E: ndarray
            Array of energy values with `float` type.

        Returns
        -------

        guess_pars: ndarray
            Array containing the initial guess of the EoS parameters with
            `float` type.

        """
        pars = np.polyfit(V, E, 4)
        first_derivative = np.polyder(pars,1)
        second_derivative = np.polyder(pars,2)
        V_roots = np.roots(first_derivative)
        Vs = V_roots * np.conj(V_roots)
        root_minimum = np.argmin(Vs)
        V0 = np.real(V_roots[root_minimum])
        E0 = np.polyval(pars, V0)
        K0 = V0 * np.polyval(second_derivative, V0)
        guess_pars = [E0, K0, 4., V0]
        return guess_pars

    def murnaghan(self, V, E0, K0, KP, V0):
        """ Volume-integrated Murnaghan EOS formulation:

        .. math::

            E = E_0 + K_0 \\frac{V}{K^{\\prime}} \\Bigg[
            \\frac{\\big(V_0/V\\big)^{K^{\\prime}}}{{K^{\\prime}} -1} + 1
            \\Bigg] - V_0 \\frac{K_0}{{K^{\\prime}} - 1}

        Reference: Fu, C. L. and Ho, K. M., *Phys. Rev. B*, 1983, **28**, 5480.

        Parameters
        ----------
        V: ndarray
            Array of unit cell volumes with `float` type.

        E0: float
            Energy of the unit cell, parameter of the EoS.

        K0: float
            Bulk modulus, parameter of the EoS.

        KP: float
            First-derivative of bulk modulus, parameter of the EoS.

        V0: float
            Unit cell volume, parameter of the EoS.

        Returns
        -------

        E: ndarray
            Array of energy values related to the provided unit cell volumes
            as obtained from the EoS.
        """
        E = E0 + K0*V/KP*(((V0/V)**KP)/(KP-1)+1) - V0*K0/(KP-1)
        return E

    def birchmurnaghan(self, V, E0, K0, KP, V0):
        """ Volume-integrated 3rd-order Birch-Murnaghan EOS formulation:

        .. math::

            E = E_0 + K_0 V_0 \\frac{9}{16} \\Big\\{ K^{\\prime} \\Big(
            \\eta^2 - 1 \\Big)^{3} + \\Big[ \\big( \\eta^{2} -1 \\big)^{2}
            \\big(6 - 4 \\eta^{2} \\big) \\Big] \\Big\\}

        .. math::

            \\eta = \\Bigg( \\frac{V_0}{V} \\Bigg)^{\\frac{1}{3}}

        Reference: Hebbache , M and Zemzemi, M., *Phys. Rev. B*, 2004, **70**,
        224107.

        Parameters
        ----------
        V: ndarray
            Array of unit cell volumes with `float` type.

        E0: float
            Energy of the unit cell, parameter of the EoS.

        K0: float
            Bulk modulus, parameter of the EoS.

        KP: float
            First-derivative of bulk modulus, parameter of the EoS.

        V0: float
            Unit cell volume, parameter of the EoS.

        Returns
        -------

        E: ndarray
            Array of energy values related to the provided unit cell volumes
            as obtained from the EoS.
        """
        eta = (V0/V)**(1./3.)
        E = E0 + (9.*K0*V0/16.) * (KP*(eta**2-1)**3 + ((eta**2-1)**2)
                                   * (6. - 4.*eta**2))
        return E

    def pouriertarantola(self, V, E0, K0, KP, V0):
        """ Volume-integrated Pourier-Tarantola EOS formulation:

        .. math::

            E = E_0 + \\frac{B_0 V_0 \\rho^2}{6} \\bigg(3 + \\rho
            \\big(K^{\\prime}-2\\big) \\bigg)

        .. math::

            \\eta = \\Bigg( \\frac{V_0}{V} \\Bigg)^{\\frac{1}{3}}

        .. math::

            \\rho = -3 ln(\\eta)

        Reference: Hebbache , M and Zemzemi, M., *Phys. Rev. B*, 2004, **70**,
        224107.

        Parameters
        ----------
        V: ndarray
            Array of unit cell volumes with `float` type.

        E0: float
            Energy of the unit cell, parameter of the EoS.

        K0: float
            Bulk modulus, parameter of the EoS.

        KP: float
            First-derivative of bulk modulus, parameter of the EoS.

        V0: float
            Unit cell volume, parameter of the EoS.

        Returns
        -------

        E: ndarray
            Array of energy values related to the provided unit cell volumes
            as obtained from the EoS.
        """
        eta = (V/V0)**(1./3.)
        rho = -3.*np.log(eta)
        E = E0 + (K0*V0*rho**2)/6.*(3. + rho*(KP - 2))
        return E

    def vinet(self, V, E0, K0, KP, V0):
        """Volume-integrated Vinet EOS formulation:

        .. math::

            E = E_0 + 2 \\frac{ K_0 V_0}{\\big(K^{\\prime} - 1)}^2 \\Big\\{
            2 - \\bigg[ 5 + 3K^{\\prime} \\big(\\eta - 1\\big) \\bigg]
            e^{-3\\big(K^{\\prime}-1\\big)\\big(\\eta -1\\big) / 2}

        .. math::

            \\eta = \\Bigg( \\frac{V_0}{V} \\Bigg)^{\\frac{1}{3}}

        Reference: Hebbache , M and Zemzemi, M., *Phys. Rev. B*, 2004, **70**,
        224107.

        Parameters
        ----------
        V: ndarray
            Array of unit cell volumes with `float` type.

        E0: float
            Energy of the unit cell, parameter of the EoS.

        K0: float
            Bulk modulus, parameter of the EoS.

        KP: float
            First-derivative of bulk modulus, parameter of the EoS.

        V0: float
            Unit cell volume, parameter of the EoS.

        Returns
        -------

        E: ndarray
            Array of energy values related to the provided unit cell volumes
            as obtained from the EoS.
        """
        eta = (V/V0)**(1./3.)
        E = (E0 + 2.*K0*V0/(KP-1.)**2
             * (2. - (5. + 3. * KP * (eta-1.)-3. * eta) *
                np.exp(-3. * (KP-1.) * (eta-1.)/2.)))
        return E

    def pressure(self, eostag, ev_pars, V):
        """ This method is called to calculate the pressure from P-V
        formulations of the equation of state.

        Parameters
        ----------
        eostag: str
            Acronym of the EOS formulation.
        ev_pars: ndarray(dtype=float, ndim=1)
            Array of the E-V EoS parameters.
        V: ndarray(dtype=float, ndim=1)
            Array of volume values.

        ..note::
            The bulk modulus parameter has to be provided in units of
            pressure. It won't be converted here!

        Returns
        -------
        pressures: ndarray(dtype=float, ndim=1)
            Array of the pressure values calculated from fitted
            EOS parameters.

        """
        # Set the EOS function
        eos = self._set_pv_eos_type(eostag)()
        #
        # Set P-V EoS parameters according to:
        #  - p[0] = K0
        #  - p[1] = K'
        #  - p[2] = K'' (usually set to zero)
        #  - p[3] = V0
        #
        pv_pars = np.zeros(4, dtype=float)
        pv_pars[0] = ev_pars[1]
        pv_pars[1] = ev_pars[2]
        pv_pars[3] = ev_pars[3]
        # Calculate pressures
        pressures = eos.eos(pv_pars, V, order=3)
        return pressures
