# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" This module contains a class that performs Equation of State fittings on 
experimental data.

.. note::

    Available P-V Equation of State formulations are: 
    
    - Murnaghan
    - Birch-Murnaghan (order 2, 3 and 4)
    - Natural Strain (order 2, 3 and 4)
    - Vinet (order 2 and 3)
    - Tait (order 2 and 3)
    
"""

import numpy as np
from scipy.odr import ODR, Model, RealData

class EOS(object):
    """ Generic EOS class
    """
    
    def __init__(self):
        return
    
    def eos(self, v, P, order=2, linear=False):
        """ Dummy method, it will be overwritten by each EOS class that inherits
        from this one.
        """
        return None
    
    
    def guess(self, V, P):
        """ Parabolic fit of the :math:`P(V)` data in order to provide an
        initial guess for the EOS fit.
        
        Attributes
        ----------
        V: ndarray(dtype=float, ndim=1)
            Array of volume values.
        P: ndarray(dtype=float, ndim=1)
            Array of pressure values.
            
        Returns
        -------
        guess_pars: (dtype=float, ndim=1)
            Array containing the estimates of the EOS parameters.
        """
        # Calculate V0
        V0 = np.polyval(np.polyfit(P,V,3), 0.)
        # Calculate K0
        dp = np.polyder(np.polyfit(V,P,3), m=1)
        K0 = -V0 * np.polyval(dp,V0)
        # Initial KP
        KP = 4.
        # Initial KPP
        KPP = 0.0
        guess_pars = [K0, KP, KPP, V0]
        return guess_pars
    
    
    def check(self, sigmaV, sigmaP):
        """ This method checks if there are weights for the input data. """ 
        check = 0
        try:
            len(sigmaV)
            check = check + 1
        except TypeError:
            None
        try:
            len(sigmaV)
            check = check + 2
        except TypeError:
            None
        return check
    
    
    def fit(self, 
            V, 
            P, 
            guess=None, 
            sigmaV=None, 
            sigmaP=None, 
            order=2, 
            linear=False,
            fixed=None):
        """ This method is called to perform the equation of state fit, using 
        one of the implemented formulations.
        
        Attributes
        ----------
        eostag: str
            Acronym of the EOS formulation.
        P: ndarray(dtype=float, ndim=1)
            Array of pressure values.
        V: ndarray(dtype=float, ndim=1)
            Array of volume values.
        guess: ndarray(dtype=float, ndim=1), optional
            Array of guessed EOS parameters.
        sigmaV: ndarray(dtype=float, ndim=1), optional
            Array of the uncertainties of volume values.
        sigmaP: ndarray(dtype=float, ndim=1), optional
            Array of the uncertainties of pressure values.
        order: int, optional
            Order of the equation of state formulation.
        linear: bool, optional
            Flag to fit the volume or axial data.
            
        Returns
        -------
        output: object
            scipy.ODR.Output object containing the results of the ODR routine.
        
        """
        #
        # Set which parameters will be kept fixed during the run
        try:
            len(fixed)
            fix = fixed
        except TypeError:    
            if order == 2:
                fix=[1,0,0,1]
            elif order == 3:
                fix=[1,1,0,1]
            else:
                fix=[1,1,1,1]
        #
        # Check if there is a guess
        try:
            len(guess)
            p0 = guess
        except TypeError:
            #
            # Perform a guess on the parameters
            p0 = self.guess(V, P)
        #
        # Check if we're going to perform a linearized fit
        if linear == True:
            p0[0] = p0[0]*3.
        model = Model(self.eos, extra_args=[order, linear, True])
        # Set the data for ODR
        data = RealData(V, P, sigmaV, sigmaP)
        odr = ODR(data, model, p0, ifixb=fix, maxit=500)
        if self.check(sigmaV, sigmaP) == 0 or self.check(sigmaV, sigmaP) == 2:
            odr.set_job(fit_type=2) # Ordinary least square
        else:
            odr.set_job(fit_type=0) # Explicit Orthogonal distance regression
        output = odr.run()
        return output


class Murnaghan(EOS):
    
    def __init__(self):
        super(Murnaghan, self).__init__()
        return

    def __repr__(self):
        return 'Murnaghan'

    def eos(self, p, V, order=2, linear=False, dofit=False):
        """ Murnaghan (1937) EOS:
     
        .. math::
             
            P_{VT} = \\frac{K_{0T}}{K^{\\prime}_{0T}} \\Bigg[ \\Big(
            \\frac{V_{0T}}{V} \\Big)^{K^{\\prime}_{0T}} - 1 \\Bigg]
     
        Attributes
        ----------
        p: ndarray(dtype=float, ndim=1)
            Array of EOS parameters (see below).
            
        V: ndarray(dtype=float, ndim=1)
            Array of volume values.
            
        order: int, optional
            Order of the EOS formulation.
            
        linear: bool, optional
            Flag to perform volumetric or axial EOS fits.
            
        dofit: bool, optional
            Flag used to fit the data or to perform calculations.
     
        Parameters
        ----------
        K0: float
            Bulk modulus value.
     
        Kp: float
            First derivative of the bulk modulus with respect to pressure.
            
        Kpp: float
            Second derivative of the bulk modulus with respect to pressure.
     
        V0: float
            Volume at zero pressure.
     
        Returns
        -------
        P: ndarray(dtype=float, ndim=1)
            Array of pressure value obtained from the EOS.
            
        
        **References:**
        
        F. Murnaghan, *American Journal of Mathematics* 1937, **49**, 235. 
        """
        p[2] = 0
        eta = (p[3]/V)**(p[1])
        P = (p[0]/p[1])*(eta - 1.)
        return P
    
    def normalize(self, p, V, P, sigmaV=None, sigmaP=None, order=2, 
                  linear=False):
        """ Since there is no normalized pressure expression for the 
        Murnaghan EOS, this is a dummy method.
        """        
        f_calc = np.zeros(len(V)) 
        NP_calc = np.zeros(len(V)) 
        f_obs = np.zeros(len(V)-1) 
        NP_obs = np.zeros(len(V)-1) 
        sigf = np.zeros(len(V)-1) 
        sigNP = np.zeros(len(V)-1)
                                           
        return f_calc, NP_calc, f_obs, NP_obs, sigf, sigNP
    

class BirchMurnaghan(EOS):
    
    def __init__(self):
        super(BirchMurnaghan, self).__init__()
        return

    def __repr__(self):
        return 'Birch-Murnaghan'
    
    def eos(self, p, V, order=2, linear=False, dofit=False):
        """ Third-order Birch-Murnaghan (1947) EOS:
      
        .. math::
              
            P_{VT} = 3K_{0T}f_E \\big(1 + 2f_E \\big)^{5/2} \\Bigg[
            1 + \\frac{3}{2}\\big(K^{\\prime}_{0T} - 4\\big)f_E + 
      
        .. math::
          
            \\frac{3}{2} \\Big(K_{0T}K^{\\prime \\prime}_{0T}+ \\big(K^{\\prime}_{0T} - 4\\big)
            \\big(K^{\\prime}_{0T} - 3\\big) + \\frac{35}{9} \\Big) f^2_E \\Bigg]
              
        with :math:`f_E` being:
      
        .. math::
      
            f_E = \\frac{1}{2}\\Bigg[
            \\bigg( \\frac{V_{0T}}{V} \\bigg)^{\\frac{2}{3}} -1
            \\Bigg]
      
        In the third-order formulation, the term :math:`K^{\\prime \\prime}_{0T}`
        is fixed and given by :
      
        .. math::
      
            K^{\\prime \\prime}_{0T} = - \\frac{1}{K_{0T}} \\Bigg[
            \\big( 3 - K^{\\prime}_{0T} \\big) \\big( 4 - K^{\\prime}_{0T} \\big) +
            \\frac{35}{9} \\Bigg]
      
        Attributes
        ----------
        p: ndarray(dtype=float, ndim=1)
            Array of EOS parameters (see below).
            
        V: ndarray(dtype=float, ndim=1)
            Array of volume values.
            
        order: int, optional
            Order of the EOS formulation.
            
        linear: bool, optional
            Flag to perform volumetric or axial EOS fits.
            
        dofit: bool, optional
            Flag used to fit the data or to perform calculations.
     
        Parameters
        ----------
        K0: float
            Bulk modulus value.
     
        Kp: float
            First derivative of the bulk modulus with respect to pressure.
            
        Kpp: float
            Second derivative of the bulk modulus with respect to pressure.
     
        V0: float
            Volume at zero pressure.
     
        Returns
        -------
        P: ndarray(dtype=float, ndim=1)
            Array of pressure value obtained from the EOS.
        
        
        **References:**
        
        F. Birch, *Physical Review* 1947, **71**, 809.
          
        F. Birch, *Journal of Geophysical Research* 1986, **91**, 4949.
        """
        if dofit == True:
            if order == 2 and linear == False:
                p[1] = 4.
            if order == 2 and linear == True:
                p[1] = 12.
            if order < 4:
                p[2] = np.float64(
                    (-1./p[0])*( (3.-p[1])*(4.-p[1])+(35./9.))
                    ) # implied value
            
        fE = np.float64(0.5 * (np.power((p[3]/V),2./3.) - 1.) )
        P = np.float64(
            3.*p[0]*fE*np.power((1.+2*fE),2.5) * (
                1.+1.5*(p[1]-4.)*fE + 1.5*(p[0]*p[2] + (p[1]-4.)*(p[1]-3)+
                                           35./9.)*np.power(fE,2.)
                )
            )
        
        return P
    
    def normalize(self, p, V, P, sigmaV=None, sigmaP=None, order=2, 
                  linear=False):
        """ Computes the normalized pressure for f-F plots.
        
        Attributes
        ----------
        p: ndarray(dtype=float, ndim=1)
            Array of EOS parameters (see below).
            
        V: ndarray(dtype=float, ndim=1)
            Array of volume values.
            
        order: int, optional
            Order of the EOS formulation.
            
        linear: bool, optional
            Flag to perform volumetric or axial EOS fits.

     
        Parameters
        ----------
        K0: float
            Bulk modulus value.
     
        Kp: float
            First derivative of the bulk modulus with respect to pressure.
            
        Kpp: float
            Second derivative of the bulk modulus with respect to pressure.
     
        V0: float
            Volume at zero pressure.
     
        Returns
        -------
        fE: ndarray(dtype=float, ndim=1)
            Array of strain values.
            
        NP: ndarray(dtype=float, ndim=1)
            Array of normalized pressure values.
        
        sigfE: ndarray(dtype=float, ndim=1)
            Array of errors on strain values.
            
        sigNP: ndarray(dtype=float, ndim=1)
            Array of errors on normalized pressure values.
            
        """
        # Store parameters locally 
        pars = p[:]

        volume = np.delete(np.copy(V), 0)
        press = np.delete(np.copy(P), 0)
        
        # Take into account linearity
        if order == 2 and linear == True:
            pars[1] = 4.
            pars[2] = np.float64(
                    (-1./pars[0])*( (3.-pars[1])*(4.-pars[1])+(35./9.))
                    ) # implied value
        #
        # Calculate the normalized pressure curve
        fE_calc = np.float64(0.5 * (np.power((pars[3]/V),2./3.) - 1.) )
        f = 3.*fE_calc*np.power((1.+2*fE_calc),2.5)
        NP_calc = self.eos(pars, V, order, linear) / f
        #
        # Calculate the normalized pressure data
        # 1. set reference
        V_0 = V[0]
        # 2. calculate fE
        fE_obs = np.float64(0.5 * (np.power((V_0/volume),2./3.) - 1.) )
        f = 3.*fE_obs*np.power((1.+2*fE_obs),2.5)
        # 3. normalize experimental pressure
        NP_obs = press / f
        #
        # Calculate errors in fE-F
        try:
            len(sigmaV)
            len(sigmaP)
            calc_eps = True
        except TypeError:
            calc_eps = False
        
        if calc_eps == True:
            # 1. store copy of sigmas
            sigV_0 = sigmaV[0]
            sigV = np.delete(np.copy(sigmaV), 0)
            sigP = np.delete(np.copy(sigmaP), 0)
            # 2. calculate eta and its sigma
            eta = volume/V_0
            sig_eta = (sigV / volume) + (sigV_0 / V_0)
            # 3. calculate sigma_fE
            sigfE = (1./3.)*np.power(eta,-5./3.)*sig_eta
            # 4. calculate sigma'
            sig_prime = (7.*np.power(eta,-2./3.)-5.)*sig_eta / \
            (3*eta*(1.-np.power(eta,-2./3.)))
            # 5. calculate sigma NP
            sigNP = NP_obs*np.sqrt(np.power(sigP/press,2.)+np.power(sig_prime,2.))
        else:
            sigfE = np.zeros(len(V)-1)
            sigNP = np.zeros(len(V)-1)                                
        
        return fE_calc, NP_calc, fE_obs, NP_obs, sigfE, sigNP


class NaturalStrain(EOS):
    
    def __init__(self):
        super(NaturalStrain, self).__init__()
        return

    def __repr__(self):
        return 'Natural Strain'
    
    def eos(self, p, V, order=2, linear=False, dofit=False):
        """ Natural strain (Poirier-Tarantola, 1998) EOS.
        
        .. math::
              
            P_{VT} = 3K_{0T}\\Big(\\frac{V_{0T}}{V_{PT}}\\Big)f_N \\Bigg[
            1 + af_N + bf^2_N \\Bigg]
            
        with:
        
        .. math::
        
            f_N = \\frac{1}{3} ln\\Big(\\frac{V_0}{V}\\Big)
            
        .. math::
        
            a = \\frac{3}{2} \\big(K^{\\prime}_{0T} - 2 \\big)
            
        .. math::
        
            b = \\frac{3}{2} \\Big[ 1 + K_{0T} K^{\\prime \\prime}_{0T} + 
            \\big(K^{\\prime}_{0T} -2\\big) + \\big(K^{\\prime}_{0T} -2\\big)^2
            \\Big]
            
        Second-order truncation implies a value of 
        
        .. math::
            
            K^{\\prime}_{0T} = 2
            
        whereas third-order truncation leads to an implied value of
        
        .. math::
        
            K^{\\prime \\prime}_{0T} = - \\frac{1}{K_{0T}} \\Big[ 1 + 
            \\big(K^{\\prime}_{0T} -2\\big) + \\big(K^{\\prime}_{0T} -2\\big)^2
            \\Big]
       
        Attributes
        ----------
        p: ndarray(dtype=float, ndim=1)
            Array of EOS parameters (see below).
            
        V: ndarray(dtype=float, ndim=1)
            Array of volume values.
            
        order: int, optional
            Order of the EOS formulation.
            
        linear: bool, optional
            Flag to perform volumetric or axial EOS fits.
            
        dofit: bool, optional
            Flag used to fit the data or to perform calculations.
     
        Parameters
        ----------
        K0: float
            Bulk modulus value.
     
        Kp: float
            First derivative of the bulk modulus with respect to pressure.
            
        Kpp: float
            Second derivative of the bulk modulus with respect to pressure.
     
        V0: float
            Volume at zero pressure.
     
        Returns
        -------
        P: ndarray(dtype=float, ndim=1)
            Array of pressure value obtained from the EOS.
        
        
        **References:**
        
        J. P. Poirier and A. Tarantola , *Physics of Earth and Planetary 
        Interiors* 1998, **109**, 335.
           
        """
        if dofit == True:
            if order == 2 and linear == False:
                p[1] = 2.
            if order == 2 and linear == True:
                p[1] = 6.
            if order < 4:
                p[2] = (-1./p[0])*(1.+(p[1]-2.)+np.power(p[1]-2.,2.)) # implied
            
        fN = (1./3.)*np.log(p[3]/V)
        a = 1.5 * (p[1] - 2.)
        b = 1.5 * (1. + p[2]*p[0] + (p[1]- 2.) + np.power(p[1]-2.,2.))
        
        P = (3.*p[0]*(p[3]/V)*fN) * (1. + a*fN + b*np.power(fN,2.))
        return P
    
    def normalize(self, p, V, P, sigmaV=None, sigmaP=None, order=2, 
                  linear=False):
        """ Computes the normalized pressure for f-F plots.
        
        Attributes
        ----------
        p: ndarray(dtype=float, ndim=1)
            Array of EOS parameters (see below).
            
        V: ndarray(dtype=float, ndim=1)
            Array of volume values.
            
        order: int, optional
            Order of the EOS formulation.
            
        linear: bool, optional
            Flag to perform volumetric or axial EOS fits.

     
        Parameters
        ----------
        K0: float
            Bulk modulus value.
     
        Kp: float
            First derivative of the bulk modulus with respect to pressure.
            
        Kpp: float
            Second derivative of the bulk modulus with respect to pressure.
     
        V0: float
            Volume at zero pressure.
     
        Returns
        -------
        fN: ndarray(dtype=float, ndim=1)
            Array of strain values.
            
        NP: ndarray(dtype=float, ndim=1)
            Array of normalized pressure values.
        
        sigfN: ndarray(dtype=float, ndim=1)
            Array of errors on strain values.
            
        sigNP: ndarray(dtype=float, ndim=1)
            Array of errors on normalized pressure values.
            
        """
        # Store parameters locally 
        pars = p[:]
        
        volume = np.delete(np.copy(V), 0)
        press = np.delete(np.copy(P), 0)
        
        # Take into account linearity
        if order == 2 and linear == True:
            pars[1] = 2.
            pars[2] = (-1./pars[0])*(1.+(pars[1]-2.)+np.power(pars[1]-2.,2.))
        
        # Calculate the normalized pressure curve from EoS
        fN_calc = (1./3.)*np.log(p[3]/V)
        f = 3.*(pars[3]/V)*fN_calc
        NP_calc = self.eos(pars, V, order, linear) / f
        #
        # Calculate the normalized pressure data
        # 1. set reference
        V_0 = V[0]
        # 2. calculate fE
        fN_obs = (1./3.)*np.log(V_0/volume)
        f = 3.*(V_0/volume)*fN_obs
        # 3. normalize experimental pressure
        NP_obs = press / f
        #
        # Calculate errors in f-F
        try:
            len(sigmaV)
            len(sigmaP)
            calc_eps = True
        except TypeError:
            calc_eps = False
        
        if calc_eps == True:
            # 1. store copy of sigmas
            sigV_0 = sigmaV[0]
            sigV = np.delete(np.copy(sigmaV), 0)
            sigP = np.delete(np.copy(sigmaP), 0)
            # 2. calculate eta and its sigma
            eta = volume/V_0
            sig_eta = (sigV / volume) + (sigV_0 / V_0)
            # 3. calculate sigma_fN
            sigfN = (1./3.)*np.power(eta,-1)*sig_eta
            # 4. calculate sigma'
            sig_prime = (1.-np.power(np.log(eta),-1.))*sig_eta/eta
            # 5. calculate sigma NP
            sigNP = NP_obs*np.sqrt(np.power(sigP/press,2.)+np.power(sig_prime,2.))
        else:
            sigfN = np.zeros(len(V)-1)
            sigNP = np.zeros(len(V)-1)                                
        
        return fN_calc, NP_calc, fN_obs, NP_obs, sigfN, sigNP                 


class Vinet(EOS):
    
    def __init__(self):
        super(Vinet, self).__init__()
        return

    def __repr__(self):
        return 'Vinet'    
    
    def eos(self, p, V, order=2, linear=False, dofit=False):
        """ Vinet (1986) EOS:
      
        .. math::
              
            P_{VT} = K_{0T} \\frac{3f_V}{\\big(1-f_V\\big)^2} 
            exp\\big(\\eta f_V\\big)
            
        with:
        
        .. math::
        
            f_V = 1- \\Big( \\frac{V_{PT}}{V_{0T}} \\Big)^{1/3}
            
        .. math::
        
            \\eta = \\frac{3}{2} \\big(K^{\\prime}_{0T} - 1 \\big)
      
        
      
        Attributes
        ----------
        p: ndarray(dtype=float, ndim=1)
            Array of EOS parameters (see below).
            
        V: ndarray(dtype=float, ndim=1)
            Array of volume values.
            
        order: int, optional
            Order of the EOS formulation.
            
        linear: bool, optional
            Flag to perform volumetric or axial EOS fits.
            
        dofit: bool, optional
            Flag used to fit the data or to perform calculations.
     
        Parameters
        ----------
        K0: float
            Bulk modulus value.
     
        Kp: float
            First derivative of the bulk modulus with respect to pressure.
            
        Kpp: float
            Second derivative of the bulk modulus with respect to pressure.
     
        V0: float
            Volume at zero pressure.
     
        Returns
        -------
        P: ndarray(dtype=float, ndim=1)
            Array of pressure value obtained from the EOS.
        
        
        **References:**
              
        P. Vinet, J. Ferrante, J. Rose, J. Smith, *Journal of Geophysical Research*
        1987, **92**, 9319.
        
        P. Vinet, J. Ferrante, J. Smith, J. Rose, *Journal of Physics C: Condensed Matter*
        1986, **19**, L467.
        """
        if dofit == True:
            if order == 2 and linear == False:
                p[1] = 1.
            if order == 2 and linear == True:
                p[1] = 3.

        # Kpp is always implied
        p[2] = (-1./p[0])*(np.power(p[1]/2.,2.)+(p[1]/2.)-(19./46.))
            
        fV = 1. - (V/p[3])**(1./3.)
        eta = (1.5) * (p[1] - 1.)
        P = 3.* p[0] * (fV / ((1.-fV)**2)) * np.exp(eta*fV)
        
        return P
    
    def normalize(self, p, V, P, sigmaV=None, sigmaP=None, order=2, 
                  linear=False):
        """ Computes the normalized pressure for f-F plots.
        
        Attributes
        ----------
        p: ndarray(dtype=float, ndim=1)
            Array of EOS parameters (see below).
            
        V: ndarray(dtype=float, ndim=1)
            Array of volume values.
            
        order: int, optional
            Order of the EOS formulation.
            
        linear: bool, optional
            Flag to perform volumetric or axial EOS fits.
     
        Parameters
        ----------
        K0: float
            Bulk modulus value.
     
        Kp: float
            First derivative of the bulk modulus with respect to pressure.
            
        Kpp: float
            Second derivative of the bulk modulus with respect to pressure.
     
        V0: float
            Volume at zero pressure.
     
        Returns
        -------
        fV: ndarray(dtype=float, ndim=1)
            Array of strain values.
            
        NP: ndarray(dtype=float, ndim=1)
            Array of normalized pressure values.
        
        sigfV: ndarray(dtype=float, ndim=1)
            Array of errors on strain values.
            
        sigNP: ndarray(dtype=float, ndim=1)
            Array of errors on normalized pressure values.
            
        """
        # Store parameters locally 
        pars = p[:]
        
        volume = np.delete(np.copy(V), 0)
        press = np.delete(np.copy(P), 0)
        
        # Take into account linearity
        if order == 2 and linear == True:
            pars[1] = 1.
        
        # Take into account linearity
        if order == 2 and linear == True:
            pars[1] = 2.
            pars[2] = (-1./pars[0])*(1.+(pars[1]-2.)+np.power(pars[1]-2.,2.))
        
        # Calculate the normalized pressure curve from EoS
        fV_calc = 1. - np.power((V/pars[3]),1./3.)
        eta_fV = (1.5) * (pars[1] - 1.)
        NP_calc = pars[0]*np.exp(eta_fV*fV_calc)
        #
        # Calculate the normalized pressure data
        # 1. set reference
        V_0 = V[0]
        # 2. calculate fE
        fV_obs = 1. - np.power((volume/V_0),1./3.)
        # 3. normalize experimental pressure
        NP_obs = (press * np.power((1. - fV_obs),2.))/(3.*fV_obs)
        #
        # Calculate errors in f-F
        try:
            len(sigmaV)
            len(sigmaP)
            calc_eps = True
        except TypeError:
            calc_eps = False
        
        if calc_eps == True:
            # 1. store copy of sigmas
            sigV_0 = sigmaV[0]
            sigV = np.delete(np.copy(sigmaV), 0)
            sigP = np.delete(np.copy(sigmaP), 0)
            # 2. calculate eta and its sigma
            eta = volume/V_0
            sig_eta = (sigV / volume) + (sigV_0 / V_0)
            # 3. calculate sigma_fN
            sigfV = (1./3.)*np.power(eta,-2./3.)*sig_eta
            # 4. calculate sigma NP
            sigNP = NP_obs*np.sqrt(np.power(sigP/press,2.)+\
                                   sig_eta/(3*np.power(eta,2./3.)-1.))
        else:
            sigfV = np.zeros(len(V)-1)
            sigNP = np.zeros(len(V)-1)                                
        
        return fV_calc, NP_calc, fV_obs, NP_obs, sigfV, sigNP


class Tait(EOS):
    
    def __init__(self):
        super(Tait, self).__init__()
        return

    def __repr__(self):
        return 'Tait'
    
    def eos(self, p, V, order=2, linear=False, dofit=False):
        """ Tait (Freund and Ingalls, 1986) EOS:
        
        .. math::
        
            P_{VT} = \\frac{1}{b} \\Bigg\\{ \\Bigg[
            \\frac{\\big(V_{PT}/V_{0T}\\big) + a - 1 }{a}
             \\Bigg]^{-1/c} - 1 \\Bigg\\}
             
        with:
        
        .. math::
        
            a = \\frac{1 + K^{\\prime}_{0T}}{1 + K^{\\prime}_{0T} + 
            K_{0T} K^{\\prime \\prime}_{0T}}
        
        .. math::
        
            b = \\frac{K^{\\prime}_{0T}}{K_{0T}} - 
            \\frac{K^{\\prime \\prime}_{0T}}{1 + K^{\\prime}_{0T}}
        
        .. math::
            
            c =  \\frac{1+K^{\\prime}_{0T}+K_{0T} K^{\\prime \\prime}_{0T}}
            {\\big(K^{\\prime}_{0T} \\big)^2 + K^{\\prime}_{0T} - 
            K_{0T} K^{\\prime \\prime}_{0T}} 
            
        
        Attributes
        ----------
        p: ndarray(dtype=float, ndim=1)
            Array of EOS parameters (see below).
            
        V: ndarray(dtype=float, ndim=1)
            Array of volume values.
            
        order: int, optional
            Order of the EOS formulation.
            
        linear: bool, optional
            Flag to perform volumetric or axial EOS fits.
            
        dofit: bool, optional
            Flag used to fit the data or to perform calculations.
     
        Parameters
        ----------
        K0: float
            Bulk modulus value.
     
        Kp: float
            First derivative of the bulk modulus with respect to pressure.
            
        Kpp: float
            Second derivative of the bulk modulus with respect to pressure.
     
        V0: float
            Volume at zero pressure.
     
        Returns
        -------
        P: ndarray(dtype=float, ndim=1)
            Array of pressure value obtained from the EOS
        
        
        
        **References:**
        
        J. Freund, R. Ingalls, *Journal of Physics & Chemistry of Solids* 1989,
        **50**, 263.
        
        Y. K. Huang, C. Y. Chow, *Journal of Physics D: Applied Physics* 1974,
        **7**, 2021.
        
        """
        if dofit == True:
            if order == 2 and linear == False:
                p[1] = 4.
            if order == 2 and linear == True:
                p[1] = 12.
            if order < 4:
                p[2] = -p[1]/p[0]

        a = (1.+p[1])/(1.+p[1]+p[0]*p[2]) 
        b = (p[1]/p[0])-(p[2]/(1.+p[1]))
        c = (1.+p[1]+p[0]*p[2])/(np.power(p[1],2.)+p[1]-p[0]*p[2])
        P = (1/b)*( np.power(((V/p[3])+a-1.)/a, -1./c) - 1.)
        
        return P
    
    def normalize(self, p, V, P, sigmaV=None, sigmaP=None, order=2, 
                  linear=False):
        """ Since there is no normalized pressure expression for the 
        Tait EOS, this is a dummy method.
        """
        f_calc = np.zeros(len(V)) 
        NP_calc = np.zeros(len(V)) 
        f_obs = np.zeros(len(V)-1) 
        NP_obs = np.zeros(len(V)-1) 
        sigf = np.zeros(len(V)-1) 
        sigNP = np.zeros(len(V)-1)
        return f_calc, NP_calc, f_obs, NP_obs, sigf, sigNP
