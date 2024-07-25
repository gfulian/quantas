# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

import math
import numpy as np
from numpy import linalg
from scipy import optimize
import itertools as it
from click import progressbar

def vector(theta, phi, chi=None):
    """
    Define a unit vector using polar coordinates.

    Parameters
    ----------

    theta: float
        First angle of the vector.

    phi: float
        Second angle of the vector.

    chi: float, optional
        Third angle of the vector (defines the measurement vector).

    Returns
    -------

    list
        List containing the vector components in Cartesian coordinates.

    """
    if chi is None:
        a = np.sin(theta)*np.cos(phi)
        b = np.sin(theta)*np.sin(phi)
        c = np.cos(theta)
    else:
        a = np.cos(theta)*np.cos(phi)*np.cos(chi) - np.sin(phi)*np.sin(chi)
        b = np.cos(theta)*np.sin(phi)*np.cos(chi) + np.cos(phi)*np.sin(chi)
        c = - np.sin(theta)*np.cos(chi)
    return [ a, b, c ]


class SOEC(object):
    """
    This class holds all the information and methods to analyze the tensor
    matrix given as input.

    This object is a forked version of the ELATE_ scripts of Gaillac and
    co-workers. [1]_

    Parameters
    ----------

    matrix: ndarray
        :math:`6 \\times 6` matrix containing the values of the symmetric
        elastic moduli in Voigt's notation with `float` type.

    density: float, optional
        Density of the material (in kg m^-3) with `float` type.


    .. _ELATE: http://progs.coudert.name/elate

    .. rubric:: References

    .. [1] Gaillac, R.; Pullumbi, P.; Coudert, F. X. Journal of
           Physics-Condensed Matter 2016, 28, 275201.

    """

    def __init__(self, matrix, density=None):
        """ Constructor method """
        # Store the symmetric 6x6 SOEC matrix
        self.C = matrix

        # Store the density
        self.density = density

        # Set the compliance tensor matrix
        try:
            self.S = linalg.inv(self.C)
        except:
            raise ValueError('matrix is singular')

        # Set the compliance matrix in its (3x3x3x3) tensorial form
        S_tmp = map(self.S_ijkl, it.product([0,1,2], repeat=4))
        self.Smat = np.array(list(S_tmp)).reshape(3,3,3,3)

        # Set the stiffness matrix in its (3x3x3x3) tensorial form
        C_tmp = map(self.C_ijkl, it.product([0,1,2], repeat=4))
        self.Cmat = np.array(list(C_tmp)).reshape(3,3,3,3)
        
        return

    def C_ijkl(self, iters):
        """
        This method calculates the stiffness tensor element from the given
        indexes.
        
        Arguments
        ---------
        iter: tuple
            Tuple containing the four indexes of the (3x3x3x3) C tensor.
            
        Returns
        -------
        C_ijkl: float
            ijkl-th element of the stiffness tensor.

        """
        i, j, k, l = iters
        # Voigt Matrix
        vm = [[0, 5, 4], [5, 1, 3], [4, 3, 2]]
        # Compliance coefficients
        return self.C[vm[i][j]][vm[k][l]]

    def S_ijkl(self, iters):
        """
        This method calculates the compliance tensor element from the given
        indexes.

        Parameters
        ----------
        iter: tuple
            Tuple containing the four indexes of the (3x3x3x3) S tensor.

        Returns
        -------
        float
            :math:`ijkl`-th element of the compliance tensor.

        """
        i, j, k, l = iters
        # Voigt Matrix
        vm = [[0, 5, 4], [5, 1, 3], [4, 3, 2]]
        # Compliance coefficients
        def sc(p,q): return 1. / ((1+p//3)*(1+q//3))
        return sc(vm[i][j], vm[k][l]) * self.S[vm[i][j]][vm[k][l]]

    @property
    def stiffness(self):
        """
        Return the elastic moduli in Voigt's notation.
        """
        return self.C

    @property
    def compliance(self):
        """
        Return the compliance tensor of the elastic constants.
        """
        return self.S

    def young_modulus(self, angles):
        """
        Calculate the Young's modulus along a specified direction.

        Parameters
        ----------

        angles: tuple
            Tuple containing the :math:`\\theta` and :math:`\\phi`
            angular values.

        Returns
        -------

        float
            Young's modulus along the direction specified by theta and phi.
        """
        a = vector(angles[0], angles[1])
        r = 0
        for i, j, k, l in it.product((0, 1, 2), repeat=4):
            r += self.Smat[i, j, k, l] * a[i] * a[j] * a[k] * a[l]
        return 1/r

    def linear_compressibility(self, angles):
        """
        Calculate the linear compressibility (LC, :math:`\\beta`) along a
        specified direction.

        Parameters
        ----------

        angles: tuple
            Tuple containing the :math:`\\theta` and :math:`\\phi`
            angular values.

        Returns
        -------

        float
             linear compressibility along the specified direction.
        """
        a = vector(angles[0], angles[1])
        r = 0
        for i, j, k in it.product((0, 1, 2), repeat=3):
            r += a[i]*a[j] * self.Smat[i][j][k][k]
        return 1000 * r

    def shear_modulus(self, angles):
        """
        Calculate the shear modulus along a specified direction. Three angles
        are necessary and two vectors (direction of the shear modulus and
        direction of measurement).

        Parameters
        ----------

        angles: tuple
            Tuple containing the :math:`\\theta`, :math:`\\phi` and
            :math:`\\chi` angular values.

        Returns
        -------

        float
             Shear modulus along the specified direction.
        """
        a = vector(angles[0], angles[1])
        b = vector(angles[0], angles[1], angles[2])
        r = 0
        for i, j, k, l in it.product((0, 1, 2), repeat=4):
            r += a[i]*b[j]*a[k]*b[l] * self.Smat[i][j][k][l]
        return 1/(4*r)

    def poisson_ratio(self, x):
        """
        Calculate the Poisson's ratio along a specified direction.Three angles
        are necessary and two vectors (direction of the shear modulus and
        direction of measurement).

        Parameters
        ----------

        angles: tuple
            Tuple containing the :math:`\\theta`, :math:`\\phi` and
            :math:`\\chi` angular values.

        Returns
        -------

        Nu: float
             Poisson's ratio along the specified direction.
        """
        a = vector(x[0], x[1])
        b = vector(x[0], x[1], x[2])
        r1 = 0
        r2 = 0
        for i, j, k, l in it.product((0, 1, 2), repeat=4):
            r1 += a[i]*a[j]*b[k]*b[l] * self.Smat[i][j][k][l]
            r2 += a[i]*a[j]*a[k]*a[l] * self.Smat[i][j][k][l]
        return -r1/r2

    def averages(self):
        """
        Calculate the elastic averages according to the Voigt-Reuss-Hill method:

        .. math::

           K_V = \\frac{1}{9} \\big[ C_{11} + C_{22} + C_{33} + 2\\big( C_{12} +
           C_{13} + C_{23} \\big) \\big]

        .. math::

           K_R = \\big[ S_{11} + S_{22} + S_{33} + 2\\big( S_{12} + S_{13} +
           S_{23} \\big) \\big]^{-1}

        .. math::

           \\mu_V = \\frac{1}{15} \\big[ C_{11} + C_{22} + C_{33} + 3\\big( C_{44} +
           C_{55} + C_{66} \\big) - \\big( C_{11} + C_{13} + C_{23} \\big) \\big]

        .. math::

           \\mu_R = \\frac{15}{4} \\Big[ S_{11} + S_{22} + S_{33} -
           \\big( S_{11} + S_{13} + S_{23} \\big) + 3\\big( S_{44} + S_{55} +
           S_{66} \\big) \\Big]^{-1}

        .. math::

           K_{VRH} = \\frac{K_V + K_R}{2}

        .. math::

           \\nu_{VRH} = \\frac{\\mu_V + \\mu_R}{2}

        .. math::

           E_{VRH} = \\frac{9 K_{VRH} \\mu_{VRH}}{3K_{VRH} + \\mu_{VRH}}


        Returns
        -------
        avg: list
             Averages of the elastic behaviour of the material.
        """
        # Bulk modulus
        KV = np.einsum('iijj', self.Cmat)/9.
        KR = 1/np.einsum('iijj', self.Smat)
        KH = 0.5 * (KV + KR)

        # Shear modulus
        GV = (self.C[0][0] + self.C[1][1] + self.C[2][2] -
              (self.C[0][1] + self.C[0][2] + self.C[1][2]) +
              3*(self.C[3][3] + self.C[4][4] + self.C[5][5])
              ) / 15.
        GR = (15.) / (4.*(self.S[0][0] + self.S[1][1] + self.S[2][2]) -
                      4.*(self.S[0][1] + self.S[0][2] + self.S[1][2]) +
                      3.*(self.S[3][3] + self.S[4][4] + self.S[5][5])
                         )
        GH = 0.5 * (GV + GR)

        return [
            [KV, 1/(1/(3*GV) + 1/(9*KV)), GV, (1 - 3*GV/(3*KV+GV))/2],
            [KR, 1/(1/(3*GR) + 1/(9*KR)), GR, (1 - 3*GR/(3*KR+GR))/2],
            [KH, 1/(1/(3*GH) + 1/(9*KH)), GH, (1 - 3*GH/(3*KH+GH))/2]
            ]

    def isotropic_velocities(self):
        """
        Seismic velocities as if the material was isotropic.
        """
        B = self.averages()[2][0] # VRH average
        G = self.averages()[2][2] # VRH average
        # Primary velocity
        Vp = np.sqrt(1000.0*(B + 4.0*G/3)/self.density)
        # Secondary velocity
        Vs = np.sqrt(1000.0*G/self.density)
        return np.array([Vs, Vs, Vp])

    def shear_2d(self, x):
        ftol = 0.001
        xtol = 0.01
        def func1(z): return self.shear_modulus([x[0], x[1], z])
        r1 = optimize.minimize(func1, np.pi/2.0, args=(), method = 'Powell', options={'xtol':xtol, 'ftol':ftol})#, bounds=[(0.0,np.pi)])
        def func2(z): return -self.shear_modulus([x[0], x[1], z])
        r2 = optimize.minimize(func2, np.pi/2.0, args=(), method = 'Powell', options={'xtol':xtol, 'ftol':ftol})#, bounds=[(0.0,np.pi)])
        return (float(r1.fun), -float(r2.fun))

    def poisson_2d(self, x):
        ftol = 0.001
        xtol = 0.01
        def func1(z): return self.poisson_ratio([x[0], x[1], z])
        r1 = optimize.minimize(func1, np.pi/2.0, args=(), method = 'Powell', options={'xtol':xtol, 'ftol':ftol})#, bounds=[(0.0,np.pi)])
        def func2(z): return -self.poisson_ratio([x[0], x[1], z])
        r2 = optimize.minimize(func2, np.pi/2.0, args=(), method = 'Powell', options={'xtol':xtol, 'ftol':ftol})#, bounds=[(0.0,np.pi)])
        return (min(0,float(r1.fun)), max(0,float(r1.fun)), -float(r2.fun))

    def phase_velocity(self, angles):
        """
        Calculate the phase velocities (longitudinal and shear) by solving
        the Christoffel's equation.

        Parameters
        ----------

        angles: tuple
            Tuple containing the :math:`\\theta` and :math:`\\phi`
            angular values.

        Returns
        -------

        velocities: ndarray
            Array containing the :math:`v_{s1}`, :math:`v_{s2}` and
            :math:`v_{l}` velocities (in ascending order).

        """
        # Convert GPa to Pa
        C = self.C * 1.e9
        # Create the divergence operator
        l = vector(angles[0], angles[1])
        divop = np.zeros((3,6), dtype=np.float64)
        divop[0][0] = l[0]  # l_x
        divop[1][1] = l[1]  # l_y
        divop[2][2] = l[2]  # l_z
        divop[0][4] = l[2]  # l_z
        divop[1][3] = l[2]  # l_z
        divop[0][5] = l[1]  # l_y
        divop[2][3] = l[1]  # l_y
        divop[2][4] = l[0]  # l_x
        divop[1][5] = l[0]  # l_x
        divop_T = divop.transpose()
        mult = np.matmul(C,divop_T)
        A = np.matmul(divop,mult)
        #
        # Diagonalize the matrix A to obtain velocities and vectors
        velocities, _ = np.linalg.eigh(A)
        velocities = np.sqrt(velocities/self.density)/1000. # km s^-1
        return velocities

    def longitudinal_wave(self, angles):
        """
        Commodity method that returns the longitudinal phase velocity
        (:math:`v_{l}`).
        """
        vel = self.phase_velocity(angles)
        return vel[2]

    def shear_wave_1(self, angles):
        """
        Commodity method that returns the slow shear phase velocity
        (:math:`v_{s1}`).
        """
        vel = self.phase_velocity(angles)
        return vel[0]

    def shear_wave_2(self, angles):
        """
        Commodity method that returns the fast shear phase velocity
        (:math:`v_{s2}`).
        """
        vel = self.phase_velocity(angles)
        return vel[1]

    def polar_young(self, theta, phi):
        result = np.zeros(len(theta),dtype=np.float64)
        with progressbar(theta, width=50) as progress:
            i = 0
            for c in progress:
                result[i] = self.young_modulus([theta[i],phi[i]])
                i += 1
        return result

    def polar_compressibility(self, theta, phi):
        result = np.zeros((len(theta),2),dtype=np.float64)
        with progressbar(theta, width=50) as progress:
            i = 0
            for c in progress:
                angles = [theta[i], phi[i]]
                result[i][0] = max(0, self.linear_compressibility(angles))
                result[i][1] = max(0, -self.linear_compressibility(angles))
                i += 1
        return result

    def polar_shear(self, theta, phi):
        result = np.zeros((len(theta),2),dtype=np.float64)
        with progressbar(theta, width=50) as progress:
            i = 0
            for c in progress:
                tmp = self.shear_2d([theta[i], phi[i]])
                result[i][0] =tmp[0]
                result[i][1] =tmp[1]
                i += 1
        return result

    def polar_poisson(self, theta, phi):
        result = np.zeros((len(theta),3),dtype=np.float64)
        with progressbar(theta, width=50) as progress:
            i = 0
            for c in progress:
                tmp = self.poisson_2d([theta[i], phi[i]])
                result[i][0] =tmp[0]
                result[i][1] =tmp[1]
                result[i][2] =tmp[2]
                i += 1
        return result

    def polar_waves(self, theta, phi):
        result = np.zeros((len(theta),3),dtype=np.float64)
        with progressbar(theta, width=50) as progress:
            i = 0
            for c in progress:
                tmp = self.phase_velocity([theta[i], phi[i]])
                result[i][0] =tmp[0]
                result[i][1] =tmp[1]
                result[i][2] =tmp[2]
                i += 1
        return result


class SOECOrtho(SOEC):
    """
    An elastic tensor, for the specific case of an orthorhombic system.

    """

    def __init__(self, arg):
        """Initialize from a matrix, or from an Elastic object"""
        if isinstance(arg, SOEC):
            self.C = arg.C
            self.S = arg.S
            self.Smat = arg.Smat
            self.Cmat = arg.Cmat
            self.density = arg.density
        else:
            raise TypeError('SOECOrtho constructor argument should be a SOEC object')

    def young_modulus(self, angles):
        ct2 = math.cos(angles[0])**2
        st2 = 1 - ct2
        cf2 = math.cos(angles[1])**2
        sf2 = 1 - cf2
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s44 = 4 * self.Smat[1][2][1][2]
        s55 = 4 * self.Smat[0][2][0][2]
        s66 = 4 * self.Smat[0][1][0][1]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
        return 1/(ct2**2*s33 + 2*cf2*ct2*s13*st2 + cf2*ct2*s55*st2 + 2*ct2*s23*sf2*st2 + ct2*s44*sf2*st2 + cf2**2*s11*st2**2 + 2*cf2*s12*sf2*st2**2 + cf2*s66*sf2*st2**2 + s22*sf2**2*st2**2)

    def linear_compressibility(self, angles):
        ct2 = math.cos(angles[0])**2
        cf2 = math.cos(angles[1])**2
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
        return 1000 * (ct2 * (s13 + s23 + s33) + (cf2 * (s11 + s12 + s13) + (s12 + s22 + s23) * (1 - cf2)) * (1 - ct2))

    def shear_modulus(self, angles):
        ct = math.cos(angles[0])
        ct2 = ct*ct
        st2 = 1 - ct2
        cf = math.cos(angles[1])
        sf = math.sin(angles[1])
        sf2 = sf*sf
        cx = math.cos(angles[2])
        cx2 = cx*cx
        sx = math.sin(angles[2])
        sx2 = 1 - cx2
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s44 = 4 * self.Smat[1][2][1][2]
        s55 = 4 * self.Smat[0][2][0][2]
        s66 = 4 * self.Smat[0][1][0][1]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
        r = (
                    ct2*ct2*cx2*s44*sf2 + cx2*s44*sf2*st2*st2 + 4*cf**3*ct*cx*(-2*s11 + 2*s12 + s66)*sf*st2*sx
                    + 2*cf*ct*cx*sf*(ct2*(s44 - s55) + (4*s13 - 4*s23 - s44 + s55 - 4*s12*sf2 + 4*s22*sf2 - 2*s66*sf2)*st2)*sx
                    + s66*sf2*sf2*st2*sx2 + cf**4*st2*(4*ct2*cx2*s11 + s66*sx2)
                    + ct2*(2*cx2*(2*s33 + sf2*(-4*s23 - s44 + 2*s22*sf2))*st2 + s55*sf2*sx2)
                    + cf**2*(ct2*ct2*cx2*s55 + ct2*(-2*cx2*(4*s13 + s55 - 2*(2*s12 + s66)*sf2)*st2 + s44*sx2)
                                     + st2*(cx2*s55*st2 + 2*(2*s11 - 4*s12 + 2*s22 - s66)*sf2*sx2))
                )
        return 1/r

    def poisson_ratio(self, angles):
        ct = math.cos(angles[0])
        ct2 = ct*ct
        st2 = 1 - ct2
        cf = math.cos(angles[1])
        sf = math.sin(angles[1])
        cx = math.cos(angles[2])
        sx = math.sin(angles[2])
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s44 = 4 * self.Smat[1][2][1][2]
        s55 = 4 * self.Smat[0][2][0][2]
        s66 = 4 * self.Smat[0][1][0][1]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
        return (
    (-(ct**2*cx**2*s33*st2) - cf**2*cx**2*s13*st2*st2 - cx**2*s23*sf**2*st2*st2 + ct*cx*s44*sf*st2*(ct*cx*sf + cf*sx) -
                    ct**2*s23*(ct*cx*sf + cf*sx)**2 - cf**2*s12*st2*(ct*cx*sf + cf*sx)**2 - s22*sf**2*st2*(ct*cx*sf + cf*sx)**2 +
                    cf*ct*cx*s55*st2*(cf*ct*cx - sf*sx) - cf*s66*sf*st2*(ct*cx*sf + cf*sx)*(cf*ct*cx - sf*sx) -
                    ct**2*s13*(cf*ct*cx - sf*sx)**2 - cf**2*s11*st2*(cf*ct*cx - sf*sx)**2 - s12*sf**2*st2*(cf*ct*cx - sf*sx)**2)/
                (ct**4*s33 + 2*cf**2*ct**2*s13*st2 + cf**2*ct**2*s55*st2 + 2*ct**2*s23*sf**2*st2 + ct**2*s44*sf**2*st2 +
                    cf**4*s11*st2*st2 + 2*cf**2*s12*sf**2*st2*st2 + cf**2*s66*sf**2*st2*st2 + s22*sf**4*st2*st2)
        )
