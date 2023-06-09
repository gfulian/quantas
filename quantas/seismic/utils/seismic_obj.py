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
This module contains the Seismic class with the methods used to calculate the
seismic wave velocities in a crystalline material. Some other auxiliary
functions are also included.
"""

import numpy as np
import itertools as it


def vector(theta, phi):
    """ Return the Cartesian components of a unit vector from polar coordinates.

    Parameters
    ----------

    theta: float
        Theta angle.

    phi: float
        Phi angle

    Returns
    -------

    ndarray(dtype=float, ndim=1)
        Array with the Cartesian components of the vector.
    """
    a = np.sin(theta)*np.cos(phi)
    b = np.sin(theta)*np.sin(phi)
    c = np.cos(theta)
    return np.asarray([ a, b, c ])


def cofactor(A):
    """
    Return the cofactor matrix of a 3x3 matrix. This function relies primarily
    on NumPy but when the matrix is singular (determinant equal to zero). In
    that case, an element-wise approach is adopted. 

    Parameters
    ----------

    A: ndarray(dtype=floar, ndim=2)
        3x3 matrix.

    Returns
    -------

    ndarray(dtype=float, ndim=2)
       3x3 cofactor matrix of matrix m.
       
    """
    cofactor = np.empty((3, 3))
    try:
        # Check if the matrix is singular
        if np.linalg.det(A) != 0:
            cofactor = np.linalg.inv(A).T * np.linalg.det(A)
            return cofactor
        else:
            # If the matrix is singular, use element-wise math
            cofactor[0][0] = A[1][1]*A[2][2] - A[1][2]*A[2][1]
            cofactor[0][1] = A[1][2]*A[2][0] - A[1][0]*A[2][2]
            cofactor[0][2] = A[1][0]*A[2][1] - A[1][1]*A[2][0]

            cofactor[1][0] = A[0][2]*A[2][1] - A[0][1]*A[2][2]
            cofactor[1][1] = A[0][0]*A[2][2] - A[0][2]*A[2][0]
            cofactor[1][2] = A[0][1]*A[2][0] - A[0][0]*A[2][1]

            cofactor[2][0] = A[0][1]*A[1][2] - A[0][2]*A[1][1]
            cofactor[2][1] = A[0][2]*A[1][0] - A[0][0]*A[1][2]
            cofactor[2][2] = A[0][0]*A[1][1] - A[0][1]*A[1][0]
        return cofactor
    except Exception as e:
        print('Could not find cofactor matrix due to ', e)


class Seismic(object):
    """ This class holds all the information and methods to analyze the tensor
    matrix given as input.

    Parameters
    ----------

    matrix: ndarray(dtype=float, ndim=2)
        2D-array (6x6) containing the values of the symmetric SOEC matrix, in
        Voigt's notation.

    density: float
        Density of the material, expressed in kg m^-3.

    """

    direction = None
    theta = None
    phi = None
    M = None
    _grad_M = None
    _eigenval = None
    _eigenvec = None
    _grad_eigenval = None
    _hessian_eig = None
    _phase_velocity = None
    _group_velocity = None
    _group_abs = None
    _group_dir = None
    _group_theta = None
    _group_phi = None
    _powerflow_angle = None
    _enhancement = None
    
    def __init__(self, matrix, density):
        """ Constructor method for the Christoffel object. """
        # Store the symmetric 6x6 elastic moduli matrix in Voigt's notation
        self.C = matrix

        # Store the density
        self.density = density

        # Set the compliance tensor matrix
        try:
            self.S = np.linalg.inv(self.C)
        except:
            raise ValueError("Matrix is singular, aborting operation")

        # Set the compliance matrix in its (3x3x3x3) tensorial form
        S_tmp = map(self.S_ijkl, it.product([0,1,2], repeat=4))
        self.Smat = np.array(list(S_tmp)).reshape(3,3,3,3)

        # Set the stiffness matrix in its (3x3x3x3) tensorial form
        C_tmp = map(self.C_ijkl, it.product([0,1,2], repeat=4))
        self.Cmat = np.array(list(C_tmp)).reshape(3,3,3,3)
        
        # Convert the values for wave velocities calculation
        self.Cmat *= 1000.0/density
        self.hessM = self.hessian_christoffel_matrix(self.Cmat)
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
        """ Read-only property
        Return the elastic constants in Voigt's notation.
        """
        return self.C

    @property
    def compliance(self):
        """ Read-only property
        Return the compliance matrix in Voigt's notation.
        """
        return self.S

    @property
    def christoffel(self):
        """ Read-only property
        Return the Christoffel matrix M.
        """
        return self.M

    @property
    def q(self):
        """ Read-only property
        Return the current direction vector q.
        """
        return self.direction

    @property
    def bulk_modulus(self):
        """ Return the Voigt-Reuss-Hill (VRH) average of the bulk modulus,
        calculated according to the following formulas:

        .. math::

           K_V = \\frac{1}{9} \\big[ C_{11} + C_{22} + C_{33} + 2\\big( C_{12} + 
           C_{13} + C_{23} \\big) \\big]
        
        .. math::
        
           K_R = \\big[ S_{11} + S_{22} + S_{33} + 2\\big( S_{12} + S_{13} + 
           S_{23} \\big) \\big]^{-1}

        .. math::

           K_{VRH} = \\frac{K_V + K_R}{2}


        Returns
        -------

        float
            The average between the Voigt and Reuss bulk moduli. 
        
        """
        # Voigt bound
        KV = np.einsum('iijj', self.Cmat*self.density/1000.)/9.
        # Reuss bound
        KR = 1/np.einsum('iijj', self.Smat)

        return (KV + KR) / 2.

    @property
    def shear_modulus(self):
        """ Return the Voigt-Reuss-Hill (VRH) average of the shear modulus, calculated
        according to the following formulas:

        .. math::
        
           \\mu_V = \\frac{1}{15} \\big[ C_{11} + C_{22} + C_{33} + 3\\big( C_{44} +
           C_{55} + C_{66} \\big) - \\big( C_{11} + C_{13} + C_{23} \\big) \\big]
        
        .. math::
        
           \\mu_R = \\frac{15}{4} \\Big[ S_{11} + S_{22} + S_{33} - 
           \\big( S_{11} + S_{13} + S_{23} \\big) + 3\\big( S_{44} + S_{55} + 
           S_{66} \\big) \\Big]^{-1}
        
        .. math::
        
           \\mu_{VRH} = \\frac{\\mu_V + \\mu_R}{2}


        Returns
        -------

        float
            The average between the Voigt and Reuss shear moduli. 
        
        """
        # Voigt bound
        GV = (self.C[0][0] + self.C[1][1] + self.C[2][2] -
              (self.C[0][1] + self.C[0][2] + self.C[1][2]) +
              3*(self.C[3][3] + self.C[4][4] + self.C[5][5])
              ) / 15.
        # Reuss bound
        GR = (15.) / (4.*(self.S[0][0] + self.S[1][1] + self.S[2][2]) -
                      4.*(self.S[0][1] + self.S[0][2] + self.S[1][2]) +
                      3.*(self.S[3][3] + self.S[4][4] + self.S[5][5])
                         )
        
        return 0.5 * (GV + GR)

    @property
    def isotropic_velocity(self):
        """ Return the sound velocities for an isotropic (polycrystalline)
        material.

        Returns
        -------

        ndarray(dtype=float, ndim=2)
            The array containing the isotropic seismic (sound) velocities.
        
        """
        K = self.bulk_modulus
        G = self.shear_modulus
        # Calculate the primary velocity
        V_P = np.sqrt(1000.0*(K + 4.0*G/3)/self.density)
        # Calculate the secondary velocity
        V_S = np.sqrt(1000.0*G/self.density)
        return np.asarray([V_S, V_S, V_P])

    @property
    def iso_S(self):
        """ Return the isotropic secondary seismic velocity.

        Returns
        -------

        float
            The isotropic secondary seismic velocity.
        """
        return self.isotropic_velocity[0]

    @property
    def iso_P(self):
        """ Return the isotropic primary seismic velocity.

        Returns
        -------

        float
            The isotropic secondary seismic velocity.
        """
        return self.isotropic_velocity[2]

    def hessian_christoffel_matrix(self, C):
        """
        Return the Hessian of the dynamical matrix, which is independent of q
        by definition.

        .. math::

           H_{ijkl} = \\frac{\\partial^2 M_{kl}}{\\partial x_i \\partial x_j}

        Parameters
        ----------

        C: ndarray
            :math:`3 \\times 3 \\times 3 \\times` tensor of the elastic moduli
            with `float` type.

        Returns
        -------
        hessianmat: ndarray
            :math:`3 \\times 3 \\times 3 \\times 3` Hessian of the dynamical
            matrix

        """
        hessM = np.empty((3, 3, 3, 3))
        for i, j, k, l in it.product(np.arange(3, dtype=int), repeat=4):
            hessM[i][j][k][l] = C[k][i][j][l] + C[k][j][i][l]
        return hessM

    def clear(self):
        """Clear all data."""
        self.direction = None
        self.theta = None
        self.phi = None
        self.M = None
        self._grad_M = None
        self._eigenval = None
        self._eigenvec = None
        self._grad_eigenval = None
        self._hessian_eig = None
        self._phase_velocity = None
        self._group_velocity = None
        self._group_abs = None
        self._group_dir = None
        self._group_theta = None
        self._group_phi = None
        self._powerflow_angle = None
        self._enhancement = None
        return

    def spherical_direction(self, theta, phi):
        """ Define a wave vector in spherical coordinates (in radians).
        Each time a new direction is set, the previously stored data are
        deleted.

        Parameters
        ----------

        theta: float
            Polar angle.

        phi: float
            Azimuth angle.

        """
        self.clear()

        self.direction = vector(theta, phi)
        self.theta = theta
        self.phi = phi

        # Set the new Christoffel matrix
        self.M = np.dot(self.direction, np.dot(self.direction, self.Cmat))
        return

    def cartesian_direction(self, vector):
        """ Define a wave vector in Cartesian coordinates. Each time a new
        direction is set, the previously stored data are deleted.

        Parameters
        ----------

        vector: ndarray
            Vector defining the direction along which the x-axis will be
            aligned.

        """
        self.clear()

        # Normalize the vector
        self.direction = vector / np.linalg.norm(vector)

        # Calculate the theta and phi angles
        x = self.direction[0]
        y = self.direction[1]
        z = self.direction[2]
        self.phi = np.arctan2(x, y)
        self.theta = np.arctan2(np.sqrt(y ** 2 + x ** 2), z)

        # Set the new Christoffel matrix
        self.M = np.dot(self.direction, np.dot(self.direction, self.Cmat))
        return

    @property
    def grad_M(self):
        """ Return the gradient of the Christoffel's matrix. """
        if type(self._grad_M) == type(None):
            self._grad_M = np.transpose(
                np.dot(self.direction,
                       self.Cmat + np.transpose(self.Cmat, (0, 2, 1, 3))),
                (1, 0, 2))
        return self._grad_M

    def solve(self):
        """ Calculate the eigenvalues and eigenvectors of the Christoffel's
        matrix M. The values are then sorted in ascending order and stored.
        """
        eigenval, eigenvec = np.linalg.eigh(self.M)
        sorting = np.argsort(eigenval)
        self._eigenval = eigenval[sorting]
        self._eigenvec = eigenvec.T[sorting]
        return

    @property
    def eigenval(self):
        """ Read-only property
        Return the eigenvalues of the Christoffel's matrix.
        """
        if type(self._eigenval) == type(None):
            self.solve()
        return self._eigenval

    @property
    def eigenvec(self):
        """ Read-only property
        Return the eigenvectors of the Christoffel's matrix.
        """
        if type(self._eigenval) == type(None):
            self.solve()
        return self._eigenvec

    @property
    def phase_velocity(self):
        """ Read-only property
        Return the phase velocity from the eigenvalues of the Christoffel's
        matrix.
        """
        if type(self._eigenval) == type(None):
            self.solve()
        return np.sign(self._eigenval)*np.sqrt(np.absolute(self._eigenval))

    @property
    def rel_phase_velocity(self):
        """ Read-only property
        Return the relative phase velocity from the eigenvalues of the
        Christoffel's matrix.
        """
        if type(self._eigenval) == type(None):
            self.solve()
        return self.phase_velocity / self.isotropic_velocity

    @property
    def hessian_eigenval(self):
        """ Read-only property
        Return the hessian of the eigenvalues of the Christoffel's matrix,
        calculated as:

        .. math::

           H_{nij} = \\frac{\\partial^2 \\lambda_n}{\\partial x_i \\partial x_j}
           
        """
        hessian_eig = np.zeros((3, 3, 3))
        for i in range(3):
            hessian_eig[i] += np.dot(
                np.dot(self.hessM, self.eigenvec[i]), self.eigenvec[i])
            pinv = np.linalg.pinv(
                self.eigenval[i]*np.identity(3)-self.M,
                rcond=1e-10)
            dvec = np.dot(self.grad_M, self.eigenvec[i])
            hessian_eig[i] += 2. * np.dot(np.dot(dvec, pinv), dvec.T)
        return hessian_eig

    def solve_group(self):
        """ Calculate the group seismic velocities as a gradient of the phase
        velocities.
        """
        # Collect data
        pvel = self.phase_velocity
        evec = self.eigenvec
        gmat = self.grad_M

        # Set arrays for results
        gvel = np.empty((3, 3))  # Group velocities
        geig = np.empty((3, 3))  # Gradient of the M eigenvectors
        gabs = np.empty((3))  # Group absolute velocities
        gdir = np.empty((3, 3))  # Group directions
        gtheta = np.empty((3))  # Group theta angle
        gphi = np.empty((3))  # Group phi angle

        # loop over slow secondary, fast secondary, primary
        for i in range(3):
            # loop over cartesian components
            for j in range(3):
                geig[i, j] = np.dot(evec[i], np.dot(gmat[j], evec[i]))
                # Calculate group velocities
                gvel[i, j] = geig[i, j] / (2.*pvel[i])
            # Calculate group absolute velocities
            gabs[i] = np.linalg.norm(gvel[i])
            # Calculate the group velocity direction
            gdir[i] = gvel[i] / gabs[i]

            x = gdir[i][0]
            y = gdir[i][1]
            z = gdir[i][2]
            gphi[i] = np.arctan2(x, y)
            gtheta[i] = np.arctan2(np.sqrt(y**2 + x**2), z)

        # Calculate the powerflow angle
        cos_powerflow_angle = np.dot(gdir, self.direction)
        
        # Store data
        self._group_velocity = gvel
        self._grad_eigenval = geig
        self._group_abs = gabs
        self._group_dir = gdir
        self._group_theta = gtheta
        self._group_phi = gphi
        self._powerflow_angle = np.arccos(np.around(cos_powerflow_angle, 10))
        return

    @property
    def group_velocity(self):
        """ Read-only property
        Return the group velocity.
        """
        if type(self._group_velocity) == type(None):
            self.solve_group()
        return self._group_velocity

    @property
    def grad_eigenval(self):
        """ Read-only property
        Return the gradient of the eigenvalues of the Christoffel's matrix.
        """
        if type(self._grad_eigenval) == type(None):
            self.solve_group()
        return self._grad_eigenval

    @property
    def group_abs(self):
        """ Read-only property
        Return the absolute values of the group velocities.
        """
        if type(self._group_abs) == type(None):
            self.solve_group()
        return self._group_abs

    @property
    def rel_group_velocity(self):
        """ Read-only property
        Return the relative group velocities.
        """
        if type(self._group_abs) == type(None):
            self.solve_group()
        return self.group_abs / self.isotropic_velocity

    @property
    def group_dir(self):
        if type(self._group_dir) == type(None):
            self.solve_group()
        return self._group_dir

    @property
    def group_theta(self):
        if type(self._group_theta) == type(None):
            self.solve_group()
        return self._group_theta

    @property
    def group_phi(self):
        if type(self._group_phi) == type(None):
            self.solve_group()
        return self._group_phi

    @property
    def powerflow_angle(self):
        """ Read-only property
        Return the powerflow angle.
        """
        if type(self._powerflow_angle) == type(None):
            self.solve_group()
        return self._powerflow_angle

    @property
    def enhancement(self):
        """ Read-only property
        Return the absolute values of the group velocities.
        """
        # Collect data
        hessian = self.hessian_eigenval
        pvel = self.phase_velocity
        gvel = self.group_velocity
        gabs = self.group_abs

        grad = np.empty((3, 3, 3))
        enhancement = np.empty(3)

        for i in range(3):
            grad[i] = hessian[i] / gabs[i]
            grad[i] -= np.outer(gvel[i], np.dot(hessian[i], gvel[i]))  \
                       / (gabs[i]**3)
            grad[i] /= 2. * pvel[i]
            
            enhancement[i] = 1. / np.linalg.norm(
                np.dot(cofactor(grad[i]), self.direction))

        self._enhancement = enhancement
        return self._enhancement

