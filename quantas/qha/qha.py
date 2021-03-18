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
This module performs the QHA calculations necessary to obtain
thermodynamic and thermomelastic data on a single solid phase at selected
pressure (P) and temperature (T) conditions.

Required input data
-------------------
volume: ndarray
    Array containing unit cell volumes with float type.
temperature: ndarray
    Array containing temperature values with float type.
pressure: ndarray
    Array containing temperature values with float type.
static_energy: ndarray
    Array containing static energy (T = 0 K) related to each volume with
    float type.
phonons: ndarray(dtype=float, ndim=2)
    Array of phonon frequencies calculated at each volume with float type.

Output data
-----------
Uzp: ndarray(dtype=float, ndim=2)
    Zero-point energy at selected P-T conditions.
Uth: ndarray(dtype=float, ndim=2)
    Thermal contribution to internal energy at selected P-T conditions.
Utot: ndarray(dtype=float, ndim=2)
    Total energy (U0 + Uzp + Uth) at selected P-T conditions.
S: ndarray(dtype=float, ndim=2)
    Entropy at selected P-T conditions.
F: ndarray(dtype=float, ndim=2)
    Helmholtz free energy at selected P-T conditions.
Fvib: ndarray(dtype=float, ndim=2)
    Vibrational contribution to free energy at selected P-T conditions.
Cv: ndarray(dtype=float, ndim=2)
    Isochoric heat capacity at selected P-T conditions.
Cp: ndarray(dtype=float, ndim=2)
    Isobaric heat capacity at selected P-T conditions.
Cp-Cv: ndarray(dtype=float, ndim=2)
    Difference between Cp and Cv at selected P-T conditions.
KT: ndarray(dtype=float, ndim=2)
    Isothermal bulk modulus at selected P-T conditions.
Kp: ndarray(dtype=float, ndim=2)
    Pressure-derivative of bulk modulus at selected P-T conditions.
KS: ndarray(dtype=float, ndim=2)
    Adiabatic bulk modulus at selected P-T conditions.
alphaV: ndarray(dtype=float, ndim=2)
    Coefficient of volumetric thermal expansion at selected P-T conditions.
"""
from time import thread_time as clock
import numpy as np
import h5py

# Physical constants
from scipy.constants import Avogadro as N
from scipy.constants import physical_constants as pc

# Basic claculator
from quantas.core.calculator import BasicCalculator

# Unit conversions
from quantas.utils.physics.units import convert_frequency
from quantas.utils.physics.units import convert_energy
from quantas.utils.physics.units import convert_temperature
from quantas.utils.physics.units import convert_pressure
from quantas.utils.physics.units import convert_volume
from quantas.utils.physics.units import energy_to_pressure
from quantas.utils.physics.units import pressure_to_energy

# Flags
from quantas.utils.flags import Flag

# File reader
from quantas.IO.qha_reader import QHAInputFileReader

# Progress bar
from click import progressbar

# Thermodynamic functions
from quantas.utils.physics.statistical_mechanics import zero_point_energy
from quantas.utils.physics.statistical_mechanics import thermal_energy
from quantas.utils.physics.statistical_mechanics import internal_energy
from quantas.utils.physics.statistical_mechanics import entropy
from quantas.utils.physics.statistical_mechanics import vibrational_free_energy
from quantas.utils.physics.statistical_mechanics import isochoric_heat_capacity
from quantas.utils.physics.statistical_mechanics import free_energy
from quantas.utils.physics.thermodynamics import enthalpy
from quantas.utils.physics.thermodynamics import gibbs
from quantas.utils.physics.thermodynamics import adiabatic_bulk_modulus
from quantas.utils.physics.thermodynamics import gruneisen_parameter

# Equation of State formulations for energy(volume) fits
from .utils.ev_eos import EnergyEOS
from .utils.kieffer import Kieffer

# Other utilities for fitting purposes
from quantas.utils.math.polynomials import polyfit
from quantas.utils.math.polynomials import R_squared
from quantas.utils.math.polynomials import interpolate
from quantas.utils.math.polynomials import find_polynomial_minimum

np.warnings.filterwarnings('ignore')


class QHACalculator(BasicCalculator):
    """ Quasi-Harmonic Approximation class derived from the
    quantas.core.calculator.BasicCalculator class.

    Parameters
    ----------

    settings: quantas.cui.settings.Settings
        Settings for the current run.

    """

    descriptions = {
        'Uzp': 'Zero-point energy',
        'Uth': 'Thermal energy',
        'Utot': 'Total (vibrational + static) internal energy',
        'S': 'Entropy',
        'Fvib': 'Vibrational Helmholtz free energy',
        'F': 'Helmholtz free energy',
        'Cv': 'Isochoric heat capacity',
        'VT': 'Unit cell volume',
        'Cp': 'Isobaric heat capacity',
        'Cp-Cv': 'Correction for isobaric heat capacity',
        'KT': 'Isothermal bulk modulus',
        'KS': 'Adiabatic bulk modulus',
        'Kp': 'First-derivative of bulk modulus',
        'alphaV': 'Volumetric coefficient of thermal expansion',
        'H': 'Enthalpy',
        'G': 'Gibbs free energy',
        'gruneisen': 'Grüneisen parameters',
        }

    def __init__(self, settings):
        """ Constructor for the Quasi-Harmonic Approximation calculator. """
        BasicCalculator.__init__(self)

        self.silent = settings['silent']
        self.debug = settings['debug']
        self.log = settings['logfile']
        
        self.temperature = settings['trange']
        self.pressure = settings['prange']
        # flags
        self.frequency_interpolation = Flag()
        self.thermodynamic_interpolation = Flag()
        self.polynomial_minimization = Flag()
        self.eos_minimization = Flag()
        self.acoustic_kieffer = Flag()

        if settings['qha_scheme'] == 'td':
            self.thermodynamic_interpolation.on()
        elif settings['qha_scheme'] == 'freq':
            self.frequency_interpolation.on()

        if settings['minimization'] == 'eos':
            self.eos_minimization.on()
        elif settings['minimization'] == 'poly':
            self.polynomial_minimization.on()
        self.eos_formulation = settings['eos_function']

        if settings['kieffer']:
            self.acoustic_kieffer.on()

        self._edeg = settings['edeg']
        self._fdeg = settings['fdeg']
        self._eunit = settings['energy_unit']
        self._funit = settings['frequency_unit']
        self._vunit = settings['lenght_unit']
        self._tunit = settings['temperature_unit']
        self._punit = settings['pressure_unit']

        return

    def init_results(self, nt, npress):
        """
        This method initiates the result dictionary with zeros.

        Parameters
        ----------

        nt: int
            Number of temperature values.

        np: int
            Number of pressure values.

        """
        self._results = {
            'Uzp': np.zeros((nt, npress), dtype=np.float64),
            'Uth': np.zeros((nt, npress), dtype=np.float64),
            'Utot': np.zeros((nt, npress), dtype=np.float64),
            'S': np.zeros((nt, npress), dtype=np.float64),
            'F': np.zeros((nt, npress), dtype=np.float64),
            'Fvib': np.zeros((nt, npress), dtype=np.float64),
            'Cv': np.zeros((nt, npress), dtype=np.float64),
            'Cp-Cv': np.zeros((nt, npress), dtype=np.float64),
            'Cp': np.zeros((nt, npress), dtype=np.float64),
            'VT': np.zeros((nt, npress), dtype=np.float64),
            'KT': np.zeros((nt, npress), dtype=np.float64),
            'KS': np.zeros((nt, npress), dtype=np.float64),
            'Kp': np.zeros((nt, npress), dtype=np.float64),
            'alphaV': np.zeros((nt, npress), dtype=np.float64),
            'H': np.zeros((nt, npress), dtype=np.float64),
            'G': np.zeros((nt, npress), dtype=np.float64),
            'gruneisen': np.zeros((nt, npress), dtype=np.float64),
            }
        return

    @property
    def volume(self):
        """ Read-only property.

        Returns
        -------
        ndarray(dtype=float, ndim=1)
            Array containing the input volume values.
        """
        return self.input.volume

    @property
    def static_energy(self):
        """ Read-only property.

        Returns
        -------
        ndarray(dtype=float, ndim=1)
            Array containing the input static (electronic) energy values.
        """
        return self.input.energy

    @property
    def phonons(self):
        """ Read-only property.

        Returns
        -------
        ndarray(dtype=float, ndim=2)
            Array containing the input frequency values at each volume.
        """
        return self.input.frequencies

    @property
    def acoustic(self):
        """ Read-only property.

        Returns
        -------
        ndarray(dtype=float, ndim=2)
            Array containing the acoustic frequency values at each volume.
        """
        return self.input.acoustic_kieffer

    @property
    def pressure(self):
        """ Read-only property.

        Returns
        -------
        ndarray(dtype=float, ndim=1)
            Array containing the target pressure values.
        """
        return self._pressure

    @pressure.setter
    def pressure(self, press):
        """ Set the pressure range for the current calculation.

        Parameters
        ----------

        press: tuple
            Tuple containing the minimum, maximum and increment of pressure
            values with float type.

        """
        self._pressure = np.arange(press[0], press[1]+press[2], press[2])
        return

    @property
    def temperature(self):
        """ Read-only property.

        Returns
        -------
        ndarray(dtype=float, ndim=1)
            Array containing the target temperature values.
        """
        return self._temperature

    @temperature.setter
    def temperature(self, temp):
        """ Set the temperature range for the current calculation.

        Parameters
        ----------

        temp: tuple
            Tuple containing the minimum, maximum and increment of temperature
            values with float type.

        """
        self._temperature = np.arange(temp[0], temp[1]+temp[2], temp[2])
        return

    def read_input(self, filename):
        """ This method read the provided input file and return a dictionary
        containing the necessary data.

        Parameters
        ----------
        filename: str
            Complete path of the inputfile.

        Returns
        -------
        dict
            Dictionary containing the input data.

        """
        self.echo('Reading input file: {0}'.format(filename))
        self.echo_debug('Instantiate file reader')
        data = QHAInputFileReader(filename)
        

        if not data.completed:
            return data.error
        if len(data.volume) < 4:
            error = 'Insufficient number of volumes for QHA, at least' +\
                    ' 4 points are required'
            return error

        self.input = data
        return

    def report_input_data(self):
        """ This method write to the output stream the main information about
        the input.
        """
        idata = self.input
        elabel = 'Energy (' + self._eunit + ')'
        vlabel = 'Volume (' + self._vunit + '^3)'
        plabel = 'Pressure (' + self._punit + ')'
        
        self.echo('')
        self.echo('Job: {0}'.format(idata.jobname))
        self.echo('')
        self.echo('System:')
        self.echo('- {0:40} {1}'.format('Number of volumes', idata.nvol))
        self.echo('- {0:40} {1}'.format('Number of atoms', idata.natoms))
        self.echo('- {0:40} {1}'.format(
            'Number of sampled k-points', idata.kpoints))
        self.echo(
            '- {0:40} {1}'.format(
                'Number of frequencies',
                int(idata.natoms * 3 * idata.total_q_points)
                )
            )
        self.echo('')
        self.echo('Volume and energy values:\n')
        if self.eos_minimization.is_on():
            eos_pars, eos_errs = self.find_eos_minimum(
                self.volume, self.static_energy, 0.)
            eos = EnergyEOS()
            eos_pars[1] = energy_to_pressure(
                eos_pars[1], self._eunit, self._vunit, self._punit)
            eos_errs[1] = energy_to_pressure(
                eos_errs[1], self._eunit, self._vunit, self._punit)
            pressures = eos.pressure(
                self.eos_formulation, eos_pars, self.volume)
        else:
            pressures = -energy_to_pressure(
                np.gradient(self.static_energy, self.volume),
                self._eunit, self._vunit, self._punit)
        self.echo(' {:15} {: ^15} {: ^30}'.format(plabel, vlabel, elabel))
        self.echo(' {:-^15} {:-^15} {:-^30}'.format('', '', ''))
        for i in range(len(self.volume)):
            self.echo(' {: 15.3f} {: ^15.6f} {: ^30.12e}'.format(
                pressures[i], self.volume[i], self.static_energy[i]))

        if self.eos_minimization.is_on():
            self.echo('')
            self.echo('EoS fitting parameters for static energy vs volume data:')
            for i in range(len(eos_pars)):
                val = eos_pars[i]
                err = eos_errs[i]
                if i == 0:
                    label = 'E0'
                if i == 1:
                    label = 'K0'
                if i == 2:
                    label = "K'"
                if i == 3:
                    label = 'V0'
                try:
                    magnitude = np.log10(err)
                    if magnitude < 0.:
                        magnitude = -int(np.floor(magnitude))
                    else:
                        magnitude = int(magnitude)
                    red_err = err*np.power(10., magnitude)
                    if np.log10(int(np.round(red_err,0))) >= 1.:
                        magnitude -= 1
                        red_err = err*np.power(10., magnitude)
                    if magnitude > 0:
                        red_val = np.round(val, magnitude)
                    else:
                        red_val = int(np.round(val, magnitude))
                    red_err = int(np.round(red_err,0))
                    self.echo('{} = {}({})'.format(label, red_val, red_err))
                except OverflowError:
                    self.echo('{} = {}({})'.format(label, val, err))
        self.echo('')
        # Check if Kieffer's approach was requested, even in the case of
        # supercell. If yes, throw a warning.
        if idata.kpoints > 1:
            if self.acoustic_kieffer.is_on():
                text = "WARNING! The Kieffer's method for acoustic "
                text += 'frequencies was requested.\n'
                text += 'However, it seems that more than Gamma-point phonons '
                text += 'are reported in \n'
                text += 'the provided input file.\n'
                text += 'Please, check the consistency of the final results.'
                self.echo_warning(text)
        self.echo('')
        return

    def run(self):
        """ Start the quasi-harmonic calculations. """
        self.init_results(self.temperature.shape[0], self.pressure.shape[0])

        t0 = clock()
        
        self.echo('{0}{1:-^78}{0}'.format(
            '#', 'Quasi-Harmonic Approximation calculation started'))
        self.echo('')
        #
        n = self.temperature.shape[0]
        m = self.pressure.shape[0]
        # Conversions
        T = convert_temperature(self.temperature, self._tunit, 'K')
        #
        grid, idx = self.fine_grid()
        self.echo_debug('Grid expansion set')
        #
        phonons = convert_frequency(self.phonons, self._funit, 'Hz')
        if self.acoustic_kieffer.is_on():
            if type(self.acoustic) == type(None):
                self.error = 'Acoustic frequency missing in input'
                return
            acoustic = convert_frequency(self.acoustic, self._funit, 'Hz')
        else:
            acoustic = None
        if self.frequency_interpolation.is_on():
            self.echo_debug('Switch to frequency interpolation mode')
            self.run_frequency_interpolation(T, phonons, grid, idx, acoustic)
        if self.thermodynamic_interpolation.is_on():
            self.echo_debug('Switch to thermodynamics interpolation mode')
            self.run_thermodynamics_interpolation(T, phonons,grid, idx,
                                                  acoustic
                                                  )
        self._finalize_run()
        
        self.echo('')
        self.echo('Total wall time: {0:6.2f} sec'.format(clock() - t0))
        self.echo('{0}{1:-^78}{0}'.format('#', 'QHA Calculation ended'))
        self.completed = True
        return

    def run_frequency_interpolation(self, T, phonons, grid, idx, acoustic):
        """
        Calculate QHA thermodynamic/thermoelastic properties using the
        interpolation of phonon frequencies scheme.

        Parameters
        ----------

        T: ndarray
            Array of temperature values with `float` type.

        phonons: ndarray
            Array of phonon frequencies with `float` type.

        grid: ndarray
            Array containing the fine grid with `float` type.

        idx: int
            Index of the central point of the grid, which represents the
            point of interest whose values will be stored.

        """
        #
        # Fit frequency values
        #
        t0_ffit = clock()
        self.echo('')
        self.echo(' - Fitting frequency using polynomial of order {0}'.format(
            self._fdeg))
        fpars = np.zeros(
            (phonons.shape[0], phonons.shape[1], self._fdeg+1),
            dtype=np.float64)
        for i in range(phonons.shape[0]):
            fpars[i] = polyfit(self.volume, phonons[i], self._fdeg)
        fR2 = np.zeros(
            (phonons.shape[0], phonons.shape[1]), dtype=np.float64
            )
        for i in range(phonons.shape[0]):
            fR2[i] = R_squared(self.volume, phonons[i], fpars[i])
        self.__report_frequency_fit(phonons, fR2)
        self.echo_time(t0_ffit, clock())
        #
        # Fit static energy values
        #
        energy_pars = polyfit(self.volume, self.static_energy, self._edeg)
        #
        # Calculate free energy F(V,T)
        #
        weights = self.input.weights / self.input.total_q_points
        self.echo(' - Calculation of harmonic Helmoltz free energy F(V,T)')
        t0_F = clock()
        F_ha = free_energy(
            self.static_energy, convert_energy(
                vibrational_free_energy(T, phonons, weights),
                'kjmol', self._eunit)
            )
        self.echo_time(t0_F, clock())
        #
        # Volume minimization
        #
        t0 = clock()
        self.echo_override(
            ' - Volume minimization using {0}:'.format(
                'polynomials' if self.polynomial_minimization.is_on()
                else 'EoS')
            )
        with progressbar(self.temperature, width=50) as progress:
            i = 0
            for c in progress:
                for j in range(self.pressure.shape[0]):
                    VT, KT, Kp = self.minimize_volume(
                        self.volume, F_ha[i], self.pressure[j])
                    self._results['VT'][i, j] = VT
                    if not self.polynomial_minimization.is_on():
                        self._results['KT'][i, j] = KT
                        self._results['Kp'][i, j] = Kp
                i += 1
        self.echo_time(t0, clock())
        #
        # Calculation of QHA data
        #
        t0 = clock()
        self.echo_override(
            ' - Calculation of quasi-harmonic thermodynamic properties')
        with progressbar(self.temperature, width=50) as progress:
            i = 0
            for c in progress:
                for j in range(self.pressure.shape[0]):
                    #
                    # Create grid of points near minimum volume
                    #
                    volume = self._results['VT'][i, j]
                    volume_grid = volume + (volume * grid)
                    #
                    # Interpolate frequency values
                    #
                    phon_grid = interpolate(volume_grid, fpars)
                    #
                    # Interpolate static energy values
                    #
                    U0_grid = interpolate(volume_grid, energy_pars)
                    #
                    # Calculate thermodynamics
                    #
                    Uzp_grid = convert_energy(
                        zero_point_energy(np.asarray([0.]), phon_grid, weights),
                        'kjmol', self._eunit)[0]
                    Uth_grid = convert_energy(
                        thermal_energy(np.asarray([T[i]]), phon_grid, weights),
                        'kjmol', self._eunit)[0]
                    S_grid = convert_energy(
                        entropy(np.asarray([T[i]]), phon_grid, weights),
                        'kjmol', self._eunit)[0]
                    Cv_grid = convert_energy(
                        isochoric_heat_capacity(np.asarray([T[i]]), phon_grid,
                                                weights),
                        'kjmol', self._eunit)[0]
                    Fvib_grid = convert_energy(
                        vibrational_free_energy(np.asarray([T[i]]), phon_grid,
                                                weights),
                        'kjmol', self._eunit)[0]
                    U_grid = U0_grid + Uzp_grid + Uth_grid
                    F_grid = U0_grid + Fvib_grid
                    #
                    # Calculate the bulk modulus and its pressure derivative if polynomial
                    # minimization was selected
                    #
                    if self.polynomial_minimization.is_on():
                        KT, Kp = self.numerical_bulk_modulus(
                            volume_grid, F_grid, idx)
                        self._results['KT'][i, j] = KT
                        self._results['Kp'][i, j] = Kp
                    #
                    # Store data
                    #
                    self._results['Uzp'][i, j] = Uzp_grid[idx]
                    self._results['Uth'][i, j] = Uth_grid[idx]
                    self._results['S'][i, j] = S_grid[idx]
                    self._results['Cv'][i, j] = Cv_grid[idx]
                    self._results['Fvib'][i, j] = Fvib_grid[idx]
                    self._results['Utot'][i, j] = U_grid[idx]
                    self._results['F'][i, j] = F_grid[idx]
                    #
                i += 1
        self.echo_time(t0, clock())
        return

    def run_thermodynamics_interpolation(self, T, phonons, grid, idx, acoustic):
        """
        Calculate QHA thermodynamic/thermoelastic properties using the
        interpolation of thermodynamics scheme.

        Parameters
        ----------

        T: ndarray
            Array of temperature values with `float` type.

        phonons: ndarray
            Array of phonon frequencies with `float` type.

        grid: ndarray
            Array containing the fine grid with `float` type.

        idx: int
            Index of the central point of the grid, which represents the
            point of interest whose values will be stored.

        """
        #
        # Calculate thermodynamics within harmonic approximation (HA)
        #
        t0 = clock()
        weights = self.input.weights / self.input.total_q_points
        self.echo(' - Calculation of harmonic thermodynamic properties:')
        # Zero-point internal energy
        t0 = clock()
        t0_ha = t0
        Uzp_ha = convert_energy(
            zero_point_energy(np.zeros(1, dtype=float), phonons, weights),
            'kjmol', self._eunit
            )
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * zero-point energy', 1000*(clock()-t0)))
        # thermal internal energy
        t0 = clock()
        Uth_ha = convert_energy(thermal_energy(T, phonons, weights),
                                'kjmol', self._eunit)
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * thermal internal energy', 1000*(clock()-t0)))
        # entropy
        t0 = clock()
        S_ha = convert_energy(entropy(T, phonons, weights),
                              'kjmol', self._eunit)
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * entropy', 1000*(clock()-t0)))
        # isochoric heat capacity
        t0 = clock()
        Cv_ha = convert_energy(isochoric_heat_capacity(T, phonons, weights),
                               'kjmol', self._eunit)
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * isochoric heat capacity', 1000*(clock()-t0)))
        # vibrational free energy
        t0 = clock()
        Fvib_ha = convert_energy(vibrational_free_energy(T, phonons, weights),
                                 'kjmol', self._eunit)
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * vibrational free energy', 1000*(clock()-t0)))
        # finalization
        t0 = clock()
        U_ha = internal_energy(self.static_energy, Uzp_ha[0], Uth_ha)
        F_ha = free_energy(self.static_energy, Fvib_ha)
        #
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * finalization (total U and F)', clock()-t0))
        self.echo('   HA calculation time {0:15.3f} msec'.format(
            1000*(clock()-t0_ha)))
        self.echo('')
        # If requested, calculate thermodynamics from acoustic frequencies
        # according to Kieffer's approach
        S_ha_optic = S_ha.copy() # to be deleted
        
        if self.acoustic_kieffer.is_on():
            t0 = clock()
            self.echo_override(
                ' - Calculation of acoustic contribution to thermodynamics:')
            with progressbar(self.temperature, width=50) as progress:
                i = 0
                for c in progress:
                    for j in range(self.volume.shape[0]):
                        f = acoustic[:, j]
                        kieffer = Kieffer(f)
                        Cv_aco = convert_energy(kieffer.heat_capacity(T[i]),
                                                'kjmol', self._eunit)
                        S_aco = convert_energy(kieffer.entropy(T[i]),
                                               'kjmol', self._eunit)
                        Cv_ha[i, j] += Cv_aco
                        S_ha[i, j] += S_aco
                        Fvib_ha[i, j] -= T[i] * S_aco * np.power(10., -3.)
                        F_ha[i, j] -= T[i] * S_aco * np.power(10., -3.)
                    i += 1
            self.echo_time(t0, clock())

##        print(S_ha_optic[0, 0], convert_energy(S_ha_optic[0, 0], self._eunit, 'kjmol'))
##        print(S_ha[0, 0], convert_energy(S_ha[0, 0], self._eunit, 'kjmol'))
##        print(100.*(S_ha[0, 0] - S_ha_optic[0, 0])/S_ha[0, 0])
        #
        # Fit harmonic thermodynamics
        #
        t0 = clock()
        t0_fit = t0
        self.echo(' - Polynomial regression of harmonic thermodynamic '
                 'properties [f(V)]:')
        Uzp_pars = polyfit(self.volume, Uzp_ha, self._edeg)
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * zero-point energy', 1000*(clock()-t0)))
        t0 = clock()
        Uth_pars = polyfit(self.volume, Uth_ha, self._edeg)
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * thermal internal energy', 1000*(clock()-t0)))
        t0 = clock()
        S_pars = polyfit(self.volume, S_ha, self._edeg)
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * entropy', 1000*(clock()-t0)))
        t0 = clock()
        Cv_pars = polyfit(self.volume, Cv_ha, self._edeg)
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * isochoric heat capacity', 1000*(clock()-t0)))
        t0 = clock()
        Fvib_pars = polyfit(self.volume, Fvib_ha, self._edeg)
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * vibrational free energy', 1000*(clock()-t0)))
        t0 = clock()
        U_pars = polyfit(self.volume, U_ha, self._edeg)
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * total internal energy', 1000*(clock()-t0)))
        t0 = clock()
        F_pars = polyfit(self.volume, F_ha, self._edeg)
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * total Helmholtz free energy', 1000*(clock()-t0)))
        self.echo_time(t0_fit, clock())
        #
        # Volume minimization
        #
        t0 = clock()
        self.echo_override(
            ' - Volume minimization using {0}:'.format(
                'polynomials' if self.polynomial_minimization.is_on()
                else 'EoS')
            )
        with progressbar(self.temperature, width=50) as progress:
            i = 0
            for c in progress:
                for j in range(self.pressure.shape[0]):
                    VT, KT, Kp = self.minimize_volume(
                    self.volume, F_ha[i], self.pressure[j], F_pars[i])
                    self._results['VT'][i, j] = VT
                    if not self.polynomial_minimization.is_on():
                        self._results['KT'][i, j] = KT
                        self._results['Kp'][i, j] = Kp
                i += 1
        self.echo_time(t0, clock())
        #
        # Thermodynamic interpolation
        #
        t0 = clock()
        t0_interp = t0
        self.echo(' - Thermodynamic interpolation at V(T,P):')
        for i in range(self.temperature.shape[0]):
            self._results['Uzp'][i] = interpolate(
                self._results['VT'][i], Uzp_pars[0])
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * zero-point energy', 1000*(clock()-t0)))
        t0 = clock()
        for i in range(self.temperature.shape[0]):
            self._results['Uth'][i] = interpolate(
                self._results['VT'][i], Uth_pars[i])
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * thermal internal energy', 1000*(clock()-t0)))
        t0 = clock()
        for i in range(self.temperature.shape[0]):
            self._results['S'][i] = interpolate(
                self._results['VT'][i], S_pars[i])
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * entropy', 1000*(clock()-t0)))
        t0 = clock()
        for i in range(self.temperature.shape[0]):
            self._results['Cv'][i] = interpolate(
                self._results['VT'][i], Cv_pars[i])
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * isochoric heat capacity', 1000*(clock()-t0)))
        t0 = clock()
        for i in range(self.temperature.shape[0]):
            self._results['Fvib'][i] = interpolate(
                self._results['VT'][i], Fvib_pars[i])
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * vibrational free energy', 1000*(clock()-t0)))
        t0 = clock()
        for i in range(self.temperature.shape[0]):
            self._results['Utot'][i] = interpolate(
                self._results['VT'][i], U_pars[i])
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * total internal energy', 1000*(clock()-t0)))
        t0 = clock()
        for i in range(self.temperature.shape[0]):
            self._results['F'][i] = interpolate(
                self._results['VT'][i], F_pars[i])
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * total Helmholtz free energy', 1000*(clock()-t0)))
        self.echo_time(t0_interp, clock())
        #
        # Bulk modulus (if minimization is based on polynomials)
        #   TODO: it could be optimized
        #
        if self.polynomial_minimization.is_on():
            t0 = clock()
            self.echo_override(' - Calculation of bulk modulus:')
            with progressbar(self.temperature, width=50) as progress:
                i = 0
                for c in progress:
                    for j in range(self.pressure.shape[0]):
                        V0 = self._results['VT'][i, j]
                        V_grid = V0 + (V0 * grid)
                        F_grid = interpolate(V_grid, F_pars[i])
                        KT, Kp = self.numerical_bulk_modulus(V_grid, F_grid, idx)
                        self._results['KT'][i][j] = KT
                        self._results['Kp'][i][j] = Kp
                    i += 1
            self.echo_time(t0, clock())
        return

    def _finalize_run(self):
        #
        
        t0 = clock()
        t0_f = t0
        self.echo(' - Calculation of enthalpy H(T,P) and Gibbs free '
                 'energy G(T,P):')
        self._results['H'] = self.set_enthalpy(
            self._results['VT'],  self.pressure, self._results['Utot'])
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * enthalpy', 1000*(clock()-t0)))
        #
        t0 = clock()
        self._results['G'] = self.set_gibbs_free_energy(
            self._results['VT'],  self.pressure, self._results['F'])
        self.echo('{0:40} - elapsed time {1:15.3f} msec'.format(
            '    * Gibbs free energy', 1000*(clock()-t0)))
        self.echo_time(t0_f, clock())
        if self.temperature.shape[0] < 50:
            self.echo_warning(' - WARNING! The calculation of:')
            self.echo_warning('   * volumetric thermal expansion coefficient')
            self.echo_warning('   * isobaric heat capacity')
            self.echo_warning('   * adiabatic bulk modulus')
            self.echo_warning('   requires at least 50 points of temperatures.')
            self.echo_warning('')
            self._results['alphaV'] = np.zeros((self.temperature.shape[0],
                                                self.pressure.shape[0]))
            self._results['Cp-Cv'] = np.zeros((self.temperature.shape[0],
                                               self.pressure.shape[0]))
            self._results['Cp'] = np.zeros((self.temperature.shape[0],
                                            self.pressure.shape[0]))
            self._results['gruneisen'] = np.zeros((self.temperature.shape[0],
                                                   self.pressure.shape[0]))
            self._results['KS'] = self._results['KT'].copy()
            
            return

        T = convert_temperature(self.temperature, self._tunit, 'K')
        t0 = clock()
        self.echo(' - Calculation of volumetric thermal expansion coefficient:')
        # Original
##        dVdT = np.gradient(self._results['VT'], axis=0)/np.gradient(T)[:, None]
##        self._results['alphaV'] = dVdT / self._results['VT']
        # Modification
        for i in range(self.pressure.shape[0]):
            V = self._results['VT'][:, i]
            dV = np.zeros(self.temperature.shape[0], dtype=np.float64)
            dV[0:-1] = np.diff(V) / np.diff(T)
            dV[-1] = (V[-1] - V[-2]) / (T[-1] - T[-2])
            self._results['alphaV'][:, i] = dV / V
        self.echo_time(t0, clock())

        t0 = clock()
        self.echo('  - Calculation of isobaric heat capacity')
        Tmat = np.repeat(T, self._results['VT'].shape[1], axis=0)
        Tmat = Tmat.reshape(len(T), self._results['VT'].shape[1])
        Cp_corr = np.power(self._results['alphaV'], 2.)
        Cp_corr *= convert_pressure(self._results['KT'], self._punit, 'Pa')
        Cp_corr *= convert_volume(self._results['VT'], self._vunit, 'm')
        Cp_corr *= N * Tmat
        self._results['Cp-Cv'] = convert_energy(Cp_corr, 'kjmol', self._eunit)
        Cp_corr = None
        self._results['Cp'] = self._results['Cv'] + self._results['Cp-Cv']
        self.echo_time(t0, clock())

        t0 = clock()
        self.echo('  - Calculation of adiabatic bulk modulus K_S(P,T)')
        KS = adiabatic_bulk_modulus(
            T, convert_volume(self._results['VT'], self._vunit, 'm'),
            convert_pressure(self._results['KT'], self._punit, 'Pa'),
            self._results['alphaV'],
            convert_energy(self._results['Cv'], self._eunit, 'kjmol'), 1.)
        self._results['KS'] = convert_pressure(KS, 'Pa', self._punit)
        self.echo_time(t0, clock())

        t0 = clock()
        self.echo('  - Calculation of Gruneisen parameters')
        self._results['gruneisen'] = gruneisen_parameter(
            convert_volume(self._results['VT'], self._vunit, 'm'),
            convert_pressure(self._results['KT'], self._punit, 'Pa'),
            self._results['alphaV'],
            convert_energy(self._results['Cv'], self._eunit, 'kjmol'))
        self.echo_time(t0, clock())
        return

    @staticmethod
    def fine_grid(npoints=5, separation=0.05):
        """ Set an array of factors to generate volumes around a central point.

        Parameters
        ----------
        n: int, optional
            Number of points, including the central point (odd number).
        separation: float, optional
            Separation factor (%) between the points.

        Returns
        -------
        ndarray(dtype=float, ndim=1)
            Array of factors.
        int
            Index of the central point.
        """
        spacing = separation / 100.
        maximum = -((npoints // 2) - npoints)
        minimum = (maximum - npoints)
        expansion = np.arange(minimum, maximum, 1) * spacing
        index = int(len(expansion)//2)  # Central point index
        return expansion, index

    def find_eos_minimum(self, V, F, P0):
        """
        Find the minimum volume from the F function fitted with a
        volume-integrated equation of state formulation.

        Parameters
        ----------

        V: ndarray(dtype=float, ndim=1)
            Array of volume values.
        F: ndarray(dtype=float, ndim=1)
            Array of Helmholtz free energy values at selected temperature.
        P0: float
            Pressure.

        Returns
        -------

        ndarray(dtype=float, ndim=1)
            EOS parameters.

        """
        P0 = self.set_target_pressure(P0)
        F_P = np.array(F) + P0*np.array(V)
        eos = EnergyEOS()
        eos_pars, eos_errs = eos.fit(self.eos_formulation, V, F_P)
        return eos_pars, eos_errs

    def set_target_pressure(self, pressure):
        """
        Set the pressure value from pressure scale to energy/volume scale.

        Parameters
        ----------

        pressure: float
            Pressure value.

        Returns
        -------

        target: float
            Converted pressure value.

        """
        target = pressure_to_energy(pressure, self._eunit, self._vunit,
                                    self._punit)
        return target

    def minimize_volume(self, volume, free_energy, pressure, parameters=None):
        """
        Minimize the :math:`F(V,T)` data using either a polynomial function
        or a volume-integrated equation of state formulation.

        Parameters
        ----------

        volume: ndarray
            Array of unit cell volumes with `float` type.

        free_energy: ndarray
            Array of harmonic Helmholtz free energy with `float` type.

        pressure: float
            Target pressure for the volume minimization.

        parameters: None or ndarray, optional
            If using polynomials for minimization procedure, the parameters
            or the :math:`F(V,T)` curve can be provided, otherwise they are
            calculated from the free_energy data.

        Returns
        -------

        volume_PT: float
            Volume at selected pressure and temperature conditions.

        KT_PT: float
            Bulk modulus at selected pressure and temperature conditions.

        Kp_PT: float
            First-derivative of bulk modulus at selected pressure and
            temperature conditions.

        .. note::

            If using polynomials, KT_PT and Kp_PT are returned as zero values.

        """
        if self.polynomial_minimization.is_on():
            if isinstance(parameters, type(None)):
                parameters = polyfit(volume, free_energy, self._edeg)
            pressure = self.set_target_pressure(pressure)
            volume_PT = find_polynomial_minimum(parameters, pressure)[0]
            return volume_PT, 0., 0.
        elif self.eos_minimization.is_on():
            eos_pars, _ = self.find_eos_minimum(volume, free_energy, pressure)
            eos_pars[1] = energy_to_pressure(eos_pars[1], self._eunit,
                                             self._vunit, self._punit)
            volume_PT = eos_pars[3]
            KT_PT = eos_pars[1]
            Kp_PT = eos_pars[2]
            return volume_PT, KT_PT, Kp_PT

    def numerical_bulk_modulus(self, V, F, idx):
        """
        Calculate the isothermal bulk modulus (KT) by second derivative of
        the Helmholtz free energy of the system with respect to volume:

        .. math::

            K_T = V_0(T) \\bigg( \\frac{\\partial^2 F(V,T)}
            {\\partial V(T)^2} \\bigg)

        The functional form of the :math:`F(V)` is obtained by a polynomial fit
        of selected order.

        Parameters
        ----------

        V0: float
            Volume at which the KT value is calculated.
        V: ndarray(dtype=float, ndim=1)
            Array of volume values.
        F: ndarray(dtype=float, ndim=1)
            Array of Helmholtz free energy values evaluated in V.

        Returns
        -------
        float
            Isothermal bulk modulus (KT) value.

        """
        n = len(V)
        KTs = np.zeros(n, dtype=np.float64)
        pressures = np.zeros(n, dtype=np.float64)

        f = np.polynomial.polynomial.polyfit(V, F, self._edeg)
        df = np.polynomial.polynomial.polyder(f, 1)
        d2f = np.polynomial.polynomial.polyder(f, 2)

        pressures = -energy_to_pressure(np.polynomial.polynomial.polyval(V, df),
                                        self._eunit, self._vunit, self._punit)
        KTs = energy_to_pressure(V * np.polynomial.polynomial.polyval(V, d2f),
                                 self._eunit, self._vunit, self._punit)

        dp = np.polynomial.polynomial.polyder(
            np.polynomial.polynomial.polyfit(pressures, KTs, 2)
            )
        Kp = np.polynomial.polynomial.polyval(pressures[idx], dp)
        return KTs[idx], Kp

    def set_enthalpy(self, volume, pressure, energy):
        """
        Calculate the enthalpy at selected pressure and temperature
        conditions.

        Parameters
        ---------

        volume: ndarray
            Array of :math:`V(P,T)` values with `float` type.

        pressure: ndarray
            Array of pressure values with `float` type.

        energy: ndarray
            Array of :math:`U(P,T)` values with `float` type.

        Returns
        -------

        ndarray
            Enthalpy values in selected units of energy.

        """
        res = enthalpy(
            convert_volume(volume, self._vunit, 'm'),
            convert_pressure(pressure, self._punit, 'Pa'),
            convert_energy(energy, self._eunit, 'kjmol')
            )
        return convert_energy(res, 'kjmol', self._eunit)

    def set_gibbs_free_energy(self, volume, pressure, energy):
        """
        Calculate the Gibbs free energy at selected pressure and temperature
        conditions.

        Parameters
        ---------

        volume: ndarray
            Array of :math:`V(P,T)` values with `float` type.

        pressure: ndarray
            Array of pressure values with `float` type.

        energy: ndarray
            Array of :math:`F(P,T)` values with `float` type.

        Returns
        -------

        ndarray
            Gibbs free energy values in selected units of energy.

        """
        res = gibbs(
            convert_volume(volume, self._vunit, 'm'),
            convert_pressure(pressure, self._punit, 'Pa'),
            convert_energy(energy, self._eunit, 'kjmol')
                )
        return convert_energy(res, 'kjmol', self._eunit)

    def __report_frequency_fit(self, f, R2):
        for i in range(f.shape[0]):
            self.echo(' - Band # {}'.format(i))
            for j in range(f.shape[1]):
                freqN = '   * Frequency {0:5d}: R^2 = {1:.6f}'.format(
                    j+1, R2[i, j])
                if R2[i, j] >= 0.90:
                    msg = freqN + '\tOK'
                elif R2[i, j] < 0.90 and R2[i, j] >= 0.80:
                    msg = freqN + '\tLow'
                else:
                    msg = freqN + '\tBAD!'
                self.echo(msg)
            self.echo('   Averaged R^2: {0:.6f}'.format(np.mean(R2[i])))
            self.echo('')
        return

    def report_results(self):
        """
        Report to the output stream the calculated data at each P-T conditions
        for a rapid check of the results.
        """
        results = self._results
        for i in range(len(self._pressure)):
            self.echo('')
            self.echo('{0}{1:-^78}{0}'.format('#', ''))
            self.echo('{0} Results at {1} {2}'.format('#', self.pressure[i],
                                                 self._punit))
            self.echo('{0}{1:-^78}{0}'.format('#', ''))
            self.echo('')
            self.echo(
                '{0: ^7} {1: ^14} {2: ^16} {3: ^9} {4: ^9} {5: ^16}'.format(
                    'T ({})'.format(self._tunit),
                    'V(T) (A^3)',
                    'a_V (K^-1)',
                    'KT ({})'.format(self._punit),
                    'KS ({})'.format(self._punit),
                    'Cp-Cv (J/mol K)')
                 )
            self.echo(
                '{0:-^7} {0:-^14} {0:-^16} {0:-^9} {0:-^9} {0:-^16}'.format('')
                 )
            for j in range(len(self._temperature)):
                msg = '{0: 7.1f} '.format(self._temperature[j])
                msg += '{0: 14.6f} '.format(results['VT'][j][i])
                msg += '{0: 16.8e} '.format(results['alphaV'][j][i])
                msg += '{0: 9.3f} '.format(results['KT'][j][i])
                msg += '{0: 9.3f} '.format(results['KS'][j][i])
                msg += '{0: 16.6f} '.format(convert_energy(
                    results['Cp-Cv'][j][i], self._eunit, 'kjmol'))
                self.echo(msg)
            self.echo('')
        

    def export_hdf5(self, filename):
        info = """
This file was created with Quantas.

Job: {0}

It contains the results of the quasi-harmonic approximation (QHA)
calculations performed with the following settings:
  - energy scale:      {1}
  - volume scale:      {2}
  - frequency scale:   {3}
  - temperature scale: {4}
  - pressure scale:    {5}
  - calculation of thermodynamics: {6} interpolation
  - volume minimizazion:           {7} fitting

Thermodynamics and mechanical properties were calculated as a function of
pressure and temperature.
""".format(self.input.job if hasattr(self.input, 'job') else 'Unknown',
           self._eunit, self._vunit, self._funit, self._tunit, self._punit,
           'frequency' if self.frequency_interpolation.is_on() else 'thermodynamic',
           'polynomial' if self.polynomial_minimization.is_on() else 'EoS'
           )

        properties = ['Uzp',
                      'Uth',
                      'Utot',
                      'S',
                      'Fvib',
                      'F',
                      'Cv',
                      'VT',
                      'Cp',
                      'Cp-Cv',
                      'KT',
                      'KS',
                      'Kp',
                      'H',
                      'G',
                      'alphaV',
                      'gruneisen'
                      ]

        with h5py.File(filename, 'w') as f:
            f.attrs['info'] = info
            dset = f.create_dataset('T', data=self.temperature)
            dset.attrs['desc'] = 'Temperature'
            dset.attrs['unit'] = self._tunit
            dset = f.create_dataset('P', data=self.pressure)
            dset.attrs['desc'] = 'Pressure'
            dset.attrs['unit'] = self._punit
            for prop in properties:
                dset = f.create_dataset(prop, data=self.results[prop])
                dset.attrs['desc'] = self.descriptions[prop]
                if prop == 'Cv' or prop == 'S':
                    dset.attrs['unit'] = 'm'+self._eunit
                elif prop == 'VT':
                    dset.attrs['unit'] = self._vunit
                elif prop == 'KT' or prop == 'KS':
                    dset.attrs['unit'] = self._punit
                elif prop == 'Kp' or prop == 'gruneisen':
                    dset.attrs['unit'] = ''
                elif prop == 'alphaV':
                    dset.attrs['unit'] = 'K^-1'
                else:
                    dset.attrs['unit'] = self._eunit
        
        self.echo('Calculated data exported to {0}'.format(filename))
        
        return
