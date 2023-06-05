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
This module contains the SeismicCalculator class and some auxiliary functions.
"""

try:
    from time import thread_time as clock
except ImportError:
    # Backward compatibility with Python < 3.7
    from time import clock

import numpy as np
import h5py

# Basic calculator
from quantas.core.calculator import BasicCalculator

# SOEC reader
from quantas.IO.soec_reader import SOECInputFileReader

# Flags
from quantas.utils.flags import Flag

# Christoffel utilities
from .utils.seismic_obj import Seismic
from .utils.seismic_plotter import SeismicPlotter

# Progress bar
from click import progressbar


def plotting_is_available():
    """ Check if plotly and matplotlib are available. """
    try:
        import plotly
        _plotly = True
    except ImportError:
        _plotly = False

    try:
        import matplotlib
        _mpl = True
    except ImportError:
        _mpl = False

    return _mpl, _plotly


class SeismicCalculator(BasicCalculator):
    """ Seismic wave velocity calculator class derived from the
    quantas.core.calculator.BasicCalculator class.

    Parameters
    ----------

    settings: quantas.cui.settings.Settings
        Settings for the current run.

    """

    _plotting = Flag()

    def __init__(self, settings):
        """ Constructor for the Seismic calculator """
        BasicCalculator.__init__(self)
        self.silent = settings['silent']
        self.debug = settings['debug']
        self.log = settings['logfile']
        # angular points
        self.ntheta = settings['ntheta']
        self.nphi = settings['nphi']
        # log settings
        self.echo_debug('Seismic calculator object instantiated')
        # units
        self._punit = settings['pressure_unit']
        # set basic data (input and output)
        self._results = self.__init_data()
        self.echo_debug('Input and output variables set')
        # plotting options
        if settings['plotting']:
            _mpl, _plotly = plotting_is_available()
            if _mpl or _plotly:
                self._plotting.on()
                self.echo_debug('Plotting enabled')
                self._plot_settings = {
                    'mpl': _mpl,
                    'plotly': _plotly,
                    'dpi': settings['dpi']
                    }
            else:
                self.echo_debug('Plotting requested, but packages not present')
                self.echo_debug('Plotting disabled')
        else:
            self.echo_debug('Plotting disabled')
        return

    def __init_data(self):
        """ Create a dictionary for the properties that will be calculated.
        """
        results_dict = {
            'primary': np.zeros(
                (self.ntheta*self.nphi, 14), dtype=float),
            'slow_secondary': np.zeros(
                (self.ntheta*self.nphi, 14), dtype=float),
            'fast_secondary': np.zeros(
                (self.ntheta*self.nphi, 14), dtype=float),
            }
        return results_dict

    @property
    def results(self):
        return self._results

    @property
    def plotting(self):
        return self._plotting.is_on()

    def read_input(self, filename):
        """ This method read an input file for Quantas.

        Parameters
        ----------
        file: str
            Path to the input file.

        Returns
        -------
        error: str or None
            Error in the input file.
        """
        self.echo('Reading input file: {0}'.format(filename))
        self.echo_debug('Instantiate file reader')
        reader = SOECInputFileReader(filename)
        if not reader.completed:
            return reader.error
        else:
            self.echo_debug('Input file correctly read')
            if reader.density == 0.:
                return 'Density not provided, calculation aborted'
            self._job = reader.jobname
            self.echo_debug('Job name set')
            self._seismic = Seismic(reader.stiffness, reader.density)
            self.echo_debug('Seismic object instantiated')
        return

    @property
    def seismic(self):
        return self._seismic

    def report_input_data(self):
        """ Write to the output stream the main information about the input.
        """
        self.echo('Analysis of the sound velocities in ' + self._job)
        self.echo('')
        self.echo('Density: {0} kg m^-3'.format(self.seismic.density))
        self.echo('')
        self.echo('Stiffness matrix (values in GPa)')
        s = ('  ' + 6*'{: 10.4f}    ')
        for i in range(6):
            self.echo(s.format(*tuple(self.seismic.stiffness[i])))
        self.echo('')
        self.echo('Compliance tensor (values in TPa^-1)')
        for i in range(6):
            s = ('  ' + 6*'{:9.6f}    ')
            self.echo(s.format(*tuple(self.seismic.compliance[i]*1000)))
        self.echo('')
        return

    def run(self):
        """ Start the calculation of the wave velocities by solving the
        Christoffel's equation.
        """
        # Start timing
        t0 = clock()
        self.echo_debug('Start Seismic analysis')

        # Get isotropic velocities
        self.echo_debug('Calculate isotropic velocities')
        iso_vel = self.seismic.isotropic_velocity
        txt = "Start calculation of velocities by solving Christoffel's"
        txt += ' equation'
        self.echo(txt)
        #
        storing = [self.results['slow_secondary'],
                   self.results['fast_secondary'],
                   self.results['primary']
                   ]

        # Loop over surface
        c = 0  # Pointer
        self.theta = np.linspace(0., 1., self.ntheta)
        with progressbar(self.theta, width=50) as progress:
            theta = 0
            for prog in progress:
                self.echo_debug(
                    'Calculation of point {0}/{1}:'.format(theta, self.ntheta)
                    )
                for phi in range(self.nphi):
                    self.seismic.spherical_direction(
                        0.5*np.pi*theta/(self.ntheta-1),
                        2.0*np.pi*phi/(self.nphi-1))

                    # Calculate everything and store in variables
                    self.seismic.solve()
                    self.seismic.solve_group()
                    
                    velocity = self.seismic.phase_velocity
                    self.echo_debug('   - phase velocity (done)')
                    polarization = self.seismic.eigenvec
                    self.echo_debug('   - polarization (done)')
                    group_vel = self.seismic.group_velocity
                    self.echo_debug('   - group velocity (done)')
                    group_abs = self.seismic.group_abs
                    self.echo_debug('   - group absolute (done)')
                    powerflow_angle = self.seismic.powerflow_angle
                    self.echo_debug('   - powerflow angle (done)')
                    ehnancement_fac = self.seismic.enhancement
                    self.echo_debug('   - enhancement factor (done)')

                    # Store the results
                    for i in range(3):
                        storing[i][c][0] = self.seismic.theta
                        storing[i][c][1] = self.seismic.phi
                        storing[i][c][2] = velocity[i]
                        storing[i][c][3] = 100*(velocity[i]/iso_vel[i]-1.0)
                        for j in range(3):
                            storing[i][c][4+j] = polarization[i][j]
                        storing[i][c][7] = group_abs[i]
                        storing[i][c][8] = 100*(group_abs[i]/iso_vel[i]-1.0)
                        for j in range(3):
                            storing[i][c][9+j] = group_vel[i][j]
                        storing[i][c][12] = np.degrees(powerflow_angle[i])
                        storing[i][c][13] = ehnancement_fac[i]
                    self.echo_debug('   - storing the results (done)')
                    # Shift pointer
                    c += 1
                theta += 1
                self.echo_debug('   - moving to the next point')
        #
        # Finish
        self.echo('')
        self.echo_time(t0, clock())
        self.completed = True
        self.echo_debug('Run completed. Exiting main routine.')
        return True

    def export_hdf5(self, filename):
        info = """
This file was created with Quantas.

It contains the results of the Christoffel equation appliet to the second-order
elastic constants (SOECs, matrix in Voigt's notation), performed with the
following settings:
  - pressure scale:    {0}
  - density: {1}

This file contains three datasets:
  - slow secondary data (ssecondary)
  - fast secondary data (fsecondary)
  - primary data (primary)

Each dataset contains the following results (measurement unit):
  - Column 1     : Theta angle (rad)
  - Column 2     : Phi angle (rad)
  - Column 3     : Phase velocity (km s^-1)
  - Column 4     : Relative phase velocity (%)
  - Columns 5-7  : Phase polarization (x,y,z)
  - Column 8     : Group absolute velocity (km s^-1)
  - Column 9     : Relative group velocity (%)
  - Columns 10-12: Group velocity (km s^-1)
  - Column 13    : Powerflow angle (degree)
  - Column 14    : Enhancement factor (log_10A)
""".format(self._punit, 'kg m^-3')
        with h5py.File(filename, 'w') as f:
            f.attrs['info'] = info

            dset = f.create_dataset(
                'ssecondary', data=self.results['slow_secondary'])
            dset = f.create_dataset(
                'fsecondary', data=self.results['fast_secondary'])
            dset = f.create_dataset(
                'primary', data=self.results['primary'])
        self.echo('')
        self.echo('Calculated data exported to {0}'.format(filename))
        self.h5File = filename
        return

    def plot(self, filename=None):
        """ Create the 3D and 2D plots of the seismic velocity data.

        Parameters
        ----------

        filename: str, optional
            Path to the hdf5 output file containing the seismic wave results.
            It is optional because, for now, the plots are made from the main
            routine.

        """
        plotter = SeismicPlotter(self._plot_settings)
        if filename is None:
            plotter.read_data(self.h5File)
        else:
            plotter.read_data(filename)
        plotter.run()
        return

