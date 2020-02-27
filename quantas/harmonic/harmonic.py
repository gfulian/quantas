''' Harmonic Approximation (HA) calculator module.

This module performs the HA calculations necessary to obtain harmonic
thermodynamic properties on a single solid phase at selected temperature
conditions.

Required input data
-------------------
volume: ndarray
    Array containing unit cell volumes with float type.
temperature: ndarray
    Array containing temperature values with float type.
static_energy: ndarray
    Array containing static energy (T = 0 K) related to each volume with
    float type.
phonons: ndarray(dtype=float, ndim=2)
    Array of phonon frequencies calculated at each volume with float type.

Output data
-----------
Uzp: ndarray(dtype=float, ndim=2)
    Zero-point energy at selected V-T conditions.
Uth: ndarray(dtype=float, ndim=2)
    Thermal contribution to internal energy at selected V-T conditions.
Utot: ndarray(dtype=float, ndim=2)
    Total energy (U0 + Uzp + Uth) at selected V-T conditions.
S: ndarray(dtype=float, ndim=2)
    Entropy at selected P-T conditions.
F: ndarray(dtype=float, ndim=2)
    Helmholtz free energy at selected V-T conditions.
Fvib: ndarray(dtype=float, ndim=2)
    Vibrational contribution to free energy at selected V-T conditions.
Cv: ndarray(dtype=float, ndim=2)
    Isochoric heat capacity at selected V-T conditions.

'''
try:
    from time import thread_time as clock
except ImportError:
    # Backward compatibility with Python < 3.7
    from time import clock

import numpy as np
import h5py

# Basic claculator
from quantas.core.calculator import BasicCalculator

# Unit conversions
from quantas.utils.physics.units import convert_frequency as cf
from quantas.utils.physics.units import convert_energy as ce
from quantas.utils.physics.units import convert_temperature as ct

# File reader
from quantas.IO.qha_reader import QHAInputFileReader

# Thermodynamic functions
from quantas.utils.physics.statistical_mechanics import zero_point_energy
from quantas.utils.physics.statistical_mechanics import thermal_energy
from quantas.utils.physics.statistical_mechanics import internal_energy
from quantas.utils.physics.statistical_mechanics import entropy
from quantas.utils.physics.statistical_mechanics import vibrational_free_energy
from quantas.utils.physics.statistical_mechanics import isochoric_heat_capacity
from quantas.utils.physics.statistical_mechanics import free_energy

np.warnings.filterwarnings('ignore')


class HACalculator(BasicCalculator):
    '''
    Harmonic calculator class derived from the
    quantas.core.calculator.BasicCalculator class.

    The calculated properties are::

      - zero-point vibrational energy;
      - thermal internal energy;
      - total (vibrational + static) energy;
      - entropy;
      - vibrational Helmholtz free energy;
      - (total) Hemlholtz free energy;
      - constant-volume (isochoric) heat capacity.

    Parameters
    ----------

    settings: dict
        Dictionary containing the settings for the current run.
    '''

    descriptions = {
        'Uzp': 'Zero-point energy',
        'Uth': 'Thermal internal energy',
        'Utot': 'Total (vibrational + static) energy',
        'S': 'Entropy',
        'Fvib': 'Vibrational Helmholtz free energy',
        'F': 'Helmholtz free energy',
        'Cv': 'Isochoric heat capacity'
        }

    def __init__(self, settings):
        ''' Constructor for the Harmonic Approximation calculator '''
        BasicCalculator.__init__(self)
        self.temperature = settings['trange']
        self._eunit = settings['energy_unit']
        self._funit = settings['frequency_unit']
        self._vunit = settings['lenght_unit']
        self._tunit = settings['temperature_unit']
        self.silent = settings['silent']
        self.log = settings['logfile']
        return

    @property
    def temperature(self):
        ''' Get the temperature range for the current calculation.

        Returns
        -------

        ndarray(dtype=float, ndim=1)
            Array containing the target temperature values.

        '''
        return self._temperature

    @temperature.setter
    def temperature(self, temp):
        ''' Set the temperature range for the current calculation.

        Parameters
        ----------

        temp: tuple
            Tuple containing the minimum, maximum and increment of temperature
            values with float type.

        '''
        self._temperature = np.arange(temp[0], temp[1]+temp[2], temp[2])
        return

    def read_input(self, filename):
        ''' This method read the provided input file and return a dictionary
        containing the necessary data.

        Parameters
        ----------
        filename: str
            Complete path of the inputfile.

        '''
        self.echo('Reading input file: {0}\n'.format(filename))
        data = QHAInputFileReader(filename)

        if not data.completed:
            return data.error

        self.input = data
        return None

    def report_input_data(self):
        ''' This method write to the output stream the main information about
        the input.
        '''
        idata = self.input
        elabel = 'Energy (' + self._eunit + ')'
        vlabel = 'Volume (' + self._vunit + '^3)'
        self.echo('')
        self.echo('Job: {0}'.format(idata.jobname))
        self.echo('')
        self.echo('System:')
        self.echo('- {0:40} {1}'.format('Number of volumes', idata.nvol))
        self.echo('- {0:40} {1}'.format('Number of atoms', idata.natoms))
        self.echo('- {0:40} {1}'.format(
            'Number of sampled k-points', idata.kpoints))
        self.echo('- {0:40} {1}'.format(
            'Number of frequencies', int(
                idata.natoms * 3 * idata.total_q_points)))
        self.echo('')
        self.echo('Volume and energy values:\n')
        self.echo(' {0: ^15} {1: ^30}'.format(vlabel, elabel))
        self.echo(' {0:-^15} {1:-^30}'.format('', ''))
        for i in range(idata.nvol):
            self.echo(' {0: ^15.6f} {1: ^30.12e}'.format(
                idata.volume[i], idata.energy[i]))
        self.echo('')
        return

    def init_results(self, nv, nt):
        '''
        This method initiates the result arrays with zeros.

        Parameters
        ----------

        nt: int
            Number of temperature values.

        nv: int
            Number of unit cell volumes.

        '''
        self._results = {
            'Uzp': np.zeros((1, nv), dtype=np.float64),
            'Uth': np.zeros((nt, nv), dtype=np.float64),
            'Utot': np.zeros((nt, nv), dtype=np.float64),
            'S': np.zeros((nt, nv), dtype=np.float64),
            'F': np.zeros((nt, nv), dtype=np.float64),
            'Fvib': np.zeros((nt, nv), dtype=np.float64),
            'Cv': np.zeros((nt, nv), dtype=np.float64),
            }
        return

    def run(self):
        ''' Start the harmonic approximation calculations. '''
        self.echo('{0}{1:-^78}{0}'.format(
            '#', 'Harmonic Approximation calculation started'))
        self.echo('')
        #
        # Start timer
        tstart = clock()
        np.seterr(over='ignore')
        #
        # Calculate weights for phonon bands
        weights = self.input.weights/self.input.total_q_points
        #
        # Convert frequency and temperature
        temperature = ct(self.temperature, self._tunit, 'K')
        freq = cf(self.input.frequencies, self._funit, 'Hz')
        #
        # Set the arrays for calculation and data storing
        self.init_results(self.input.nvol, self._temperature.shape[0])
        #
        # Calculate zero-point energy
        t0 = clock()
        self.echo(' - Start calculation of zero-point energy')
        self._results['Uzp'] = ce(
            zero_point_energy(np.zeros(1, dtype=float), freq, weights),
            'kjmol', self._eunit)
        self.echo('   Finished, elapsed time {0:8.2f} sec'.format(clock()-t0))
        self.echo('')
        #
        t0 = clock()
        self.echo(' - Start calculation of thermal internal energy')
        self._results['Uth'] = ce(
            thermal_energy(temperature, freq, weights),
            'kjmol', self._eunit)
        self.echo('   Finished, elapsed time {0:8.2f} sec'.format(clock()-t0))
        self.echo('')
        #
        t0 = clock()
        self.echo(' - Start calculation of entropy')
        self._results['S'] = ce(
            entropy(temperature, freq, weights), 'kjmol', self._eunit)
        self.echo('   Finished, elapsed time {0:8.2f} sec'.format(clock()-t0))
        self.echo('')
        #
        t0 = clock()
        self.echo(' - Start calculation of isochoric heat capacity')
        self._results['Cv'] = ce(
            isochoric_heat_capacity(temperature, freq, weights),
            'kjmol', self._eunit)
        self.echo('   Finished, elapsed time {0:8.2f} sec'.format(clock()-t0))
        self.echo('')
        #
        t0 = clock()
        self.echo(' - Start calculation of vibrational Helmholtz free energy')
        self._results['Fvib'] = ce(
            vibrational_free_energy(temperature, freq, weights),
            'kjmol', self._eunit)
        self.echo('   Finished, elapsed time {0:8.2f} sec'.format(clock()-t0))
        self.echo('')
        #
        t0 = clock()
        self.echo(' - Calculate total internal internal energy and '
                 'Helmholtz free energy')
        self._results['Utot'] = internal_energy(
            self.input.energy, self._results['Uzp'][0], self._results['Uth'])
        self._results['F'] = free_energy(
            self.input.energy, self._results['Fvib'])
        #
        self.echo('   Finished, elapsed time {0:8.2f} sec'.format(clock()-t0))
        self.echo('')
        #
        self.echo('All done!')
        self.echo('')
        #
        # Stop timer
        tfinish = clock()
        elapse = tfinish-tstart
        self.echo('Total calculation time: {0:6.2f} sec'.format(elapse))
        self.echo('{0}{1:-^78}{0}'.format('#', 'HA calculations finished'))
        self.echo('')
        #
        # Calculation finished
        self.completed = True
        return

    def export_hdf5(self, filename):
        info = '''
This file was created with Quantas.

Job: {0}

It contains the results of the harmonic approximation (HA) calculations
performed with the following settings:
  - energy scale: {1}
  - volume scale: {2}
  - frequency scale: {3}
  - temperature scale: {4}

Thermodynamics properties were calculated as a function of volume and
temperature.
'''.format(self.input.jobname, self._eunit, self._vunit, self._funit,
           self._tunit)
        properties = ['Uzp', 'Uth', 'Utot', 'S', 'Fvib', 'F', 'Cv']

        with h5py.File(filename, 'w') as f:
            f.attrs['info'] = info
            dset = f.create_dataset('V', data=self.input.volume)
            dset.attrs['desc'] = 'Unit cell volume'
            dset.attrs['unit'] = self._vunit

            dset = f.create_dataset('T', data=self.temperature)
            dset.attrs['desc'] = 'Temperature'
            dset.attrs['unit'] = self._tunit

            dset = f.create_dataset('U0', data=self.input.energy)
            dset.attrs['desc'] = 'Static (electronic) energies'
            dset.attrs['unit'] = self._eunit
            for prop in properties:
                dset = f.create_dataset(prop, data=self.results[prop])
                dset.attrs['desc'] = self.descriptions[prop]
                if prop == 'Cv' or prop == 'S':
                    dset.attrs['unit'] = 'm'+self._eunit
                else:
                    dset.attrs['unit'] = self._eunit
        self.echo('Calculated data exported to {0}'.format(filename))
        return
