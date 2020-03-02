# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

try:
    from time import thread_time as clock
except ImportError:
    # Backward compatibility with Python < 3.7
    from time import clock

import numpy as np
from scipy import optimize
import h5py
try:
    import matplotlib.pyplot as plt
    mpl = True
except ImportError:
    mpl = False

# Basic claculator
from quantas.core.calculator import BasicCalculator

# SOEC
from .utils.soec_obj import SOEC, SOECOrtho, vector
from .utils.symmetry import check_symmetry
from .utils.plotting import young_polar_plot
from .utils.plotting import compressibility_polar_plot
from .utils.plotting import shear_polar_plot
from .utils.plotting import poisson_polar_plot
from .utils.plotting import seismic_waves_polar_plot

# File reader
from quantas.IO.soec_reader import SOECInputFileReader

# Output stream formatting constants
tl = '    {0:-^10}---{0:-^32}---{0:-^32}'           # Top/bottom line
ml = '    {0:-^10}-|-{0:-^32}-|-{0:-^32}'           # Middle line
dh = '    {:10}' + 2*' | {: ^32}'                   # Double header
dv = '    {:10}' + 2*' | {: ^32.4f}'                # Double value
qh = '    {:10}' + 2*' | {: ^15}  {: ^15}'          # Quad header
qv = '    {:10}' + 2*' | {: ^15.4f}  {: ^15.4f}'    # Quad value
ttl = '    {0:-^10}'+ 3*'---{0:-^22}'               # Top/bottom line
tml = '    {0:-^10}'+ 3*'-|-{0:-^22}'               # Middle line
th = '    {:10}' + 3*' | {: ^22}'                   # Triple header
tv = '    {:10}' + 3*' | {: ^22.4f}'                # Triple value
hh = '    {:10}' + 3*' | {: ^10}  {: ^10}'          # Hexa header
hv = '    {:10}' + 3*' | {: ^10.4f}  {: ^10.4f}'    # Hexa values

def minimize(func, dim, method='fmin'):
    if dim == 2:
        r = ((0, np.pi), (0, np.pi))
        n = 25
        x0 = [1., 1.]
    elif dim == 3:
        r = ((0, np.pi), (0, np.pi), (0, np.pi))
        n = 10
        x0 = [1., 1., 1.]
    if method == 'fmin':
        return optimize.brute(func, r, Ns = n, full_output = True,
                              finish = optimize.fmin)[0:2]
    elif method == 'bh':
        res = optimize.basinhopping(func, x0, niter=200)
        return [ res.x, res.fun ]


def maximize(func, dim):
    res = minimize(lambda x: -func(x), dim)
    return (res[0], -res[1])


class SOECCalculator(BasicCalculator):
    """

    """

    _plotting = False
    _plot_basename = ''
    _plot_dpi = 80

    _plot_functions = {
        'E': young_polar_plot,
        'LC': compressibility_polar_plot,
        'G': shear_polar_plot,
        'Nu': poisson_polar_plot,
        'waves': seismic_waves_polar_plot,
        }

    def __init__(self, settings):
        """ Constructor for the second-order elastic constants calculator """
        BasicCalculator.__init__(self)
        self.polar = settings['polar']
        self.silent = settings['silent']
        self.debug = settings['debug']
        self.log = settings['logfile']
        # log settings
        self.echo_debug('SOEC calculator object instantiated')
        # units
        self._punit = settings['pressure_unit']
        # set basic data (input and output)
        self._results = self.__init_data()
        self.echo_debug('Input and output variables set')
        # plotting options
        if settings['plotting']:
            self._plotting = True
            self._plot_basename = settings['outfig']
            self._plot_dpi = settings['dpi']
            self.echo_debug('Plotting enabled')
        else:
            self.echo_debug('Plotting disabled')
        return

    def __init_data(self):
        """
        Create a dictionary for the elastic properties that will be
        calculated.
        """
        self._soec = None
        self._job = None
        self._system = None
        results_dict = {
            'avgs': None,
            'eigv': None,
            'polar': {
                'xy': {
                    'E': np.zeros((1, 2)),
                    'LC': np.zeros((1, 3)),
                    'G': np.zeros((1, 3)),
                    'Nu': np.zeros((1, 4)),
                    'waves': np.zeros((1, 4)),
                    },
                'xz': {
                    'E': np.zeros((1, 2)),
                    'LC': np.zeros((1, 3)),
                    'G': np.zeros((1, 3)),
                    'Nu': np.zeros((1, 4)),
                    'waves': np.zeros((1, 4)),
                    },
                'yz': {
                    'E': np.zeros((1, 2)),
                    'LC': np.zeros((1, 3)),
                    'G': np.zeros((1, 3)),
                    'Nu': np.zeros((1, 4)),
                    'waves': np.zeros((1, 4)),
                    }
                },
            }
        return results_dict

    @property
    def results(self):
        return self._results

    @property
    def soec(self):
        return self._soec

    def read_input(self, filename):
        """ Read an input file for Quantas and store the elastic data.

        Parameters
        ----------
        file: str
            Path to the input file.
        """
        
        self.echo('Reading input file: {0}'.format(filename))
        self.echo_debug('Instantiate file reader')
        reader = SOECInputFileReader(filename)
        if not reader.completed:
            return reader.error
        else:
            
            self.echo_debug('Input file correctly read')
            self._job = reader.jobname
            self.echo_debug('Job name set')
            self._soec = SOEC(reader.stiffness, reader.density)
            self.echo_debug('SOEC object instantiated')
            self._system = check_symmetry(self._soec.stiffness)
            self.echo_debug('Crystal system: {}'.format(self._system))
            if self._system == 'orthorhombic':
                self.echo_debug('Switching to orthorombic SOEC object')
                self._soec = SOECOrtho(self._soec)
        return

    def report_input_data(self):
        """ This method write to the output stream the main information about
        the input.
        """
        self.echo('\nElastic analysis of ' + self._job)
        self.echo('')
        self.echo('System is ' + self._system)
        if self.soec.density == 0.:
            self.echo('Density: not provided')
        else:
            self.echo('Density: {0} kg m^-3'.format(self.soec.density))
        self.echo('')
        self.echo('Stiffness matrix (values in GPa)')
        s = ('  ' + 6*'{: 10.4f}    ')
        for i in range(6):
            self.echo(s.format(*tuple(self.soec.stiffness[i])))
        self.echo('')
        self.echo('Compliance tensor (values in TPa^-1)')
        for i in range(6):
            s = ('  ' + 6*'{:9.6f}    ')
            self.echo(s.format(*tuple(self.soec.compliance[i]*1000)))
        self.echo('')
        
        return

    def run(self):
        """
        Start SOEC analysis.
        """
        #
        # Start timing
        tstart = clock()
        self.echo_debug('Start SOEC analysis')
        #
        # Calculate averages
        self.echo_debug(' - Calculation of VRH average values')
        avg = self.soec.averages()
        #
        # Calculate stiffness eigenvalues
        self.echo_debug(' - Calculation stiffness eigenvalues')
        eigv = sorted(np.linalg.eig(self.soec.stiffness)[0])
        #
        # Report current data
        
        self.report_initial_results(avg, eigv)
        #
        # Store results
        
        self._results['avgs'] = avg
        self._results['eigv'] = eigv
        #
        # Check SOEC eigenvalues consistency
        if any(eigv) <= 0:
            msg = 'SOEC matrix is not definite positive, '+\
            'crystal is mechanically unstable'
            self.echo(msg)
            self.echo('No further analysis will be performed.')
            return
        #
        # Calculate variations of the mechanical properties
        self.echo_debug(
            ' - Calculation Young modulus and Linear Compressibility variation')
        E, LC = self.single_vector_variation()
        self.report_variation_single(E, LC)
        self.echo_debug(
            ' - Calculation shear modulus and Poisson ratio variation')
        G, Nu = self.dual_vector_variation()
        self.report_variation_dual(G, Nu)
        if self.soec.density > 0.:
            self.echo_debug(' - Calculation wave velocities variation')
            Vp, Vs1, Vs2 = self.phase_velocity_variation()
            self.report_variation_seismic(Vp, Vs1, Vs2)
        #
        # Polar properties
        if self.polar:
            self.calculate_polar_properties()
        #
        # Finish
        self.echo('')
        msg = 'Calculation time: {0:8.1f} sec'.format(clock()-tstart)
        self.echo(msg)
        #
        # Plotting
        if self.polar and self._plotting and mpl:
            self.echo('')
            self.echo('Plotting results as requested:')
            if self.soec.density > 0.:
                prop_list = ['E', 'LC', 'G', 'Nu', 'waves']
            else:
                prop_list = ['E', 'LC', 'G', 'Nu']

            for prop in prop_list:
                figname = self._plot_basename + '_' + prop + '.png'
                figure = self._plot_functions[prop](
                    polar['xy'][prop], polar['xz'][prop], polar['yz'][prop])
                figure.savefig(figname, dpi=self._plot_dpi)
                self.echo(' - figure {} generated'.format(figname))

        elif self._plotting and not mpl:
            self.echo('')
            self.echo('Plotting requested, but Matplotlib not found')
        self.completed = True
        return

    def calculate_polar_properties(self, nrad=360):
        """
        Calculate elastic properties on the (xy), (xz) and (yz) planes.

        Parameters
        ----------
        
        nrad: int, optional
            Number of radians to be considered
        """
        polar = {}
        self.echo(' - Calculation of polar (2D) properties:')
        phi = np.linspace(0, 2*np.pi, nrad, dtype=np.float64)
        self.echo('     * along (xy)')
        polar['xy'] = {}
        theta = (np.pi/2)*np.ones(nrad, dtype=np.float64) # set XY
        self.echo("         a. Young's modulus")
        young_xy = self.soec.polar_young(theta, phi)
        self.echo('         b. Linear compressibility')
        lc_xy = self.soec.polar_compressibility(theta, phi)
        self.echo('         c. Shear modulus')
        shear_xy = self.soec.polar_shear(theta, phi)
        self.echo("         d. Poisson's ratio")
        poisson_xy = self.soec.polar_poisson(theta, phi)
        if self.soec.density > 0.:
            self.echo('         e. Wave velocities')
            waves_xy = self.soec.polar_waves(theta, phi)

        self.echo('     * along (xz)')
        polar['xz'] = {}
        theta = np.zeros(nrad, dtype=np.float64) # set XZ
        self.echo("         a. Young's modulus")
        young_xz = self.soec.polar_young(phi, theta)
        self.echo('         b. Linear compressibility')
        lc_xz = self.soec.polar_compressibility(phi, theta)
        self.echo('         c. Shear modulus')
        shear_xz = self.soec.polar_shear(phi, theta)
        self.echo("         d. Poisson's ratio")
        poisson_xz = self.soec.polar_poisson(phi, theta)
        if self.soec.density > 0.:
            self.echo('         e. Wave velocities')
            waves_xz = self.soec.polar_waves(phi, theta)

        self.echo('     * along (yz)')
        polar['yz'] = {}
        theta = (np.pi/2)*np.ones(nrad, dtype=np.float64) # set XY
        self.echo("         a. Young's modulus")
        young_yz = self.soec.polar_young(phi, theta)
        self.echo('         b. Linear compressibility')
        lc_yz = self.soec.polar_compressibility(phi, theta)
        self.echo('         c. Shear modulus')
        shear_yz = self.soec.polar_shear(phi, theta)
        self.echo("         d. Poisson's ratio")
        poisson_yz = self.soec.polar_poisson(phi, theta)
        if self.soec.density > 0.:
            self.echo('         e. Wave velocities')
            waves_yz = self.soec.polar_waves(phi, theta)
        #
        # Store polar properties
        polar['xy']['E'] = np.column_stack((np.degrees(phi), young_xy))
        polar['xy']['LC'] = np.column_stack((np.degrees(phi).T, lc_xy))
        polar['xy']['G'] = np.column_stack((np.degrees(phi).T, shear_xy))
        polar['xy']['Nu'] = np.column_stack((np.degrees(phi).T, poisson_xy))
        if self.soec.density > 0.:
            polar['xy']['waves'] = np.column_stack((np.degrees(phi).T, waves_xy))

        polar['xz']['E'] = np.column_stack((np.degrees(phi), young_xz))
        polar['xz']['LC'] = np.column_stack((np.degrees(phi).T, lc_xz))
        polar['xz']['G'] = np.column_stack((np.degrees(phi).T, shear_xz))
        polar['xz']['Nu'] = np.column_stack((np.degrees(phi).T, poisson_xz))
        if self.soec.density > 0.:
            polar['xz']['waves'] = np.column_stack((np.degrees(phi).T, waves_xz))

        polar['yz']['E'] = np.column_stack((np.degrees(phi), young_yz))
        polar['yz']['LC'] = np.column_stack((np.degrees(phi).T, lc_yz))
        polar['yz']['G'] = np.column_stack((np.degrees(phi).T, shear_yz))
        polar['yz']['Nu'] = np.column_stack((np.degrees(phi).T, poisson_yz))
        if self.soec.density > 0.:
            polar['yz']['waves'] = np.column_stack((np.degrees(phi).T, waves_yz))
        self._results['polar'] = polar
        return

    def report_initial_results(self, avg, eigenval):
        """
        """
        self.echo('Average properties')
        avgnames = ['Voigt', 'Reuss', 'Hill']
        head1 = ('', 'Bulk', "Young's", 'Shear', "Poisson's")
        head2 = ('', 'modulus', 'modulus', 'modulus', 'ratio')
        head3 = ('', '(GPa)', '(GPa)', '(GPa)', '')
        sh = '{:5}' + 4*'  {: ^10}'
        self.echo(sh.format(*head1)+'\n'+sh.format(*head2)+'\n'+sh.format(*head3))
        sv = '{:5}' + 4*'  {: ^10.5f}'
        for i in range(3):
            self.echo(sv.format(avgnames[i], *tuple(avg[i])))
        self.echo('')
        self.echo('Eigenvalues of the stiffness matrix:')
        for i in range(6):
            self.echo('    lambda_{0}: {1: ^7.5f}'.format(i+1, eigenval[i]))
        self.echo('')
        return

    def directional_variation(self, func, dim):
        """
        """
        #
        # Minimum and maximum value
        minval = minimize(func, dim)
        maxval = maximize(func, dim)
        #
        # Anisotropy
        if minval[1] > 0.:
            anisotropy = maxval[1]/minval[1]
        else:
            anisotropy = np.inf
        #
        # First (or Only) axis
        minax1 = vector(*minval[0])
        maxax1 = vector(*maxval[0])
        if dim == 2:
            minax2 = None
            maxax2 = None
        else:
            #
            # Second axis
            minax2 = vector(*minval[0])
            maxax2 = vector(*maxval[0])
        return [ (minval[1], maxval[1]), anisotropy, minax1, maxax1, minax2, maxax2 ]

    def single_vector_variation(self):
        """
        """
        E_var = self.directional_variation(self.soec.young_modulus, 2)
        LC_var = self.directional_variation(self.soec.linear_compressibility, 2)
        return E_var, LC_var

    def dual_vector_variation(self):
        """
        """
        G_var = self.directional_variation(self.soec.shear_modulus, 3)
        Nu_var = self.directional_variation(self.soec.poisson_ratio, 3)
        return G_var, Nu_var

    def phase_velocity_variation(self):
        """
        """
        Vp_var = self.directional_variation(self.soec.longitudinal_wave, 2)
        Vs1_var = self.directional_variation(self.soec.shear_wave_1, 2)
        Vs2_var = self.directional_variation(self.soec.shear_wave_2, 2)
        return Vp_var, Vs1_var, Vs2_var

    def report_variation_single(self, E_var, LC_var):
        """
        """
        self.echo('Variations of the elastic moduli:')
        self.echo('')
        self.echo(tl.format(''))
        self.echo(dh.format('', *("Young's modulus", 'Linear compressibility')))
        self.echo(ml.format(''))
        self.echo(qh.format('',*('E_min', 'E_max', 'beta_min', 'beta_max')))
        self.echo(qv.format('Values', *tuple(E_var[0]), *tuple(LC_var[0])))
        self.echo(ml.format(''))
        self.echo(dv.format('Anisotropy', E_var[1], LC_var[1]))
        self.echo(ml.format(''))
        label = ['', 'Axis', '']
        for i in range(3):
            self.echo(qv.format(label[i], E_var[2][i], E_var[3][i],
                           LC_var[2][i], LC_var[3][i]))
        self.echo(tl.format(''))
        self.echo('Notes: E min/max values in GPa, beta min/max values in TPa^-1')
        self.echo('')
        return

    def report_variation_dual(self, G_var, Nu_var):
        """
        """
        self.echo(tl.format(''))
        self.echo(dh.format('', *('Shear modulus', "Poisson's ratio")))
        self.echo(ml.format(''))
        self.echo(qh.format('', *('G_min', 'G_max', 'nu_min', 'nu_max')))
        self.echo(qv.format('Values', *tuple(G_var[0]), *tuple(Nu_var[0])))
        self.echo(ml.format(''))
        self.echo(dv.format('Anisotropy', G_var[1], Nu_var[1]))
        self.echo(ml.format(''))
        label = ['', '1st Axis', '']
        for i in range(3):
            self.echo(qv.format(label[i], G_var[2][i], G_var[3][i],
                            Nu_var[2][i], Nu_var[3][i]))
        self.echo(ml.format(''))
        label = ['', '2nd Axis', '']
        for i in range(3):
            self.echo(qv.format(label[i], G_var[4][i], G_var[5][i],
                            Nu_var[4][i], Nu_var[5][i]))
        self.echo(tl.format(''))
        self.echo('Notes: G min/max values in GPa')
        self.echo('')
        return

    def report_variation_seismic(self, Vp_var, Vs1_var, Vs2_var):
        self.echo('')
        self.echo('Variations of the seismic velocities:')
        self.echo('')
        self.echo(ttl.format(''))
        self.echo(th.format('', *('V_s1', 'V_s2', 'V_p')))
        self.echo(tml.format(''))
        self.echo(hh.format('', *('min', 'max', 'min', 'max', 'min', 'max')))
        self.echo(hv.format('Values', *tuple(Vs1_var[0]), *tuple(Vs2_var[0]),
                        *tuple(Vp_var[0])))
        self.echo(tml.format(''))
        self.echo(tv.format('Anisotropy', Vs1_var[1], Vs2_var[1], Vp_var[1]))
        self.echo(tml.format(''))
        label = ['', 'Axis', '']
        for i in range(3):
            self.echo(hv.format(label[i],
                           Vs1_var[2][i], Vs1_var[3][i],
                           Vs2_var[2][i], Vs2_var[3][i],
                           Vp_var[2][i], Vp_var[3][i],
                           )
                 )
        self.echo(ttl.format(''))
        self.echo('Notes: min/max values in km s^-1')
        self.echo('')
        return

    def export_hdf5(self, filename):
        info = '''
This file was created with Quantas.

It contains the results of the analysis of the second-order elastic constants
(SOECs) matrix in Voigt's notation, performed with the following settings:
  - pressure scale:    {0}
  - density: {1}

'''.format(self._punit, 'kg m^-3')
        with h5py.File(filename, 'w') as f:
            f.attrs['info'] = info
            dset = f.create_dataset('stiffness', data=self.soec.stiffness)
            dset.attrs['unit'] = self._punit
            dset = f.create_dataset('compliance', data=self.soec.compliance)
            dset.attrs['unit'] = 'TPa^-1'
            dset = f.create_dataset('averages', data=self.results['avgs'])
            dset.attrs['unit'] = self._punit
            dset = f.create_dataset('eigenvals', data=self.results['eigv'])
            dset.attrs['unit'] = self._punit
            for key in self.results['polar']:
                for var in self.results['polar'][key]:
                    dsetname = 'polar_'+var+'_'+key
                    dset = f.create_dataset(
                        dsetname, data=self.results['polar'][key][var])
                    if var in ['E', 'G']:
                        dset.attrs['unit'] = self._punit
                    # TO-DO: add support for more units
                    elif var == 'LC':
                        dset.attrs['unit'] = 'TPa^-1'
                    elif var == 'waves':
                        dset.attrs['unit'] = 'km s^-1'
                    else:
                        dset.attrs['unit'] = ''
        
        self.echo('')
        self.echo('Calculated data exported to {0}'.format(filename))
        return
