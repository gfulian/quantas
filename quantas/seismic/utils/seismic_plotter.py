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
Module containing the class used to make 2D and 3D plots of the seismic
wave velocities.
"""

import h5py
import numpy as np

# Basic plotter
from quantas.core.plotter import BasicPlotter

# Other Quantas utilities
from quantas.IO.hdf5_reader import QuantasHDF5Reader
from quantas.utils.flags import Flag

# Check if plotly is available
try:
    import plotly
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    _plotly = Flag(activated=True)
except ImportError:
    _plotly = Flag(activated=False)

# Check if matplotlib is available
try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as tick
    from matplotlib.colors import LinearSegmentedColormap
    _mpl = Flag(activated=True)
except ImportError:
    _mpl = Flag(activated=False)

# --------------------
# Colorscales
# --------------------
default_colorscale = ['blue', 'white', 'red']
powerflow_colorscale = ['black', 'cyan']
enhancement_colorscale = ['blue', 'white', 'red', 'black']

class SeismicPlotter(BasicPlotter):
    """
    The SeismicPlotter class that is used to make 2D and 3D plots of the
    seismic wave velocities.

    Parameters
    ----------
    
    settings: dict
        Dictionary containing the general settings used to plot the
        results.

    Attributes
    ----------

    titles: dict
        Dictionary containing the names of ...

    datasets: dict
        Dictionary containing the information on each datasets calculated
        with the SeismicCalculator.

    layout_3D: dict
        Dictionary with the layout settings for a 3D plot made with plotly.

    axis_x_3D: dict
        Dictionary with the necessary information to draw the x axis on a 3D
        plot with plotly.

    axis_y_3D: dict
        Dictionary with the necessary information to draw the y axis on a 3D
        plot with plotly.

    axis_y_3D: dict
        Dictionary with the necessary information to draw the z axis on a 3D
        plot with plotly.

    axis_x_3D: dict
        Dictionary with the settings used to draw a colorbar on the 3D plot
        made with plotly.

    """
    

    titles = {
        'ssecondary': 'Slow Secondary',
        'fsecondary': 'Fast Secondary',
        'primary': 'Primary'
        }

    datasets = {
        2: {
            'name': 'phase velocity',
            'title': r'$\huge{v_p \: (km \: s^{-1})}$',
            'title2D': r'$v_p \: (km \: s^{-1})$',
            'suffix': 'Vp',
            'colorscale': default_colorscale,
            'cmin': {
                'ssecondary': None,
                'fsecondary': None,
                'primary': None
                },
            'cmax': {
                'ssecondary': None,
                'fsecondary': None,
                'primary': None
                },
            'levels': [10, 10, 10],
            'anisotropy': True
            },
        3: {
            'name': 'relative phase velocity',
            'title': r'$\huge{v_p - v_{iso} \: (\%)}$',
            'title2D': r'$v_p - v_{iso} \: (\%)$',
            'suffix': 'Vp_rel',
            'colorscale': default_colorscale,
            'cmin': {
                'ssecondary': None,
                'fsecondary': None,
                'primary': None
                },
            'cmax': {
                'ssecondary': None,
                'fsecondary': None,
                'primary': None
                },
            'levels': [10, 10, 10],
            'anisotropy': False
            },
        7: {
            'name': 'group velocity',
            'title': r'$\huge{v_g \: (km \: s^{-1})}$',
            'title2D': r'$v_g \: (km \: s^{-1})$',
            'suffix': 'Vg',
            'colorscale': default_colorscale,
            'cmin': {
                'ssecondary': None,
                'fsecondary': None,
                'primary': None
                },
            'cmax': {
                'ssecondary': None,
                'fsecondary': None,
                'primary': None
                },
            'levels': [10, 10, 10],
            'anisotropy': True
            },
        8: {
            'name': 'relative group velocity',
            'title': r'$\huge{v_g - v_{iso} \: (\%)}$',
            'title2D': r'$v_g - v_{iso} \: (\%)$',
            'suffix': 'Vg_rel',
            'colorscale': default_colorscale,
            'cmin': {
                'ssecondary': None,
                'fsecondary': None,
                'primary': None
                },
            'cmax': {
                'ssecondary': None,
                'fsecondary': None,
                'primary': None
                },
            'levels': [10, 10, 10],
            'anisotropy': False
            },
        12: {
            'name': 'powerflow angle',
            'title': r'$\huge{\psi \: (°)}$',
            'title2D': r'$\psi \: (°)}$',
            'suffix': 'pf_angle',
            'colorscale': powerflow_colorscale,
            'cmin': {
                'ssecondary': None,
                'fsecondary': None,
                'primary': None
                },
            'cmax': {
                'ssecondary': None,
                'fsecondary': None,
                'primary': None
                },
            'levels': [5, 5, 5],
            'anisotropy': False
            },
        13: {
            'name': 'enhancement factor',
            'title': r'$\huge{\log_{10}(A)}$',
            'title2D': r'$\log_{10}(A)$',
            'suffix': 'enhancement_factor',
            'colorscale': enhancement_colorscale,
            'cmin': {
                'ssecondary': -2.5,
                'fsecondary': -2.5,
                'primary': -0.6
                },
            'cmax': {
                'ssecondary': 5,
                'fsecondary': 5,
                'primary': 1.2
                }
            ,
            'levels': [0, 0, 0],
            'anisotropy': False
            },
        
        }

    layout_3D = {
        'scene': {
            'xaxis': {
                'visible': False,
                'range': [-1, 1.5]
                },
            'yaxis': {
                'visible': False,
                'range': [-1, 1.5]
                },
            'zaxis': {
                'visible': False,
                'range': [-1, 1.5]
                },
            'camera': {
                'center': {
                    'x': 0,
                    'y': 0,
                    'z': 0
                    },
                'eye': {
                    'x': 1.1,
                    'y': 0.7,
                    'z': 0.7
                    },
                'projection': {
                    'type': 'perspective'
                    },
                'up': {
                    'x': 0,
                    'y': 0,
                    'z': 1
                    }
                },
        'aspectratio': {
            'x': 1,
            'y': 1,
            'z': 1
            },
        'annotations': [
                {
                    'showarrow': False,
                    'x': 1.5,
                    'y': 0,
                    'z': -0.04,
                    'text': 'x',
                    'font': {
                        'color': 'black',
                        'size': 40
                        }
                    },
                {
                    'showarrow': False,
                    'x': 0,
                    'y': 1.5,
                    'z': -0.04,
                    'text': 'y',
                    'font': {
                        'color': 'black',
                        'size': 40
                        }
                    },
                {
                    'showarrow': False,
                    'x': 0,
                    'y': 0.05,
                    'z': 1.5,
                    'text': 'z',
                    'font': {
                        'color': 'black',
                        'size': 40
                        }
                    }
                ]
            },
        'width': 1000,
        'height': 1000,
        'autosize': False,
        'margin': {
            'l': 0,
            'r': 0,
            'b': 0,
            't': 0,
            'pad': 1
            }
        }

    axis_x_3D = {
        'line': {
            'color': 'rgb(0, 0, 0)',
            'width': 3
            },
        'mode': 'lines',
        'type': 'scatter3d',
        'x': [1, 1.5],
        'y': [0, 0],
        'z': [0, 0],
        'showlegend': False
        }

    axis_y_3D = {
        'line': {
            'color': 'rgb(0, 0, 0)',
            'width': 3
            },
        'mode': 'lines',
        'type': 'scatter3d',
        'x': [0, 0],
        'y': [1, 1.5],
        'z': [0, 0],
        'showlegend': False
        }

    axis_z_3D = {
        'line': {
            'color': 'rgb(0, 0, 0)',
            'width': 3
            },
        'mode': 'lines',
        'type': 'scatter3d',
        'x': [0, 0],
        'y': [0, 0],
        'z': [1, 1.5],
        'showlegend': False
        }
    
    colorbar = {
        'outlinecolor': 'black',
        'outlinewidth': 1.5,
        'len': 0.8,
        'tickfont': {
            'color': 'black',
            'size': 20
            }
        }
    
    def __init__(self, settings):
        """ Construction method of the class.

        Parameters
        ----------

        settings: dict
            Dictionary containing the general settings used to plot the
            results.

        """
        BasicPlotter.__init__(self)
        self.__mpl = settings['mpl']
        self.__plotly = settings['plotly']
        self.__dpi = settings['dpi']
        self.__2D = settings['2D']
        self.__3D = settings['3D']
        self.__projection = settings['projection']
        return

    def read_data(self, filename):
        """
        """
        error = None
        # Read the data
        self.data = QuantasHDF5Reader(filename)

        try:
            self.data.read()
        except OSError:
            error = 'This is not a QUANTAS hdf5 file'
            return error

        if not 'Phase velocity' in self.data.info:
            error = 'This file does not contain the required results'
            return error

        # Save a basic name for plotting
        self._basename = filename.replace('.hdf5', '')
        self._basename = self._basename.replace('_SEISMIC', '')
        # Get the number of theta and phi angular values
        self._ntheta = 0
        self._nphi = 0
        theta = self.data.get_data_column('primary', 0)
        # Find nphi
        stop = False
        while not stop:
            if theta[self._nphi] != 0.:
                stop = True
                continue
            self._nphi += 1
        # Find ntheta
        self._ntheta = int(theta.shape[0] / self.nphi)
        return

    @property
    def ntheta(self):
        return self._ntheta

    @property
    def nphi(self):
        return self._nphi

    def run(self):
        """ Start making the various 3D and 2D plots in an automated way. """
        for key in self.datasets:
            # Make 3D plots
            if self.__plotly:
                if self.__3D:
                    msg = '- making 3D plots of {}'
                    self.echo(msg.format(self.datasets[key]['name']))
                    self.make_3D_plot(key)
            # Make 2D plots
            if self.__mpl:
                if self.__2D:
                    msg = '- making 2D plots of {}'
                    self.echo(msg.format(self.datasets[key]['name']))
                    self.make_2D_plot(key)

        self.echo('- making 2D plots of different ratios:')
        self.echo('   * S-wave anisotropy = 200*(v_s1-v_s2)/(v_s1+v_s2)')
        self.echo('   * v_P/v_s1')
        self.echo('   * v_P/v_s2')
        self.make_2D_plot_ratios()
        
        return

    def anisotropy(self, xmin, xmax):
        """ Return the percentage of anisotropy A of a certain quantity x,
        defined as:

        .. math::

           A = 200 \\frac{x_{max}-x_{min}}{x_{max}+x_{min}}

        Parameters
        ----------

        xmin: floar or ndarray
            Minimum value of a certain quantity.

        xmax: floar or ndarray
            Maximum value of a certain quantity.

        Returns
        -------

        float or ndarray
            Percentage of anisotropy.

        """
        return 200*(xmax-xmin)/(xmax+xmin)

    def make_3D_plot(self, quantity):
        """ Make and save a set of three figures containind a 3D plot of a
        specific result. The plotly package is employed in this method.

        Parameters
        ----------

        quantity: int
            Index of the column of the results dataset that will be plotted.

        """

        # Read the dataset and create the meshgrids
        u, v = np.mgrid[0:np.pi/2:self.ntheta*1j, 0:2*np.pi:self.nphi*1j]
        x = np.sin(u) * np.cos(v)
        y = np.sin(u) * np.sin(v)
        z = np.cos(u)
        
        for p in ['ssecondary', 'fsecondary', 'primary']:
            c = self.data.get_data_column(p, quantity).reshape(self.ntheta,
                                                               self.nphi)
            if quantity == 13:
                c = np.log10(c)

            fig = go.Figure(layout=self.layout_3D)
            # Plot the results on a sphere
            fig.add_traces(
                [
                    go.Surface(x=x, y=y, z=z, surfacecolor=c,
                               colorscale=self.datasets[quantity]['colorscale'],
                               colorbar=self.colorbar,
                               cmin=self.datasets[quantity]['cmin'][p],
                               cmax=self.datasets[quantity]['cmax'][p],
                            ),
                    go.Surface(x=-x, y=-y, z=-z, surfacecolor=c,
                               colorscale=self.datasets[quantity]['colorscale'],
                               cmin=self.datasets[quantity]['cmin'][p],
                               cmax=self.datasets[quantity]['cmax'][p],
                               showscale=False)
                    ]
                )

            # Add x, y and z axes
            fig.add_traces([
                self.axis_x_3D,
                self.axis_y_3D,
                self.axis_z_3D
                ])

            # Add the title to the colorbar
            fig = self.add_colorbar_title(
                fig, self.datasets[quantity]['title']
                )

            # Add a title at the bottom
            fig = self.add_figure_title_3D(fig, self.titles[p])

            # Save the figure
            fig.write_image(
                '{}_{}_3D_{}.png'.format(
                    self._basename, self.datasets[quantity]['suffix'], p)
                )
        return

    def make_2D_plot(self, quantity):
        """ Make and save a figure containing three 2D (polar) plots of a
        specific result, subdivided into slow secondary, fast secondary and
        primary. The matplotlib package is employed in this method.

        Two type of 2D projection on the xy plane (upper part of the
        sphere) are implemented, i.e., 'equal_area' (default) and 'stereo'
        projections.

        Parameters
        ----------

        quantity: int
            Index of the column of the results dataset that will be plotted.

        """
        
        # Read the dataset and create the meshgrids for stereo picture
        u, v = np.mgrid[0:np.pi/2:self.ntheta*1j, 0:2*np.pi:self.nphi*1j]
        if self.__projection == 'stereo':
            x = np.tan(0.5*u) * np.cos(v)
            y = np.tan(0.5*u) * np.sin(v)
        elif self.__projection == 'eqar':
            x = np.sqrt(2)*np.sin(0.5*u)*np.cos(v)
            y = np.sqrt(2)*np.sin(0.5*u)*np.sin(v)
        
        fig, axes = plt.subplots(1, 3, figsize=(16, 4))

        c = 0  # Pointer
        for p in ['ssecondary', 'fsecondary', 'primary']:
            # Collect the data
            z = self.data.get_data_column(p, quantity).reshape(self.ntheta,
                                                               self.nphi)
            levels=256
            if quantity == 13:
                z = np.log10(z)
                levels = np.linspace(
                    self.datasets[quantity]['cmin'][p],
                    self.datasets[quantity]['cmax'][p],
                    256
                    )
            
            cmap = LinearSegmentedColormap.from_list(
                'Custom', self.datasets[quantity]['colorscale'], N=256)

            # Select the ax
            ax = axes[c]

            # Plot the filled contour map
            cs = ax.contourf(x, y, z,
                             levels=levels,
                             cmap=cmap, extend='both',
                             vmin=self.datasets[quantity]['cmin'][p],
                             vmax=self.datasets[quantity]['cmax'][p]
                             )

            # Plot isolines 
            if self.datasets[quantity]['levels'][c] != 0:
                clevels = self.datasets[quantity]['levels'][c]
                cl = ax.contour(x, y, z, colors='white', linewidths=1,
                                levels=clevels)

            # Plot a symbol for the minimum and maximum values of the quantity
            i, j = np.unravel_index(z.argmax(), z.shape)
            ax.scatter(x[i, j], y[i, j], s=50, c='k', marker='s',
                       edgecolor='white')
            i, j = np.unravel_index(z.argmin(), z.shape)
            ax.scatter(x[i, j], y[i, j], s=40, c='white', edgecolor='k')

            if self.datasets[quantity]['anisotropy']:
                msg = '   * {: >15}: anisotropy = {: >5.1f} %'
                self.echo(msg.format(
                    self.titles[p], self.anisotropy(z.min(), z.max()))
                          )
            
            cbar = fig.colorbar(cs, ax=ax,
                                format=tick.FormatStrFormatter('%.2f'))
            cbar.set_label(self.datasets[quantity]['title2D'], size=16)

            ax.arrow(-0.95, -0.95, 0.35, 0, width=0.01, color='black')
            ax.text(-0.55, -0.98, 'x', size=14)
            ax.arrow(-0.95, -0.95, 0, 0.35, width=0.01, color='black')
            ax.text(-0.98, -0.51, 'y', size=14)

            # Remove the axis
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            plt.setp(ax.spines.values(), visible=False)

            # Add the calculated quantity
            ax.set_xlabel('{}'.format(self.titles[p]), size=16)

            # Set the graph limit
            ax.set_xlim(-1, 1)
            ax.set_ylim(-1, 1)
            ax.set_aspect('equal')

            c += 1

        fig.tight_layout()
        
        fig.savefig('{}_{}_2D_{}.png'.format(self._basename,
                                             self.datasets[quantity]['suffix'],
                                             self.__projection),
                    dpi=self.__dpi)
        return

    def make_2D_plot_ratios(self):
        """ Make and save a figure containing three 2D (polar) plots of three
        specific results:

        - S-wave anisotropy, calculated as AVs = 200*(Vs1-Vs2)/(Vs1+Vs2);
        - v_P / v_s1 ratio (fast secondary);
        - v_P / v_s2 ratio (slow secondary).

        """
        
        # Read the dataset and create the meshgrids for stereo picture
        u, v = np.mgrid[0:np.pi/2:self.ntheta*1j, 0:2*np.pi:self.nphi*1j]
        if self.__projection == 'stereo':
            x = np.tan(0.5*u) * np.cos(v)
            y = np.tan(0.5*u) * np.sin(v)
        elif self.__projection == 'eqar':
            x = np.sqrt(2)*np.sin(0.5*u)*np.cos(v)
            y =  np.sqrt(2)*np.sin(0.5*u)*np.sin(v)
        
        fig, axes = plt.subplots(1, 3, figsize=(16,4))

        # Take the phase velocities
        vp = self.data.get_data_column('primary', 2).reshape(self.ntheta,
                                                            self.nphi)
        vs1 = self.data.get_data_column('fsecondary', 2).reshape(self.ntheta,
                                                                self.nphi)
        vs2 = self.data.get_data_column('ssecondary', 2).reshape(self.ntheta,
                                                                self.nphi)

        # General settings
        levels=256
        
        # Plot 1: S-wave anisotropy
        ax = axes[0]
        A = self.anisotropy(vs2, vs1)
        #AVs = 200 * (vs1 - vs2) / (vs1 + vs2)

        cmap = LinearSegmentedColormap.from_list('Custom', default_colorscale,
                                                 N=levels)        

        # Filled contours
        cs = ax.contourf(x, y, A, levels=levels, cmap=cmap, extend='both')

        # Levels
        cl = ax.contour(x, y, A, colors='white', linewidths=1, levels=10)

        # Min/Max values
        i, j = np.unravel_index(A.argmax(), A.shape)
        ax.scatter(x[i, j], y[i, j], s=50, c='k', marker='s',
                   edgecolor='white')
        i, j = np.unravel_index(A.argmin(), A.shape)
        ax.scatter(x[i, j], y[i, j], s=40, c='white', edgecolor='k')
        

        cbar = fig.colorbar(cs, ax=ax,
                            format=tick.FormatStrFormatter('%.0f'))
        cbar.set_label('S-wave anisotropy (%)', size=16)

        # Plot 2: V_P / V_S1
        ax = axes[1]
        ratio = vp / vs1

        # Filled contours
        cs = ax.contourf(x, y, ratio, levels=levels, cmap=cmap, extend='both')

        # Levels
        cl = ax.contour(x, y, ratio, colors='white', linewidths=1, levels=10)

        # Min/Max values
        i, j = np.unravel_index(ratio.argmax(), ratio.shape)
        ax.scatter(x[i, j], y[i, j], s=50, c='k', marker='s',
                   edgecolor='white')
        i, j = np.unravel_index(ratio.argmin(), ratio.shape)
        ax.scatter(x[i, j], y[i, j], s=40, c='white', edgecolor='k')
        

        cbar = fig.colorbar(cs, ax=ax,
                            format=tick.FormatStrFormatter('%.1f'))
        cbar.set_label('$v_P / v_{s1}$', size=16)

        # Plot 3: V_P / V_S2
        ax = axes[2]
        ratio = vp / vs2

        # Filled contours
        cs = ax.contourf(x, y, ratio, levels=levels, cmap=cmap, extend='both')

        # Levels
        cl = ax.contour(x, y, ratio, colors='white', linewidths=1, levels=10)

        # Min/Max values
        i, j = np.unravel_index(ratio.argmax(), ratio.shape)
        ax.scatter(x[i, j], y[i, j], s=50, c='k', marker='s',
                   edgecolor='white')
        i, j = np.unravel_index(ratio.argmin(), ratio.shape)
        ax.scatter(x[i, j], y[i, j], s=40, c='white', edgecolor='k')

        cbar = fig.colorbar(cs, ax=ax,
                            format=tick.FormatStrFormatter('%.1f'))
        cbar.set_label('$v_P / v_{s2}$', size=16)
        
        for ax in axes:
            ax.arrow(-0.95, -0.95, 0.35, 0, width=0.01, color='black')
            ax.text(-0.55, -0.98, 'x', size=14)
            ax.arrow(-0.95, -0.95, 0, 0.35, width=0.01, color='black')
            ax.text(-0.98, -0.51, 'y', size=14)

            ax.set_axis_off()
            
            ax.set_xlim(-1, 1)
            ax.set_ylim(-1, 1)
            ax.set_aspect('equal')


        fig.tight_layout()
        
        fig.savefig('{}_{}_2D.png'.format(self._basename, 'ratios'),
                    dpi=self.__dpi)
        return

    def make_2D_plot_polarization(self):
        """ Make and save a figure containing three 2D (polar) plots of the
        phase velocities and the acoustic wave polarization.
        """
        quantity = 2
        
        # Read the dataset and create the meshgrids for stereo picture
        u, v = np.mgrid[0:np.pi/2:self.ntheta*1j, 0:2*np.pi:self.nphi*1j]
        if self.__projection == 'stereo':
            x = np.tan(0.5*u) * np.cos(v)
            y = np.tan(0.5*u) * np.sin(v)
        elif self.__projection == 'eqar':
            x = np.sqrt(2)*np.sin(0.5*u)*np.cos(v)
            y =  np.sqrt(2)*np.sin(0.5*u)*np.sin(v)
        
        fig, axes = plt.subplots(1, 3, figsize=(16,4))

        c = 0  # Pointer
        for p in ['ssecondary', 'fsecondary', 'primary']:

            # Collect the phase velocity
            z = self.data.get_data_column(p, quantity).reshape(self.ntheta,
                                                               self.nphi)

            # Collect the phase polarization
            u = self.data.get_data_column(p, 9).reshape(self.ntheta,
                                                        self.nphi)
            v = self.data.get_data_column(p, 10).reshape(self.ntheta,
                                                        self.nphi)

            levels=256
            if quantity == 13:
                z = np.log10(z)
                levels = np.linspace(
                    self.datasets[quantity]['cmin'][p],
                    self.datasets[quantity]['cmax'][p],
                    256
                    )
            
            cmap = LinearSegmentedColormap.from_list(
                'Custom', self.datasets[quantity]['colorscale'], N=256)

            # Select the ax
            ax = axes[c]

            # Plot the filled contour map
            cs = ax.contourf(x, y, z,
                             levels=levels,
                             cmap=cmap, extend='both',
                             vmin=self.datasets[quantity]['cmin'][p],
                             vmax=self.datasets[quantity]['cmax'][p]
                             )

            # Plot isolines 
            if self.datasets[quantity]['levels'][c] != 0:
                clevels = self.datasets[quantity]['levels'][c]
                cl = ax.contour(x, y, z, colors='white', linewidths=1,
                                levels=clevels)
            
            cbar = fig.colorbar(cs, ax=ax,
                                format=tick.FormatStrFormatter('%.2f'))
            cbar.set_label(self.datasets[quantity]['title2D'], size=16)

            # Plot the quivers
            ax.quiver(x[::20, ::20], y[::20, ::20],  # position of the arrow
                      u[::20, ::20], v[::20, ::20],  # direction of the arrow
                      pivot='mid')

            ax.arrow(-0.95, -0.95, 0.35, 0, width=0.01, color='black')
            ax.text(-0.55, -0.98, 'x', size=14)
            ax.arrow(-0.95, -0.95, 0, 0.35, width=0.01, color='black')
            ax.text(-0.98, -0.51, 'y', size=14)

            # Remove the axis
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            plt.setp(ax.spines.values(), visible=False)

            # Add the calculated quantity
            ax.set_xlabel('{}'.format(self.titles[p]), size=16)

            # Set the graph limit
            ax.set_xlim(-1, 1)
            ax.set_ylim(-1, 1)
            ax.set_aspect('equal')

            c += 1

        fig.tight_layout()
        
        fig.savefig('{}_{}_2D_{}_polarization.png'.format(
            self._basename, self.datasets[quantity]['suffix'],
            self.__projection), dpi=self.__dpi)
        return

    def add_figure_title_3D(self, figure, text):
        """ Add a title to the 3D figure.

        Parameters
        ----------

        figure: plotly.graphical_object.Figure
            The figure on which the title will be added.

        text: str
            The title that will be written.

        Returns
        -------
        figure: plotly.graphical_object.Figure
            Updated figure with the title .

        """
        figure.update_layout(
            title={
                'text': text,
                'x': 0.5,
                'y': 0.05,
                'xanchor': 'center',
                'yanchor': 'bottom',
                'font': {
                    'color': 'black',
                    'size': 50
                    }
                }
            )
        return figure

    def add_colorbar_title(self, figure, text):
        """ Add a title to the colorbar of the 3D figure. Due to the not so
        nice default positioning of title in the plotly colorbar, it was
        implemented this method that includes an annotation near the
        colorbar.

        Parameters
        ----------

        figure: plotly.graphical_object.Figure
            The figure on which the colorbar title will be added.

        text: str
            The title of the colorbar that will be written on the figure, near
            the colorbar.

        Returns
        -------

        figure: plotly.graphical_object.Figure
            Updated figure with the colorbar title.

        """
        figure.add_annotation(
            text=text,
            font=dict(
                color='black',
                size=20,
                ),
            showarrow=False,
            xref='paper', yref='paper',
            x=1.06, y=0.94
            )
        return figure

