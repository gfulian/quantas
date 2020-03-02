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
try:
    import matplotlib.pyplot as plt
    plotting = True
except ImportError:
    plotting = False
        

def make_polar_figure():
    """
    Create a generic polar figure with Matplotlib.

    Returns
    -------
    fig: matplotlib.pyplot.figure
        Figure object.

    ax1: matplotlib.pyplot.axis
        Axis for the (xy) plane plot.

    ax2: matplotlib.pyplot.axis
        Axis for the (xz) plane plot.

    ax3: matplotlib.pyplot.axis
        Axis for the (yz) plane plot.
    """
    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(18,6), subplot_kw={'projection':'polar'})
    ax1.grid(linestyle='dotted')
    ax1.set_theta_zero_location("N")
    ax2.grid(linestyle='dotted')
    ax2.set_theta_zero_location("N")
    ax3.grid(linestyle='dotted')
    ax3.set_theta_zero_location("N")
    return fig, ax1, ax2, ax3


def add_description(ax, text):
    """
    Add a description in the axis plot.

    Parameters
    ----------

    ax: matplotlib.pyplot.axis
        Axis where the description will be written.

    text: str
        Textto be inserted in the axis.
    """
    ax.text(0.2, -0.05, text, ha='center', va='center',
            transform=ax.transAxes, fontsize=12)
    return


def young_polar_plot(xy, xz, yz):
    """
    Plot the Young's modulus on the (xy), (xz) and (yz) planes.

    Parameters
    ----------

    xy: ndarray
        Array containing the polar coordinates and the values of Young's
        modulus on the (xy) plane.

    xz: ndarray
        Array containing the polar coordinates and the values of Young's
        modulus on the (xz) plane.

    yz: ndarray
        Array containing the polar coordinates and the values of Young's
        modulus on the (yz) plane.
        
    Returns
    -------

    fig: matplotlib.pyplot.figure
        Figure object containing the polar plot of Young's modulus.
    """
    title = "Young's modulus"
    fig, ax1, ax2, ax3 = make_polar_figure()
    ax1.plot(np.radians(xy[:, 0]), xy[:, 1], color='#009010', lw=2)
    ax1.set_rmax(1.1 * np.max(xy[:, 1]))
    add_description(ax1, title + " (xy)")

    ax2.plot(np.radians(xz[:, 0]), xz[:, 1], color='#009010', lw=2)
    ax2.set_rmax(1.1 * np.max(xz[:, 1]))
    add_description(ax2, title + " (xz)")

    ax3.plot(np.radians(yz[:, 0]), yz[:, 1], color='#009010', lw=2)
    ax3.set_rmax(1.1 * np.max(yz[:, 1]))
    add_description(ax3, title + " (yz)")

    fig.tight_layout()
    return fig


def compressibility_polar_plot(xy, xz, yz):
    """
    Plot the linear compressibility (LC) on the (xy), (xz) and (yz) planes.

    Parameters
    ----------

    xy: ndarray
        Array containing the polar coordinates and the values of linear
        compressibility on the (xy) plane.

    xz: ndarray
        Array containing the polar coordinates and the values of linear
        compressibility on the (xz) plane.

    yz: ndarray
        Array containing the polar coordinates and the values of linear
        compressibility on the (yz) plane.
        
    Returns
    -------

    fig: matplotlib.pyplot.figure
        Figure object containing the polar plot of linear compressibility.
    """
    title = "Linear compressibility"
    fig, ax1, ax2, ax3 = make_polar_figure()
    ax1.plot(np.radians(xy[:, 0]), xy[:, 1], color='#009010', lw=2)
    ax1.plot(np.radians(xy[:, 0]), xy[:, 2], color='red', lw=2)
    ax1.set_rmax(1.1 * np.max(xy[:, 1]))
    add_description(ax1, title + " (xy)")

    ax2.plot(np.radians(xz[:, 0]), xz[:, 1], color='#009010', lw=2)
    ax2.plot(np.radians(xz[:, 0]), xz[:, 2], color='red', lw=2)
    ax2.set_rmax(1.1 * np.max(xz[:, 1]))
    add_description(ax2, title + " (xz)")

    ax3.plot(np.radians(yz[:, 0]), yz[:, 1], color='#009010', lw=2)
    ax3.plot(np.radians(yz[:, 0]), yz[:, 2], color='red', lw=2)
    ax3.set_rmax(1.1 * np.max(yz[:, 1]))
    add_description(ax3, title + " (yz)")

    fig.legend(("Positive", "Negative"), fontsize=12)
    fig.tight_layout()
    return fig


def shear_polar_plot(xy, xz, yz):
    """
    Plot the shear modulus on the (xy), (xz) and (yz) planes.

    Parameters
    ----------

    xy: ndarray
        Array containing the polar coordinates and the values of shear
        modulus on the (xy) plane.

    xz: ndarray
        Array containing the polar coordinates and the values of shear
        modulus on the (xz) plane.

    yz: ndarray
        Array containing the polar coordinates and the values of shear
        modulus on the (yz) plane.
        
    Returns
    -------

    fig: matplotlib.pyplot.figure
        Figure object containing the polar plot of shear modulus.
    """
    title = "Shear modulus"
    fig, ax1, ax2, ax3 = make_polar_figure()
    ax1.plot(np.radians(xy[:, 0]), xy[:, 1], "--", color='#009010', lw=2)
    ax1.plot(np.radians(xy[:, 0]), xy[:, 2], color='red', lw=2)
    ax1.set_rmax(1.1 * np.max(xy[:, 2]))
    add_description(ax1, title + " (xy)")

    ax2.plot(np.radians(xz[:, 0]), xz[:, 1], "--", color='#009010', lw=2)
    ax2.plot(np.radians(xz[:, 0]), xz[:, 2], color='red', lw=2)
    ax2.set_rmax(1.1 * np.max(xz[:, 2]))
    add_description(ax2, title + " (xz)")

    ax3.plot(np.radians(yz[:, 0]), yz[:, 1], "--", color='#009010', lw=2)
    ax3.plot(np.radians(yz[:, 0]), yz[:, 2], color='red', lw=2)
    ax3.set_rmax(1.1 * np.max(yz[:, 2]))
    add_description(ax3, title + " (yz)")

    fig.legend(("Minimum", "Maximum"), fontsize=12)
    fig.tight_layout()
    return fig


def poisson_polar_plot(xy, xz, yz):
    """
    Plot the Poisson's ratio on the (xy), (xz) and (yz) planes.

    Parameters
    ----------

    xy: ndarray
        Array containing the polar coordinates and the values of Poisson's
        ratio on the (xy) plane.

    xz: ndarray
        Array containing the polar coordinates and the values of Poisson's
        ratio on the (xz) plane.

    yz: ndarray
        Array containing the polar coordinates and the values of Poisson's
        ratio on the (yz) plane.
        
    Returns
    -------

    fig: matplotlib.pyplot.figure
        Figure object containing the polar plot of Poisson's ratio.
    """
    title = "Poisson's ratio"
    fig, ax1, ax2, ax3 = make_polar_figure()
    ax1.plot(np.radians(xy[:, 0]), xy[:, 1], "--", color='blue', lw=2)
    ax1.plot(np.radians(xy[:, 0]), xy[:, 2], color='#009010', lw=2)
    ax1.plot(np.radians(xy[:, 0]), xy[:, 3], color='red', lw=2)
    ax1.set_rmax(1.1 * np.max(xy[:, 3]))
    add_description(ax1, title + " (xy)")

    ax2.plot(np.radians(xz[:, 0]), xz[:, 1], "--", color='blue', lw=2)
    ax2.plot(np.radians(xz[:, 0]), xz[:, 2], color='#009010', lw=2)
    ax2.plot(np.radians(xz[:, 0]), xz[:, 3], color='red', lw=2)
    ax2.set_rmax(1.1 * np.max(xz[:, 3]))
    add_description(ax2, title + " (xz)")

    ax3.plot(np.radians(yz[:, 0]), yz[:, 1], "--", color='blue', lw=2)
    ax3.plot(np.radians(yz[:, 0]), yz[:, 2], color='#009010', lw=2)
    ax3.plot(np.radians(yz[:, 0]), yz[:, 3], color='red', lw=2)
    ax3.set_rmax(1.1 * np.max(yz[:, 3]))
    add_description(ax3, title + " (yz)")

    fig.legend(("Negative", "Minimum", "Maximum"), fontsize=12)
    fig.tight_layout()
    return fig


def seismic_waves_polar_plot(xy, xz, yz):
    """
    Plot the seismic wave velocities on the (xy), (xz) and (yz) planes.

    Parameters
    ----------

    xy: ndarray
        Array containing the polar coordinates and the values of seismic wave
        velocities on the (xy) plane.

    xz: ndarray
        Array containing the polar coordinates and the values of seismic wave
        velocities on the (xz) plane.

    yz: ndarray
        Array containing the polar coordinates and the values of seismic wave
        velocities on the (yz) plane.
        
    Returns
    -------

    fig: matplotlib.pyplot.figure
        Figure object containing the polar plot of seismic wave velocities.
    """
    title = "Seismic waves"
    fig, ax1, ax2, ax3 = make_polar_figure()
    ax1.plot(np.radians(xy[:, 0]), xy[:, 1], "--", color='blue', lw=2)
    ax1.plot(np.radians(xy[:, 0]), xy[:, 2], color='#009010', lw=2)
    ax1.plot(np.radians(xy[:, 0]), xy[:, 3], color='red', lw=2)
    ax1.set_rmax(1.1 * np.max(xy[:, 3]))
    add_description(ax1, title + " (xy)")

    ax2.plot(np.radians(xz[:, 0]), xz[:, 1], "--", color='blue', lw=2)
    ax2.plot(np.radians(xz[:, 0]), xz[:, 2], color='#009010', lw=2)
    ax2.plot(np.radians(xz[:, 0]), xz[:, 3], color='red', lw=2)
    ax2.set_rmax(1.1 * np.max(xz[:, 3]))
    add_description(ax2, title + " (xz)")

    ax3.plot(np.radians(yz[:, 0]), yz[:, 1], "--", color='blue', lw=2)
    ax3.plot(np.radians(yz[:, 0]), yz[:, 2], color='#009010', lw=2)
    ax3.plot(np.radians(yz[:, 0]), yz[:, 3], color='red', lw=2)
    ax3.set_rmax(1.1 * np.max(yz[:, 3]))
    add_description(ax3, title + " (yz)")

    fig.legend(("$v_{s1}$", "$v_{s2}$", "$v_{p}$"), fontsize=12)
    fig.tight_layout()
    return fig
