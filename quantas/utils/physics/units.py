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
This module contains functions used to convert values between different
scales.
"""
from scipy import constants as cs
import numpy as np
from numpy import power

# Physical Constants
AMU = cs.atomic_mass        # (kg)
N = cs.Avogadro             # Avogadro's number
h = cs.Planck               # (J s)
Hbar = cs.hbar              # (J s)
kB = cs.Boltzmann           # (J/K)
c = cs.c                    # (m/s)

# Metric conversions
nano = power(10.,-9.)
micro = power(10.,-6.)
milli = power(10.,-3.)
kilo = power(10.,3.)
mega = power(10.,6.)
giga = power(10.,9.)
tera = power(10.,12.)

# Temperature
zero_Celsius = 273.15

# Units typically used in QM software
#
# 1. Distances
Angstrom = cs.value(u'Angstrom star')   # (m)
Bohr = cs.value(u'Bohr radius')         # (m)
#
# 2. Electronic Energy
eV = cs.value(u'electron volt')                     # (J)
Ha = cs.value(u'Hartree energy')                    # (J)
Ry = cs.value(u'Rydberg constant times hc in J')    # (J)
#
def convert_temperature(val, old_scale, new_scale):
    """
    Convert from a temperature scale to another one among Celsius, Kelvin,
    Fahrenheit and Rankine scales. Supported scales are::

        Celsius ('Celsius', 'celsius', 'C' or 'c')
        Kelvin ('Kelvin', 'kelvin', 'K', 'k')
        Fahrenheit ('Fahrenheit', 'fahrenheit', 'F' or 'f')
        Rankine ('Rankine', 'rankine', 'R', 'r')

    Parameters
    ----------
    val : array_like
        Value(s) of the temperature(s) to be converted expressed in the
        original scale.

    old_scale: str
        Specifies as a string the original scale from which the temperature
        value(s) will be converted.

    new_scale: str
        Specifies as a string the new scale to which the temperature
        value(s) will be converted.

    Returns
    -------
    res : float or array of floats
        Value(s) of the converted temperature(s) expressed in the new scale.


    Examples
    --------
    >>> from quantascli.utils.units import convert_temperature
    >>> convert_temperature(np.array([-40, 40.0]), 'Celsius', 'Kelvin')
    array([ 233.15,  313.15])

    """
    # Convert from `old_scale` to Kelvin
    if old_scale.lower() in ['celsius', '°c', 'c']:
        tempo = np.asanyarray(val) + zero_Celsius
    elif old_scale.lower() in ['kelvin', 'k']:
        tempo = np.asanyarray(val)
    elif old_scale.lower() in ['fahrenheit', 'f']:
        tempo = (np.asanyarray(val) - 32.) * 5. / 9. + zero_Celsius
    elif old_scale.lower() in ['rankine', 'r']:
        tempo = np.asanyarray(val) * 5. / 9.
    else:
        raise NotImplementedError('%s scale is unsupported: supported scales '
                                  'are Celsius, Kelvin, Fahrenheit and '
                                  'Rankine' % old_scale)

    # and from Kelvin to `new_scale`.
    if new_scale.lower() in ['celsius', '°c', 'c']:
        res = tempo - zero_Celsius
    elif new_scale.lower() in ['kelvin', 'k']:
        res = tempo
    elif new_scale.lower() in ['fahrenheit', 'f']:
        res = (tempo - zero_Celsius) * 9. / 5. + 32.
    elif new_scale.lower() in ['rankine', 'r']:
        res = tempo * 9. / 5.
    else:
        raise NotImplementedError('%s scale is unsupported: supported scales '
                                  'are Celsius, Kelvin, Fahrenheit and '
                                  'Rankine' % new_scale)

    return res

def convert_energy(val, old_scale, new_scale):
    """
    Convert from an energy scale to another one among::

        Hartree ('hartree', 'ha')
        electron Volt ('electron volt', 'ev')
        Rydberg ('rydberg', 'ry')
        kilo Joule / mol ('kj/mol', 'kjmol')

    Parameters
    ----------
    val : array_like
        Energy value(s) to be converted expressed in the original scale.

    old_scale: str
        Specifies as a string the original scale from which the energy
        value(s) will be converted.

    new_scale: str
        Specifies as a string the new scale to which the energy
        value(s) will be converted.

    Returns
    -------
    res : float or array of floats
        Converted energy value(s) expressed in the new scale.

    Examples
    --------
    >>> from quantascli.utils.units import convert_energy
    >>> convert_energy(np.array([-1.0, 1.0]), 'Ha', 'eV')
    array([ -27.21138602,  27.21138602])

    """
    # Convert from `old_scale` to Hartree
    if old_scale.lower() in ['Hartree', 'Ha', 'ha']:
        energy = np.asanyarray(val)
    elif old_scale.lower() in ['electron volt', 'eV', 'ev']:
        energy = np.asanyarray(val) * (eV / Ha)
    elif old_scale.lower() in ['Rydberg', 'Ry', 'ry']:
        energy = np.asanyarray(val) * (Ry / Ha)
    elif old_scale.lower() in ['kJ/mol', 'kj/mol', 'kjmol']:
        energy = np.asanyarray(val) * (1000. / Ha / N)
    else:
        raise NotImplementedError('%s scale is unsupported: supported scales '
                                  'are Hartree, eV, Rydberg and kJ/mol'
                                  % old_scale)
    # and from Hartree to `new_scale`.
    if new_scale.lower() in ['Hartree', 'Ha', 'ha']:
        res = energy
    elif new_scale.lower() in ['electron volt', 'eV', 'ev']:
        res = energy * (Ha / eV)
    elif new_scale.lower() in ['Rydberg', 'Ry', 'ry']:
        res = energy * (Ha / Ry)
    elif new_scale.lower() in ['kJ/mol', 'kj/mol', 'kjmol']:
        res = energy * (Ha * N / 1000.)
    else:
        raise NotImplementedError('%s scale is unsupported: supported scales '
                                  'are Hartree, eV, Rydberg and kJ/mol'
                                  % new_scale)
    return res


def convert_pressure(val, old_scale, new_scale):
    """
    Convert from a pressure scale to another one among the following
    supported scales::

        Pascal ('pa')
        kiloPascal ('kpa')
        megaPascal ('mpa')
        gigaPascal ('gpa')
        bar ('bar')
        kilobar ('kbar')
        megabar ('mbar')
        gigabar ('gbar').

    Parameters
    ----------
    val : array_like
        Pressure value(s) to be converted with `float` type.

    old_scale: str
        Specifies as a string the original scale from which the energy
        value(s) will be converted.

    new_scale: str
        Specifies as a string the new scale to which the energy
        value(s) will be converted.

    Returns
    -------
    res : float or ndarray
        Converted pressure value(s) expressed in the new scale.

    Examples
    --------
    >>> from quantascli.utils.units import convert_pressure
    >>> convert_pressure([20., 25], 'GPa', 'Mbar')
    array([0.2 , 0.25])

    """
    # Convert from `old_scale` to Pascals
    if old_scale.lower() in ['gpa']:
        pressure = np.asanyarray(val) * giga
    elif old_scale.lower() in ['mpa']:
        pressure = np.asanyarray(val) * mega
    elif old_scale.lower() in ['kpa']:
        pressure = np.asanyarray(val) * kilo
    elif old_scale.lower() in ['pa']:
        pressure = np.asanyarray(val)
    elif old_scale.lower() in ['gbar']:
        pressure = np.asanyarray(val) * np.power(10., 14.)
    elif old_scale.lower() in ['mbar']:
        pressure = np.asanyarray(val) * np.power(10., 11.)
    elif old_scale.lower() in ['kbar']:
        pressure = np.asanyarray(val) * np.power(10., 8.)
    elif old_scale.lower() in ['bar']:
        pressure = np.asanyarray(val) * np.power(10., 5.)
    else:
        raise NotImplementedError('%s scale is unsupported' % old_scale)
    # and from Pascals to `new_scale`.
    if new_scale.lower() in ['gpa']:
        res = np.asanyarray(pressure) / giga
    elif new_scale.lower() in ['mpa']:
        res = np.asanyarray(pressure) / mega
    elif new_scale.lower() in ['kpa']:
        res = np.asanyarray(pressure) / kilo
    elif new_scale.lower() in ['pa']:
        res = np.asanyarray(pressure)
    elif new_scale.lower() in ['gbar']:
        res = np.asanyarray(pressure) / np.power(10., 14.)
    elif new_scale.lower() in ['mbar']:
        res = np.asanyarray(pressure) / np.power(10., 11.)
    elif new_scale.lower() in ['kbar']:
        res = np.asanyarray(pressure) / np.power(10., 8.)
    elif new_scale.lower() in ['bar']:
        res = np.asanyarray(pressure) / np.power(10., 5.)
    else:
        raise NotImplementedError('%s scale is unsupported' % new_scale)
    return res


def energy_to_pressure(val, energy_scale, volume_scale, pressure_scale):
    """
    Convert values from an energy/volume scale to a pressure one. Supported
    scales are, for energy::

        Hartree ('hartree', 'ha')
        electron Volt ('electron volt', 'ev')
        Rydberg ('rydberg', 'ry')

    for volume (as lenghts)::

        Angstrom ('angstrom', 'a')
        Bohr ('bohr')

    and for pressure::

        Pascal ('pa')
        kiloPascal ('kpa')
        megaPascal ('mpa')
        gigaPascal ('gpa')
        bar ('bar')
        kilobar ('kbar')
        megabar ('mbar')
        gigabar ('gbar')

    Parameters
    ----------
    val : array_like
        Energy value(s) to be converted expressed in the original scale.

    energy_scale: str
        Specifies as a string the energy scale from which the energy/volume
        value(s) will be converted.

    volume_scale: str
        Specifies as a string the volume scale from which the energy/volume
        value(s) will be converted.

    pressure_scale: str
        Specifies as a string the pressure scale to which the energy/volume
        value(s) will be converted.

    Returns
    -------
    res : float or array of floats
        Converted energy/volume value(s) expressed in the pressure scale.

    Examples
    --------
    >>> from quantascli.utils.units import energy_to_pressure
    >>> energy_to_pressure([0.001, 0.2], 'eV', 'angstrom', 'gpa')
    [  1.60210478 320.42095571]

    """
    res = np.asarray(val)
    # Convert energy in J
    if energy_scale.lower() in ['hartree', 'ha']:
        res *= Ha
    elif energy_scale.lower() in ['electron volt', 'ev']:
        res *= eV
    elif energy_scale.lower() in ['rydberg', 'ry']:
        res *= Ry
    else:
        raise NotImplementedError('%s scale is unsupported' % energy_scale)

    # Divide energy by volume, in m^3
    if volume_scale.lower() in ['angstrom', 'a', 'a^3']:
        res /= np.power(Angstrom, 3.)
    elif volume_scale.lower() in ['bohr', 'bohr^3']:
        res /= np.power(Bohr, 3.)
    else:
        raise NotImplementedError('%s scale is unsupported' % volume_scale)

    # Scale for pressure unit
    if pressure_scale.lower() in ['gpa']:
        res /= giga
    elif pressure_scale.lower() in ['mpa']:
        res /= mega
    elif pressure_scale.lower() in ['kpa']:
        res /= kilo
    elif pressure_scale.lower() in ['pa']:
        res /= 1.
    elif pressure_scale.lower() in ['gbar']:
        res /= (np.power(10., 14.))
    elif pressure_scale.lower() in ['mbar']:
        res /= (np.power(10., 11.))
    elif pressure_scale.lower() in ['kbar']:
        res /= (np.power(10., 8.))
    elif pressure_scale.lower() in ['bar']:
        res /= (np.power(10., 5.))
    else:
        raise NotImplementedError('%s scale is unsupported' % pressure_scale)

    return res

def pressure_to_energy(val, energy_scale, volume_scale, pressure_scale):
    """
    Convert values from a pressure scale to an energy/volume one Supported
    scales are, for energy::

        Hartree ('hartree', 'ha')
        electron Volt ('electron volt', 'ev')
        Rydberg ('rydberg', 'ry')

    for volume (as lenghts)::

        Angstrom ('angstrom', 'a')
        Bohr ('bohr')

    and for pressure::

        Pascal ('pa')
        kiloPascal ('kpa')
        megaPascal ('mpa')
        gigaPascal ('gpa')
        bar ('bar')
        kilobar ('kbar')
        megabar ('mbar')
        gigabar ('gbar').

    Parameters
    ----------
    val : array_like
        Pressure value(s) to be converted expressed in the original scale.

    energy_scale: str
        Specifies as a string the energy scale to which the pressure
        value(s) will be converted.

    volume_scale: str
        Specifies as a string the volume scale to which the pressure
        value(s) will be converted.

    pressure_scale: str
        Specifies as a string the pressure scale from which the value(s)
        will be converted.

    Returns
    -------
    res : float or ndarray
        Converted pressure value(s) expressed in the energy/volume scale.

    Examples
    --------
    >>> from quantascli.utils.units import pressure_to_energy
    >>> pressure_to_energy([1., 2.], 'ha', 'angstrom', 'gpa')
    [0.00022938 0.00045876]

    """
    res = np.asarray(val)
    # Convert pressures in Pascals
    if pressure_scale.lower() in ['gpa']:
        res *= giga
    elif pressure_scale.lower() in ['mpa']:
        res *= mega
    elif pressure_scale.lower() in ['kpa']:
        res *= kilo
    elif pressure_scale.lower() in ['pa']:
        res *= 1.
    elif pressure_scale.lower() in ['gbar']:
        res *= (np.power(10., 14.))
    elif pressure_scale.lower() in ['mbar']:
        res *= (np.power(10., 11.))
    elif pressure_scale.lower() in ['kbar']:
        res *= (np.power(10., 8.))
    elif pressure_scale.lower() in ['bar']:
        res *= (np.power(10., 5.))
    else:
        raise NotImplementedError('%s scale is unsupported' % pressure_scale)

    # Multiply by volume scale, in m^3
    if volume_scale.lower() in ['angstrom', 'a', 'a^3']:
        res *= np.power(Angstrom, 3.)
    elif volume_scale.lower() in ['bohr', 'bohr^3']:
        res *= np.power(Bohr, 3.)
    else:
        raise NotImplementedError('%s scale is unsupported' % volume_scale)

    # Convert energy in J
    if energy_scale.lower() in ['hartree', 'ha']:
        res /= Ha
    elif energy_scale.lower() in ['electron volt', 'ev']:
        res /= eV
    elif energy_scale.lower() in ['rydberg', 'ry']:
        res /= Ry
    else:
        raise NotImplementedError('%s scale is unsupported' % energy_scale)

    return res


def convert_frequency(val, old_scale, new_scale):
    """
    Convert from an energy scale to another one among::

        wavenumbers ('cm^-1')
        Hertz ('Hz', 'hz')
        TeraHertz ('THz', 'thz')

    Parameters
    ----------
    val : array_like
        Frequency value(s) to be converted.

    old_scale: str
        Specifies as a string the frequency scale from which the frequency
        value(s) will be converted.

    new_scale: str
        Specifies as a string the frequency scale to which the frequency
        value(s) will be converted.

    Returns
    -------
    res : float or ndarray
        Converted frequency value(s).

    Examples
    --------
    >>> from quantascli.utils.units import convert_frequency
    >>> convert_frequency(1500., 'cm^-1', 'thz')
    44.9688687

    """
    # Convert from `old_scale` to cm^-1
    if old_scale.lower() in ['cm^-1']:
        freq = np.asanyarray(val)
    elif old_scale.lower() in ['Hz', 'hz']:
        freq = np.asanyarray(val) / (c * 100.)
    elif old_scale.lower() in ['THz', 'thz']:
        freq = tera * np.asanyarray(val) / (c * 100.)
    else:
        raise NotImplementedError('%s scale is unsupported: supported scales '
                                  'are cm^-1, Hz, and THz'
                                  % old_scale)
    # and from cm^-1 to `new_scale`.
    if new_scale.lower() in ['cm^-1']:
        res = freq
    elif new_scale.lower() in ['Hz', 'hz']:
        res = freq * (c * 100.)
    elif new_scale.lower() in ['THz', 'thz']:
        res = freq * ( c * 100.) / tera
    else:
        raise NotImplementedError('%s scale is unsupported: supported scales '
                                  'are cm^-1, Hz, and THz'
                                  % new_scale)
    return res

def convert_lenght(val, old_scale, new_scale):
    """
    Convert from an lenght scale to another one. Supported scales are::

        metre ('metre', 'm')
        Angstrom ('angstrom', 'a')
        Bohr ('bohr')

    Parameters
    ----------
    val : array_like
        Lenght value(s) to be converted.

    old_scale: str
        Specifies as a string the lenght scale from which the value(s) will
        be converted.

    new_scale: str
        Specifies as a string the lenght scale to which the value(s) will
        be converted.

    Returns
    -------
    res : float or array of floats
        Converted lenght value(s) expressed in the new scale.

    Examples
    --------
    >>> from quantascli.utils.units import convert_lenght
    >>> convert_lenght(20., 'Bohr', 'm')
    1.058354421806e-09
    """
    # Convert from `old_scale` to metres
    if old_scale.lower() in ['metre', 'm']:
        lenght = np.asanyarray(val)
    elif old_scale.lower() in ['angstrom', 'a']:
        lenght = np.asanyarray(val) * Angstrom
    elif old_scale.lower() in ['bohr']:
        lenght = np.asanyarray(val) * Bohr
    else:
        raise NotImplementedError('%s scale is unsupported' % old_scale)
    # and from metres to `new_scale`
    if new_scale.lower() in ['metre', 'm']:
        res = lenght
    elif new_scale.lower() in ['angstrom', 'a']:
        res = lenght / Angstrom
    elif new_scale.lower() in ['bohr']:
        res = lenght / Bohr
    else:
        raise NotImplementedError('%s scale is unsupported' % new_scale)

    return res

def convert_volume(val, old_scale, new_scale):
    """
    Convert from a volume scale to another one. Supported scales are::

        metre ('metre', 'm')
        Angstrom ('angstrom', 'a')
        Bohr ('bohr')

    Parameters
    ----------
    val : array_like
        Volume value(s) to be converted.

    old_scale: str
        Specifies as a string the length scale from which the volume value(s)
        will be converted.

    new_scale: str
        Specifies as a string the lenght scale to which the volume value(s)
        will be converted.

    Returns
    -------
    res : float or array of floats
        Converted volume value(s) expressed in the new scale.

    Examples
    --------
    >>> from quantascli.utils.units import convert_volume
    >>> convert_volume(10.9897, 'Angstrom', 'bohr')
    74.1654978184323

    """
    # Convert from `old_scale` to metres
    if old_scale.lower() in ['metre', 'm']:
        volume = np.asanyarray(val)
    elif old_scale.lower() in ['angstrom', 'a']:
        volume = np.asanyarray(val) * np.power(Angstrom, 3.)
    elif old_scale.lower() in ['bohr']:
        volume = np.asanyarray(val) * np.power(Bohr, 3.)
    else:
        raise NotImplementedError('%s scale is unsupported' % old_scale)
    # and from metres to `new_scale`
    if new_scale.lower() in ['metre', 'm']:
        res = volume
    elif new_scale.lower() in ['angstrom', 'a']:
        res = volume / np.power(Angstrom, 3.)
    elif new_scale.lower() in ['bohr']:
        res = volume / np.power(Bohr, 3.)
    else:
        raise NotImplementedError('%s scale is unsupported' % new_scale)

    return res
