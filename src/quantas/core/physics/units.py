# -*- coding: utf-8 -*-

"""Unit conversions for scalar and array-valued scientific quantities.

The public functions use explicit ``from_scale`` and ``to_scale`` arguments and
preserve array shapes across thermodynamic, elastic, and spectroscopic units.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import constants as cs


__all__ = [
    "AMU",
    "N",
    "h",
    "Hbar",
    "kB",
    "c",
    "nano",
    "micro",
    "milli",
    "kilo",
    "mega",
    "giga",
    "tera",
    "zero_Celsius",
    "Angstrom",
    "Bohr",
    "eV",
    "Ha",
    "Ry",
    "convert_temperature",
    "convert_energy",
    "convert_energy_per_temperature",
    "convert_pressure",
    "energy_to_pressure",
    "pressure_to_energy",
    "convert_frequency",
    "convert_length",
    "convert_volume",
]


# Physical constants
AMU = cs.atomic_mass  # kg
N = cs.Avogadro  # mol^-1
h = cs.Planck  # J s
Hbar = cs.hbar  # J s
kB = cs.Boltzmann  # J K^-1
c = cs.c  # m s^-1


# Metric prefixes
nano = 1.0e-9
micro = 1.0e-6
milli = 1.0e-3
kilo = 1.0e3
mega = 1.0e6
giga = 1.0e9
tera = 1.0e12


# Temperature
zero_Celsius = 273.15


# Units typically used in quantum-mechanical codes
Angstrom = 1.0e-10  # m
Bohr = cs.value("Bohr radius")  # m
eV = cs.value("electron volt")  # J
Ha = cs.value("Hartree energy")  # J
Ry = cs.value("Rydberg constant times hc in J")  # J


_ENERGY_TO_J = {
    "j": 1.0,
    "joule": 1.0,
    "joules": 1.0,
    "hartree": Ha,
    "ha": Ha,
    "electron volt": eV,
    "electronvolt": eV,
    "ev": eV,
    "rydberg": Ry,
    "ry": Ry,
    "mha": milli * Ha,
    "millihartree": milli * Ha,
    "mev": milli * eV,
    "millielectronvolt": milli * eV,
    "mry": milli * Ry,
    "millirydberg": milli * Ry,
    "j/mol": 1.0 / N,
    "jmol": 1.0 / N,
    "joule/mol": 1.0 / N,
    "joule per mol": 1.0 / N,
    "kj/mol": 1000.0 / N,
    "kjmol": 1000.0 / N,
    "kilojoule/mol": 1000.0 / N,
    "kilojoule per mol": 1000.0 / N,
}

_PRESSURE_TO_PA = {
    "pa": 1.0,
    "pascal": 1.0,
    "kpa": kilo,
    "kilopascal": kilo,
    "mpa": mega,
    "megapascal": mega,
    "gpa": giga,
    "gigapascal": giga,
    "bar": 1.0e5,
    "kbar": 1.0e8,
    "kilobar": 1.0e8,
    "mbar": 1.0e11,
    "megabar": 1.0e11,
    "gbar": 1.0e14,
    "gigabar": 1.0e14,
}

_LENGTH_TO_M = {
    "metre": 1.0,
    "meter": 1.0,
    "m": 1.0,
    "centimetre": 1.0e-2,
    "centimeter": 1.0e-2,
    "cm": 1.0e-2,
    "millimetre": milli,
    "millimeter": milli,
    "mm": milli,
    "micrometre": micro,
    "micrometer": micro,
    "um": micro,
    "µm": micro,
    "nanometre": nano,
    "nanometer": nano,
    "nm": nano,
    "picometre": 1.0e-12,
    "picometer": 1.0e-12,
    "pm": 1.0e-12,
    "angstrom": Angstrom,
    "ang": Angstrom,
    "a": Angstrom,
    "å": Angstrom,
    "bohr": Bohr,
    "bohr radius": Bohr,
}

_FREQUENCY_TO_CM = {
    "cm^-1": 1.0,
    "cm-1": 1.0,
    "wavenumber": 1.0,
    "wavenumbers": 1.0,
    "hz": 1.0 / (c * 100.0),
    "hertz": 1.0 / (c * 100.0),
    "thz": tera / (c * 100.0),
    "terahertz": tera / (c * 100.0),
}


def _normalize_scale(scale: str) -> str:
    """
    Normalize a unit label.

    Parameters
    ----------
    scale : str
        Unit label to normalize.

    Returns
    -------
    str
        Normalized unit label.
    """
    return scale.strip().lower()


def _asarray(value: Any) -> np.ndarray:
    """
    Convert a value to a NumPy array of floats.

    Parameters
    ----------
    value : scalar or array_like
        Value to convert.

    Returns
    -------
    ndarray
        NumPy representation of the input value.
    """
    return np.asarray(value, dtype=float)


def _return(value: np.ndarray, original: Any) -> float | np.ndarray:
    """
    Return a scalar if the original input was scalar.

    Parameters
    ----------
    value : ndarray
        Converted value.
    original : scalar or array_like
        Original input value.

    Returns
    -------
    float or ndarray
        Converted value with a convenient shape.
    """
    if np.isscalar(original):
        return float(value)
    return value


def _linear_conversion(
    value: Any,
    from_scale: str,
    to_scale: str,
    factors: dict[str, float],
    reference_unit: str,
) -> float | np.ndarray:
    """
    Convert values using multiplicative factors to a reference unit.

    Parameters
    ----------
    value : scalar or array_like
        Value to convert.
    from_scale : str
        Unit of the input value.
    to_scale : str
        Unit of the output value.
    factors : dict
        Dictionary containing conversion factors to the reference unit.
    reference_unit : str
        Name of the reference unit, used only in error messages.

    Returns
    -------
    float or ndarray
        Converted value.

    Raises
    ------
    NotImplementedError
        If one of the requested scales is not available.
    """
    from_key = _normalize_scale(from_scale)
    to_key = _normalize_scale(to_scale)

    if from_key not in factors:
        raise NotImplementedError(
            f"{from_scale} scale is unsupported for this conversion."
        )
    if to_key not in factors:
        raise NotImplementedError(
            f"{to_scale} scale is unsupported for this conversion."
        )

    converted = _asarray(value) * factors[from_key] / factors[to_key]
    return _return(converted, value)


def convert_temperature(
    value: Any,
    from_scale: str,
    to_scale: str,
) -> float | np.ndarray:
    """
    Convert temperatures between Celsius, Kelvin, Fahrenheit, and Rankine.

    Parameters
    ----------
    value : scalar or array_like
        Temperature value or values to convert.
    from_scale : str
        Original temperature scale.
    to_scale : str
        Target temperature scale.

    Returns
    -------
    float or ndarray
        Converted temperature value or values.

    Raises
    ------
    NotImplementedError
        If one of the requested scales is not supported.

    Examples
    --------
    >>> convert_temperature([-40.0, 40.0], "Celsius", "Kelvin")
    array([233.15, 313.15])
    """
    from_key = _normalize_scale(from_scale)
    to_key = _normalize_scale(to_scale)
    array = _asarray(value)

    if from_key in {"celsius", "°c", "c"}:
        kelvin = array + zero_Celsius
    elif from_key in {"kelvin", "k"}:
        kelvin = array
    elif from_key in {"fahrenheit", "f"}:
        kelvin = (array - 32.0) * 5.0 / 9.0 + zero_Celsius
    elif from_key in {"rankine", "r"}:
        kelvin = array * 5.0 / 9.0
    else:
        raise NotImplementedError(
            f"{from_scale} scale is unsupported: supported scales are "
            "Celsius, Kelvin, Fahrenheit, and Rankine."
        )

    if to_key in {"celsius", "°c", "c"}:
        result = kelvin - zero_Celsius
    elif to_key in {"kelvin", "k"}:
        result = kelvin
    elif to_key in {"fahrenheit", "f"}:
        result = (kelvin - zero_Celsius) * 9.0 / 5.0 + 32.0
    elif to_key in {"rankine", "r"}:
        result = kelvin * 9.0 / 5.0
    else:
        raise NotImplementedError(
            f"{to_scale} scale is unsupported: supported scales are "
            "Celsius, Kelvin, Fahrenheit, and Rankine."
        )

    return _return(result, value)


def convert_energy(
    value: Any,
    from_scale: str,
    to_scale: str,
) -> float | np.ndarray:
    """
    Convert energies between common atomistic simulation units.

    Supported scales include Hartree, electronvolt, Rydberg, and kJ/mol.

    Parameters
    ----------
    value : scalar or array_like
        Energy value or values to convert.
    from_scale : str
        Original energy scale.
    to_scale : str
        Target energy scale.

    Returns
    -------
    float or ndarray
        Converted energy value or values.

    Examples
    --------
    >>> convert_energy([-1.0, 1.0], "Ha", "eV")
    array([-27.211386...,  27.211386...])
    """
    return _linear_conversion(value, from_scale, to_scale, _ENERGY_TO_J, "J")


def convert_energy_per_temperature(
    value: Any,
    from_scale: str,
    to_scale: str,
) -> float | np.ndarray:
    """
    Convert quantities with dimensions of energy per temperature.

    The conversion is suitable for entropy and heat capacities. Energy units
    may be expressed per cell, for example ``Ha cell^-1 K^-1``, or per mole,
    for example ``J mol^-1 K^-1``. The function converts the energy numerator
    and the width of the temperature interval independently. It does not
    change the normalization between a cell and a formula unit.

    Parameters
    ----------
    value : scalar or array_like
        Entropy or heat-capacity value or values to convert.
    from_scale : str
        Original energy-per-temperature scale.
    to_scale : str
        Target energy-per-temperature scale.

    Returns
    -------
    float or ndarray
        Converted value or values.

    Raises
    ------
    NotImplementedError
        If an energy or temperature scale is unsupported.

    Examples
    --------
    >>> convert_energy_per_temperature(1.0, "J/mol/K", "Ha/K")
    3.808798...e-07
    """
    from_energy, from_temperature = _split_energy_temperature_scale(from_scale)
    to_energy, to_temperature = _split_energy_temperature_scale(to_scale)

    converted = np.asarray(
        convert_energy(value, from_energy, to_energy),
        dtype=float,
    )
    converted *= _temperature_interval_to_kelvin(
        to_temperature
    ) / _temperature_interval_to_kelvin(from_temperature)
    return _return(converted, value)


def _split_energy_temperature_scale(scale: str) -> tuple[str, str]:
    """
    Split an energy-per-temperature label into numerator and denominator.

    Parameters
    ----------
    scale : str
        Unit label to split.

    Returns
    -------
    tuple of str
        Energy scale accepted by :func:`convert_energy` and temperature scale.

    Raises
    ------
    NotImplementedError
        If the label does not contain a supported energy scale.
    """
    normalized = scale.strip().lower()
    normalized = normalized.replace("−", "-").replace("⁻¹", "^-1")
    normalized = normalized.replace("degrees", "degree")
    normalized = normalized.replace("kilojoule", "kj").replace("joule", "j")
    normalized = normalized.replace("electron volt", "electronvolt")
    normalized = normalized.replace(" ", "").replace("per", "/")

    suffixes = (
        (
            "Celsius",
            (
                "/degreecelsius",
                "/celsius",
                "/degc",
                "/°c",
                "degreecelsius^-1",
                "degreecelsius-1",
                "celsius^-1",
                "celsius-1",
                "°c^-1",
                "°c-1",
            ),
        ),
        (
            "Fahrenheit",
            (
                "/degreefahrenheit",
                "/fahrenheit",
                "/degf",
                "/°f",
                "degreefahrenheit^-1",
                "degreefahrenheit-1",
                "fahrenheit^-1",
                "fahrenheit-1",
                "°f^-1",
                "°f-1",
            ),
        ),
        (
            "Rankine",
            (
                "/degreerankine",
                "/rankine",
                "/degr",
                "/°r",
                "degreerankine^-1",
                "degreerankine-1",
                "rankine^-1",
                "rankine-1",
                "°r^-1",
                "°r-1",
            ),
        ),
        ("K", ("/kelvin", "/k", "kelvin^-1", "kelvin-1", "k^-1", "k-1")),
    )

    temperature = "K"
    numerator = normalized
    for candidate, endings in suffixes:
        ending = next((item for item in endings if normalized.endswith(item)), None)
        if ending is not None:
            numerator = normalized[: -len(ending)]
            temperature = candidate
            break

    numerator = numerator.rstrip("/")
    for token in (
        "cell^-1",
        "cell-1",
        "/cell",
        "formulaunit^-1",
        "formulaunit-1",
        "/formulaunit",
        "f.u.^-1",
        "f.u.-1",
        "/f.u.",
    ):
        numerator = numerator.replace(token, "")
    numerator = numerator.rstrip("/")

    aliases = {
        "j": "J",
        "joule": "J",
        "j/mol": "J/mol",
        "jmol^-1": "J/mol",
        "jmol-1": "J/mol",
        "jmol": "J/mol",
        "kj/mol": "kJ/mol",
        "kjmol^-1": "kJ/mol",
        "kjmol-1": "kJ/mol",
        "kjmol": "kJ/mol",
        "ha": "Ha",
        "hartree": "Ha",
        "mha": "mHa",
        "millihartree": "mHa",
        "ev": "eV",
        "electronvolt": "eV",
        "mev": "meV",
        "millielectronvolt": "meV",
        "ry": "Ry",
        "rydberg": "Ry",
        "mry": "mRy",
        "millirydberg": "mRy",
    }
    if numerator not in aliases:
        raise NotImplementedError(
            f"{scale} scale is unsupported for energy-per-temperature conversion."
        )
    return aliases[numerator], temperature


def _temperature_interval_to_kelvin(scale: str) -> float:
    """
    Return the kelvin width represented by one temperature-unit interval.

    Parameters
    ----------
    scale : str
        Temperature scale.

    Returns
    -------
    float
        Width of one interval expressed in kelvin.

    Raises
    ------
    NotImplementedError
        If the scale is unsupported.
    """
    key = _normalize_scale(scale).replace("°", "")
    if key in {"k", "kelvin", "c", "celsius", "degc", "degreecelsius"}:
        return 1.0
    if key in {
        "f",
        "fahrenheit",
        "degf",
        "degreefahrenheit",
        "r",
        "rankine",
        "degr",
        "degreerankine",
    }:
        return 5.0 / 9.0
    raise NotImplementedError(
        f"{scale} scale is unsupported for temperature-interval conversion."
    )


def convert_pressure(
    value: Any,
    from_scale: str,
    to_scale: str,
) -> float | np.ndarray:
    """
    Convert pressures between Pascal, bar, and their multiples.

    Parameters
    ----------
    value : scalar or array_like
        Pressure value or values to convert.
    from_scale : str
        Original pressure scale.
    to_scale : str
        Target pressure scale.

    Returns
    -------
    float or ndarray
        Converted pressure value or values.

    Examples
    --------
    >>> convert_pressure([20.0, 25.0], "GPa", "Mbar")
    array([0.2 , 0.25])
    """
    return _linear_conversion(value, from_scale, to_scale, _PRESSURE_TO_PA, "Pa")


def energy_to_pressure(
    value: Any,
    energy_scale: str,
    volume_scale: str,
    pressure_scale: str,
) -> float | np.ndarray:
    """
    Convert an energy density to a pressure.

    The input is interpreted as an energy divided by a volume, where the energy
    and volume units are specified separately. Volume units are specified by
    the corresponding length scale, for example ``angstrom`` means
    Angstrom^3.

    Parameters
    ----------
    value : scalar or array_like
        Energy density value or values to convert.
    energy_scale : str
        Energy unit of the numerator.
    volume_scale : str
        Length unit used to define the volume unit.
    pressure_scale : str
        Target pressure unit.

    Returns
    -------
    float or ndarray
        Converted pressure value or values.

    Examples
    --------
    >>> energy_to_pressure([0.001, 0.2], "eV", "angstrom", "GPa")
    array([  1.60217663, 320.435326...])
    """
    energy_key = _normalize_scale(energy_scale)
    volume_key = _normalize_scale(volume_scale)
    pressure_key = _normalize_scale(pressure_scale)

    if energy_key not in _ENERGY_TO_J:
        raise NotImplementedError(f"{energy_scale} scale is unsupported.")
    if volume_key not in _LENGTH_TO_M:
        raise NotImplementedError(f"{volume_scale} scale is unsupported.")
    if pressure_key not in _PRESSURE_TO_PA:
        raise NotImplementedError(f"{pressure_scale} scale is unsupported.")

    result = (
        _asarray(value)
        * _ENERGY_TO_J[energy_key]
        / (_LENGTH_TO_M[volume_key] ** 3)
        / _PRESSURE_TO_PA[pressure_key]
    )
    return _return(result, value)


def pressure_to_energy(
    value: Any,
    energy_scale: str,
    volume_scale: str,
    pressure_scale: str,
) -> float | np.ndarray:
    """
    Convert a pressure to an energy density.

    The output is expressed as an energy divided by a volume, where the energy
    and volume units are specified separately. Volume units are specified by
    the corresponding length scale, for example ``angstrom`` means
    Angstrom^3.

    Parameters
    ----------
    value : scalar or array_like
        Pressure value or values to convert.
    energy_scale : str
        Target energy unit of the numerator.
    volume_scale : str
        Length unit used to define the volume unit.
    pressure_scale : str
        Original pressure unit.

    Returns
    -------
    float or ndarray
        Converted energy density value or values.

    Examples
    --------
    >>> pressure_to_energy([1.0, 2.0], "Ha", "angstrom", "GPa")
    array([0.00022937, 0.00045875])
    """
    energy_key = _normalize_scale(energy_scale)
    volume_key = _normalize_scale(volume_scale)
    pressure_key = _normalize_scale(pressure_scale)

    if energy_key not in _ENERGY_TO_J:
        raise NotImplementedError(f"{energy_scale} scale is unsupported.")
    if volume_key not in _LENGTH_TO_M:
        raise NotImplementedError(f"{volume_scale} scale is unsupported.")
    if pressure_key not in _PRESSURE_TO_PA:
        raise NotImplementedError(f"{pressure_scale} scale is unsupported.")

    result = (
        _asarray(value)
        * _PRESSURE_TO_PA[pressure_key]
        * (_LENGTH_TO_M[volume_key] ** 3)
        / _ENERGY_TO_J[energy_key]
    )
    return _return(result, value)


def convert_frequency(
    value: Any,
    from_scale: str,
    to_scale: str,
) -> float | np.ndarray:
    """
    Convert frequencies between wavenumbers, Hertz, and Terahertz.

    Parameters
    ----------
    value : scalar or array_like
        Frequency value or values to convert.
    from_scale : str
        Original frequency scale.
    to_scale : str
        Target frequency scale.

    Returns
    -------
    float or ndarray
        Converted frequency value or values.

    Examples
    --------
    >>> convert_frequency(1500.0, "cm^-1", "THz")
    44.9688687...
    """
    return _linear_conversion(value, from_scale, to_scale, _FREQUENCY_TO_CM, "cm^-1")


def convert_length(
    value: Any,
    from_scale: str,
    to_scale: str,
) -> float | np.ndarray:
    """
    Convert lengths between supported metric units, Angstrom, and Bohr.

    Parameters
    ----------
    value : scalar or array_like
        Length value or values to convert.
    from_scale : str
        Original length scale.
    to_scale : str
        Target length scale.

    Returns
    -------
    float or ndarray
        Converted length value or values.

    Examples
    --------
    >>> convert_length(20.0, "Bohr", "m")
    1.0583544218e-09
    """
    return _linear_conversion(value, from_scale, to_scale, _LENGTH_TO_M, "m")


def convert_volume(
    value: Any,
    from_scale: str,
    to_scale: str,
) -> float | np.ndarray:
    """
    Convert volumes between cubes of supported metric units, Angstrom, and Bohr.

    The ``from_scale`` and ``to_scale`` parameters refer to the length unit
    used to define the volume unit. For example, ``angstrom`` means
    Angstrom^3.

    Parameters
    ----------
    value : scalar or array_like
        Volume value or values to convert.
    from_scale : str
        Original length scale defining the volume unit.
    to_scale : str
        Target length scale defining the volume unit.

    Returns
    -------
    float or ndarray
        Converted volume value or values.

    Examples
    --------
    >>> convert_volume(10.9897, "Angstrom", "bohr")
    74.165497...
    """
    from_key = _normalize_scale(from_scale)
    to_key = _normalize_scale(to_scale)

    if from_key not in _LENGTH_TO_M:
        raise NotImplementedError(f"{from_scale} scale is unsupported.")
    if to_key not in _LENGTH_TO_M:
        raise NotImplementedError(f"{to_scale} scale is unsupported.")

    result = (
        _asarray(value) * (_LENGTH_TO_M[from_key] ** 3) / (_LENGTH_TO_M[to_key] ** 3)
    )
    return _return(result, value)
