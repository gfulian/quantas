import numpy as np

from quantas.core.physics.units import (
    convert_energy,
    convert_frequency,
    convert_length,
    convert_pressure,
    convert_temperature,
    convert_volume,
    energy_to_pressure,
    pressure_to_energy,
)


def test_temperature_celsius_to_kelvin():
    result = convert_temperature([-40.0, 40.0], "Celsius", "Kelvin")
    np.testing.assert_allclose(result, [233.15, 313.15])


def test_energy_hartree_to_ev():
    result = convert_energy(1.0, "Ha", "eV")
    np.testing.assert_allclose(result, 27.211386, rtol=1e-6)


def test_pressure_gpa_to_mbar():
    result = convert_pressure([20.0, 25.0], "GPa", "Mbar")
    np.testing.assert_allclose(result, [0.2, 0.25])


def test_length_bohr_to_angstrom():
    result = convert_length(1.0, "bohr", "angstrom")
    np.testing.assert_allclose(result, 0.529177, rtol=1e-6)


def test_volume_angstrom_to_bohr():
    result = convert_volume(1.0, "angstrom", "bohr")
    np.testing.assert_allclose(result, 6.748334, rtol=1e-6)


def test_frequency_wavenumber_to_thz():
    result = convert_frequency(1500.0, "cm^-1", "THz")
    np.testing.assert_allclose(result, 44.9688687, rtol=1e-7)


def test_energy_pressure_roundtrip():
    pressure = energy_to_pressure(0.001, "eV", "angstrom", "GPa")
    energy = pressure_to_energy(pressure, "eV", "angstrom", "GPa")

    np.testing.assert_allclose(energy, 0.001)


def test_energy_per_temperature_jmol_to_hartree_per_cell():
    from scipy import constants as cs
    from quantas.core.physics.units import convert_energy_per_temperature

    result = convert_energy_per_temperature(
        1.0,
        "J mol^-1 K^-1",
        "Ha cell^-1 K^-1",
    )
    expected = 1.0 / (cs.Avogadro * cs.physical_constants["Hartree energy"][0])
    np.testing.assert_allclose(result, expected, rtol=1.0e-12)


def test_energy_per_temperature_interval_conversion():
    from quantas.core.physics.units import convert_energy_per_temperature

    result = convert_energy_per_temperature(
        1.0,
        "J/mol/K",
        "J/mol/Fahrenheit",
    )
    np.testing.assert_allclose(result, 5.0 / 9.0)
