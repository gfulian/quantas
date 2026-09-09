"""Tests for pressure reconstructed from static energy--volume data."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.eos import (
    EnergyEOS,
    pressure_from_energy_eos,
    pressure_from_energy_polynomial,
)
from quantas.core.physics.units import energy_to_pressure


def test_polynomial_pressure_is_negative_energy_derivative() -> None:
    """A quadratic E(V) recovers its analytic pressure in canonical units."""
    volume = np.linspace(8.0, 12.0, 5)
    energy = -10.0 + 0.01 * (volume - 10.0) ** 2

    estimate = pressure_from_energy_polynomial(
        volume,
        energy,
        degree=2,
        energy_unit="Ha",
        volume_unit="angstrom",
    )

    expected = energy_to_pressure(
        -0.02 * (volume - 10.0),
        "Ha",
        "angstrom",
        "GPa",
    )
    assert estimate.success
    assert estimate.method == "polynomial"
    assert estimate.metadata == {"degree": 2}
    np.testing.assert_allclose(estimate.pressure, expected, rtol=2.0e-12, atol=1e-11)


def test_eos_pressure_retains_canonical_model_and_fit() -> None:
    """An exact BM3 dataset returns its pressure and model provenance."""
    model = EnergyEOS()
    volume = np.linspace(66.0, 78.0, 9)
    parameters = np.array([-100.0, 0.55, 4.2, 72.0])
    energy = model.evaluate("BM3", volume, parameters)

    estimate = pressure_from_energy_eos(
        volume,
        energy,
        eos="BM3",
        energy_unit="Ha",
        volume_unit="angstrom",
    )

    expected = energy_to_pressure(
        model.pressure("BM3", parameters, volume),
        "Ha",
        "angstrom",
        "GPa",
    )
    assert estimate.success
    assert estimate.metadata == {
        "eos": "BM3",
        "eos_family": "birchmurnaghan",
        "eos_order": 3,
    }
    assert estimate.fit.parameter_names
    np.testing.assert_allclose(estimate.pressure, expected, rtol=1.0e-8, atol=1.0e-7)


def test_energy_pressure_requires_positive_volume() -> None:
    """Nonphysical volumes fail before any fit or unit conversion."""
    with pytest.raises(ValueError, match="positive volumes"):
        pressure_from_energy_polynomial(
            [0.0, 1.0, 2.0],
            [1.0, 0.0, 1.0],
            degree=2,
            energy_unit="Ha",
            volume_unit="angstrom",
        )


def test_tait_energy_pressure_reconstruction_uses_matching_integrated_form() -> None:
    """The Tait E(V) fit reconstructs its analytical pressure counterpart."""
    model = EnergyEOS()
    volume = np.linspace(66.0, 78.0, 11)
    parameters = np.array([-100.0, 0.55, 4.2, 72.0])
    energy = model.evaluate("T3", volume, parameters)

    estimate = pressure_from_energy_eos(
        volume,
        energy,
        eos="T3",
        energy_unit="Ha",
        volume_unit="angstrom",
    )

    expected = energy_to_pressure(
        model.pressure("T3", parameters, volume),
        "Ha",
        "angstrom",
        "GPa",
    )
    assert estimate.success
    assert estimate.metadata == {
        "eos": "T3",
        "eos_family": "tait",
        "eos_order": 3,
    }
    np.testing.assert_allclose(estimate.pressure, expected, rtol=2.0e-8, atol=2.0e-7)
