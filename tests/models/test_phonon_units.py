"""Tests for self-describing HA/QHA measurement-unit resolution."""

from __future__ import annotations

import pytest

from quantas.models import resolve_phonon_measurement_units


def test_phonon_measurement_units_prefer_metadata_and_accept_overrides() -> None:
    units = {
        "energy": "eV",
        "volume": "angstrom^3",
        "frequency": "cm^-1",
        "length": "angstrom",
    }

    resolved = resolve_phonon_measurement_units(units)
    overridden = resolve_phonon_measurement_units(
        units,
        energy_unit="Ry",
        length_unit="bohr",
        frequency_unit="THz",
    )

    assert resolved.energy == "eV"
    assert resolved.length == "A"
    assert resolved.volume == "angstrom^3"
    assert resolved.frequency == "cm^-1"
    assert overridden.energy == "Ry"
    assert overridden.length == "bohr"
    assert overridden.volume == "bohr^3"
    assert overridden.frequency == "THz"


def test_phonon_measurement_units_reject_inconsistent_length_and_volume() -> None:
    with pytest.raises(ValueError, match="length and volume units are inconsistent"):
        resolve_phonon_measurement_units(
            {
                "energy": "Ha",
                "volume": "bohr^3",
                "frequency": "cm^-1",
                "length": "angstrom",
            }
        )
