# -*- coding: utf-8 -*-

"""Measurement-unit resolution for Quantas phonon input contracts.

HA and QHA input files are self-describing.  This module resolves the stored
measurement units and optional caller overrides into the compact unit labels
used by the numerical workflows.  Pressure and temperature are intentionally
outside this contract because they describe calculation/output domains rather
than measurements stored in the phonon input file.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class PhononMeasurementUnits:
    """Resolved measurement units carried by one phonon input dataset.

    Parameters
    ----------
    energy : str
        Energy unit used to interpret stored static energies.
    length : str
        Length unit whose cube defines stored volumes.
    frequency : str
        Unit used to interpret stored ordinary phonon frequencies.
    volume : str
        Human-readable cubic unit consistent with ``length``.
    """

    energy: str
    length: str
    frequency: str
    volume: str


_ENERGY_ALIASES = {
    "ha": "Ha",
    "hartree": "Ha",
    "hartrees": "Ha",
    "ev": "eV",
    "electronvolt": "eV",
    "electronvolts": "eV",
    "electron volt": "eV",
    "electron volts": "eV",
    "ry": "Ry",
    "rydberg": "Ry",
    "rydbergs": "Ry",
}
_LENGTH_ALIASES = {
    "a": "A",
    "ang": "A",
    "angstrom": "A",
    "angstroms": "A",
    "å": "A",
    "bohr": "bohr",
    "bohrs": "bohr",
    "bohr radius": "bohr",
}
_FREQUENCY_ALIASES = {
    "cm-1": "cm^-1",
    "cm^-1": "cm^-1",
    "wavenumber": "cm^-1",
    "wavenumbers": "cm^-1",
    "thz": "THz",
    "terahertz": "THz",
    "hz": "Hz",
    "hertz": "Hz",
}
_VOLUME_ALIASES = {
    "a^3": "A",
    "a3": "A",
    "angstrom^3": "A",
    "angstrom3": "A",
    "angstroms^3": "A",
    "angstroms3": "A",
    "å^3": "A",
    "å3": "A",
    "bohr^3": "bohr",
    "bohr3": "bohr",
    "bohrs^3": "bohr",
    "bohrs3": "bohr",
}
_VOLUME_LABELS = {"A": "angstrom^3", "bohr": "bohr^3"}


def _normalize(label: object) -> str:
    """Return one case-insensitive normalized unit token."""
    return str(label).strip().lower().replace("³", "^3")


def _canonical(label: object, aliases: Mapping[str, str], quantity: str) -> str:
    """Return one supported canonical unit label.

    Parameters
    ----------
    label : object
        Unit label to normalize.
    aliases : mapping
        Supported normalized aliases and their canonical values.
    quantity : str
        Quantity name used in validation messages.

    Returns
    -------
    str
        Canonical unit label.

    Raises
    ------
    ValueError
        If the requested unit is not supported by the phonon contract.
    """
    key = _normalize(label)
    try:
        return aliases[key]
    except KeyError as exc:
        raise ValueError(f"unsupported phonon {quantity} unit: {label}") from exc


def resolve_phonon_measurement_units(
    units: Mapping[str, str],
    *,
    energy_unit: str | None = None,
    length_unit: str | None = None,
    frequency_unit: str | None = None,
) -> PhononMeasurementUnits:
    """Resolve phonon-input measurement units and optional explicit overrides.

    Stored input metadata are authoritative unless the caller explicitly
    supplies an override.  Historical inputs should already carry the legacy
    defaults injected by :class:`~quantas.io.phonons.PhononInputFileReader`.
    The stored ``length`` and ``volume`` declarations are required to describe
    the same length basis.

    Parameters
    ----------
    units : mapping
        Input unit metadata containing ``energy``, ``length``, ``volume``, and
        ``frequency`` entries.
    energy_unit : str or None, optional
        Explicit energy-unit override.
    length_unit : str or None, optional
        Explicit length-unit override.  Stored volumes are interpreted in the
        cube of this unit.
    frequency_unit : str or None, optional
        Explicit ordinary-frequency-unit override.

    Returns
    -------
    PhononMeasurementUnits
        Canonical units ready for HA/QHA options.

    Raises
    ------
    ValueError
        If required metadata are missing, unsupported, or internally
        inconsistent.
    """
    required = ("energy", "length", "volume", "frequency")
    missing = [key for key in required if not str(units.get(key, "")).strip()]
    if missing:
        raise ValueError(
            "phonon unit metadata missing: " + ", ".join(missing)
        )

    stored_length = _canonical(units["length"], _LENGTH_ALIASES, "length")
    volume_length = _canonical(units["volume"], _VOLUME_ALIASES, "volume")
    if stored_length != volume_length:
        raise ValueError(
            "phonon length and volume units are inconsistent: "
            f"{units['length']} vs {units['volume']}"
        )

    resolved_energy = _canonical(
        units["energy"] if energy_unit is None else energy_unit,
        _ENERGY_ALIASES,
        "energy",
    )
    resolved_length = _canonical(
        units["length"] if length_unit is None else length_unit,
        _LENGTH_ALIASES,
        "length",
    )
    resolved_frequency = _canonical(
        units["frequency"] if frequency_unit is None else frequency_unit,
        _FREQUENCY_ALIASES,
        "frequency",
    )
    return PhononMeasurementUnits(
        energy=resolved_energy,
        length=resolved_length,
        frequency=resolved_frequency,
        volume=_VOLUME_LABELS[resolved_length],
    )


__all__ = ["PhononMeasurementUnits", "resolve_phonon_measurement_units"]
