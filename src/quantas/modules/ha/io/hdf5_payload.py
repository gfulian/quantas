# -*- coding: utf-8 -*-

"""HDF5 payload readers and writers for harmonic workflow results.

This module owns only the HA-specific scientific payload under ``/results``.
The shared Quantas HDF5 envelope is provided by :mod:`quantas.io.hdf5`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np

from quantas.io.hdf5 import (
    decode_scalar,
    read_array_dataset,
    read_group_mapping,
    write_array_dataset,
    write_mapping,
)
from quantas.models.kieffer import KiefferThermodynamicContribution
from quantas.modules.ha.models import HAResult


HA_DATASETS = {
    "volume": ("V", "Unit cell volume"),
    "temperature": ("T", "Temperature"),
    "static_energy": ("U0", "Static (electronic) energies"),
    "zero_point_energy": ("Uzp", "Zero-point energy"),
    "thermal_energy": ("Uth", "Thermal energy"),
    "internal_energy": ("Utot", "Internal energy"),
    "entropy": ("S", "Entropy"),
    "vibrational_free_energy": ("Fvib", "Vibrational Helmholtz free energy"),
    "free_energy": ("F", "Helmholtz free energy"),
    "isochoric_heat_capacity": ("Cv", "Isochoric heat capacity"),
}


def write_ha_payload(h5: h5py.File, result: HAResult) -> h5py.Group:
    """Write the HA-specific result payload.

    Parameters
    ----------
    h5 : h5py.File
        Open destination file.
    result : HAResult
        Harmonic-approximation result payload.

    Returns
    -------
    h5py.Group
        Created ``results`` group.
    """
    group = h5.create_group("results")
    group.attrs["jobname"] = result.jobname
    metadata_group = group.create_group("metadata")
    write_mapping(metadata_group, result.metadata)
    for attribute in HA_DATASETS:
        value = getattr(result, attribute)
        if value is not None:
            write_array_dataset(group, attribute, value)
    if result.kieffer_contribution is not None:
        _write_kieffer_contribution(group, result.kieffer_contribution)
    return group


def read_current_ha_payload(group: h5py.Group) -> HAResult:
    """Read the current HA payload from the ``results`` group.

    Parameters
    ----------
    group : h5py.Group
        HDF5 group containing HA result arrays.

    Returns
    -------
    HAResult
        Reconstructed HA payload.
    """
    metadata_group = group.get("metadata") or group.get("ha_metadata")
    return HAResult(
        jobname=str(group.attrs.get("jobname", "Unknown")),
        temperature=_read_optional_array(group, "temperature"),
        volume=_read_optional_array(group, "volume"),
        static_energy=_read_optional_array(group, "static_energy"),
        zero_point_energy=_read_optional_array(group, "zero_point_energy"),
        thermal_energy=_read_optional_array(group, "thermal_energy"),
        internal_energy=_read_optional_array(group, "internal_energy"),
        entropy=_read_optional_array(group, "entropy"),
        vibrational_free_energy=_read_optional_array(group, "vibrational_free_energy"),
        free_energy=_read_optional_array(group, "free_energy"),
        isochoric_heat_capacity=_read_optional_array(group, "isochoric_heat_capacity"),
        kieffer_contribution=_read_kieffer_contribution(group),
        metadata=_read_metadata_group(metadata_group),
    )


def _write_kieffer_contribution(
    group: h5py.Group,
    contribution: KiefferThermodynamicContribution,
) -> None:
    """Write the separately traceable acoustic contribution."""
    acoustic = group.create_group("kieffer_contribution")
    for name in (
        "cutoff_frequencies_hz",
        "effective_velocities_km_s",
        "zero_point_energy",
        "thermal_energy",
        "entropy",
        "vibrational_free_energy",
        "isochoric_heat_capacity",
    ):
        write_array_dataset(acoustic, name, getattr(contribution, name))
    metadata = acoustic.create_group("metadata")
    write_mapping(metadata, contribution.metadata)


def _read_kieffer_contribution(
    group: h5py.Group,
) -> KiefferThermodynamicContribution | None:
    """Read an optional acoustic contribution from a native HA result."""
    acoustic = group.get("kieffer_contribution")
    if not isinstance(acoustic, h5py.Group):
        return None
    required = (
        "cutoff_frequencies_hz",
        "effective_velocities_km_s",
        "zero_point_energy",
        "thermal_energy",
        "entropy",
        "vibrational_free_energy",
        "isochoric_heat_capacity",
    )
    values = {name: read_array_dataset(acoustic, name) for name in required}
    return KiefferThermodynamicContribution(
        cutoff_frequencies_hz=np.asarray(values["cutoff_frequencies_hz"]),
        effective_velocities_km_s=np.asarray(values["effective_velocities_km_s"]),
        zero_point_energy=np.asarray(values["zero_point_energy"]),
        thermal_energy=np.asarray(values["thermal_energy"]),
        entropy=np.asarray(values["entropy"]),
        vibrational_free_energy=np.asarray(values["vibrational_free_energy"]),
        isochoric_heat_capacity=np.asarray(values["isochoric_heat_capacity"]),
        metadata=_read_metadata_group(acoustic.get("metadata")),
    )


def read_historical_ha_payload(h5: h5py.File, filename: Path) -> tuple[HAResult, str]:
    """Read an earlier Quantas HA layout for controlled migration.

    Parameters
    ----------
    h5 : h5py.File
        Open historical HDF5 file.
    filename : Path
        Source filename stored in migrated metadata.

    Returns
    -------
    tuple of HAResult and str
        Reconstructed payload and historical schema label.

    Raises
    ------
    ValueError
        If the file does not contain a recognizable HA layout.
    """
    if "results" in h5 and "volume" in h5["results"]:
        schema_version = "1.0"
        if "metadata" in h5:
            schema_version = str(
                decode_scalar(h5["metadata"].attrs.get("schema_version", "1.0"))
            )
        return read_current_ha_payload(h5["results"]), schema_version

    required = {"V", "T", "U0"}
    if not required.issubset(set(h5.keys())):
        raise ValueError("This does not appear to be a Quantas HA HDF5 file")

    metadata: dict[str, Any] = {
        "source": str(filename),
        "historical_layout": True,
        "info": decode_scalar(h5.attrs.get("info", "")),
        "units": {
            "volume": _dataset_unit(h5, "V"),
            "temperature": _dataset_unit(h5, "T"),
            "energy": _dataset_unit(h5, "U0"),
            "entropy": _dataset_unit(h5, "S"),
            "heat_capacity": _dataset_unit(h5, "Cv"),
        },
    }
    result = HAResult(
        jobname="Unknown",
        temperature=h5["T"][:],
        volume=h5["V"][:],
        static_energy=h5["U0"][:],
        zero_point_energy=_read_optional_array(h5, "Uzp"),
        thermal_energy=_read_optional_array(h5, "Uth"),
        internal_energy=_read_optional_array(h5, "Utot"),
        entropy=_read_optional_array(h5, "S"),
        vibrational_free_energy=_read_optional_array(h5, "Fvib"),
        free_energy=_read_optional_array(h5, "F"),
        isochoric_heat_capacity=_read_optional_array(h5, "Cv"),
        metadata=metadata,
    )
    return result, "historical"


def migrate_schema_1_0_thermodynamic_units(
    result: HAResult,
    schema_version: str,
) -> bool:
    """Convert schema 1.0 HA thermal derivatives to native energy units.

    Parameters
    ----------
    result : HAResult
        HA result reconstructed from the modern HDF5 layout.
    schema_version : str
        HDF5 schema version.

    Returns
    -------
    bool
        ``True`` when a migration was applied.
    """
    metadata = result.metadata if isinstance(result.metadata, dict) else {}
    if metadata.get("thermodynamic_unit_convention"):
        return False
    if str(schema_version).strip() not in {"1", "1.0"}:
        return False

    converted = False
    for name in ("entropy", "isochoric_heat_capacity"):
        value = getattr(result, name, None)
        if value is None:
            continue
        setattr(result, name, np.asarray(value, dtype=np.float64) * 1.0e-3)
        converted = True
    if not converted:
        return False

    units = metadata.setdefault("units", {})
    energy_unit = str(units.get("energy", "Ha"))
    units["entropy"] = f"{energy_unit} cell^-1 K^-1"
    units["heat_capacity"] = f"{energy_unit} cell^-1 K^-1"
    metadata["thermodynamic_unit_convention"] = "native_energy_per_cell_per_kelvin"
    input_metadata = metadata.get("input", {})
    try:
        formula_units = int(input_metadata.get("formula_units", 1))
    except (AttributeError, TypeError, ValueError):
        formula_units = 1
    if formula_units <= 0:
        formula_units = 1
    try:
        natoms = int(input_metadata.get("natoms", 0))
    except (AttributeError, TypeError, ValueError):
        natoms = 0
    metadata.setdefault(
        "normalization",
        {
            "native_basis": "cell",
            "formula_units_per_cell": formula_units,
            "natoms_per_cell": natoms,
            "natoms_per_formula_unit": (
                float(natoms) / float(formula_units) if natoms > 0 else 0.0
            ),
            "molar_basis": "formula_unit",
        },
    )
    metadata["unit_migration"] = {
        "from_schema": str(schema_version),
        "scale_factor": 1.0e-3,
        "assumed_formula_units_per_cell": formula_units,
    }
    return True


def _read_optional_array(group: h5py.Group | h5py.File | None, name: str):
    """Read a HDF5 dataset if present."""
    if group is None:
        return None
    return read_array_dataset(group, name, required=False)


def _read_metadata_group(group: h5py.Group | None) -> dict[str, Any]:
    """Read simple metadata from a HDF5 group."""
    if group is None:
        return {}
    return read_group_mapping(group)


def _dataset_unit(group: h5py.Group | h5py.File, name: str) -> str | None:
    """Read the unit attribute of a HDF5 dataset."""
    if name not in group:
        return None
    return decode_scalar(group[name].attrs.get("unit"))
