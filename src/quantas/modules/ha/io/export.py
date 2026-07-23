# -*- coding: utf-8 -*-

"""
Exporters for harmonic-approximation results.

The HDF5 exporter stores HA results in the native Quantas 2.0 schema.
Generic metadata, input, options, diagnostics, and workflow events are written
through the shared HDF5 infrastructure.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, TextIO

import h5py
import numpy as np

from quantas.core.physics.units import (
    convert_energy,
    convert_energy_per_temperature,
)
from quantas.io.hdf5 import (
    write_diagnostics,
    write_events,
    write_input_data,
    write_options,
    write_precision_metadata,
    write_result_metadata,
)
from quantas.io.path import ensure_suffix
from quantas.models import BasicExport, BasicHDF5Export, ResultData
from quantas.modules.ha.io.hdf5_payload import HA_DATASETS, write_ha_payload
from quantas.modules.ha.models import HAResult
from quantas.modules.ha.report import THERMODYNAMIC_LABELS


class HAHDF5Export(BasicHDF5Export):
    """Export harmonic-approximation results to native Quantas HDF5."""

    def export(
        self,
        result: ResultData,
        filename: str | Path,
        report_text: str | None = None,
    ) -> None:
        """Export one harmonic result using the shared result schema.

        Parameters
        ----------
        result : ResultData
            Generic Quantas result containing a :class:`HAResult` under
            ``result.results["ha"]``.
        filename : str or Path
            Destination HDF5 file.
        report_text : str or None, optional
            Optional frontend-generated report text.

        Raises
        ------
        TypeError
            If the result payload is not a :class:`HAResult`.
        """
        filename = ensure_suffix(filename, ".hdf5")
        ha_result = result.results["ha"]
        if not isinstance(ha_result, HAResult):
            raise TypeError("result.results['ha'] must contain a HAResult object")

        with h5py.File(filename, "w") as h5:
            write_result_metadata(h5, result)
            write_precision_metadata(h5)
            write_input_data(h5, result.input_data)
            write_options(h5, result.options)
            write_ha_payload(h5, ha_result)
            write_diagnostics(h5, result, report_text=report_text)
            write_events(h5, result.events)
        self.completed = True


class HATableExport(BasicExport):
    """
    Text exporter for selected HA thermodynamic arrays.
    """

    def export(
        self,
        result: ResultData,
        filename: str | Path,
        property_name: str = "F",
        unit: str | None = None,
    ) -> None:
        """
        Export a selected HA property to a text table.

        Parameters
        ----------
        result : ResultData
            Generic Quantas result containing HA results.
        filename : str or Path
            Output text file.
        property_name : str, optional
            Historical property key or HAResult attribute name to export.
        unit : str or None, optional
            Requested output unit. If ``None``, values are exported in the
            units stored in the HA result metadata. Energy-like quantities use
            energy units, while entropy and heat capacity use energy per kelvin
            labels such as ``"kJ mol^-1 K^-1"``.

        Raises
        ------
        KeyError
            If the requested property is unknown.
        ValueError
            If the requested property is unavailable.
        """
        ha_result = result.results["ha"]
        if not isinstance(ha_result, HAResult):
            raise TypeError("result.results['ha'] must contain a HAResult object")

        attr = _resolve_property_name(property_name)
        value = getattr(ha_result, attr)
        if value is None:
            raise ValueError(f"HA property '{property_name}' is not available")

        filename = ensure_suffix(filename, ".dat")
        units = ha_result.metadata.get("units", {})
        values, property_unit = convert_property_values(
            value=np.asarray(value),
            attr=attr,
            units=units,
            target_unit=unit,
            normalization=ha_result.metadata.get("normalization"),
        )
        with filename.open("w", encoding="utf8") as stream:
            _write_table(
                stream=stream,
                temperature=np.asarray(ha_result.temperature),
                volume=np.asarray(ha_result.volume),
                values=np.asarray(values),
                property_key=HA_DATASETS[attr][0],
                property_description=HA_DATASETS[attr][1],
                temperature_unit=str(units.get("temperature", "K")),
                volume_unit=str(units.get("volume", "unknown")),
                property_unit=property_unit,
            )

        self.completed = True


def convert_property_values(
    value: np.ndarray,
    attr: str,
    units: dict[str, Any],
    target_unit: str | None = None,
    normalization: dict[str, Any] | None = None,
) -> tuple[np.ndarray, str]:
    """
    Convert HA property values for tabular export.

    Parameters
    ----------
    value : ndarray
        Property values as stored in the HA result.
    attr : str
        HAResult attribute name.
    units : dict
        Unit metadata stored in the HA result.
    target_unit : str or None, optional
        Requested output unit. If omitted, the stored unit is retained.
    normalization : dict or None, optional
        Result normalization metadata. When converting from native per-cell
        units to molar units, ``formula_units_per_cell`` is used so that one
        mole refers to one mole of formula units.

    Returns
    -------
    tuple of ndarray and str
        Converted values and output unit label.

    Raises
    ------
    NotImplementedError
        If the requested energy unit conversion is unsupported.
    """
    source_unit = _property_unit(attr, units)
    if target_unit is None:
        return np.asarray(value, dtype=np.float64), source_unit

    target_unit = _canonical_output_unit(target_unit, attr)
    if attr in {"volume", "temperature"}:
        if target_unit != source_unit:
            raise NotImplementedError(
                "HA table export currently supports unit conversion only for "
                "energy-like, entropy, and heat-capacity quantities."
            )
        return np.asarray(value, dtype=np.float64), source_unit

    formula_units = _formula_units_per_cell(normalization)
    if attr in {"entropy", "isochoric_heat_capacity"}:
        converted = convert_energy_per_temperature(value, source_unit, target_unit)
    else:
        source_energy = _energy_unit_from_property_unit(source_unit)
        target_energy = _energy_unit_from_property_unit(target_unit)
        converted = convert_energy(value, source_energy, target_energy)

    converted = np.asarray(converted, dtype=np.float64)
    source_is_molar = _is_molar_unit(source_unit)
    target_is_molar = _is_molar_unit(target_unit)
    if source_is_molar and not target_is_molar:
        converted *= formula_units
    elif target_is_molar and not source_is_molar:
        converted /= formula_units
    return converted, target_unit


def _canonical_output_unit(unit: str, attr: str) -> str:
    """
    Return a canonical unit label for HA table export.

    Parameters
    ----------
    unit : str
        User-provided unit label.
    attr : str
        HAResult attribute name.

    Returns
    -------
    str
        Canonical unit label.
    """
    compact = _compact_unit(unit)
    per_temperature = attr in {"entropy", "isochoric_heat_capacity"}
    if compact.startswith("kjmol"):
        return "kJ mol^-1 K^-1" if per_temperature else "kJ mol^-1"
    if compact.startswith("jmol"):
        return "J mol^-1 K^-1" if per_temperature else "J mol^-1"
    if compact.startswith(("ha", "hartree")):
        return "Ha cell^-1 K^-1" if per_temperature else "Ha"
    if compact.startswith(("ev", "electronvolt")):
        return "eV cell^-1 K^-1" if per_temperature else "eV"
    if compact.startswith(("ry", "rydberg")):
        return "Ry cell^-1 K^-1" if per_temperature else "Ry"
    return unit


def _energy_unit_from_property_unit(unit: str) -> str:
    """
    Extract the energy part from an energy-per-temperature unit label.

    Parameters
    ----------
    unit : str
        Energy or energy-per-temperature unit label.

    Returns
    -------
    str
        Unit label accepted by :func:`convert_energy`.
    """
    compact = _compact_unit(unit)
    if compact.startswith("kj"):
        return "kjmol"
    if compact.startswith("j"):
        return "j/mol"
    if compact.startswith("ha") or compact.startswith("hartree"):
        return "Ha"
    if compact.startswith("ev") or compact.startswith("electronvolt"):
        return "eV"
    if compact.startswith("ry") or compact.startswith("rydberg"):
        return "Ry"
    return unit


def _compact_unit(unit: str) -> str:
    """
    Normalize a unit label for permissive command-line matching.

    Parameters
    ----------
    unit : str
        Unit label.

    Returns
    -------
    str
        Compact lowercase unit label.
    """
    return (
        unit.strip()
        .lower()
        .replace(" ", "")
        .replace("_", "")
        .replace("/", "")
        .replace("−", "-")
    )


def _formula_units_per_cell(normalization: dict[str, Any] | None) -> float:
    """Return the number of formula units per cell from result metadata."""
    if not isinstance(normalization, dict):
        return 1.0
    value = normalization.get("formula_units_per_cell", 1)
    try:
        formula_units = float(value)
    except (TypeError, ValueError):
        return 1.0
    if formula_units <= 0.0:
        return 1.0
    return formula_units


def _is_molar_unit(unit: str) -> bool:
    """Return whether a unit label is normalized per mole."""
    compact = unit.strip().lower().replace(" ", "")
    return "mol" in compact


def _write_table(
    stream: TextIO,
    temperature: np.ndarray,
    volume: np.ndarray,
    values: np.ndarray,
    property_key: str,
    property_description: str,
    temperature_unit: str,
    volume_unit: str,
    property_unit: str,
) -> None:
    """
    Write a simple rectangular HA data table.

    Parameters
    ----------
    stream : file-like object
        Open text stream.
    temperature : ndarray
        Temperature grid.
    volume : ndarray
        Volume grid.
    values : ndarray
        Property values with shape ``(nt, nv)`` or ``(nv,)``.
    property_key : str
        Historical property key.
    property_description : str
        Human-readable property description.
    temperature_unit : str
        Temperature unit label.
    volume_unit : str
        Volume unit label.
    property_unit : str
        Property unit label.

    Raises
    ------
    ValueError
        If the array dimensions are incompatible.
    """
    array = values
    if array.ndim == 1:
        array = array[np.newaxis, :]
    if array.ndim != 2:
        raise ValueError("HA table export requires a one- or two-dimensional array")
    if array.shape[1] != volume.shape[0]:
        raise ValueError("property volume dimension is incompatible with volumes")
    if array.shape[0] not in {1, temperature.shape[0]}:
        raise ValueError(
            "property temperature dimension is incompatible with temperatures"
        )

    stream.write("# Quantas HA table export\n")
    stream.write(f"# Property: {property_key} - {property_description}\n")
    stream.write(f"# Units: {property_unit}\n")
    stream.write("#\n")
    headers = [f"T / {temperature_unit}"]
    headers.extend(f"V={value:.8g} / {volume_unit}" for value in volume)
    stream.write("\t".join(headers))
    stream.write("\n")

    nrows = temperature.shape[0] if array.shape[0] != 1 else 1
    for irow in range(nrows):
        t_value = temperature[irow] if array.shape[0] != 1 else 0.0
        row = [f"{float(t_value):.8f}"]
        row.extend(f"{float(value):.12E}" for value in array[irow])
        stream.write("\t".join(row))
        stream.write("\n")


def _resolve_property_name(name: str) -> str:
    """
    Resolve a result attribute or symbol property key.

    Parameters
    ----------
    name : str
        HAResult attribute name or historical key.

    Returns
    -------
    str
        HAResult attribute name.

    Raises
    ------
    KeyError
        If the property name is unknown.
    """
    if name in HA_DATASETS:
        return name
    for attr, (symbol, _) in HA_DATASETS.items():
        if name == symbol:
            return attr
    # Also accept names defined in the report mapping.
    for attr, (symbol, _) in THERMODYNAMIC_LABELS.items():
        if name == symbol:
            return attr
    raise KeyError(f"unknown HA property: {name}")


def _property_unit(attr: str, units: dict[str, Any]) -> str:
    """
    Return the unit label for a HA result property.

    Parameters
    ----------
    attr : str
        HAResult attribute name.
    units : dict
        Unit metadata stored in HAResult.

    Returns
    -------
    str
        Unit label.
    """
    if attr == "entropy":
        return str(units.get("entropy", "unknown"))
    if attr == "isochoric_heat_capacity":
        return str(units.get("heat_capacity", "unknown"))
    if attr == "temperature":
        return str(units.get("temperature", "K"))
    if attr == "volume":
        return str(units.get("volume", "unknown"))
    return str(units.get("energy", "unknown"))
