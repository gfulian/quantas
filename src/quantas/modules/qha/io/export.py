# -*- coding: utf-8 -*-

"""Exporters for quasi-harmonic approximation results.

The exporters defined here write QHA results either to the native Quantas HDF5
container or to plain text tables suitable for inspection, plotting and external
post-processing. Generic result metadata and workflow records use the shared
Quantas 2.0 HDF5 schema.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, TextIO

import h5py
import numpy as np

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
from quantas.modules.qha.io.hdf5_payload import (
    write_qha_fit_diagnostics,
    write_qha_payload,
)
from quantas.modules.qha.models import QHAResult
from quantas.modules.qha.io.export_delimited import _write_delimited_table
from quantas.modules.qha.io.export_human import _write_human_table
from quantas.modules.qha.io.export_schema import (
    EXPORT_PROPERTY_ORDER,
    STRUCTURAL_PROPERTY_ALIASES,
    _StructuralExportData,
)
from quantas.modules.qha.report import QHA_PROPERTY_LABELS, resolve_property_name
from quantas.modules.qha.formatting import (
    QHATableFormat,
)











@dataclass(slots=True)


class QHAHDF5Export(BasicHDF5Export):
    """Export QHA results to the native Quantas HDF5 schema."""

    def export(
        self,
        result: ResultData,
        filename: str | Path,
        report_text: str | None = None,
    ) -> None:
        """Export one QHA result and its complete workflow context.

        Parameters
        ----------
        result : ResultData
            Generic Quantas result containing a :class:`QHAResult`.
        filename : str or Path
            Destination HDF5 file.
        report_text : str or None, optional
            Optional frontend-generated report text.

        Raises
        ------
        TypeError
            If the stored module payload is not a :class:`QHAResult`.
        """
        filename = ensure_suffix(filename, ".hdf5")
        qha_result = _extract_qha_result(result)

        with h5py.File(filename, "w") as h5:
            write_result_metadata(h5, result)
            write_precision_metadata(h5)
            write_input_data(h5, result.input_data)
            write_options(h5, result.options)
            write_qha_payload(h5, qha_result)
            diagnostics = write_diagnostics(h5, result, report_text=report_text)
            write_qha_fit_diagnostics(diagnostics, qha_result)
            write_events(h5, result.events)
        self.completed = True


class QHATableExport(BasicExport):
    """Text exporter for QHA pressure-temperature result tables."""

    def export(
        self,
        result: ResultData,
        filename: str | Path,
        property_name: str | None = None,
        include_uncertainty: bool = True,
        delimiter: str | None = None,
        table_format: QHATableFormat | None = None,
    ) -> None:
        """Export QHA results to a long pressure-temperature table.

        Parameters
        ----------
        result : ResultData
            Generic Quantas result containing QHA results.
        filename : str or Path
            Output text file.
        property_name : str or None, optional
            Historical property key or QHAResult attribute name to export. If
            ``None`` or ``"all"``, all available pressure-temperature
            properties are exported. The aliases ``"structure"``, ``"cell"``,
            and ``"lattice"`` export only equilibrium cell parameters, axial
            thermal expansion, and the Cartesian thermal-expansion tensor.
        include_uncertainty : bool, optional
            If ``True``, export matching uncertainty columns when available.
        delimiter : str or None, optional
            Column delimiter. If omitted, ``.dat`` and ``.txt`` files use a
            fixed-width layout with a two-line header, while ``.csv`` files
            use commas.
        table_format : QHATableFormat or None, optional
            Formatting rules for numerical values.

        Raises
        ------
        KeyError
            If the requested property is unknown.
        ValueError
            If no exportable property is available.
        """
        qha_result = _extract_qha_result(result)
        filename = Path(filename)
        if filename.suffix.lower() not in {".dat", ".txt", ".csv"}:
            filename = ensure_suffix(filename, ".dat")
        if delimiter is None and filename.suffix.lower() == ".csv":
            delimiter = ","
        fmt = table_format or QHATableFormat()

        structural_only = _is_structural_property_request(property_name)
        include_structural = (
            structural_only
            or property_name is None
            or (property_name is not None and property_name.lower() == "all")
        )
        properties = (
            [] if structural_only else _export_property_specs(qha_result, property_name)
        )
        structural = (
            _prepare_structural_export_data(
                qha_result,
                include_uncertainty=include_uncertainty,
            )
            if include_structural
            else None
        )
        if not properties and structural is None:
            raise ValueError("no QHA pressure-temperature properties are available")

        with filename.open("w", encoding="utf8", newline="") as stream:
            _write_long_property_table(
                stream=stream,
                result=qha_result,
                properties=properties,
                structural=structural,
                include_uncertainty=include_uncertainty,
                delimiter=delimiter,
                table_format=fmt,
            )
        self.completed = True


def _export_property_specs(
    result: QHAResult,
    property_name: str | None,
) -> list[tuple[str, str, str, np.ndarray]]:
    """Return exportable QHA property specifications.

    Parameters
    ----------
    result : QHAResult
        QHA result object.
    property_name : str or None
        Selected property or ``None``/``"all"`` for all properties.

    Returns
    -------
    list of tuple
        Tuples containing attribute name, symbol key, description and values.
    """
    if property_name is not None and property_name.lower() != "all":
        attr = resolve_property_name(property_name)
        value = getattr(result, attr)
        if value is None:
            raise ValueError(f"QHA property '{property_name}' is not available")
        symbol, description = QHA_PROPERTY_LABELS[attr]
        return [(attr, symbol, description, np.asarray(value, dtype=np.float64))]

    specs: list[tuple[str, str, str, np.ndarray]] = []
    for attr in EXPORT_PROPERTY_ORDER:
        symbol, description = QHA_PROPERTY_LABELS[attr]
        value = getattr(result, attr, None)
        if value is None:
            continue
        array = np.asarray(value, dtype=np.float64)
        if _is_pressure_temperature_property(result, array):
            specs.append((attr, symbol, description, array))
    return specs


def _is_structural_property_request(property_name: str | None) -> bool:
    """Return whether a table request selects structural QHA properties."""
    if property_name is None:
        return False
    return property_name.strip().lower() in STRUCTURAL_PROPERTY_ALIASES


def _prepare_structural_export_data(
    result: QHAResult,
    *,
    include_uncertainty: bool,
) -> _StructuralExportData | None:
    """Validate and collect structural pressure-temperature arrays.

    Parameters
    ----------
    result : QHAResult
        QHA result containing structural-path properties.
    include_uncertainty : bool
        Include matching one-standard-deviation arrays when available.

    Returns
    -------
    _StructuralExportData or None
        Prepared structural data, or ``None`` when no lattice parameters are
        stored.

    Raises
    ------
    ValueError
        If a stored structural array is incompatible with the pressure-
        temperature grid.
    """
    if result.lattice_parameters is None:
        return None
    temperature = _required_grid(result.temperature, "temperature")
    pressure = _required_grid(result.pressure, "pressure")
    shape = (temperature.size, pressure.size)

    parameters = _required_structural_array(
        result.lattice_parameters,
        shape + (6,),
        "lattice_parameters",
    )
    sigma_parameters = None
    if include_uncertainty:
        value = result.uncertainties.get("sigma_cell_parameters")
        if value is not None:
            sigma_parameters = _required_structural_array(
                value,
                shape + (6,),
                "sigma_cell_parameters",
            )

    axial = _optional_structural_array(
        result.axial_thermal_expansion,
        shape + (3,),
        "axial_thermal_expansion",
    )
    sigma_axial = None
    if include_uncertainty:
        for key in (
            "sigma_axial_thermal_expansion",
            "sigma_alphaABC",
        ):
            value = result.uncertainties.get(key)
            if value is not None:
                sigma_axial = _required_structural_array(
                    value,
                    shape + (3,),
                    key,
                )
                break

    tensor = _optional_structural_array(
        result.thermal_expansion_tensor,
        shape + (3, 3),
        "thermal_expansion_tensor",
    )
    alpha_v = None
    if result.thermal_expansion is not None:
        alpha_v = _as_pressure_temperature_array(
            result,
            np.asarray(result.thermal_expansion, dtype=np.float64),
            "alphaV",
        )

    trace = None
    residual = None
    if tensor is not None:
        trace = np.trace(tensor, axis1=-2, axis2=-1)
        if alpha_v is not None:
            residual = trace - alpha_v

    extrapolation = None
    if result.structural_extrapolation_mask is not None:
        extrapolation = np.asarray(result.structural_extrapolation_mask, dtype=bool)
        if extrapolation.shape != shape:
            raise ValueError(
                "QHA structural array 'structural_extrapolation_mask' must "
                f"have shape {shape}"
            )

    return _StructuralExportData(
        lattice_parameters=parameters,
        lattice_parameter_uncertainties=sigma_parameters,
        axial_thermal_expansion=axial,
        axial_thermal_expansion_uncertainties=sigma_axial,
        thermal_expansion_tensor=tensor,
        volumetric_thermal_expansion=alpha_v,
        tensor_trace=trace,
        trace_residual=residual,
        extrapolation_mask=extrapolation,
    )


def _required_structural_array(
    value: Any,
    expected_shape: tuple[int, ...],
    name: str,
) -> np.ndarray:
    """Return one structural array after strict shape validation."""
    array = np.asarray(value, dtype=np.float64)
    if array.shape != expected_shape:
        raise ValueError(
            f"QHA structural array '{name}' must have shape {expected_shape}"
        )
    return array


def _optional_structural_array(
    value: Any | None,
    expected_shape: tuple[int, ...],
    name: str,
) -> np.ndarray | None:
    """Return a validated optional structural array."""
    if value is None:
        return None
    return _required_structural_array(value, expected_shape, name)


def _write_long_property_table(
    stream: TextIO,
    result: QHAResult,
    properties: list[tuple[str, str, str, np.ndarray]],
    structural: _StructuralExportData | None,
    include_uncertainty: bool,
    delimiter: str | None,
    table_format: QHATableFormat,
) -> None:
    """Write a pressure-temperature table.

    Parameters
    ----------
    stream : file-like object
        Open text stream.
    result : QHAResult
        QHA result object.
    properties : list of tuple
        Property specifications returned by :func:`_export_property_specs`.
    structural : _StructuralExportData or None
        Optional structural pressure-temperature arrays.
    include_uncertainty : bool
        If ``True``, matching uncertainty columns are written.
    delimiter : str or None
        Column delimiter. If ``None``, a compact fixed-width human-readable
        layout is used. Delimited output is written as a flat machine-readable
        table without decorative metadata rows.
    table_format : QHATableFormat
        Numerical formatting rules.
    """
    temperature = _required_grid(result.temperature, "temperature")
    pressure = _required_grid(result.pressure, "pressure")
    units = result.metadata.get("units", {})
    temperature_unit = str(units.get("temperature", "K"))
    pressure_unit = str(units.get("pressure", "unknown"))
    thermal_expansion_source = _thermal_expansion_source_names(result)

    arrays: list[tuple[str, str, str, np.ndarray, np.ndarray | None]] = []
    for attr, symbol, description, values in properties:
        array = _as_pressure_temperature_array(result, values, symbol)
        sigma = None
        if include_uncertainty:
            sigma_value = _uncertainty_for(result, attr, symbol)
            if sigma_value is not None:
                sigma = _as_pressure_temperature_array(
                    result, np.asarray(sigma_value, dtype=np.float64), f"sigma_{symbol}"
                )
        arrays.append((attr, symbol, description, array, sigma))

    if delimiter is None:
        _write_human_table(
            stream,
            result,
            arrays,
            structural,
            temperature,
            pressure,
            temperature_unit,
            pressure_unit,
            table_format,
            thermal_expansion_source,
        )
        return

    _write_delimited_table(
        stream,
        arrays,
        structural,
        temperature,
        pressure,
        temperature_unit,
        pressure_unit,
        units,
        delimiter,
        table_format,
        thermal_expansion_source,
    )


























def _is_pressure_temperature_property(result: QHAResult, array: np.ndarray) -> bool:
    """Return whether an array can be exported on the pressure-temperature grid."""
    if result.temperature is None or result.pressure is None:
        return False
    temperature = np.asarray(result.temperature)
    pressure = np.asarray(result.pressure)
    if array.ndim == 2 and array.shape == (temperature.size, pressure.size):
        return True
    if array.ndim == 1 and array.shape[0] == temperature.size:
        return True
    return False


def _as_pressure_temperature_array(
    result: QHAResult,
    values: np.ndarray,
    property_key: str,
) -> np.ndarray:
    """Return values as an array with shape ``(nT, nP)``.

    Parameters
    ----------
    result : QHAResult
        QHA result object.
    values : ndarray
        Property values.
    property_key : str
        Name used in error messages.

    Returns
    -------
    ndarray
        Pressure-temperature array.

    Raises
    ------
    ValueError
        If the values are incompatible with the pressure-temperature grid.
    """
    temperature = _required_grid(result.temperature, "temperature")
    pressure = _required_grid(result.pressure, "pressure")
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 1 and array.shape[0] == temperature.shape[0]:
        return np.repeat(array[:, np.newaxis], pressure.shape[0], axis=1)
    if array.ndim == 2 and array.shape == (temperature.shape[0], pressure.shape[0]):
        return array
    raise ValueError(f"QHA property '{property_key}' must have shape (T,) or (T, P)")






def _extract_qha_result(result: ResultData) -> QHAResult:
    """Return the QHA result stored in a generic result container.

    Parameters
    ----------
    result : ResultData
        Generic Quantas result.

    Returns
    -------
    QHAResult
        Stored QHA result.

    Raises
    ------
    KeyError
        If the result does not contain a ``"qha"`` entry.
    TypeError
        If the entry is not a :class:`QHAResult` instance.
    """
    qha_result = result.results["qha"]
    if not isinstance(qha_result, QHAResult):
        raise TypeError("result.results['qha'] must contain a QHAResult object")
    return qha_result


def _write_property_table(
    stream: TextIO,
    result: QHAResult,
    values: np.ndarray,
    sigma: np.ndarray | None,
    property_key: str,
    property_description: str,
    property_unit: str,
    pressure_unit: str,
    temperature_unit: str,
) -> None:
    """Write a rectangular pressure-temperature QHA data table.

    Parameters
    ----------
    stream : file-like object
        Open text stream.
    result : QHAResult
        QHA result object.
    values : ndarray
        Property values.
    sigma : ndarray or None
        Optional uncertainty array with the same shape as ``values``.
    property_key : str
        Historical property key.
    property_description : str
        Human-readable property description.
    property_unit : str
        Property unit label.
    pressure_unit : str
        Pressure unit label.
    temperature_unit : str
        Temperature unit label.

    Raises
    ------
    ValueError
        If the array shape cannot be represented as a table.
    """
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 1:
        temperature = _required_grid(result.temperature, "temperature")
        if array.shape[0] != temperature.shape[0]:
            raise ValueError("one-dimensional QHA arrays must match temperature")
        pressure = np.asarray([0.0], dtype=np.float64)
        array = array[:, np.newaxis]
    elif array.ndim == 2:
        pressure = _required_grid(result.pressure, "pressure")
        temperature = _required_grid(result.temperature, "temperature")
        if array.shape != (temperature.shape[0], pressure.shape[0]):
            raise ValueError("two-dimensional QHA arrays must have shape (T, P)")
    else:
        raise ValueError("QHA table export requires a one- or two-dimensional array")

    sigma_array: np.ndarray | None = None
    if sigma is not None:
        sigma_array = np.asarray(sigma, dtype=np.float64)
        if sigma_array.ndim == 1:
            sigma_array = sigma_array[:, np.newaxis]
        if sigma_array.shape != array.shape:
            raise ValueError("uncertainty array shape is incompatible with values")

    stream.write("# Quantas QHA table export\n")
    stream.write(f"# Property: {property_key} - {property_description}\n")
    stream.write(f"# Units: {property_unit}\n")
    stream.write("#\n")
    headers = [
        f"T / {temperature_unit}",
        f"P / {pressure_unit}",
        f"{property_key} / {property_unit}",
    ]
    if sigma_array is not None:
        headers.append(f"sigma_{property_key} / {property_unit}")
    stream.write("\t".join(headers))
    stream.write("\n")

    for it, temperature_value in enumerate(temperature):
        for ip, pressure_value in enumerate(pressure):
            row = [
                f"{float(temperature_value):.8f}",
                f"{float(pressure_value):.8f}",
                f"{float(array[it, ip]):.12E}",
            ]
            if sigma_array is not None:
                row.append(f"{float(sigma_array[it, ip]):.12E}")
            stream.write("\t".join(row))
            stream.write("\n")




def _thermal_expansion_source_names(result: QHAResult) -> np.ndarray | None:
    """Return per-point thermal-expansion source names."""
    if result.thermal_expansion_source is None:
        return None
    source = np.asarray(result.thermal_expansion_source)
    metadata = result.metadata.get("thermal_expansion", {})
    codes = metadata.get("source_codes", {}) if isinstance(metadata, dict) else {}
    if not isinstance(codes, dict) or not codes:
        codes = {
            "invalid": 0,
            "mixed_derivative": 1,
            "mode_gruneisen": 2,
            "numerical": 3,
            "numerical_fallback": 4,
        }
    names = np.full(source.shape, "unknown", dtype=object)
    for name, code in codes.items():
        names[source == int(code)] = str(name)
    return names


def _required_grid(value: np.ndarray | None, name: str) -> np.ndarray:
    """Return a one-dimensional grid array.

    Parameters
    ----------
    value : ndarray or None
        Grid value.
    name : str
        Name used in error messages.

    Returns
    -------
    ndarray
        One-dimensional grid.

    Raises
    ------
    ValueError
        If the grid is missing or not one-dimensional.
    """
    if value is None:
        raise ValueError(f"QHA result is missing the {name} grid")
    array = np.asarray(value, dtype=np.float64)
    if array.ndim != 1:
        raise ValueError(f"{name} grid must be one-dimensional")
    return array


def _uncertainty_for(result: QHAResult, attr: str, symbol: str) -> np.ndarray | None:
    """Return the uncertainty array associated with a QHA property.

    Parameters
    ----------
    result : QHAResult
        QHA result object.
    attr : str
        QHAResult attribute name.
    symbol : str
        Historical property key.

    Returns
    -------
    ndarray or None
        Matching uncertainty array if available.
    """
    candidates = (f"sigma_{attr}", f"sigma_{symbol}", attr, symbol)
    for key in candidates:
        value = result.uncertainties.get(key)
        if value is not None:
            return np.asarray(value)
    return None




def _is_scalar(value: Any) -> bool:
    """Return whether a value can be stored in a simple HDF5 dataset.

    Parameters
    ----------
    value : Any
        Value to inspect.

    Returns
    -------
    bool
        ``True`` for scalar-like values.
    """
    return isinstance(value, (str, bytes, bool, int, float, np.number))
