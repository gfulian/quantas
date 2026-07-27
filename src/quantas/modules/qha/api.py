# -*- coding: utf-8 -*-

"""Internal QHA workflow facade behind :mod:`quantas.api.qha`.

The functions read normalized quasi-harmonic inputs, execute the workflow, and
build reports and persistent results without creating frontend objects.
"""

from __future__ import annotations

from pathlib import Path

from quantas.io.path import ensure_suffix
from quantas.core.events import Observer
from quantas.models import (
    ModuleContract,
    PlotCollection,
    PlotInventory,
    ReportTable,
    ResultData,
    input_data_table,
    mapping_table,
)
from quantas.models.phonons import PhononInputData
from quantas.modules.qha.calculator import QHACalculator
from quantas.modules.qha.io.hdf5 import read_qha_hdf5 as _read_qha_hdf5
from quantas.modules.qha.io.reader import phonon_to_qha_input, read_qha_input
from quantas.modules.qha.inspect import PressureVolumePreview, pressure_volume_preview
from quantas.modules.qha.io.export import QHAHDF5Export
from quantas.modules.qha.models import QHAInput, QHAOptions, QHAResult
from quantas.modules.qha.plot import (
    QHAPlotOptions,
    build_qha_plot_collection,
    describe_qha_plots,
)
from quantas.modules.qha.report import (
    failed_points_table,
    result_summary_table,
    thermal_expansion_provenance_table,
    structural_property_tables,
)


def run_qha(
    input_data: QHAInput | PhononInputData | str | Path,
    options: QHAOptions | None = None,
    observer: Observer | None = None,
) -> ResultData:
    """Run a quasi-harmonic approximation calculation.

    Parameters
    ----------
    input_data : QHAInput, PhononInputData, str, or Path
        Normalized QHA input object, normalized phonon input object, or path to a
        Quantas phonon YAML input file.
    options : QHAOptions or None, optional
        Options controlling the QHA calculation. If ``None``, default options
        are used.
    observer : Observer or None, optional
        Observer receiving workflow events. If ``None``, the calculation runs
        silently through the default null observer used by the calculator.

    Returns
    -------
    ResultData
        Generic Quantas result object containing metadata, input summary,
        options, and the calculated QHA result.

    Raises
    ------
    TypeError
        If ``input_data`` is not a supported input type.
    ValueError
        If the input file or input object is invalid.
    """
    qha_input = normalize_qha_input(input_data)
    calculator = QHACalculator(
        qha_input=qha_input,
        options=options,
        observer=observer,
    )
    return calculator.execute()


def inspect_qha_input(
    input_data: QHAInput | PhononInputData | str | Path,
    options: QHAOptions | None = None,
    *,
    include_polynomial: bool = True,
    include_eos: bool = True,
    polynomial_degree: int | None = None,
    eos: str | None = None,
    maxfev: int | None = None,
) -> PressureVolumePreview:
    """Inspect the pressure-volume range covered by QHA input data.

    Parameters
    ----------
    input_data : QHAInput, PhononInputData, str, or Path
        Normalized QHA input object, normalized phonon input object, or path to a
        Quantas phonon YAML input file.
    options : QHAOptions or None, optional
        Options providing units, polynomial degree and default EOS.
    include_polynomial : bool, optional
        If ``True``, include pressure estimates from the polynomial fit.
    include_eos : bool, optional
        If ``True``, include pressure estimates from the EOS fit.
    polynomial_degree : int or None, optional
        Polynomial degree used for the static energy fit.
    eos : str or None, optional
        EOS used for the energy-volume fit.
    maxfev : int or None, optional
        Maximum number of optimizer evaluations for the EOS fit.

    Returns
    -------
    PressureVolumePreview
        Structured pressure-volume preview with fit diagnostics.

    Raises
    ------
    TypeError
        If ``input_data`` is not a supported input type.
    ValueError
        If input volume and energy arrays are missing or inconsistent.
    """
    qha_input = normalize_qha_input(input_data)
    return pressure_volume_preview(
        qha_input,
        options,
        include_polynomial=include_polynomial,
        include_eos=include_eos,
        polynomial_degree=polynomial_degree,
        eos=eos,
        maxfev=maxfev,
    )


def normalize_qha_input(
    input_data: QHAInput | PhononInputData | str | Path,
) -> QHAInput:
    """Return a normalized QHA input object.

    Parameters
    ----------
    input_data : QHAInput, PhononInputData, str, or Path
        Input object or path to a Quantas phonon YAML input file.

    Returns
    -------
    QHAInput
        Normalized QHA input object.

    Raises
    ------
    TypeError
        If the input value cannot be converted to :class:`QHAInput`.
    """
    if isinstance(input_data, QHAInput):
        return input_data
    if isinstance(input_data, PhononInputData):
        return phonon_to_qha_input(input_data)
    if isinstance(input_data, (str, Path)):
        return read_qha_input(input_data)
    raise TypeError(
        "input_data must be a QHAInput instance, a PhononInputData instance, "
        "or a path to a YAML file"
    )


def read_qha_hdf5(filename: str | Path) -> ResultData:
    """Read a complete native Quantas QHA result."""
    return _read_qha_hdf5(filename)


def write_qha_hdf5(
    result: ResultData,
    filename: str | Path,
    *,
    report_text: str | None = None,
) -> Path:
    """Write a complete native Quantas qha result.

    Parameters
    ----------
    result : ResultData
        Complete Quantas qha result.
    filename : str or Path
        Destination HDF5 file.
    report_text : str or None, optional
        Optional frontend-rendered report stored in diagnostics.

    Returns
    -------
    Path
        Final output path with the ``.hdf5`` suffix.
    """
    output = ensure_suffix(filename, ".hdf5")
    QHAHDF5Export().export(result, output, report_text=report_text)
    return output


def build_qha_report(result: ResultData) -> list[ReportTable]:
    """Build neutral report tables from a complete QHA result."""
    payload = result.results.get("qha")
    if result.metadata.module != "qha" or not isinstance(payload, QHAResult):
        raise ValueError("ResultData does not contain a valid QHA result.")
    tables: list[ReportTable] = []
    if result.input_data is not None:
        tables.append(
            input_data_table(
                "Input data",
                result.input_data.data,
                source=result.input_data.source,
            )
        )
    tables.append(mapping_table("Selected options", result.options))
    tables.append(result_summary_table(payload))
    provenance = thermal_expansion_provenance_table(payload)
    if provenance is not None:
        tables.append(provenance)
    tables.extend(structural_property_tables(payload))
    if payload.failed_points:
        tables.append(failed_points_table(payload))
    return tables


def build_qha_plots(
    result: ResultData,
    property_names: list[str] | tuple[str, ...] | None = None,
    options: QHAPlotOptions | None = None,
) -> PlotCollection:
    """Build neutral QHA plot specifications from a complete result."""
    payload = result.results.get("qha")
    if result.metadata.module != "qha" or not isinstance(payload, QHAResult):
        raise ValueError("ResultData does not contain a valid QHA result.")
    return build_qha_plot_collection(
        result,
        property_names=property_names,
        options=options,
    )


def describe_qha_plot_inventory(result: ResultData) -> PlotInventory:
    """Describe exact-grid QHA plot properties and representations."""
    payload = result.results.get("qha")
    if result.metadata.module != "qha" or not isinstance(payload, QHAResult):
        raise ValueError("ResultData does not contain a valid QHA result.")
    return describe_qha_plots(payload)


MODULE_CONTRACT = ModuleContract(
    name="qha",
    result_key="qha",
    read_input=read_qha_input,
    run=run_qha,
    read_hdf5=read_qha_hdf5,
    write_hdf5=write_qha_hdf5,
    build_report=build_qha_report,
    build_plots=build_qha_plots,
)
