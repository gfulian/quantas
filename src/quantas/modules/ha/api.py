# -*- coding: utf-8 -*-

"""Internal HA workflow facade behind :mod:`quantas.api.ha`.

The functions read normalized phonon inputs, execute harmonic thermodynamics,
and build persistent results without printing or creating frontend objects.
"""

from __future__ import annotations

from pathlib import Path

from quantas.io.path import ensure_suffix
from quantas.core.events import Observer
from quantas.models import (
    ModuleContract,
    PhononInputData,
    PlotCollection,
    ReportTable,
    ResultData,
    input_data_table,
    mapping_table,
)
from quantas.modules.ha.calculator import HACalculator
from quantas.modules.ha.io.reader import (
    HAInputFileReader,
    phonon_to_ha_input,
    read_ha_hdf5,
)
from quantas.modules.ha.io.inpgen import create_ha_input as create_ha_input_file
from quantas.modules.ha.io.export import HAHDF5Export
from quantas.modules.ha.models import HAInput, HAOptions, HAResult
from quantas.modules.ha.plot import build_ha_plot_collection
from quantas.modules.ha.report import (
    static_data_table,
    thermodynamic_property_tables,
    thermodynamic_summary_table,
)


def read_ha_input(filename: str | Path) -> HAInput:
    """
    Read a Quantas phonon YAML input file.

    Parameters
    ----------
    filename : str or Path
        Path to the Quantas YAML input file containing volumes, static
        energies, q-point weights, and phonon frequencies.

    Returns
    -------
    HAInput
        Normalized harmonic-approximation input data.

    Raises
    ------
    ValueError
        If the input file cannot be read, parsed, or validated by the reader.
    """
    reader = HAInputFileReader(filename)

    if not reader.completed:
        raise ValueError(reader.error or "Unable to read HA input file")

    return reader.to_input(source=Path(filename))


def run_ha(
    input_data: HAInput | PhononInputData | str | Path,
    options: HAOptions | None = None,
    observer: Observer | None = None,
) -> ResultData:
    """
    Run a harmonic-approximation calculation.

    Parameters
    ----------
    input_data : HAInput, PhononInputData, str, or Path
        Normalized HA input object or path to a Quantas phonon YAML input file.
    options : HAOptions or None, optional
        Options controlling the harmonic calculation. If ``None``, default
        options are used.
    observer : Observer or None, optional
        Observer receiving workflow events. If ``None``, the calculation runs
        silently through the default null observer used by the calculator.

    Returns
    -------
    ResultData
        Generic Quantas result object containing metadata, input summary,
        options, and the calculated HA result.

    Raises
    ------
    ValueError
        If the input file or input object is invalid.
    """
    ha_input = normalize_ha_input(input_data)

    calculator = HACalculator(
        ha_input=ha_input,
        options=options,
        observer=observer,
    )

    return calculator.execute()


def normalize_ha_input(
    input_data: HAInput | PhononInputData | str | Path,
) -> HAInput:
    """Return a normalized harmonic input object.

    Parameters
    ----------
    input_data : HAInput, PhononInputData, str, or Path
        Harmonic input, neutral phonon input, or YAML file path.

    Returns
    -------
    HAInput
        Normalized harmonic-approximation input data.

    Raises
    ------
    TypeError
        If the input value cannot be converted to :class:`HAInput`.
    """
    if isinstance(input_data, HAInput):
        return input_data
    if isinstance(input_data, PhononInputData):
        return phonon_to_ha_input(input_data)
    if isinstance(input_data, (str, Path)):
        return read_ha_input(input_data)
    raise TypeError(
        "input_data must be a HAInput instance, a PhononInputData instance, "
        "or a path to a YAML file"
    )


def create_ha_input(
    filename: str | Path,
    outfile: str | Path,
    *,
    interface: str = "crystal",
    is_list: bool = False,
    reference: int = 0,
    jobname: str = "Quantas HA input",
    formula_units: int = 1,
) -> Path:
    """
    Generate a Quantas phonon YAML input file from QM output data.

    Parameters
    ----------
    filename : str or Path
        QM output file or text file containing a list of QM output files.
    outfile : str or Path
        Destination YAML file.
    interface : str, optional
        Interface used to read the QM output data.
    is_list : bool, optional
        If ``True``, ``filename`` is interpreted as a file list.
    reference : int, optional
        Reference input index for multiple-file HA/QHA generation.
    jobname : str, optional
        Description written to the YAML ``job`` field.
    formula_units : int, optional
        Number of chemical formula units in the crystallographic cell.

    Returns
    -------
    Path
        Path of the generated YAML file.

    Raises
    ------
    ValueError
        If the selected interface cannot read the provided files.
    OSError
        If the output file cannot be written.
    """
    return create_ha_input_file(
        filename,
        outfile,
        interface=interface,
        is_list=is_list,
        reference=reference,
        jobname=jobname,
        formula_units=formula_units,
    )


def write_ha_hdf5(
    result: ResultData,
    filename: str | Path,
    *,
    report_text: str | None = None,
) -> Path:
    """Write a complete native Quantas ha result.

    Parameters
    ----------
    result : ResultData
        Complete Quantas ha result.
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
    HAHDF5Export().export(result, output, report_text=report_text)
    return output


def build_ha_report(result: ResultData) -> list[ReportTable]:
    """Build neutral report tables from a complete HA result.

    Parameters
    ----------
    result : ResultData
        Generic Quantas HA result.

    Returns
    -------
    list of ReportTable
        Frontend-neutral report tables.

    Raises
    ------
    ValueError
        If the result does not contain a HA payload.
    """
    payload = result.results.get("ha")
    if result.metadata.module != "ha" or not isinstance(payload, HAResult):
        raise ValueError("ResultData does not contain a valid HA result.")
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
    tables.append(static_data_table(payload))
    tables.append(thermodynamic_summary_table(payload))
    tables.extend(
        thermodynamic_property_tables(
            payload,
            row_indices=(0, -1),
            include_zero_point=False,
        )
    )
    return tables


def build_ha_plots(
    result: ResultData,
    properties: str | list[str] | tuple[str, ...] | None = None,
    *,
    unit: str | None = None,
) -> PlotCollection:
    """Build neutral HA plot specifications from a complete result.

    Parameters
    ----------
    result : ResultData
        Generic Quantas HA result.
    properties : str, list of str, tuple of str, or None, optional
        Requested thermodynamic properties.
    unit : str or None, optional
        Requested display energy unit.

    Returns
    -------
    PlotCollection
        Neutral HA plot specifications.
    """
    payload = result.results.get("ha")
    if result.metadata.module != "ha" or not isinstance(payload, HAResult):
        raise ValueError("ResultData does not contain a valid HA result.")
    return build_ha_plot_collection(payload, properties=properties, unit=unit)


MODULE_CONTRACT = ModuleContract(
    name="ha",
    result_key="ha",
    read_input=read_ha_input,
    run=run_ha,
    read_hdf5=read_ha_hdf5,
    write_hdf5=write_ha_hdf5,
    build_report=build_ha_report,
    build_plots=build_ha_plots,
)
