# -*- coding: utf-8 -*-

"""Public API for harmonic thermodynamic calculations."""

from __future__ import annotations

from pathlib import Path

from quantas.core.events import Observer
from quantas.models import PlotCollection, ReportTable, ResultData
from quantas.models.phonons import PhononInputData
from quantas.modules.ha.api import (
    build_ha_plots as _build_plots,
    build_ha_report as _build_report,
    create_ha_input as _create_input,
    normalize_ha_input as _normalize_input,
    read_ha_hdf5 as _read_result,
    read_ha_input as _read_input,
    run_ha as _run,
    write_ha_hdf5 as _write_result,
)
from quantas.modules.ha.models import (
    HAInput as Input,
    HAOptions as Options,
    HAResult as Result,
)

from .common import _public_dir, get_result_payload


def create_input(
    source: str | Path,
    destination: str | Path,
    *,
    interface: str = "crystal",
    is_list: bool = False,
    reference: int = 0,
    jobname: str = "Quantas HA input",
    formula_units: int = 1,
) -> Path:
    """Create a normalized HA YAML input from an interface output.

    Parameters
    ----------
    source : str or Path
        Quantum-mechanical output file or file list.
    destination : str or Path
        Destination YAML path.
    interface : str, optional
        Code-specific reader identifier.
    is_list : bool, optional
        Interpret ``source`` as a text file containing multiple output paths.
    reference : int, optional
        Reference structure index for multi-structure inputs.
    jobname : str, optional
        Human-readable workflow title.
    formula_units : int, optional
        Formula units represented by the crystallographic cell.

    Returns
    -------
    Path
        Written YAML input path.

    Raises
    ------
    ValueError
        If source data cannot be parsed or normalized.
    """
    return _create_input(
        source,
        destination,
        interface=interface,
        is_list=is_list,
        reference=reference,
        jobname=jobname,
        formula_units=formula_units,
    )


def read_input(source: str | Path) -> Input:
    """Read one Quantas HA input file.

    Parameters
    ----------
    source : str or Path
        Quantas phonon YAML path.

    Returns
    -------
    Input
        Validated harmonic input contract.

    Raises
    ------
    ValueError
        If the input is malformed or incomplete.
    """
    return _read_input(source)


def normalize_input(source: Input | PhononInputData | str | Path) -> Input:
    """Return a normalized HA input contract.

    Parameters
    ----------
    source : Input, PhononInputData, str, or Path
        Existing harmonic/phonon contract or YAML path.

    Returns
    -------
    Input
        Validated harmonic input suitable for :func:`run`.

    Raises
    ------
    TypeError
        If the input type is unsupported.
    ValueError
        If a supplied file cannot be parsed.
    """
    return _normalize_input(source)


def run(
    input_data: Input | PhononInputData | str | Path,
    options: Options | None = None,
    observer: Observer | None = None,
) -> ResultData:
    """Run a harmonic thermodynamic workflow.

    Parameters
    ----------
    input_data : Input, PhononInputData, str, or Path
        Harmonic input contract, neutral phonon data, or YAML path.
    options : Options or None, optional
        Temperature grid, units, and scientific calculation controls.
    observer : Observer or None, optional
        Frontend-neutral event observer.

    Returns
    -------
    ResultData
        Complete result envelope containing a harmonic payload.

    Raises
    ------
    ValueError
        If the input or selected temperature domain is invalid.
    """
    return _run(input_data, options=options, observer=observer)


def get_result(result: ResultData) -> Result:
    """Return the typed HA payload from a result envelope.

    Parameters
    ----------
    result : ResultData
        Complete Quantas result envelope.

    Returns
    -------
    Result
        Module-specific harmonic result.

    Raises
    ------
    ValueError
        If the envelope is not a valid HA result.
    """
    return get_result_payload(result, module="ha", key="ha", expected_type=Result)


def read_result(source: str | Path) -> ResultData:
    """Read a native Quantas HA HDF5 result.

    Parameters
    ----------
    source : str or Path
        Native Quantas HDF5 path.

    Returns
    -------
    ResultData
        Restored result envelope.

    Raises
    ------
    ValueError
        If the file is not a supported HA result.
    """
    return _read_result(source)


def write_result(
    result: ResultData,
    destination: str | Path,
    *,
    report_text: str | None = None,
) -> Path:
    """Write a native Quantas HA HDF5 result.

    Parameters
    ----------
    result : ResultData
        Complete HA result envelope.
    destination : str or Path
        Destination path.
    report_text : str or None, optional
        Deterministic report text to embed in diagnostics.

    Returns
    -------
    Path
        Final native HDF5 path.

    Raises
    ------
    ValueError
        If the result envelope is invalid.
    """
    get_result(result)
    return _write_result(result, destination, report_text=report_text)


def build_report(result: ResultData) -> list[ReportTable]:
    """Build frontend-neutral HA report tables.

    Parameters
    ----------
    result : ResultData
        Complete harmonic result envelope.

    Returns
    -------
    list of ReportTable
        Ordered raw-value report tables.

    Raises
    ------
    ValueError
        If the result envelope is invalid.
    """
    get_result(result)
    return _build_report(result)


def build_plots(
    result: ResultData,
    properties: str | list[str] | tuple[str, ...] | None = None,
    *,
    unit: str | None = None,
) -> PlotCollection:
    """Build frontend-neutral HA plots.

    Parameters
    ----------
    result : ResultData
        Complete harmonic result envelope.

    Returns
    -------
    PlotCollection
        Neutral thermodynamic plot specifications.

    Raises
    ------
    ValueError
        If the result envelope is invalid.
    """
    get_result(result)
    return _build_plots(result, properties=properties, unit=unit)


def __dir__() -> list[str]:
    """Return only names in the supported public contract."""
    return _public_dir(__all__)


__all__ = [
    "Input",
    "Options",
    "Result",
    "build_plots",
    "build_report",
    "create_input",
    "get_result",
    "normalize_input",
    "read_input",
    "read_result",
    "run",
    "write_result",
]
