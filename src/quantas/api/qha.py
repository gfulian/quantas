# -*- coding: utf-8 -*-

"""Public API for quasi-harmonic thermodynamic calculations."""

from __future__ import annotations

from pathlib import Path

from quantas.core.events import Observer
from quantas.modules.qha.api import (
    build_qha_plots as _build_plots,
    build_qha_report as _build_report,
    describe_qha_plot_inventory as _describe_plots,
    inspect_qha_input as _inspect,
    normalize_qha_input as _normalize_input,
    read_qha_hdf5 as _read_result,
    read_qha_input as _read_input,
    run_qha as _run,
    write_qha_hdf5 as _write_result,
)
from quantas.modules.qha.inspect import PressureVolumePreview as Preview
from quantas.modules.qha.models import (
    QHAFitFailurePolicy as FitFailurePolicy,
    QHAInput as Input,
    QHAMinimization as Minimization,
    QHAModeContinuity as ModeContinuity,
    QHAOptions as Options,
    QHAPolynomialDerivativeMethod as PolynomialDerivativeMethod,
    QHAResult as Result,
    QHAScheme as Scheme,
    QHAThermalExpansionMethod as ThermalExpansionMethod,
)
from quantas.modules.qha.plot import (
    QHACurveAxis as CurveAxis,
    QHAPlotOptions as PlotOptions,
    list_available_plot_properties,
)
from quantas.modules.qha.validation import (
    PropertyDifference,
    QHAValidationSummary as ValidationSummary,
    compare_qha_results as compare_results,
    validate_qha_result as validate_result,
)

from .common import (
    PhononInputData,
    ReportTable,
    ResultData,
    _public_dir,
    get_result_payload,
)
from .plotting import PlotCollection, PlotInventory


def read_input(source: str | Path) -> Input:
    """Read one Quantas QHA input file.

    Parameters
    ----------
    source : str or Path
        Quantas multi-volume phonon YAML path.

    Returns
    -------
    Input
        Validated quasi-harmonic input contract.

    Raises
    ------
    ValueError
        If the file is malformed or lacks required volume/phonon data.
    """
    return _read_input(source)


def normalize_input(source: Input | PhononInputData | str | Path) -> Input:
    """Return a normalized QHA input contract.

    Parameters
    ----------
    source : Input, PhononInputData, str, or Path
        Existing QHA/phonon contract or YAML path.

    Returns
    -------
    Input
        Validated QHA input suitable for :func:`inspect` and :func:`run`.

    Raises
    ------
    TypeError
        If the input type is unsupported.
    ValueError
        If a supplied file cannot be parsed.
    """
    return _normalize_input(source)


def inspect(
    input_data: Input | PhononInputData | str | Path,
    options: Options | None = None,
    *,
    include_polynomial: bool = True,
    include_eos: bool = True,
    polynomial_degree: int | None = None,
    eos: str | None = None,
    maxfev: int | None = None,
) -> Preview:
    """Inspect pressure-volume behavior without running the full workflow.

    Parameters
    ----------
    input_data : Input, PhononInputData, str, or Path
        QHA input contract, neutral phonon data, or YAML path.
    options : Options or None, optional
        QHA defaults used by the preview.
    include_polynomial : bool, optional
        Include the static polynomial preview.
    include_eos : bool, optional
        Include the energy-EOS preview.
    polynomial_degree : int or None, optional
        Explicit static-energy polynomial degree.
    eos : str or None, optional
        Explicit integrated energy-EOS tag.
    maxfev : int or None, optional
        Maximum optimizer evaluations for the EOS preview.

    Returns
    -------
    Preview
        Structured pressure-volume range and fit diagnostics.

    Raises
    ------
    ValueError
        If required volume/energy data are unavailable or a fit is invalid.
    """
    return _inspect(
        input_data,
        options=options,
        include_polynomial=include_polynomial,
        include_eos=include_eos,
        polynomial_degree=polynomial_degree,
        eos=eos,
        maxfev=maxfev,
    )


def run(
    input_data: Input | PhononInputData | str | Path,
    options: Options | None = None,
    observer: Observer | None = None,
) -> ResultData:
    """Run a quasi-harmonic thermodynamic workflow.

    Parameters
    ----------
    input_data : Input, PhononInputData, str, or Path
        QHA input contract, neutral phonon data, or YAML path.
    options : Options or None, optional
        Thermodynamic domain, fitting, minimization, and unit controls.
    observer : Observer or None, optional
        Frontend-neutral event observer.

    Returns
    -------
    ResultData
        Complete result envelope containing a QHA payload.

    Raises
    ------
    ValueError
        If the input, fits, or requested thermodynamic domain are invalid.
    """
    return _run(input_data, options=options, observer=observer)


def get_result(result: ResultData) -> Result:
    """Return the typed QHA payload from a result envelope.

    Parameters
    ----------
    result : ResultData
        Complete Quantas result envelope.

    Returns
    -------
    Result
        Module-specific quasi-harmonic result.

    Raises
    ------
    ValueError
        If the envelope is not a valid QHA result.
    """
    return get_result_payload(result, module="qha", key="qha", expected_type=Result)


def read_result(source: str | Path) -> ResultData:
    """Read a native Quantas QHA HDF5 result.

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
        If the file is not a supported QHA result.
    """
    return _read_result(source)


def write_result(
    result: ResultData,
    destination: str | Path,
    *,
    report_text: str | None = None,
) -> Path:
    """Write a native Quantas QHA HDF5 result.

    Parameters
    ----------
    result : ResultData
        Complete QHA result envelope.
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
    """Build frontend-neutral QHA report tables.

    Parameters
    ----------
    result : ResultData
        Complete QHA result envelope.

    Returns
    -------
    list of ReportTable
        Ordered raw-value thermodynamic and structural tables.

    Raises
    ------
    ValueError
        If the result envelope is invalid.
    """
    get_result(result)
    return _build_report(result)


def build_plots(
    result: ResultData,
    properties: list[str] | tuple[str, ...] | None = None,
    options: PlotOptions | None = None,
) -> PlotCollection:
    """Build frontend-neutral QHA plots.

    Parameters
    ----------
    result : ResultData
        Complete QHA result envelope.
    properties : list of str, tuple of str, or None, optional
        Public property identifiers. Defaults select the standard plot family.
    options : PlotOptions or None, optional
        Units, contour, and styling metadata for neutral plot construction.

    Returns
    -------
    PlotCollection
        Neutral one- and two-dimensional plot specifications.

    Raises
    ------
    ValueError
        If a property is unknown or required result data are unavailable.
    """
    get_result(result)
    return _build_plots(result, property_names=properties, options=options)


def describe_plots(result: ResultData) -> PlotInventory:
    """Return exact-grid QHA properties, sections, maps, and coordinates.

    Parameters
    ----------
    result : ResultData
        Complete quasi-harmonic result envelope.

    Returns
    -------
    PlotInventory
        Result-aware property metadata, temperature and pressure grids, and
        supported line-section and contour representations.
    """
    get_result(result)
    return _describe_plots(result)


def list_plot_properties(
    result: ResultData | Result,
) -> list[tuple[str, str, str]]:
    """List plot properties available in one QHA result.

    Parameters
    ----------
    result : ResultData or Result
        Complete envelope or module-specific QHA result.

    Returns
    -------
    list of tuple of str
        ``(key, attribute, description)`` records in stable display order.
    """
    return list_available_plot_properties(result)


def __dir__() -> list[str]:
    """Return only names in the supported public contract."""
    return _public_dir(__all__)


__all__ = [
    "CurveAxis",
    "FitFailurePolicy",
    "Input",
    "Minimization",
    "ModeContinuity",
    "Options",
    "PlotOptions",
    "PolynomialDerivativeMethod",
    "Preview",
    "PropertyDifference",
    "Result",
    "Scheme",
    "ThermalExpansionMethod",
    "ValidationSummary",
    "list_plot_properties",
    "build_plots",
    "build_report",
    "compare_results",
    "describe_plots",
    "get_result",
    "inspect",
    "normalize_input",
    "read_input",
    "read_result",
    "run",
    "validate_result",
    "write_result",
]
