# -*- coding: utf-8 -*-

"""Public API for equation-of-state fitting and post-fit analysis.

EOS archives can contain several datasets and immutable fit records, so this
namespace intentionally exposes a richer lifecycle than single-shot modules.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np

from quantas.core.events import Observer
from quantas.core.math.fitting import (
    CovarianceScaling,
    EffectiveVarianceOptions,
    FitMethod,
    ODRDifferenceScheme,
    ODROptions,
    OLSOptions,
    SolverOptions,
    WLSOptions,
    default_solver_options,
)
from quantas.core.physics.eos import MGDNormalization, PVTModel
from quantas.models import PlotCollection

from quantas.modules.eos import (
    EOSArchive as Archive,
    EOSBatchFailurePolicy as BatchFailurePolicy,
    EOSBatchJob as BatchJob,
    EOSBatchPlan as BatchPlan,
    EOSBatchResult as BatchResult,
    EOSCalculationResult as CalculationResult,
    EOSCapabilityStatus as CapabilityStatus,
    EOSDataset as Dataset,
    EOSDiagnosticResult as DiagnosticResult,
    EOSDomainCapability as DomainCapability,
    EOSFitDomain as FitDomain,
    EOSFitOptions as FitOptions,
    EOSFitRequest as FitRequest,
    EOSFitResult as FitResult,
    EOSPlotOptions as PlotOptions,
    EOSReportDetail as ReportDetail,
    EOSReportOptions as ReportOptions,
    EOSResultSlot as ResultSlot,
    EOSSession as Session,
    EOS_SPEC_TEMPLATE_FILENAME as SPEC_TEMPLATE_FILENAME,
    ParameterConstraint,
    build_eos_batch_preamble as build_batch_preamble,
    build_eos_batch_result_tables as build_batch_report,
    eos_calculation_summary_table as calculation_summary_table,
    eos_calculation_table as calculation_table,
    eos_diagnostic_summary_table as diagnostic_summary_table,
    eos_diagnostic_table as diagnostic_table,
    eos_domain_capability as domain_capability,
    read_eos_spec as read_spec,
    resolve_eos_spec as resolve_spec,
    write_eos_calculation_csv as write_calculation_csv,
    write_eos_diagnostic_csv as write_diagnostic_csv,
    write_eos_spec_template as write_spec_template,
)
from quantas.modules.eos.api import EOSFitter
from quantas.modules.eos.batch import EOSBatchWorkflow
from quantas.modules.eos.calculator import EOSCalculator
from quantas.modules.eos.plot import EOSPlotter
from quantas.modules.eos.facade import (
    build_eos_plots as _build_plots,
    calculate_eos as _calculate,
    diagnose_eos as _diagnose,
    fit_eos as _fit,
    open_eos_archive as _open_archive,
)
from quantas.modules.eos.io import read_eos_input as _read_input
from quantas.modules.eos.contracts import EOS_DOMAIN_CAPABILITIES as DOMAIN_CAPABILITIES

from .common import _public_dir


def read_input(
    source: str | Path,
    *,
    pressure_unit: str | None = None,
    length_unit: str | None = None,
    temperature_unit: str | None = None,
) -> Dataset:
    """Read and normalize one EOS dataset.

    Parameters
    ----------
    source : str or Path
        Keyword-directed EOS input table.
    pressure_unit, length_unit, temperature_unit : str or None, optional
        Explicit unit overrides. File declarations and documented defaults are
        used when omitted.

    Returns
    -------
    Dataset
        Validated float64 dataset with normalized units and selection metadata.

    Raises
    ------
    ValueError
        If the file format, units, or scientific columns are invalid.
    """
    return _read_input(
        source,
        pressure_unit=pressure_unit,
        length_unit=length_unit,
        temperature_unit=temperature_unit,
    )


def normalize_input(
    source: Dataset | str | Path,
    *,
    pressure_unit: str | None = None,
    length_unit: str | None = None,
    temperature_unit: str | None = None,
) -> Dataset:
    """Return a normalized EOS dataset.

    Parameters
    ----------
    source : Dataset, str, or Path
        Existing dataset contract or input path.
    pressure_unit, length_unit, temperature_unit : str or None, optional
        Explicit unit overrides used only when reading a path.

    Returns
    -------
    Dataset
        Validated EOS dataset.

    Raises
    ------
    TypeError
        If ``source`` is neither a dataset nor a path.
    """
    if isinstance(source, Dataset):
        return source
    if isinstance(source, (str, Path)):
        return read_input(
            source,
            pressure_unit=pressure_unit,
            length_unit=length_unit,
            temperature_unit=temperature_unit,
        )
    raise TypeError("source must be an EOS Dataset object or path")


def fit(
    input_data: Dataset | str | Path,
    request: FitRequest,
    *,
    pressure_unit: str | None = None,
    length_unit: str | None = None,
    temperature_unit: str | None = None,
) -> FitResult:
    """Fit one EOS request.

    Parameters
    ----------
    input_data : Dataset, str, or Path
        Normalized EOS dataset or input table path.
    request : FitRequest
        Model, scientific domain, target, parameter constraints, and solver
        options.
    pressure_unit, length_unit, temperature_unit : str or None, optional
        Unit overrides used when ``input_data`` is a path.

    Returns
    -------
    FitResult
        Immutable fit result including parameters, uncertainties, residuals,
        warnings, and numerical diagnostics.

    Raises
    ------
    ValueError
        If the request is incompatible with the dataset or the selected model.
    """
    return _fit(
        input_data,
        request,
        pressure_unit=pressure_unit,
        length_unit=length_unit,
        temperature_unit=temperature_unit,
    )


def run_batch(
    input_data: Dataset | str | Path,
    plan: BatchPlan,
    archive: str | Path,
    *,
    observer: Observer | None = None,
    overwrite: bool = False,
    creator: str = "quantas eos run",
    pressure_unit: str | None = None,
    length_unit: str | None = None,
    temperature_unit: str | None = None,
) -> BatchResult:
    """Run and persist an EOS batch plan.

    Parameters
    ----------
    input_data : Dataset, str, or Path
        Normalized EOS dataset or input path.
    plan : BatchPlan
        Ordered fitting jobs and acceptance policy.
    archive : str or Path
        Destination native EOS HDF5 archive.
    observer : Observer or None, optional
        Frontend-neutral workflow observer.
    overwrite : bool, optional
        Replace an existing archive.
    creator : str, optional
        Provenance identifier written to archive metadata.
    pressure_unit, length_unit, temperature_unit : str or None, optional
        Unit overrides used when reading a path.

    Returns
    -------
    BatchResult
        Structured summary of attempted, accepted, and failed fits.

    Raises
    ------
    FileExistsError
        If the archive exists and ``overwrite`` is false.
    ValueError
        If the plan or input dataset is invalid.
    """
    workflow = EOSBatchWorkflow(observer=observer)
    return workflow.run(
        input_data,
        plan,
        archive,
        overwrite=overwrite,
        creator=creator,
        pressure_unit=pressure_unit,
        length_unit=length_unit,
        temperature_unit=temperature_unit,
    )


def validate_request(input_data: Dataset | str | Path, request: FitRequest) -> None:
    """Validate one EOS request against a normalized dataset.

    Parameters
    ----------
    input_data : Dataset, str, or Path
        EOS dataset or input path.
    request : FitRequest
        Fit request to validate without executing the solver.

    Raises
    ------
    ValueError
        If the domain, target, model, data selection, or constraints are
        incompatible.
    """
    dataset = normalize_input(input_data)
    EOSFitter().validate_request(dataset, request)


def available_plot_types(
    archive: str | Path,
    *,
    slot: str | ResultSlot | None = None,
    record_id: int | None = None,
) -> tuple[str, ...]:
    """Return plots available for one EOS archive record.

    Parameters
    ----------
    archive : str or Path
        Native EOS HDF5 archive.
    slot : str, ResultSlot, or None, optional
        Accepted ``domain/target`` slot.
    record_id : int or None, optional
        Explicit immutable record identifier.

    Returns
    -------
    tuple of str
        Stable plot-type identifiers supported by the selected record.

    Raises
    ------
    ValueError
        If record selection is ambiguous or invalid.
    """
    plotter = EOSPlotter.from_archive(archive, slot=slot, record_id=record_id)
    return tuple(plotter.available_plot_types())


def record_domain(
    archive: str | Path,
    *,
    slot: str | ResultSlot | None = None,
    record_id: int | None = None,
) -> FitDomain:
    """Return the scientific domain of one EOS archive record.

    Parameters
    ----------
    archive : str or Path
        Native EOS HDF5 archive.
    slot : str, ResultSlot, or None, optional
        Accepted result slot.
    record_id : int or None, optional
        Explicit immutable record identifier.

    Returns
    -------
    FitDomain
        ``pv``, ``vt``, or ``pvt`` domain of the selected record.

    Raises
    ------
    ValueError
        If record selection is ambiguous or invalid.
    """
    calculator = EOSCalculator.from_archive(archive, slot=slot, record_id=record_id)
    return calculator.slot.domain


def open_archive(path: str | Path, *, writable: bool = False) -> Archive:
    """Open a native EOS archive.

    Parameters
    ----------
    path : str or Path
        Existing native EOS HDF5 archive.
    writable : bool, optional
        Open in append/update mode when true.

    Returns
    -------
    Archive
        Active archive handle. Use it as a context manager or call ``close``.

    Raises
    ------
    ValueError
        If the file is not a supported EOS archive.
    """
    return _open_archive(path, writable=writable)


def calculate(
    archive: str | Path,
    *,
    slot: str | ResultSlot | None = None,
    record_id: int | None = None,
    pressure: np.ndarray | Sequence[float] | float | None = None,
    volume: np.ndarray | Sequence[float] | float | None = None,
    temperature: np.ndarray | Sequence[float] | float | None = None,
    propagate_uncertainty: bool = True,
    relative_step: float = 1.0e-5,
) -> CalculationResult:
    """Evaluate fitted EOS properties from a native archive.

    Parameters
    ----------
    archive : str or Path
        Native EOS HDF5 archive.
    slot : str, ResultSlot, or None, optional
        Accepted result slot.
    record_id : int or None, optional
        Explicit immutable record identifier.
    pressure, volume, temperature : array-like, float, or None, optional
        Coordinates for domain-appropriate property evaluation.
    propagate_uncertainty : bool, optional
        Propagate the stored parameter covariance to calculated properties.
    relative_step : float, optional
        Relative finite-difference step used for uncertainty transformations
        without an analytical Jacobian.

    Returns
    -------
    CalculationResult
        Calculated state table, uncertainties, metadata, and warnings.

    Raises
    ------
    ValueError
        If coordinates are incomplete, incompatible, or outside mathematical
        model bounds.
    """
    return _calculate(
        archive,
        slot=slot,
        record_id=record_id,
        pressure=pressure,
        volume=volume,
        temperature=temperature,
        propagate_uncertainty=propagate_uncertainty,
        relative_step=relative_step,
    )


def diagnose(
    archive: str | Path,
    *,
    slot: str | ResultSlot | None = None,
    record_id: int | None = None,
    include_normalized_pressure: bool = True,
) -> DiagnosticResult:
    """Build residual and finite-strain diagnostics from an EOS archive.

    Parameters
    ----------
    archive : str or Path
        Native EOS HDF5 archive.
    slot : str, ResultSlot, or None, optional
        Accepted result slot.
    record_id : int or None, optional
        Explicit immutable record identifier.
    include_normalized_pressure : bool, optional
        Include finite-strain and normalized-pressure columns where defined.

    Returns
    -------
    DiagnosticResult
        Observation-level residuals, summary statistics, and warnings.

    Raises
    ------
    ValueError
        If record selection is invalid or the diagnostic transformation is not
        defined for the selected model.
    """
    return _diagnose(
        archive,
        slot=slot,
        record_id=record_id,
        include_normalized_pressure=include_normalized_pressure,
    )


def build_plots(
    archive: str | Path,
    plot_types: Sequence[str] | str | None = None,
    *,
    slot: str | ResultSlot | None = None,
    record_id: int | None = None,
    options: PlotOptions | None = None,
) -> PlotCollection:
    """Build frontend-neutral EOS plot specifications.

    Parameters
    ----------
    archive : str or Path
        Native EOS HDF5 archive.
    plot_types : sequence of str, str, or None, optional
        Requested plot identifiers. Standard plots are selected when omitted.
    slot : str, ResultSlot, or None, optional
        Accepted result slot.
    record_id : int or None, optional
        Explicit immutable record identifier.
    options : PlotOptions or None, optional
        Uncertainty, curve, residual, and surface presentation metadata.

    Returns
    -------
    PlotCollection
        Neutral plot specifications ready for a frontend renderer.

    Raises
    ------
    ValueError
        If a plot type is unsupported or record selection is invalid.
    """
    return _build_plots(
        archive,
        plot_types,
        slot=slot,
        record_id=record_id,
        options=options,
    )


def __dir__() -> list[str]:
    """Return only names in the supported public contract."""
    return _public_dir(__all__)


__all__ = [
    "Archive",
    "BatchFailurePolicy",
    "BatchJob",
    "BatchPlan",
    "BatchResult",
    "CalculationResult",
    "CovarianceScaling",
    "CapabilityStatus",
    "DOMAIN_CAPABILITIES",
    "Dataset",
    "DiagnosticResult",
    "DomainCapability",
    "EffectiveVarianceOptions",
    "FitDomain",
    "FitMethod",
    "FitOptions",
    "FitRequest",
    "FitResult",
    "MGDNormalization",
    "ODRDifferenceScheme",
    "ODROptions",
    "OLSOptions",
    "ParameterConstraint",
    "PVTModel",
    "PlotOptions",
    "ReportDetail",
    "ReportOptions",
    "ResultSlot",
    "SPEC_TEMPLATE_FILENAME",
    "Session",
    "SolverOptions",
    "available_plot_types",
    "build_batch_preamble",
    "build_batch_report",
    "build_plots",
    "calculate",
    "diagnose",
    "domain_capability",
    "calculation_summary_table",
    "calculation_table",
    "diagnostic_summary_table",
    "diagnostic_table",
    "fit",
    "normalize_input",
    "open_archive",
    "read_input",
    "read_spec",
    "record_domain",
    "resolve_spec",
    "run_batch",
    "validate_request",
    "WLSOptions",
    "default_solver_options",
    "write_calculation_csv",
    "write_diagnostic_csv",
    "write_spec_template",
]
