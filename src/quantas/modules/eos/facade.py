# -*- coding: utf-8 -*-

"""Internal EOS workflow facade used by :mod:`quantas.api.eos`.

The public namespace owns the supported application contract.  This module
keeps implementation functions focused and movable without exposing concrete
calculator, diagnostic, or plotter classes to frontends.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np

from .api import EOSFitter
from .archive import EOSArchive
from .batch import run_eos_batch
from .calculator import EOSCalculationResult, EOSCalculator
from .contracts import (
    EOS_ARCHIVE_SCHEMA_VERSION,
    EOS_DOMAIN_CAPABILITIES,
    EOS_SUPPORTED_ARCHIVE_SCHEMA_VERSIONS,
    EOSModuleContract,
)
from .diagnostics import EOSDiagnosticResult, EOSDiagnostics
from .history import EOSResultSlot
from .io import read_eos_input
from .models import EOSDataset, EOSFitRequest, EOSFitResult
from .plot import EOSPlotOptions, EOSPlotter, describe_eos_plots


def fit_eos(
    input_data: EOSDataset | str | Path,
    request: EOSFitRequest,
    *,
    pressure_unit: str | None = None,
    length_unit: str | None = None,
    temperature_unit: str | None = None,
    fitter: EOSFitter | None = None,
) -> EOSFitResult:
    """Fit one EOS request through the stable Python facade.

    Parameters
    ----------
    input_data : EOSDataset or str or Path
        Normalized dataset or keyword-directed EOS input file.
    request : EOSFitRequest
        Complete model, target, constraints, and solver request.
    pressure_unit, length_unit, temperature_unit : str or None, optional
        Unit overrides applied only when ``input_data`` is a path.
    fitter : EOSFitter or None, optional
        Optional reusable fitting service.

    Returns
    -------
    EOSFitResult
        Complete frontend-neutral result.
    """
    dataset = (
        read_eos_input(
            input_data,
            pressure_unit=pressure_unit,
            length_unit=length_unit,
            temperature_unit=temperature_unit,
        )
        if isinstance(input_data, (str, Path))
        else input_data
    )
    service = EOSFitter() if fitter is None else fitter
    return service.fit(dataset, request)


def open_eos_archive(path: str | Path, *, writable: bool = False) -> EOSArchive:
    """Open one native EOS archive through the stable facade.

    Parameters
    ----------
    path : str or Path
        Native EOS HDF5 archive.
    writable : bool, optional
        Open in update mode when ``True``; read-only otherwise.

    Returns
    -------
    EOSArchive
        Active context-manager-compatible archive.
    """
    return EOSArchive(path, mode="r+" if writable else "r")


def calculate_eos(
    archive: str | Path,
    *,
    slot: str | EOSResultSlot | None = None,
    record_id: int | None = None,
    pressure: np.ndarray | Sequence[float] | float | None = None,
    volume: np.ndarray | Sequence[float] | float | None = None,
    temperature: np.ndarray | Sequence[float] | float | None = None,
    propagate_uncertainty: bool = True,
    relative_step: float = 1.0e-5,
) -> EOSCalculationResult:
    """Evaluate fitted EOS properties from a native archive."""
    calculator = EOSCalculator.from_archive(
        archive,
        slot=slot,
        record_id=record_id,
    )
    return calculator.calculate(
        pressure=pressure,
        volume=volume,
        temperature=temperature,
        propagate_uncertainty=propagate_uncertainty,
        relative_step=relative_step,
    )


def diagnose_eos(
    archive: str | Path,
    *,
    slot: str | EOSResultSlot | None = None,
    record_id: int | None = None,
    include_normalized_pressure: bool = True,
) -> EOSDiagnosticResult:
    """Build residual and finite-strain diagnostics from a native archive."""
    diagnostics = EOSDiagnostics.from_archive(
        archive,
        slot=slot,
        record_id=record_id,
    )
    return diagnostics.build(include_normalized_pressure=include_normalized_pressure)


def build_eos_plots(
    archive: str | Path,
    plot_types: Sequence[str] | str | None = None,
    *,
    slot: str | EOSResultSlot | None = None,
    record_id: int | None = None,
    options: EOSPlotOptions | None = None,
):
    """Build frontend-neutral EOS plot specifications from an archive.

    Returns
    -------
    PlotCollection
        Neutral plots suitable for the Matplotlib renderer or another frontend adapter.
    """
    plotter = EOSPlotter.from_archive(
        archive,
        slot=slot,
        record_id=record_id,
    )
    return plotter.build(plot_types, options=options)


MODULE_CONTRACT = EOSModuleContract(
    name="eos",
    archive_schema_version=EOS_ARCHIVE_SCHEMA_VERSION,
    supported_archive_schema_versions=tuple(
        sorted(EOS_SUPPORTED_ARCHIVE_SCHEMA_VERSIONS)
    ),
    capabilities=EOS_DOMAIN_CAPABILITIES,
    read_input=read_eos_input,
    fit=fit_eos,
    run_batch=run_eos_batch,
    open_archive=open_eos_archive,
    calculate=calculate_eos,
    diagnose=diagnose_eos,
    describe_plots=describe_eos_plots,
    build_plots=build_eos_plots,
)


__all__ = [
    "MODULE_CONTRACT",
    "build_eos_plots",
    "calculate_eos",
    "describe_eos_plots",
    "diagnose_eos",
    "fit_eos",
    "open_eos_archive",
]
