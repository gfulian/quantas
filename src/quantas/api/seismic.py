# -*- coding: utf-8 -*-

"""Public API for directional seismic-wave calculations."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from quantas.core.events import Observer
from quantas.core.physics.seismic import ElasticMedium
from quantas.models import PlotCollection, ReportTable, ResultData, SphericalSummarySpec
from quantas.modules.seismic.api import (
    build_seismic_plots as _build_plots,
    build_seismic_report as _build_report,
    build_seismic_summary as _build_summary,
    build_seismic_surfaces as _build_surfaces,
    normalize_seismic_input as _normalize_input,
    read_seismic_hdf5 as _read_result,
    read_seismic_input as _read_input,
    run_seismic as _run,
    write_seismic_csv as _write_csv,
    write_seismic_hdf5 as _write_result,
)
from quantas.modules.seismic.models import (
    SeismicInput as Input,
    SeismicOptions as Options,
    SeismicResult as Result,
)
from quantas.modules.seismic.plot import (
    SeismicPlotOptions as PlotOptions,
    SeismicSurfaceOptions as SurfaceOptions,
)
from quantas.modules.seismic.plot.spec import (
    SurfaceGeometry,
    SurfaceType,
)

from .common import _public_dir, get_result_payload


def read_input(source: str | Path) -> Input:
    """Read one Quantas seismic input file.

    Parameters
    ----------
    source : str or Path
        Text input containing stiffness and density.

    Returns
    -------
    Input
        Validated passive seismic input contract.

    Raises
    ------
    ValueError
        If the file is malformed or lacks valid density/stiffness data.
    """
    return _read_input(source)


def normalize_input(source: Input | ElasticMedium | str | Path) -> Input:
    """Return a normalized seismic input contract.

    Parameters
    ----------
    source : Input, ElasticMedium, str, or Path
        Existing input, physical elastic medium, or text input path.

    Returns
    -------
    Input
        Validated seismic input suitable for :func:`run`.

    Raises
    ------
    TypeError
        If the source type is unsupported.
    ValueError
        If physical input data are invalid.
    """
    return _normalize_input(source)


def run(
    input_data: Input | ElasticMedium | str | Path,
    options: Options | None = None,
    observer: Observer | None = None,
) -> ResultData:
    """Run a directional seismic-wave workflow.

    Parameters
    ----------
    input_data : Input, ElasticMedium, str, or Path
        Seismic input, physical medium, or text input path.
    options : Options or None, optional
        Sphere sampling, degeneracy, tracking, and field controls.
    observer : Observer or None, optional
        Frontend-neutral event observer.

    Returns
    -------
    ResultData
        Complete result envelope containing a seismic payload.

    Raises
    ------
    ValueError
        If stiffness, density, or numerical options are invalid.
    """
    return _run(input_data, options=options, observer=observer)


def get_result(result: ResultData) -> Result:
    """Return the typed seismic payload from a result envelope.

    Parameters
    ----------
    result : ResultData
        Complete Quantas result envelope.

    Returns
    -------
    Result
        Module-specific seismic result.

    Raises
    ------
    ValueError
        If the envelope is not a valid seismic result.
    """
    return get_result_payload(
        result,
        module="seismic",
        key="seismic",
        expected_type=Result,
    )


def read_result(source: str | Path) -> ResultData:
    """Read a native Quantas seismic HDF5 result.

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
        If the file is not a supported seismic result.
    """
    return _read_result(source)


def write_result(
    result: ResultData,
    destination: str | Path,
    *,
    report_text: str | None = None,
) -> Path:
    """Write a native Quantas seismic HDF5 result.

    Parameters
    ----------
    result : ResultData
        Complete seismic result envelope.
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


def write_csv(result: ResultData, destination: str | Path) -> Path:
    """Write sampled seismic data in CSV form.

    Parameters
    ----------
    result : ResultData
        Complete seismic result envelope.
    destination : str or Path
        Destination CSV path.

    Returns
    -------
    Path
        Final CSV path.

    Raises
    ------
    ValueError
        If the result envelope is invalid.
    """
    get_result(result)
    return _write_csv(result, destination)


def build_report(
    result: ResultData,
    *,
    level: Literal["standard", "extended", "debug"] = "extended",
) -> list[ReportTable]:
    """Build frontend-neutral seismic report tables.

    Parameters
    ----------
    result : ResultData
        Complete seismic result envelope.
    level : {"standard", "extended", "debug"}, optional
        Scientific report detail.

    Returns
    -------
    list of ReportTable
        Ordered extrema, anisotropy, wave-property, and diagnostic tables.

    Raises
    ------
    ValueError
        If the result envelope is invalid.
    """
    get_result(result)
    return _build_report(result, level=level)


def build_plots(
    result: ResultData,
    options: PlotOptions | None = None,
) -> PlotCollection:
    """Build frontend-neutral seismic plots.

    Parameters
    ----------
    result : ResultData
        Complete seismic result envelope.
    options : PlotOptions or None, optional
        Selected quantities, sections, and plot metadata.

    Returns
    -------
    PlotCollection
        Neutral two-dimensional and summary plot specifications.

    Raises
    ------
    ValueError
        If requested quantities are unavailable.
    """
    get_result(result)
    return _build_plots(result, options=options)


def build_surfaces(
    result: ResultData,
    options: SurfaceOptions | None = None,
) -> PlotCollection:
    """Build frontend-neutral seismic surface plots.

    Parameters
    ----------
    result : ResultData
        Complete seismic result envelope.
    options : SurfaceOptions or None, optional
        Surface geometry, quantities, color fields, and mesh controls.

    Returns
    -------
    PlotCollection
        Neutral three-dimensional surface specifications.

    Raises
    ------
    ValueError
        If requested fields are unavailable or incompatible.
    """
    get_result(result)
    return _build_surfaces(result, options=options)


def build_summary(
    result: ResultData,
    options: PlotOptions | None = None,
) -> SphericalSummarySpec:
    """Build a frontend-neutral summary of extrema and anisotropy.

    Parameters
    ----------
    result : ResultData
        Complete seismic result envelope.
    options : PlotOptions or None, optional
        Quantity and branch selections used in the summary.

    Returns
    -------
    SphericalSummarySpec
        Structured extrema, directions, and anisotropy metadata.

    Raises
    ------
    ValueError
        If required spherical fields are unavailable.
    """
    get_result(result)
    return _build_summary(result, options=options)


def __dir__() -> list[str]:
    """Return only names in the supported public contract."""
    return _public_dir(__all__)


__all__ = [
    "ElasticMedium",
    "Input",
    "Options",
    "PlotOptions",
    "Result",
    "SurfaceGeometry",
    "SurfaceOptions",
    "SurfaceType",
    "build_plots",
    "build_report",
    "build_summary",
    "build_surfaces",
    "get_result",
    "normalize_input",
    "read_input",
    "read_result",
    "run",
    "write_csv",
    "write_result",
]
