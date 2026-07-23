# -*- coding: utf-8 -*-

"""Internal seismic workflow facade used by :mod:`quantas.api.seismic`."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np

from quantas.core.events import Observer
from quantas.core.physics.seismic import ElasticMedium
from quantas.io.path import ensure_suffix
from quantas.models import (
    ModuleContract,
    PlotCollection,
    ReportTable,
    ResultData,
    SphericalSummarySpec,
    mapping_table,
)
from quantas.modules.seismic.calculator import SeismicCalculator
from quantas.modules.seismic.io.export import SeismicHDF5Export, SeismicTableExport
from quantas.modules.seismic.io.reader import (
    SeismicInputFileReader,
    read_seismic_hdf5 as _read_seismic_hdf5,
)
from quantas.modules.seismic.models import SeismicInput, SeismicOptions, SeismicResult
from quantas.modules.seismic.plot import (
    SeismicPlotOptions,
    SeismicSurfaceOptions,
    build_seismic_plot_collection,
    build_seismic_summary_spec,
    build_seismic_surface_collection,
)
from quantas.modules.seismic.report import (
    build_seismic_report_tables,
    stiffness_matrix_table,
    tensor_rotation_metadata_table,
)


def read_seismic_input(filename: str | Path) -> SeismicInput:
    """Read a Quantas seismic text input file.

    Parameters
    ----------
    filename : str or Path
        Text file containing stiffness in GPa and density in kg m^-3.

    Returns
    -------
    SeismicInput
        Parsed seismic input.

    Raises
    ------
    ValueError
        If the file content is incomplete or invalid.
    """
    reader = SeismicInputFileReader(filename)
    if not reader.completed or reader.stiffness is None:
        raise ValueError(reader.error or "The seismic input could not be read.")
    return SeismicInput(
        jobname=reader.jobname,
        stiffness=reader.stiffness,
        density=reader.density,
        source=filename,
        raw=reader.raw,
    )


def read_seismic_hdf5(filename: str | Path) -> ResultData:
    """Read a complete native Quantas seismic HDF5 result.

    Parameters
    ----------
    filename : str or Path
        Native Quantas seismic HDF5 file.

    Returns
    -------
    ResultData
        Generic result containing a seismic payload.
    """
    return _read_seismic_hdf5(filename)


def write_seismic_hdf5(
    result: ResultData,
    filename: str | Path,
    *,
    report_text: str | None = None,
) -> Path:
    """Write a complete native Quantas seismic HDF5 result.

    Parameters
    ----------
    result : ResultData
        Complete Quantas seismic result.
    filename : str or Path
        Destination filename.
    report_text : str or None, optional
        Optional frontend-rendered report text.

    Returns
    -------
    Path
        Final output path with the ``.hdf5`` suffix.
    """
    output = ensure_suffix(filename, ".hdf5")
    SeismicHDF5Export().export(result, output, report_text=report_text)
    return output


def run_seismic(
    input_data: SeismicInput | ElasticMedium | str | Path,
    options: SeismicOptions | None = None,
    observer: Observer | None = None,
) -> ResultData:
    """Run a sampled acoustic-wave propagation calculation.

    Parameters
    ----------
    input_data : SeismicInput, ElasticTensor, str or Path
        Normalized seismic input, elastic medium, or a Quantas seismic text input
        file.
    options : SeismicOptions or None, optional
        Scientific sampling and numerical options.
    observer : Observer or None, optional
        Observer receiving workflow events.

    Returns
    -------
    ResultData
        Generic Quantas result containing a :class:`SeismicResult` payload.

    Raises
    ------
    TypeError
        If ``input_data`` has an unsupported type.
    ValueError
        If an elastic medium or text input is invalid.
    """
    normalized = normalize_seismic_input(input_data)
    return SeismicCalculator(
        seismic_input=normalized,
        options=options,
        observer=observer,
    ).execute()


def normalize_seismic_input(
    input_data: SeismicInput | ElasticMedium | str | Path,
) -> SeismicInput:
    """Return a normalized seismic input object.

    Parameters
    ----------
    input_data : SeismicInput, ElasticTensor, str or Path
        Seismic input, elastic medium, or text input file.

    Returns
    -------
    SeismicInput
        Normalized stiffness matrix and density.

    Raises
    ------
    TypeError
        If ``input_data`` has an unsupported type.
    ValueError
        If an elastic medium or file is invalid.
    """
    if isinstance(input_data, SeismicInput):
        return input_data
    if isinstance(input_data, (str, Path)):
        return read_seismic_input(input_data)
    if isinstance(input_data, ElasticMedium):
        return SeismicInput(
            jobname="Elastic medium",
            stiffness=input_data.elastic_tensor.stiffness.copy(),
            density=input_data.density,
        )
    raise TypeError(
        "input_data must be a SeismicInput, ElasticMedium, string, or Path."
    )


def write_seismic_csv(
    result: ResultData,
    filename: str | Path,
) -> Path:
    """Write sampled seismic fields as a named long-form CSV table.

    Parameters
    ----------
    result : ResultData
        Complete Quantas seismic result.
    filename : str or Path
        Destination CSV file.

    Returns
    -------
    Path
        Final output path with the ``.csv`` suffix.
    """
    output = ensure_suffix(filename, ".csv")
    SeismicTableExport().export(result, output)
    return output


def build_seismic_report(
    result: ResultData,
    *,
    level: Literal["standard", "extended", "debug"] = "extended",
) -> list[ReportTable]:
    """Build detailed neutral report tables from a seismic result.

    Parameters
    ----------
    result : ResultData
        Generic Quantas seismic result.
    level : {"standard", "extended", "debug"}, optional
        Scientific report detail.

    Returns
    -------
    list of ReportTable
        Ordered tables with extrema, anisotropy, derived wave properties, and
        numerical diagnostics.

    Raises
    ------
    ValueError
        If the result does not contain a seismic payload.
    """
    payload = result.results.get("seismic")
    if result.metadata.module != "seismic" or not isinstance(payload, SeismicResult):
        raise ValueError("ResultData does not contain a valid seismic result.")
    tables = build_seismic_report_tables(payload, level=level)
    insertion = 1
    if result.options:
        tables.insert(
            insertion, mapping_table("Selected numerical options", result.options)
        )
        insertion += 1
    frame = payload.metadata.get("tensor_frame", {})
    if bool(frame.get("transformed", False)):
        tables.insert(insertion, tensor_rotation_metadata_table(frame))
        insertion += 1
        if level in {"extended", "debug"} and result.input_data is not None:
            source_stiffness = result.input_data.data.get("stiffness")
            if source_stiffness is not None:
                tables.insert(
                    insertion,
                    stiffness_matrix_table(
                        np.asarray(source_stiffness, dtype=float),
                        title="Stiffness matrix before rotation / GPa",
                    ),
                )
                insertion += 1
        if level in {"extended", "debug"}:
            tables.insert(
                insertion,
                stiffness_matrix_table(
                    np.asarray(payload.stiffness, dtype=float),
                    title="Stiffness matrix after rotation / GPa",
                ),
            )
    return tables


def build_seismic_plots(
    result: ResultData,
    options: SeismicPlotOptions | None = None,
) -> PlotCollection:
    """Build neutral spherical-map specifications from a seismic result.

    Parameters
    ----------
    result : ResultData
        Generic Quantas seismic result.
    options : SeismicPlotOptions or None, optional
        Plot-property selection and spherical-map preparation options.

    Returns
    -------
    PlotCollection
        Frontend-neutral spherical-map specifications.

    Raises
    ------
    ValueError
        If the result does not contain a seismic payload.
    """
    payload = result.results.get("seismic")
    if result.metadata.module != "seismic" or not isinstance(payload, SeismicResult):
        raise ValueError("ResultData does not contain a valid seismic result.")
    return build_seismic_plot_collection(payload, options=options)


def build_seismic_surfaces(
    result: ResultData,
    options: SeismicSurfaceOptions | None = None,
) -> PlotCollection:
    """Build neutral three-dimensional acoustic surfaces.

    Parameters
    ----------
    result : ResultData
        Generic Quantas seismic result.
    options : SeismicSurfaceOptions or None, optional
        Surface type, acoustic modes, coloring, and polarization options.

    Returns
    -------
    PlotCollection
        Frontend-neutral phase, slowness, or group-velocity surfaces.

    Raises
    ------
    ValueError
        If the result does not contain a seismic payload.
    """
    payload = result.results.get("seismic")
    if result.metadata.module != "seismic" or not isinstance(payload, SeismicResult):
        raise ValueError("ResultData does not contain a valid seismic result.")
    return build_seismic_surface_collection(payload, options=options)


def build_seismic_summary(
    result: ResultData,
    options: SeismicPlotOptions | None = None,
) -> SphericalSummarySpec:
    """Build a compact six-panel seismic velocity summary.

    Parameters
    ----------
    result : ResultData
        Generic Quantas seismic result.
    options : SeismicPlotOptions or None, optional
        Map projection, colormap, extrema, and polarization options.

    Returns
    -------
    SphericalSummarySpec
        Frontend-neutral multi-panel summary.

    Raises
    ------
    ValueError
        If the result does not contain a seismic payload.
    """
    payload = result.results.get("seismic")
    if result.metadata.module != "seismic" or not isinstance(payload, SeismicResult):
        raise ValueError("ResultData does not contain a valid seismic result.")
    return build_seismic_summary_spec(payload, options=options)


MODULE_CONTRACT = ModuleContract(
    name="seismic",
    result_key="seismic",
    read_input=read_seismic_input,
    run=run_seismic,
    read_hdf5=read_seismic_hdf5,
    write_hdf5=write_seismic_hdf5,
    build_report=build_seismic_report,
    build_plots=build_seismic_plots,
)
