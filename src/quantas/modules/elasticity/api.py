# -*- coding: utf-8 -*-

"""Internal elasticity workflow facade used by :mod:`quantas.api.elasticity`."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np

from quantas.io.path import ensure_suffix
from quantas.core.events import Observer
from quantas.models import (
    ModuleContract,
    PlotCollection,
    PlotInventory,
    ReportTable,
    ResultData,
    mapping_table,
)
from quantas.modules.elasticity.calculator import ElasticityCalculator
from quantas.modules.elasticity.io.reader import (
    ElasticityInputFileReader,
    read_elasticity_hdf5 as _read_elasticity_hdf5,
)
from quantas.modules.elasticity.io.export import ElasticityHDF5Export
from quantas.core.physics.elasticity import ElasticSurfaceProperty
from quantas.modules.elasticity.models import (
    ElasticityInput,
    ElasticityOptions,
    ElasticityResult,
    ElasticitySurfaceOptions,
)
from quantas.modules.elasticity.plot import (
    ElasticityPlotProperty,
    describe_elasticity_plots,
    ElasticitySurfaceGeometry,
    build_elasticity_2d_plot_collection,
    build_elasticity_plot_collection,
    build_elasticity_surface_plot_collection,
)
from quantas.modules.elasticity.surface import resolve_elasticity_surfaces
from quantas.modules.elasticity.report import (
    averages_table,
    compliance_table,
    stability_table,
    source_stiffness_table,
    stiffness_table,
    tensor_rotation_metadata_table,
    variations_tables,
)


def read_elasticity_hdf5(filename: str | Path) -> ResultData:
    """Read a complete Quantas elasticity HDF5 result.

    Parameters
    ----------
    filename : str or Path
        Native Quantas elasticity HDF5 file.

    Returns
    -------
    ResultData
        Generic result containing an ``ElasticityResult`` payload.
    """
    return _read_elasticity_hdf5(filename)


def read_elasticity_input(filename: str | Path) -> ElasticityInput:
    """Read a Quantas elasticity input file.

    Parameters
    ----------
    filename : str or Path
        Input filename.

    Returns
    -------
    ElasticityInput
        Parsed input data.

    Raises
    ------
    ValueError
        If the file cannot be read.
    """
    reader = ElasticityInputFileReader(filename)
    if not reader.completed:
        raise ValueError(reader.error)
    return ElasticityInput(
        jobname=reader.jobname,
        stiffness=reader.stiffness,
        source=filename,
    )


def run_elasticity(
    input_data: ElasticityInput | str | Path,
    options: ElasticityOptions | None = None,
    observer: Observer | None = None,
) -> ResultData:
    """Run a second-order elasticity calculation.

    Parameters
    ----------
    input_data : ElasticityInput or str or Path
        Normalized input data or a Quantas elasticity input file.
    options : ElasticityOptions or None, optional
        Scientific workflow options.
    observer : Observer or None, optional
        Observer receiving Quantas events.

    Returns
    -------
    ResultData
        Generic Quantas result data.
    """
    elasticity_input = (
        read_elasticity_input(input_data)
        if isinstance(input_data, (str, Path))
        else input_data
    )
    calculator = ElasticityCalculator(
        elasticity_input=elasticity_input,
        options=options,
        observer=observer,
    )
    return calculator.execute()


def write_elasticity_hdf5(
    result: ResultData,
    filename: str | Path,
    *,
    report_text: str | None = None,
) -> Path:
    """Write a complete native Quantas elasticity result.

    Parameters
    ----------
    result : ResultData
        Complete Quantas elasticity result.
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
    ElasticityHDF5Export().export(result, output, report_text=report_text)
    return output


def build_elasticity_report(result: ResultData) -> list[ReportTable]:
    """Build frontend-neutral report tables from a complete elasticity result.

    Parameters
    ----------
    result : ResultData
        Generic Quantas result containing an :class:`ElasticityResult` payload.

    Returns
    -------
    list of ReportTable
        Ordered input, option, tensor, average, stability and directional-extrema
        tables. Source- and analysis-frame tensor tables are both included when a
        component transformation was applied.

    Raises
    ------
    ValueError
        If ``result`` is not an elasticity result.
    """
    payload = result.results.get("elasticity")
    if result.metadata.module != "elasticity" or not isinstance(
        payload, ElasticityResult
    ):
        raise ValueError("ResultData does not contain a valid elasticity result.")
    tables: list[ReportTable] = []
    if result.input_data is not None:
        tables.append(mapping_table("Input data", result.input_data.data))
    tables.append(mapping_table("Selected options", result.options))
    frame = payload.metadata.get("tensor_frame", {})
    transformed = bool(frame.get("transformed", False))
    if transformed and result.input_data is not None:
        source_stiffness = result.input_data.data.get("stiffness")
        if source_stiffness is not None:
            tables.append(
                source_stiffness_table(np.asarray(source_stiffness, dtype=float))
            )
        tables.append(tensor_rotation_metadata_table(frame))
        tables.append(
            stiffness_table(
                payload,
                title="Stiffness matrix after rotation (GPa)",
            )
        )
    else:
        tables.append(stiffness_table(payload))
    tables.extend(
        [
            compliance_table(payload),
            averages_table(payload),
            stability_table(payload),
        ]
    )
    if payload.variations:
        tables.extend(variations_tables(payload))
    return tables


def describe_elasticity_plot_inventory(result: ResultData) -> PlotInventory:
    """Describe plot representations available for an elasticity result.

    Parameters
    ----------
    result : ResultData
        Generic Quantas result containing an :class:`ElasticityResult` payload.

    Returns
    -------
    PlotInventory
        Result-aware inventory of available two- and three-dimensional elastic
        properties and plotting contexts.

    Raises
    ------
    ValueError
        If ``result`` is not an elasticity result.
    """
    payload = result.results.get("elasticity")
    if result.metadata.module != "elasticity" or not isinstance(
        payload, ElasticityResult
    ):
        raise ValueError("ResultData does not contain a valid elasticity result.")
    return describe_elasticity_plots(payload)


def build_elasticity_plots(result: ResultData) -> PlotCollection:
    """Build the default frontend-neutral elasticity plots.

    Parameters
    ----------
    result : ResultData
        Generic Quantas result containing persisted elasticity data.

    Returns
    -------
    PlotCollection
        Polar plot specifications built from persisted principal-plane data.

    Raises
    ------
    ValueError
        If ``result`` is not an elasticity result.
    """
    payload = result.results.get("elasticity")
    if result.metadata.module != "elasticity" or not isinstance(
        payload, ElasticityResult
    ):
        raise ValueError("ResultData does not contain a valid elasticity result.")
    return build_elasticity_plot_collection(payload)


def build_elasticity_2d_plots(
    result: ResultData,
    properties: tuple[ElasticityPlotProperty, ...] | None = None,
) -> PlotCollection:
    """Build selected principal-plane elasticity plot specifications.

    Parameters
    ----------
    result : ResultData
        Generic Quantas result containing persisted two-dimensional elasticity
        data.
    properties : tuple of str or None, optional
        Requested plot-property keys. ``None`` selects all available properties.

    Returns
    -------
    PlotCollection
        Frontend-neutral three-plane polar specifications and any preparation
        warnings.

    Raises
    ------
    ValueError
        If ``result`` is not an elasticity result.
    """
    payload = result.results.get("elasticity")
    if result.metadata.module != "elasticity" or not isinstance(
        payload, ElasticityResult
    ):
        raise ValueError("ResultData does not contain a valid elasticity result.")
    return build_elasticity_2d_plot_collection(payload, properties=properties)


def build_elasticity_3d_plots(
    result: ResultData,
    options: ElasticitySurfaceOptions | None = None,
    *,
    properties: tuple[ElasticSurfaceProperty, ...] | None = None,
    geometry: ElasticitySurfaceGeometry = "physical",
    color_mode: Literal["solid", "property"] = "property",
    colormap: str = "viridis",
    show_mesh: bool = False,
    mesh_color: str = "black",
    mesh_line_width: float = 0.5,
) -> PlotCollection:
    """Build frontend-neutral three-dimensional elasticity surfaces.

    Persisted surfaces are reused when ``options`` is ``None``. Supplying explicit
    sampling options requests a fresh in-memory calculation and does not modify the
    source result.

    Parameters
    ----------
    result : ResultData
        Generic Quantas elasticity result providing stiffness and, when available,
        persisted surfaces.
    options : ElasticitySurfaceOptions or None, optional
        Sampling options for a transient surface calculation.
    properties : tuple of ElasticSurfaceProperty or None, optional
        Selected directional property families.
    geometry : {"physical", "unit_sphere", "normalized"}, optional
        Radial representation. ``"normalized"`` is accepted as a compatibility
        alias of ``"physical"``.
    color_mode : {"solid", "property"}, optional
        Whether surfaces use semantic solid colors or property-value coloring.
    colormap : str, optional
        Matplotlib-compatible colormap name used for property coloring.
    show_mesh : bool, optional
        Whether a surface mesh overlay is requested.
    mesh_color : str, optional
        Mesh-line color passed to frontend renderers.
    mesh_line_width : float, optional
        Mesh-line width passed to frontend renderers.

    Returns
    -------
    PlotCollection
        One frontend-neutral surface specification for each selected branch.
    """
    surfaces = resolve_elasticity_surfaces(
        result,
        options=options,
        properties=properties,
    )
    return build_elasticity_surface_plot_collection(
        surfaces,
        geometry=geometry,
        color_mode=color_mode,
        colormap=colormap,
        show_mesh=show_mesh,
        mesh_color=mesh_color,
        mesh_line_width=mesh_line_width,
    )


MODULE_CONTRACT = ModuleContract(
    name="elasticity",
    result_key="elasticity",
    read_input=read_elasticity_input,
    run=run_elasticity,
    read_hdf5=read_elasticity_hdf5,
    write_hdf5=write_elasticity_hdf5,
    build_report=build_elasticity_report,
    build_plots=build_elasticity_plots,
)
