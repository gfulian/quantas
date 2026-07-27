# -*- coding: utf-8 -*-

"""Public API for second-order elasticity calculations."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from quantas.core.events import Observer
from quantas.core.geometry import TensorRotation, TensorRotationKind
from quantas.core.physics.elasticity import ElasticSurfaceProperty as SurfaceProperty
from quantas.io.path import ensure_suffix
from quantas.modules.elasticity.api import (
    build_elasticity_2d_plots as _build_2d_plots,
    describe_elasticity_plot_inventory as _describe_plots,
    build_elasticity_3d_plots as _build_3d_plots,
    build_elasticity_plots as _build_plots,
    build_elasticity_report as _build_report,
    read_elasticity_hdf5 as _read_result,
    read_elasticity_input as _read_input,
    run_elasticity as _run,
    write_elasticity_hdf5 as _write_result,
)
from quantas.modules.elasticity.models import (
    ElasticityInput as Input,
    ElasticityOptions as Options,
    ElasticityResult as Result,
    ElasticitySurfaceOptions as SurfaceOptions,
)
from quantas.modules.elasticity.io import ElasticityInputCreator, ElasticityTableExport
from quantas.modules.elasticity.plot import (
    ElasticityPlotProperty as PlotProperty,
    ElasticitySurfaceGeometry as SurfaceGeometry,
)

from .common import ReportTable, ResultData, _public_dir, get_result_payload
from .plotting import PlotCollection, PlotInventory


InputInterface = Literal["crystal", "vasp"]


def create_input(
    source: str | Path,
    destination: str | Path,
    *,
    interface: InputInterface = "crystal",
    jobname: str = "Unknown",
) -> Path:
    """Create a Quantas elasticity input from an external-code output.

    Parameters
    ----------
    source : str or Path
        CRYSTAL or VASP output containing an elastic stiffness tensor and,
        when available, density metadata.
    destination : str or Path
        Destination text path. The ``.dat`` suffix is applied when absent.
    interface : {"crystal", "vasp"}, optional
        External-code reader used to interpret ``source``.
    jobname : str, optional
        Human-readable title written to the generated input.

    Returns
    -------
    Path
        Written Quantas elasticity input path.

    Raises
    ------
    ValueError
        If the interface is unsupported or the source cannot be parsed.
    OSError
        If the destination cannot be written.
    """
    output = ensure_suffix(destination, ".dat")
    try:
        creator = ElasticityInputCreator(interface)
    except KeyError as exc:
        raise ValueError(f"Unsupported elasticity interface: {interface}") from exc
    completed, error = creator.read(source)
    if not completed:
        raise ValueError(error or "Unable to create elasticity input")
    creator.write(output, jobname=jobname)
    return output


def read_input(source: str | Path) -> Input:
    """Read one Quantas elasticity input file.

    Parameters
    ----------
    source : str or Path
        Text input containing the job name, elastic stiffness tensor, and
        optional density metadata.

    Returns
    -------
    Input
        Validated passive elasticity input contract.

    Raises
    ------
    ValueError
        If the file is malformed or scientifically inconsistent.
    """
    return _read_input(source)


def normalize_input(source: Input | str | Path) -> Input:
    """Return a normalized elasticity input contract.

    Parameters
    ----------
    source : Input, str, or Path
        Existing passive input contract or text input path.

    Returns
    -------
    Input
        Validated elasticity input suitable for :func:`run`.

    Raises
    ------
    TypeError
        If ``source`` is neither an :class:`Input` object nor a path.
    ValueError
        If a supplied file cannot be parsed.
    """
    if isinstance(source, Input):
        return source
    if isinstance(source, (str, Path)):
        return read_input(source)
    raise TypeError("source must be an elasticity Input object or path")


def run(
    input_data: Input | str | Path,
    options: Options | None = None,
    observer: Observer | None = None,
) -> ResultData:
    """Run a second-order elasticity workflow.

    Parameters
    ----------
    input_data : Input, str, or Path
        Elasticity input contract or text input path.
    options : Options or None, optional
        Scientific and numerical calculation options. Defaults are used when
        omitted.
    observer : Observer or None, optional
        Frontend-neutral event observer.

    Returns
    -------
    ResultData
        Complete result envelope containing an elasticity payload.

    Raises
    ------
    ValueError
        If the input or selected calculation is invalid.
    """
    return _run(input_data, options=options, observer=observer)


def get_result(result: ResultData) -> Result:
    """Return the typed elasticity payload from a result envelope.

    Parameters
    ----------
    result : ResultData
        Complete Quantas result envelope.

    Returns
    -------
    Result
        Module-specific elasticity result.

    Raises
    ------
    ValueError
        If the envelope belongs to another module or lacks a valid payload.
    """
    return get_result_payload(
        result,
        module="elasticity",
        key="elasticity",
        expected_type=Result,
    )


def read_result(source: str | Path) -> ResultData:
    """Read a native Quantas elasticity HDF5 result.

    Parameters
    ----------
    source : str or Path
        Native Quantas HDF5 file.

    Returns
    -------
    ResultData
        Restored result envelope.

    Raises
    ------
    ValueError
        If the file is not a supported elasticity result.
    """
    return _read_result(source)


def write_result(
    result: ResultData,
    destination: str | Path,
    *,
    report_text: str | None = None,
) -> Path:
    """Write a native Quantas elasticity HDF5 result.

    Parameters
    ----------
    result : ResultData
        Complete elasticity result envelope.
    destination : str or Path
        Destination path. The native HDF5 suffix is applied when required.
    report_text : str or None, optional
        Deterministic plain-text report to embed in diagnostics.

    Returns
    -------
    Path
        Final HDF5 path.

    Raises
    ------
    ValueError
        If ``result`` does not contain an elasticity payload.
    """
    get_result(result)
    return _write_result(result, destination, report_text=report_text)


def write_table(result: ResultData, destination: str | Path) -> Path:
    """Write principal-plane elasticity data as a neutral text table.

    Parameters
    ----------
    result : ResultData
        Complete elasticity result envelope.
    destination : str or Path
        Destination path. The ``.dat`` suffix is applied when absent.

    Returns
    -------
    Path
        Written table path.

    Raises
    ------
    ValueError
        If ``result`` is not a valid elasticity result.
    """
    get_result(result)
    output = ensure_suffix(destination, ".dat")
    ElasticityTableExport().export(result, output)
    return output


def build_report(result: ResultData) -> list[ReportTable]:
    """Build frontend-neutral elasticity report tables.

    Parameters
    ----------
    result : ResultData
        Complete elasticity result envelope.

    Returns
    -------
    list of ReportTable
        Ordered raw-value tables ready for terminal, text, CSV, or GUI
        rendering.

    Raises
    ------
    ValueError
        If the result envelope is invalid.
    """
    get_result(result)
    return _build_report(result)


def describe_plots(result: ResultData) -> PlotInventory:
    """Return result-aware elasticity plot properties and representations.

    Parameters
    ----------
    result : ResultData
        Complete elasticity result envelope.

    Returns
    -------
    PlotInventory
        Available directional properties, branch metadata, principal-plane
        context, and three-dimensional geometries.
    """
    get_result(result)
    return _describe_plots(result)


def build_plots(result: ResultData) -> PlotCollection:
    """Build the default frontend-neutral elasticity plots.

    Parameters
    ----------
    result : ResultData
        Complete elasticity result envelope.

    Returns
    -------
    PlotCollection
        Neutral plot specifications for properties present in the result.

    Raises
    ------
    ValueError
        If the result envelope is invalid.
    """
    get_result(result)
    return _build_plots(result)


def build_2d_plots(
    result: ResultData,
    properties: tuple[PlotProperty, ...] | None = None,
) -> PlotCollection:
    """Build selected two-dimensional elasticity plots.

    Parameters
    ----------
    result : ResultData
        Complete elasticity result envelope.
    properties : tuple of PlotProperty or None, optional
        Directional properties to include. Module defaults are used when
        omitted.

    Returns
    -------
    PlotCollection
        Neutral two-dimensional plot specifications.

    Raises
    ------
    ValueError
        If required directional data are unavailable.
    """
    get_result(result)
    return _build_2d_plots(result, properties=properties)


def build_3d_plots(
    result: ResultData,
    options: SurfaceOptions | None = None,
    *,
    properties: tuple[SurfaceProperty, ...] | None = None,
    geometry: SurfaceGeometry = "physical",
    color_mode: Literal["solid", "property"] = "property",
    colormap: str = "viridis",
    show_mesh: bool = False,
    mesh_color: str = "black",
    mesh_line_width: float = 0.5,
) -> PlotCollection:
    """Build selected three-dimensional elasticity surfaces.

    Parameters
    ----------
    result : ResultData
        Complete elasticity result envelope.
    options : SurfaceOptions or None, optional
        Directional sampling controls.
    properties : tuple of SurfaceProperty or None, optional
        Surface properties to include.
    geometry : SurfaceGeometry, optional
        Unit-sphere or physical-radius representation.
    color_mode : {"solid", "property"}, optional
        Surface coloring strategy.
    colormap : str, optional
        Renderer-neutral colormap identifier.
    show_mesh : bool, optional
        Include mesh-line metadata.
    mesh_color : str, optional
        Mesh-line color identifier.
    mesh_line_width : float, optional
        Mesh-line width metadata.

    Returns
    -------
    PlotCollection
        Neutral three-dimensional surface specifications.

    Raises
    ------
    ValueError
        If options are invalid or required surface data are unavailable.
    """
    get_result(result)
    return _build_3d_plots(
        result,
        options=options,
        properties=properties,
        geometry=geometry,
        color_mode=color_mode,
        colormap=colormap,
        show_mesh=show_mesh,
        mesh_color=mesh_color,
        mesh_line_width=mesh_line_width,
    )


def __dir__() -> list[str]:
    """Return only names in the supported public contract."""
    return _public_dir(__all__)


__all__ = [
    "Input",
    "InputInterface",
    "Options",
    "PlotProperty",
    "Result",
    "SurfaceGeometry",
    "SurfaceOptions",
    "SurfaceProperty",
    "TensorRotation",
    "TensorRotationKind",
    "build_2d_plots",
    "build_3d_plots",
    "build_plots",
    "build_report",
    "create_input",
    "describe_plots",
    "get_result",
    "normalize_input",
    "read_input",
    "read_result",
    "run",
    "write_result",
    "write_table",
]
