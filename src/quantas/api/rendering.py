# -*- coding: utf-8 -*-

"""Public rendering entry points for frontend-neutral Quantas output.

Scientific workflows return neutral tables and plot specifications.  This
namespace provides a small supported bridge to deterministic plain text and to
the bundled Matplotlib backend without exposing concrete renderer classes as
part of the public API.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from quantas.renderers.tables.text import (
    render_table as _render_table,
    render_tables as _render_tables,
)

from .common import ReportTable, _public_dir
from .plotting import PlotCollection


@dataclass(frozen=True, slots=True)
class RenderedPlot:
    """Describe one concrete plot produced by a public rendering call.

    Parameters
    ----------
    key : str
        Stable machine-readable plot identifier.
    kind : str
        Neutral plot family such as ``"curve"``, ``"surface"``, or
        ``"spherical_map"``.
    figure : object
        Backend figure object.  The bundled backend returns a Matplotlib
        figure without making Matplotlib renderer classes public contracts.
    path : Path or None
        Saved file path when ``output_dir`` was supplied.
    """

    key: str
    kind: str
    figure: object
    path: Path | None = None


@dataclass(frozen=True, slots=True)
class PlotRenderResult:
    """Concrete artifacts and warnings returned by :func:`render_plots`.

    Parameters
    ----------
    artifacts : tuple of RenderedPlot
        Rendered figures in the same order as the neutral collection.
    warnings : tuple of str
        Preparation warnings propagated by the scientific plot builder.
    """

    artifacts: tuple[RenderedPlot, ...]
    warnings: tuple[str, ...] = ()

    @property
    def figures(self) -> tuple[object, ...]:
        """Return all backend figure objects.

        Returns
        -------
        tuple of object
            Figure objects in rendering order.
        """
        return tuple(artifact.figure for artifact in self.artifacts)

    @property
    def paths(self) -> tuple[Path, ...]:
        """Return paths of all figures written to disk.

        Returns
        -------
        tuple of Path
            Saved paths, excluding artifacts that were not written.
        """
        return tuple(
            artifact.path for artifact in self.artifacts if artifact.path is not None
        )


def render_table(table: ReportTable) -> str:
    """Render one neutral report table as deterministic plain text.

    Parameters
    ----------
    table : ReportTable
        Frontend-neutral table returned by a Quantas report builder.

    Returns
    -------
    str
        Aligned plain text without ANSI sequences or box-drawing characters.

    Raises
    ------
    TypeError
        If ``table`` is not a :class:`~quantas.api.common.ReportTable`.
    """
    if not isinstance(table, ReportTable):
        raise TypeError("table must be a ReportTable object")
    return _render_table(table)


def render_tables(tables: Sequence[ReportTable]) -> str:
    """Render an ordered collection of neutral report tables.

    Parameters
    ----------
    tables : sequence of ReportTable
        Tables returned by a Quantas report builder.

    Returns
    -------
    str
        Concatenated deterministic plain-text report.

    Raises
    ------
    TypeError
        If any item is not a :class:`~quantas.api.common.ReportTable`.
    """
    values = tuple(tables)
    if not all(isinstance(table, ReportTable) for table in values):
        raise TypeError("tables must contain only ReportTable objects")
    return _render_tables(values)


def render_plots(
    collection: PlotCollection,
    *,
    output_dir: str | Path | None = None,
    filename_prefix: str = "",
    image_format: str = "png",
    preset: Literal["screen", "publication", "monochrome"] = "screen",
    dpi: int | None = None,
    figure_size: tuple[float, float] | None = None,
    show: bool = False,
    close: bool = False,
    axis_label_mode: Literal["cartesian", "crystal"] = "cartesian",
    annotate_extrema: bool = True,
    show_title: bool = False,
) -> PlotRenderResult:
    """Render neutral Quantas plot specifications with the bundled backend.

    Parameters
    ----------
    collection : PlotCollection
        Frontend-neutral plot collection returned by ``build_plots`` or a
        module-specific plot builder.
    output_dir : str, Path, or None, optional
        Directory in which figures are written.  Figures are only returned in
        memory when omitted.
    filename_prefix : str, optional
        Prefix prepended to generated filenames.
    image_format : str, optional
        Matplotlib-supported output format such as ``"png"``, ``"pdf"``, or
        ``"svg"``.
    preset : {"screen", "publication", "monochrome"}, optional
        Standard Quantas figure preset.
    dpi : int or None, optional
        Explicit raster resolution.  ``None`` preserves the preset value.
    figure_size : tuple of float or None, optional
        Explicit width and height in inches.  ``None`` preserves the preset or
        plot-specific dimensions.
    show : bool, optional
        Display figures through the active Matplotlib backend.
    close : bool, optional
        Close figures after rendering.  Use this for batch file generation;
        closed figures remain listed in the returned artifact metadata.
    axis_label_mode : {"cartesian", "crystal"}, optional
        Labels used for directional axes.
    annotate_extrema : bool, optional
        Annotate extrema when supported by the selected plot type.
    show_title : bool, optional
        Show titles on three-dimensional surfaces.

    Returns
    -------
    PlotRenderResult
        Public artifact metadata, backend figures, saved paths, and propagated
        warnings.

    Raises
    ------
    TypeError
        If ``collection`` is not a
        :class:`~quantas.api.plotting.PlotCollection`.
    ValueError
        If rendering options or a plot specification are invalid.
    """
    if not isinstance(collection, PlotCollection):
        raise TypeError("collection must be a PlotCollection object")
    try:
        from quantas.renderers.plots.matplotlib import (
            MatplotlibOptions,
            render_plot_collection as _render_plot_collection,
        )
    except ModuleNotFoundError as exc:
        if exc.name == "matplotlib":
            raise ModuleNotFoundError(
                "Matplotlib rendering requires the optional 'plot' extra: "
                "python -m pip install 'quantas[plot]'"
            ) from exc
        raise
    options = MatplotlibOptions.from_preset(
        preset,
        output_dir=output_dir,
        filename_prefix=filename_prefix,
        image_format=image_format,
        dpi=dpi,
        figure_size=figure_size,
        show=show,
        close=close,
        axis_label_mode=axis_label_mode,
        annotate_extrema=annotate_extrema,
        show_title=show_title,
    )
    rendered = _render_plot_collection(collection, options)
    return PlotRenderResult(
        artifacts=tuple(
            RenderedPlot(
                key=artifact.property_name,
                kind=artifact.kind,
                figure=artifact.figure,
                path=artifact.path,
            )
            for artifact in rendered.artifacts
        ),
        warnings=tuple(rendered.warnings),
    )


def __dir__() -> list[str]:
    """Return only names in the supported public contract."""
    return _public_dir(__all__)


__all__ = [
    "PlotRenderResult",
    "RenderedPlot",
    "render_plots",
    "render_table",
    "render_tables",
]
