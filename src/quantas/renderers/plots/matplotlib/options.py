# -*- coding: utf-8 -*-

"""Public data containers for the Matplotlib plotting backend."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

from quantas.models.plot import (
    ContourPlotSpec,
    LinePlotSpec,
    PanelPlotSpec,
    PlotSpec,
    PolarPlotSpec,
    SphericalMapSpec,
    SphericalSummarySpec,
    SurfacePlotSpec,
)

from quantas.renderers.plots.presets import FigurePreset, figure_preset


@dataclass(slots=True)
class MatplotlibOptions:
    """Options controlling Matplotlib figure rendering and file output.

    Parameters
    ----------
    output_dir : str or Path or None, optional
        Directory used to save figures. If ``None``, figures are returned but
        not written.
    filename_prefix : str, optional
        Prefix prepended to each specification filename stem.
    image_format : str, optional
        Matplotlib-supported image format.
    dpi : int, optional
        Resolution used when saving raster formats.
    figure_size : tuple of float or None, optional
        Requested figure width and height in inches. ``None`` preserves the
        backend or plot-specific default.
    tight_layout : bool, optional
        Whether the dispatcher applies ``Figure.tight_layout`` after the
        renderer has created, resized, and styled the figure.
    axis_label_font_size : float or None, optional
        Font size for Cartesian axis labels. ``None`` preserves the renderer
        or Matplotlib default.
    legend_font_size : float or None, optional
        Font size for legend entries and legend titles. ``None`` preserves the
        renderer or Matplotlib default.
    title_font_size : float or None, optional
        Font size for axis and figure titles. ``None`` preserves the renderer
        or Matplotlib default.
    tick_label_font_size : float or None, optional
        Font size for numeric and categorical tick labels. ``None`` preserves
        the renderer or Matplotlib default.
    show : bool, optional
        Whether to call :func:`matplotlib.pyplot.show` after rendering.
    close : bool, optional
        Whether to close figures after rendering.
    savefig_kwargs : dict, optional
        Additional keyword arguments forwarded to ``Figure.savefig``.
    axis_label_mode : {"cartesian", "crystal"}, optional
        Direction labels used on spherical maps and Cartesian surfaces.
        Crystal labels assume the tensor frame follows ``[100]``, ``[010]``
        and ``[001]``.
    polarization_color : str or None, optional
        Optional renderer-level override for polarization-axis color.
    polarization_line_width : float or None, optional
        Optional renderer-level override for polarization-axis line width.
    polarization_scale : float or None, optional
        Optional renderer-level override for polarization-axis length.
    annotate_extrema : bool, optional
        Whether minimum and maximum markers receive text annotations.
    show_title : bool, optional
        Whether three-dimensional surface figures display their title.
    preset : {"screen", "publication", "monochrome"}, optional
        Identifier of the resolved preset. Use :meth:`from_preset` to apply
        its defaults.
    monochrome : bool, optional
        Whether the rendered artifact is converted to a grayscale-safe style.
    """

    output_dir: str | Path | None = None
    filename_prefix: str = ""
    image_format: str = "png"
    dpi: int = 150
    figure_size: tuple[float, float] | None = None
    tight_layout: bool = False
    axis_label_font_size: float | None = None
    legend_font_size: float | None = None
    title_font_size: float | None = None
    tick_label_font_size: float | None = None
    show: bool = False
    close: bool = False
    savefig_kwargs: dict[str, Any] = field(default_factory=dict)
    axis_label_mode: Literal["cartesian", "crystal"] = "cartesian"
    polarization_color: str | None = None
    polarization_line_width: float | None = None
    polarization_scale: float | None = None
    annotate_extrema: bool = True
    show_title: bool = False
    preset: FigurePreset = "screen"
    monochrome: bool = False

    @classmethod
    def from_preset(
        cls,
        preset: str = "screen",
        **overrides: Any,
    ) -> MatplotlibOptions:
        """Build renderer options from one standard preset.

        Parameters
        ----------
        preset : {"screen", "publication", "monochrome"}, optional
            Standard Quantas figure preset.
        **overrides
            Explicit :class:`MatplotlibOptions` field values. Values set to
            ``None`` leave the preset default unchanged.

        Returns
        -------
        MatplotlibOptions
            Validated renderer options.
        """

        resolved = figure_preset(preset)
        values: dict[str, Any] = {
            "dpi": resolved.dpi,
            "figure_size": resolved.figure_size,
            "tight_layout": resolved.tight_layout,
            "axis_label_font_size": resolved.axis_label_font_size,
            "legend_font_size": resolved.legend_font_size,
            "title_font_size": resolved.title_font_size,
            "tick_label_font_size": resolved.tick_label_font_size,
            "savefig_kwargs": dict(resolved.savefig_kwargs),
            "preset": resolved.name,
            "monochrome": resolved.monochrome,
        }
        explicit = {key: value for key, value in overrides.items() if value is not None}
        savefig_overrides = explicit.pop("savefig_kwargs", None)
        values.update(explicit)
        if savefig_overrides is not None:
            values["savefig_kwargs"].update(dict(savefig_overrides))
        return cls(**values)

    def __post_init__(self) -> None:
        """Validate renderer options.

        Returns
        -------
        None
            The instance is validated in place.

        Raises
        ------
        ValueError
            If a numeric option is outside its accepted range.
        """
        figure_preset(self.preset)
        if self.dpi < 1:
            raise ValueError("dpi must be positive.")
        if self.figure_size is not None:
            if len(self.figure_size) != 2:
                raise ValueError("figure_size must contain width and height.")
            width, height = (float(value) for value in self.figure_size)
            if width <= 0.0 or height <= 0.0:
                raise ValueError("figure_size values must be positive.")
            self.figure_size = (width, height)
        for name in (
            "axis_label_font_size",
            "legend_font_size",
            "title_font_size",
            "tick_label_font_size",
        ):
            value = getattr(self, name)
            if value is None:
                continue
            size = float(value)
            if size <= 0.0:
                raise ValueError(f"{name} must be positive.")
            setattr(self, name, size)
        if (
            self.polarization_line_width is not None
            and self.polarization_line_width <= 0
        ):
            raise ValueError("polarization_line_width must be positive.")
        if self.polarization_scale is not None and self.polarization_scale <= 0:
            raise ValueError("polarization_scale must be positive.")


@dataclass(slots=True)
class PlotArtifact:
    """Concrete Matplotlib artifact generated from one neutral plot spec.

    Parameters
    ----------
    spec : PlotSpec
        Source neutral specification.
    figure : object
        Matplotlib figure instance.
    path : Path or None, optional
        Saved file path when file output was requested.
    """

    spec: PlotSpec
    figure: object
    path: Path | None = None

    @property
    def property_name(self) -> str:
        """Return the stable plot key.

        Returns
        -------
        str
            Machine-readable key stored in the neutral specification.
        """
        return self.spec.key

    @property
    def kind(self) -> str:
        """Return the neutral plot kind represented by the artifact.

        Returns
        -------
        str
            Stable renderer kind.
        """
        if isinstance(self.spec, LinePlotSpec):
            return "curve"
        if isinstance(self.spec, ContourPlotSpec):
            return "contour"
        if isinstance(self.spec, SurfacePlotSpec):
            return "surface"
        if isinstance(self.spec, SphericalMapSpec):
            return "spherical_map"
        if isinstance(self.spec, SphericalSummarySpec):
            return "spherical_summary"
        if isinstance(self.spec, PolarPlotSpec):
            return "polar"
        if isinstance(self.spec, PanelPlotSpec):
            return "panel"
        return "plot"


@dataclass(slots=True)
class PlotRenderResult:
    """Rendered artifacts together with preparation warnings.

    Parameters
    ----------
    artifacts : list of PlotArtifact, optional
        Ordered concrete figure artifacts.
    warnings : list of str, optional
        Non-fatal warnings propagated from the neutral plot collection.
    """

    artifacts: list[PlotArtifact] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


__all__ = ["MatplotlibOptions", "PlotArtifact", "PlotRenderResult"]
