# -*- coding: utf-8 -*-

"""Frontend-neutral plotting models used by Quantas workflows.

The data containers in this module describe prepared plotting data without
binding it to Matplotlib, Dash, or another graphical frontend. Module
specific builders create these objects from scientific results, while renderer
packages convert them into concrete figures.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, TypeAlias

import numpy as np


LineStyle: TypeAlias = Literal["solid", "dashed", "dotted", "dashdot", "none"]
BandOrientation: TypeAlias = Literal["vertical", "horizontal"]
AxisOrientation: TypeAlias = Literal["x", "y"]
AxisLocation: TypeAlias = Literal["top", "bottom", "left", "right"]


@dataclass(slots=True)
class PlotAxis:
    """Description of one plot axis.

    Parameters
    ----------
    key : str
        Stable machine-readable name of the represented quantity.
    label : str
        Complete human-readable axis label, including symbols or units when
        required.
    unit : str or None, optional
        Physical unit represented by the axis.
    limits : tuple or None, optional
        Optional lower and upper axis limits. Either bound may be ``None``.
    metadata : dict, optional
        Additional frontend-neutral axis information.
    """

    key: str
    label: str
    unit: str | None = None
    limits: tuple[float | None, float | None] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class PlotSeriesStyle:
    """Frontend-neutral style hints for a plotted series.

    Parameters
    ----------
    color : str or None, optional
        Portable color specification. ``None`` delegates color selection to
        the renderer.
    line_style : {"solid", "dashed", "dotted", "dashdot"}, optional
        Semantic line style.
    line_width : float, optional
        Requested line width.
    marker : str or None, optional
        Optional portable marker name or symbol.
    marker_size : float or None, optional
        Marker size in renderer-independent typographic points.
    marker_edge_color : str or None, optional
        Portable marker-edge color. ``None`` delegates to the renderer.
    marker_edge_width : float or None, optional
        Marker-edge width in typographic points.
    alpha : float, optional
        Portable opacity hint between zero and one.
    errorbar_line_width : float, optional
        Requested line width for error bars.
    errorbar_capsize : float, optional
        Requested error-bar cap size in typographic points.
    """

    color: str | None = None
    line_style: LineStyle = "solid"
    line_width: float = 1.5
    marker: str | None = None
    marker_size: float | None = None
    marker_edge_color: str | None = None
    marker_edge_width: float | None = None
    alpha: float = 1.0
    errorbar_line_width: float = 1.0
    errorbar_capsize: float = 2.0


@dataclass(slots=True)
class PlotBandStyle:
    """Frontend-neutral style hints for a confidence or uncertainty band.

    Parameters
    ----------
    color : str or None, optional
        Portable fill color. ``None`` delegates selection to the renderer.
    alpha : float, optional
        Fill opacity between zero and one.
    edge_color : str or None, optional
        Optional portable edge color.
    edge_width : float, optional
        Width of the band boundary.
    line_style : LineStyle, optional
        Boundary line style.
    """

    color: str | None = None
    alpha: float = 0.2
    edge_color: str | None = None
    edge_width: float = 0.0
    line_style: LineStyle = "solid"


@dataclass(slots=True)
class PlotBand:
    """One symmetric or asymmetric interval around a Cartesian curve.

    Parameters
    ----------
    key : str
        Stable machine-readable band name.
    label : str
        Human-readable legend label.
    coordinates : ndarray
        Coordinates along the varying axis.
    lower, upper : ndarray
        Lower and upper bounds perpendicular to ``coordinates``.
    orientation : {"vertical", "horizontal"}, optional
        ``"vertical"`` means ``coordinates`` are x values and the bounds are
        y values. ``"horizontal"`` means ``coordinates`` are y values and
        the bounds are x values.
    style : PlotBandStyle, optional
        Portable display hints.
    metadata : dict, optional
        Additional frontend-neutral information.
    """

    key: str
    label: str
    coordinates: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    orientation: BandOrientation = "vertical"
    style: PlotBandStyle = field(default_factory=PlotBandStyle)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize aligned arrays and validate interval ordering."""
        self.coordinates = np.asarray(self.coordinates, dtype=np.float64).copy()
        self.lower = np.asarray(self.lower, dtype=np.float64).copy()
        self.upper = np.asarray(self.upper, dtype=np.float64).copy()
        if self.coordinates.ndim != 1:
            raise ValueError("band coordinates must be one-dimensional")
        if (
            self.lower.shape != self.coordinates.shape
            or self.upper.shape != self.coordinates.shape
        ):
            raise ValueError("band coordinates and bounds must have identical shapes")
        if np.any(self.lower > self.upper):
            raise ValueError("band lower bounds cannot exceed upper bounds")
        if self.orientation not in {"vertical", "horizontal"}:
            raise ValueError("invalid band orientation")
        self.metadata = dict(self.metadata)


@dataclass(slots=True)
class ColoredPathStyle:
    """Frontend-neutral style hints for a scalar-colored path.

    Parameters
    ----------
    colormap : str, optional
        Portable colormap name.
    line_style : LineStyle, optional
        Semantic path line style.
    line_width : float, optional
        Path width.
    marker : str or None, optional
        Optional marker symbol.
    marker_size : float or None, optional
        Marker size in typographic points.
    marker_edge_color : str or None, optional
        Portable marker-edge color. ``None`` delegates to the renderer.
    marker_edge_width : float or None, optional
        Marker-edge width in typographic points.
    alpha : float, optional
        Path opacity.
    show_colorbar : bool, optional
        Whether the mapped scalar should receive a colorbar.
    value_limits : tuple or None, optional
        Optional lower and upper scalar limits shared by line and markers.
    """

    colormap: str = "viridis"
    line_style: LineStyle = "solid"
    line_width: float = 1.8
    marker: str | None = None
    marker_size: float | None = None
    marker_edge_color: str | None = None
    marker_edge_width: float | None = None
    alpha: float = 1.0
    show_colorbar: bool = True
    value_limits: tuple[float, float] | None = None


@dataclass(slots=True)
class ColoredPathSeries:
    """Cartesian path whose color is controlled by a third scalar variable.

    Parameters
    ----------
    key : str
        Stable machine-readable path name.
    label : str
        Human-readable legend label.
    x, y : ndarray
        Aligned Cartesian path coordinates.
    values : ndarray
        Scalar values mapped to the colormap.
    value_axis : PlotAxis
        Description of the color-mapped scalar quantity.
    style : ColoredPathStyle, optional
        Portable display hints.
    metadata : dict, optional
        Additional frontend-neutral information.
    """

    key: str
    label: str
    x: np.ndarray
    y: np.ndarray
    values: np.ndarray
    value_axis: PlotAxis
    style: ColoredPathStyle = field(default_factory=ColoredPathStyle)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize and validate aligned path arrays."""
        self.x = np.asarray(self.x, dtype=np.float64).copy()
        self.y = np.asarray(self.y, dtype=np.float64).copy()
        self.values = np.asarray(self.values, dtype=np.float64).copy()
        if (
            self.x.ndim != 1
            or self.y.shape != self.x.shape
            or self.values.shape != self.x.shape
        ):
            raise ValueError(
                "colored path x, y, and values must be aligned one-dimensional arrays"
            )
        if self.style.value_limits is not None:
            lower, upper = (float(value) for value in self.style.value_limits)
            if not lower < upper:
                raise ValueError("colored-path value limits must be increasing")
            self.style.value_limits = (lower, upper)
        self.metadata = dict(self.metadata)


@dataclass(slots=True)
class SecondaryAxis:
    """Prepared secondary Cartesian axis with explicit tick mapping.

    Parameters
    ----------
    key : str
        Stable machine-readable axis name.
    label : str
        Complete human-readable axis label.
    orientation : {"x", "y"}
        Axis orientation.
    location : {"top", "bottom", "left", "right"}
        Side on which the secondary axis is displayed.
    positions : ndarray
        Tick positions expressed in primary-axis coordinates.
    labels : tuple of str
        Labels corresponding one-to-one with ``positions``.
    metadata : dict, optional
        Additional frontend-neutral information.
    """

    key: str
    label: str
    orientation: AxisOrientation
    location: AxisLocation
    positions: np.ndarray
    labels: tuple[str, ...]
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize ticks and validate orientation and location."""
        self.positions = np.asarray(self.positions, dtype=np.float64).copy()
        self.labels = tuple(str(value) for value in self.labels)
        if self.positions.ndim != 1 or self.positions.size != len(self.labels):
            raise ValueError("secondary-axis positions and labels must be aligned")
        valid = (self.orientation == "x" and self.location in {"top", "bottom"}) or (
            self.orientation == "y" and self.location in {"left", "right"}
        )
        if not valid:
            raise ValueError("secondary-axis orientation and location are incompatible")
        self.metadata = dict(self.metadata)


@dataclass(slots=True)
class PlotSpan:
    """One highlighted interval parallel to a Cartesian axis.

    Parameters
    ----------
    key : str
        Stable interval identifier.
    label : str
        Human-readable legend label.
    axis : {"x", "y"}
        Axis along which ``start`` and ``end`` are defined.
    start, end : float
        Interval bounds in data coordinates.
    color : str or None, optional
        Portable fill color.
    alpha : float, optional
        Fill opacity.
    hatch : str or None, optional
        Portable hatch hint.
    metadata : dict, optional
        Additional frontend-neutral information.
    """

    key: str
    label: str
    axis: AxisOrientation
    start: float
    end: float
    color: str | None = "0.85"
    alpha: float = 0.25
    hatch: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate interval bounds."""
        self.start = float(self.start)
        self.end = float(self.end)
        if (
            not np.isfinite(self.start)
            or not np.isfinite(self.end)
            or self.start > self.end
        ):
            raise ValueError("plot-span bounds must be finite and increasing")
        if self.axis not in {"x", "y"}:
            raise ValueError("invalid plot-span axis")
        self.metadata = dict(self.metadata)


@dataclass(slots=True)
class PlotMask:
    """Two-dimensional boolean overlay for contours and domain maps.

    Parameters
    ----------
    key : str
        Stable mask identifier.
    label : str
        Human-readable legend label.
    x, y : ndarray
        One-dimensional Cartesian coordinates.
    mask : ndarray
        Boolean values with shape ``(len(y), len(x))``.
    hatch : str, optional
        Portable hatch pattern.
    color : str or None, optional
        Optional portable face color.
    alpha : float, optional
        Face opacity.
    metadata : dict, optional
        Additional frontend-neutral information.
    """

    key: str
    label: str
    x: np.ndarray
    y: np.ndarray
    mask: np.ndarray
    hatch: str = "///"
    color: str | None = "none"
    alpha: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize coordinates and validate mask shape."""
        self.x = np.asarray(self.x, dtype=np.float64).copy()
        self.y = np.asarray(self.y, dtype=np.float64).copy()
        self.mask = np.asarray(self.mask, dtype=np.bool_).copy()
        if self.x.ndim != 1 or self.y.ndim != 1:
            raise ValueError("plot-mask coordinates must be one-dimensional")
        if self.mask.shape != (self.y.size, self.x.size):
            raise ValueError("plot mask must have shape (len(y), len(x))")
        self.metadata = dict(self.metadata)


@dataclass(slots=True)
class ScalarBackground:
    """Scalar field varying along one axis and painted behind a line plot.

    Parameters
    ----------
    key : str
        Stable background identifier.
    coordinates, values : ndarray
        Aligned one-dimensional data coordinates and mapped scalar values.
    value_axis : PlotAxis
        Description of the mapped scalar.
    axis : {"x", "y"}, optional
        Axis along which the scalar varies.
    colormap : str, optional
        Portable colormap name.
    alpha : float, optional
        Background opacity.
    show_colorbar : bool, optional
        Whether a colorbar should be displayed.
    value_limits : tuple or None, optional
        Optional mapped-value limits.
    metadata : dict, optional
        Additional frontend-neutral information.
    """

    key: str
    coordinates: np.ndarray
    values: np.ndarray
    value_axis: PlotAxis
    axis: AxisOrientation = "y"
    colormap: str = "viridis"
    alpha: float = 0.2
    show_colorbar: bool = False
    value_limits: tuple[float, float] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize and validate the one-dimensional scalar field."""
        self.coordinates = np.asarray(self.coordinates, dtype=np.float64).copy()
        self.values = np.asarray(self.values, dtype=np.float64).copy()
        if self.coordinates.ndim != 1 or self.values.shape != self.coordinates.shape:
            raise ValueError("background coordinates and values must be aligned")
        if self.axis not in {"x", "y"}:
            raise ValueError("invalid scalar-background axis")
        if self.value_limits is not None:
            lower, upper = (float(value) for value in self.value_limits)
            if not lower < upper:
                raise ValueError("background value limits must be increasing")
            self.value_limits = (lower, upper)
        self.metadata = dict(self.metadata)


@dataclass(slots=True)
class PlotSeries:
    """One numerical series in a neutral plot specification.

    Parameters
    ----------
    key : str
        Stable machine-readable series name.
    label : str
        Human-readable legend label.
    x : ndarray
        One-dimensional horizontal or angular coordinates.
    y : ndarray
        One-dimensional dependent values.
    x_error, y_error : ndarray or None, optional
        Optional symmetric one-sigma uncertainties associated with the
        horizontal and vertical coordinates.
    style : PlotSeriesStyle, optional
        Portable style hints.
    metadata : dict, optional
        Additional frontend-neutral series information.
    """

    key: str
    label: str
    x: np.ndarray
    y: np.ndarray
    x_error: np.ndarray | None = None
    y_error: np.ndarray | None = None
    style: PlotSeriesStyle = field(default_factory=PlotSeriesStyle)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class LinePlotSpec:
    """Neutral specification for a Cartesian line plot.

    Parameters
    ----------
    key : str
        Stable plot identifier.
    title : str
        Figure title.
    filename_stem : str
        Default filename stem used by file renderers.
    x_axis, y_axis : PlotAxis
        Axis descriptions.
    series : list of PlotSeries
        Ordered ordinary line and error-bar series.
    bands : list of PlotBand, optional
        Confidence or uncertainty intervals drawn around Cartesian curves.
    colored_paths : list of ColoredPathSeries, optional
        Paths whose color is controlled by a third scalar quantity.
    secondary_axes : list of SecondaryAxis, optional
        Prepared secondary tick mappings.
    spans : list of PlotSpan, optional
        Highlighted intervals parallel to either Cartesian axis.
    backgrounds : list of ScalarBackground, optional
        Scalar fields painted behind the primary line data.
    legend_title : str or None, optional
        Optional legend title.
    legend_columns : int, optional
        Preferred number of legend columns.
    show_legend : bool, optional
        Whether renderers should display a legend.
    grid : bool, optional
        Whether renderers should display a grid.
    invert_x_axis, invert_y_axis : bool, optional
        Whether renderers should invert the corresponding primary axis.
    metadata : dict, optional
        Additional frontend-neutral plot information.
    """

    key: str
    title: str
    filename_stem: str
    x_axis: PlotAxis
    y_axis: PlotAxis
    series: list[PlotSeries]
    bands: list[PlotBand] = field(default_factory=list)
    colored_paths: list[ColoredPathSeries] = field(default_factory=list)
    secondary_axes: list[SecondaryAxis] = field(default_factory=list)
    spans: list[PlotSpan] = field(default_factory=list)
    backgrounds: list[ScalarBackground] = field(default_factory=list)
    legend_title: str | None = None
    legend_columns: int = 1
    show_legend: bool = True
    grid: bool = True
    invert_x_axis: bool = False
    invert_y_axis: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class ContourPlotSpec:
    """Neutral specification for a filled Cartesian contour plot.

    Parameters
    ----------
    key : str
        Stable plot identifier.
    title : str
        Figure title.
    filename_stem : str
        Default filename stem used by file renderers.
    x_axis, y_axis, value_axis : PlotAxis
        Horizontal, vertical, and mapped-value descriptions.
    x, y : ndarray
        One-dimensional coordinate grids.
    z : ndarray
        Two-dimensional mapped values with shape ``(len(y), len(x))``.
    colormap : str, optional
        Portable colormap name.
    mode : {"discrete", "smooth"}, optional
        Preferred filled-contour mode.
    levels : int, optional
        Number of principal contour levels.
    isolines : bool, optional
        Whether contour lines should be drawn.
    isoline_labels : bool, optional
        Whether isolines should be labelled.
    value_limits : tuple or None, optional
        Optional lower and upper mapped-value limits.
    center : float or None, optional
        Optional center for a diverging color normalization.
    masks : list of PlotMask, optional
        Boolean overlays used to mark diagnostic regions.
    series : list of PlotSeries, optional
        Ordinary Cartesian paths overlaid on the contour field.
    colored_paths : list of ColoredPathSeries, optional
        Scalar-colored Cartesian paths overlaid on the contour field.
    metadata : dict, optional
        Additional frontend-neutral plot information.
    """

    key: str
    title: str
    filename_stem: str
    x_axis: PlotAxis
    y_axis: PlotAxis
    value_axis: PlotAxis
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    colormap: str = "viridis"
    mode: Literal["discrete", "smooth"] = "smooth"
    levels: int = 12
    isolines: bool = True
    isoline_labels: bool = True
    value_limits: tuple[float, float] | None = None
    center: float | None = None
    masks: list[PlotMask] = field(default_factory=list)
    series: list[PlotSeries] = field(default_factory=list)
    colored_paths: list[ColoredPathSeries] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class PolarPlotPanel:
    """One polar panel in a multi-plane plot.

    Parameters
    ----------
    key : str
        Stable panel identifier, for example ``"xy"``.
    title : str
        Human-readable panel title.
    series : list of PlotSeries
        Ordered radial series. Their ``x`` values are angular coordinates.
    angle_unit : {"degree", "radian"}, optional
        Unit used by the angular series coordinates.
    radial_limit : float or None, optional
        Optional upper radial limit.
    theta_zero_location : str, optional
        Cardinal location used as angular zero.
    theta_direction : int, optional
        Angular direction, ``1`` for counter-clockwise and ``-1`` for
        clockwise.
    grid : bool, optional
        Whether the polar grid should be displayed.
    metadata : dict, optional
        Additional frontend-neutral panel information.
    """

    key: str
    title: str
    series: list[PlotSeries]
    angle_unit: Literal["degree", "radian"] = "degree"
    radial_limit: float | None = None
    theta_zero_location: str = "N"
    theta_direction: int = 1
    grid: bool = True
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class PolarPlotSpec:
    """Neutral specification for a multi-panel polar plot.

    Parameters
    ----------
    key : str
        Stable plot identifier.
    title : str
        Figure title.
    filename_stem : str
        Default filename stem used by file renderers.
    panels : list of PolarPlotPanel
        Ordered polar panels.
    show_legend : bool, optional
        Whether a figure-level legend should be displayed.
    legend_columns : int, optional
        Preferred number of legend columns.
    metadata : dict, optional
        Additional frontend-neutral plot information.
    """

    key: str
    title: str
    filename_stem: str
    panels: list[PolarPlotPanel]
    show_legend: bool = True
    legend_columns: int = 1
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class SurfaceStyle:
    """Frontend-neutral style hints for a three-dimensional surface.

    Parameters
    ----------
    color : str or None, optional
        Portable solid color used when ``colormap`` is ``None``.
    colormap : str or None, optional
        Portable colormap name used to map the physical values.
    opacity : float, optional
        Surface opacity between zero and one.
    show_colorbar : bool, optional
        Whether renderers should display a color scale.
    value_limits : tuple or None, optional
        Optional common lower and upper color limits.
    show_mesh : bool, optional
        Whether the surface mesh is outlined by the renderer.
    mesh_color : str or None, optional
        Optional edge color used when ``show_mesh`` is enabled.
    mesh_line_width : float, optional
        Width of the mesh edges when displayed.
    """

    color: str | None = None
    colormap: str | None = None
    opacity: float = 1.0
    show_colorbar: bool = True
    value_limits: tuple[float, float] | None = None
    show_mesh: bool = False
    mesh_color: str | None = None
    mesh_line_width: float = 0.5


@dataclass(slots=True)
class SurfaceLayer:
    """One Cartesian mesh in a neutral three-dimensional plot.

    Parameters
    ----------
    key : str
        Stable branch identifier.
    label : str
        Human-readable surface label.
    x, y, z : ndarray
        Two-dimensional Cartesian coordinates.
    values : ndarray
        Physical values used for color mapping and tooltips.
    theta, phi : ndarray or None, optional
        Optional angular grids in radians.
    style : SurfaceStyle, optional
        Portable style hints.
    metadata : dict, optional
        Additional frontend-neutral layer information.
    """

    key: str
    label: str
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    values: np.ndarray
    theta: np.ndarray | None = None
    phi: np.ndarray | None = None
    style: SurfaceStyle = field(default_factory=SurfaceStyle)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class VectorFieldStyle:
    """Portable style hints for a Cartesian vector or axial field.

    Parameters
    ----------
    color : str or None, optional
        Portable color specification. ``None`` delegates color selection to
        the renderer.
    line_width : float, optional
        Requested line width for arrows or axial line segments.
    scale : float, optional
        Multiplicative display scale applied by the renderer.
    opacity : float, optional
        Layer opacity between zero and one.
    """

    color: str | None = None
    line_width: float = 1.0
    scale: float = 1.0
    opacity: float = 1.0


@dataclass(slots=True)
class VectorFieldLayer:
    """Cartesian vector or axial field attached to a surface plot.

    Parameters
    ----------
    key : str
        Stable machine-readable layer identifier.
    label : str
        Human-readable layer label.
    origins : ndarray
        Cartesian origins with shape ``(n, 3)``.
    vectors : ndarray
        Cartesian directions with shape ``(n, 3)``.
    axial : bool, optional
        Whether opposite vector signs represent the same physical axis. Axial
        layers are rendered as centred line segments rather than arrows.
    resolved_mask : ndarray or None, optional
        Optional mask identifying uniquely resolved vectors or axes.
    style : VectorFieldStyle, optional
        Portable display hints.
    metadata : dict, optional
        Additional frontend-neutral layer information.
    """

    key: str
    label: str
    origins: np.ndarray
    vectors: np.ndarray
    axial: bool = False
    resolved_mask: np.ndarray | None = None
    style: VectorFieldStyle = field(default_factory=VectorFieldStyle)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class SurfacePlotSpec:
    """Neutral specification for one three-dimensional surface figure.

    Parameters
    ----------
    key : str
        Stable plot identifier.
    title : str
        Figure title.
    filename_stem : str
        Default filename stem used by file renderers.
    value_axis : PlotAxis
        Description of the radial and color-mapped physical quantity.
    layers : list of SurfaceLayer
        Ordered Cartesian surface meshes.
    vector_layers : list of VectorFieldLayer, optional
        Optional Cartesian vector or axial overlays.
    equal_aspect : bool, optional
        Whether all Cartesian axes should use the same scale.
    show_axes : bool, optional
        Whether renderers should display axes and labels.
    metadata : dict, optional
        Additional frontend-neutral plot information.
    """

    key: str
    title: str
    filename_stem: str
    value_axis: PlotAxis
    layers: list[SurfaceLayer]
    vector_layers: list[VectorFieldLayer] = field(default_factory=list)
    equal_aspect: bool = True
    show_axes: bool = True
    metadata: dict[str, Any] = field(default_factory=dict)


SphericalProjection: TypeAlias = Literal["equal_area", "stereographic"]


@dataclass(slots=True)
class SphericalMarker:
    """Marker attached to one or more directions on a spherical map.

    Parameters
    ----------
    key : str
        Stable machine-readable marker identifier.
    label : str
        Human-readable marker label.
    directions : ndarray
        Cartesian unit directions with shape ``(n, 3)``.
    marker : str, optional
        Portable marker symbol.
    metadata : dict, optional
        Additional frontend-neutral marker information.
    """

    key: str
    label: str
    directions: np.ndarray
    marker: str = "circle"
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class AxisFieldLayer:
    """Axial-vector field sampled on directions of a spherical map.

    Parameters
    ----------
    key : str
        Stable machine-readable layer identifier.
    label : str
        Human-readable layer label.
    directions : ndarray
        Cartesian unit directions at which the axes are drawn.
    axes : ndarray
        Cartesian unit axes with shape ``(n, 3)``. Opposite signs represent
        the same physical axis.
    resolved_mask : ndarray or None, optional
        Optional mask identifying uniquely resolved axes.
    metadata : dict, optional
        Additional frontend-neutral layer information.
    """

    key: str
    label: str
    directions: np.ndarray
    axes: np.ndarray
    resolved_mask: np.ndarray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class SphericalMapSpec:
    """Neutral specification for a scalar field on a spherical domain.

    Parameters
    ----------
    key : str
        Stable plot identifier.
    title : str
        Figure title.
    filename_stem : str
        Default filename stem used by file renderers.
    theta, phi : ndarray
        One-dimensional polar and azimuthal coordinates in radians.
    values : ndarray
        Scalar field with shape ``(len(theta), len(phi))``.
    value_axis : PlotAxis
        Description of the mapped physical quantity.
    hemisphere : {"upper", "lower", "full"}
        Spherical domain represented by the data.
    projection : {"equal_area", "stereographic"}, optional
        Preferred map projection.
    colormap : str, optional
        Portable colormap name.
    levels : int, optional
        Preferred number of contour levels.
    isolines : bool, optional
        Whether isolines should be displayed.
    markers : list of SphericalMarker, optional
        Directional extrema or other annotations.
    axis_layers : list of AxisFieldLayer, optional
        Axial-vector overlays such as polarization axes.
    metadata : dict, optional
        Additional frontend-neutral plot information.
    """

    key: str
    title: str
    filename_stem: str
    theta: np.ndarray
    phi: np.ndarray
    values: np.ndarray
    value_axis: PlotAxis
    hemisphere: Literal["upper", "lower", "full"]
    projection: SphericalProjection = "equal_area"
    colormap: str = "viridis"
    levels: int = 12
    isolines: bool = True
    markers: list[SphericalMarker] = field(default_factory=list)
    axis_layers: list[AxisFieldLayer] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class SphericalSummarySpec:
    """Multi-panel summary assembled from spherical scalar maps.

    Parameters
    ----------
    key : str
        Stable plot identifier.
    title : str
        Figure title.
    filename_stem : str
        Default filename stem used by file renderers.
    maps : list of SphericalMapSpec
        Ordered spherical maps shown in the summary.
    columns : int, optional
        Preferred number of panel columns.
    metadata : dict, optional
        Additional frontend-neutral layout information.
    """

    key: str
    title: str
    filename_stem: str
    maps: list[SphericalMapSpec]
    columns: int = 3
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class PanelPlotSpec:
    """Neutral multi-panel layout for Cartesian line and contour plots.

    Parameters
    ----------
    key : str
        Stable plot identifier.
    title : str
        Figure-level title.
    filename_stem : str
        Default filename stem used by file renderers.
    panels : list of LinePlotSpec or ContourPlotSpec
        Ordered Cartesian panels.
    columns : int, optional
        Preferred number of panel columns.
    share_x, share_y : bool, optional
        Whether compatible panels should share their respective axes.
    metadata : dict, optional
        Additional frontend-neutral layout information.

    Raises
    ------
    ValueError
        If no panels are provided or ``columns`` is not positive.
    """

    key: str
    title: str
    filename_stem: str
    panels: list[LinePlotSpec | ContourPlotSpec]
    columns: int = 2
    share_x: bool = False
    share_y: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate the Cartesian panel collection."""
        self.panels = list(self.panels)
        if not self.panels:
            raise ValueError("a panel plot requires at least one panel")
        if self.columns < 1:
            raise ValueError("panel-plot columns must be positive")
        self.metadata = dict(self.metadata)


PlotSpec: TypeAlias = (
    LinePlotSpec
    | ContourPlotSpec
    | PolarPlotSpec
    | SurfacePlotSpec
    | SphericalMapSpec
    | SphericalSummarySpec
    | PanelPlotSpec
)


@dataclass(slots=True)
class PlotCollection:
    """Collection of neutral plot specifications and non-fatal warnings.

    Parameters
    ----------
    plots : list of PlotSpec, optional
        Ordered specifications ready for rendering.
    warnings : list of str, optional
        Non-fatal conditions encountered while preparing plot data.
    """

    plots: list[PlotSpec] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
