# -*- coding: utf-8 -*-

"""Passive options for frontend-neutral thermoelastic plot builders."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, TypeAlias


ThermoelasticPlotPreset: TypeAlias = Literal[
    "analysis",
    "publication",
    "monochrome",
]
ThermoelasticPlotLayout: TypeAlias = Literal[
    "auto",
    "overlay",
    "facets",
    "separate",
]
ThermoelasticUncertaintyMode: TypeAlias = Literal["auto", "none", "band", "bars"]
ThermoelasticPTQuantity: TypeAlias = Literal[
    "value",
    "uncertainty",
    "relative-uncertainty",
]
ThermoelasticProfileMode: TypeAlias = Literal["absolute", "relative"]
ThermoelasticTensorCondition: TypeAlias = Literal["isothermal", "adiabatic"]
ThermoelasticProfileColor: TypeAlias = Literal[
    "component",
    "temperature",
    "none",
]
ThermoelasticProfileBackground: TypeAlias = Literal["none", "temperature"]


@dataclass(slots=True)
class ThermoelasticPlotStyleOptions:
    """Common presentation options shared by thermoelastic plot builders.

    Parameters
    ----------
    preset : {"analysis", "publication", "monochrome"}, optional
        High-level style profile. Explicit fields below override only those
        choices that are represented directly in this passive contract.
    show_title : bool or None, optional
        Whether plot specifications retain titles. ``None`` selects the preset
        default.
    grid : bool or None, optional
        Whether line plots request Cartesian grids. ``None`` selects the preset
        default.
    line_width : float or None, optional
        Width of component curves. ``None`` selects the preset default.
    marker_size : float or None, optional
        Marker size in typographic points. ``None`` selects the preset default.
    marker_edge_color : str or None, optional
        Marker border color. The default black edge keeps colored markers
        legible on scalar backgrounds and dense contour maps.
    marker_edge_width : float, optional
        Marker border width in typographic points.
    errorbar_width : float, optional
        Width of error-bar strokes.
    errorbar_capsize : float, optional
        Error-bar cap size in typographic points.
    """

    preset: ThermoelasticPlotPreset = "analysis"
    show_title: bool | None = None
    grid: bool | None = None
    line_width: float | None = None
    marker_size: float | None = None
    marker_edge_color: str | None = "black"
    marker_edge_width: float = 0.7
    errorbar_width: float = 1.0
    errorbar_capsize: float = 2.5

    def __post_init__(self) -> None:
        """Validate common style controls."""
        if self.preset not in {"analysis", "publication", "monochrome"}:
            raise ValueError("invalid thermoelastic plot preset")
        for name in ("line_width", "marker_size"):
            value = getattr(self, name)
            if value is not None and float(value) <= 0.0:
                raise ValueError(f"{name} must be positive")
        if self.marker_edge_width < 0.0:
            raise ValueError("marker_edge_width must be non-negative")
        if self.errorbar_width <= 0.0:
            raise ValueError("errorbar_width must be positive")
        if self.errorbar_capsize < 0.0:
            raise ValueError("errorbar_capsize must be non-negative")


@dataclass(slots=True)
class ThermoelasticFitPlotOptions:
    """Options controlling observed-versus-fitted ``C_IJ(V)`` plots.

    Parameters
    ----------
    style : ThermoelasticPlotStyleOptions, optional
        Common presentation controls.
    residuals : bool, optional
        Include a residual panel beneath each fit.
    confidence : float or None, optional
        Central confidence probability used for the propagated prediction
        band. ``None`` disables the band.
    curve_points : int, optional
        Number of volumes used for the smooth fitted curve.
    show_symmetry_spread : bool, optional
        Display the point-level spread among symmetry-equivalent tensor
        entries as vertical diagnostic error bars. This spread is not a
        statistical measurement uncertainty.
    """

    style: ThermoelasticPlotStyleOptions = field(
        default_factory=ThermoelasticPlotStyleOptions
    )
    residuals: bool = True
    confidence: float | None = 0.95
    curve_points: int = 300
    show_symmetry_spread: bool = False

    def __post_init__(self) -> None:
        """Validate fit-plot controls."""
        if self.curve_points < 20:
            raise ValueError("curve_points must be at least 20")
        if self.confidence is not None and not 0.0 < self.confidence < 1.0:
            raise ValueError("confidence must lie strictly between zero and one")


@dataclass(slots=True)
class ThermoelasticPTPlotOptions:
    """Options controlling pressure-temperature component maps.

    Parameters
    ----------
    style : ThermoelasticPlotStyleOptions, optional
        Common presentation controls.
    tensor_condition : {"isothermal", "adiabatic"}, optional
        Thermodynamic stiffness field represented by the map.
    quantity : {"value", "uncertainty", "relative-uncertainty"}, optional
        Quantity mapped by color.
    layout : {"auto", "facets", "separate"}, optional
        Multi-component layout. ``"auto"`` uses facets for up to four maps and
        separate figures otherwise.
    colormap : str or None, optional
        Explicit portable colormap. ``None`` chooses a scientifically suitable
        sequential or diverging map from the requested quantity.
    levels : int, optional
        Number of principal contour levels.
    isolines, isoline_labels : bool, optional
        Control contour-line overlays and labels.
    show_extrapolation : bool, optional
        Overlay distinct QHA-coordinate and elastic-volume extrapolation masks.
    panel_columns : int, optional
        Preferred number of columns for faceted maps.
    """

    style: ThermoelasticPlotStyleOptions = field(
        default_factory=ThermoelasticPlotStyleOptions
    )
    tensor_condition: ThermoelasticTensorCondition = "isothermal"
    quantity: ThermoelasticPTQuantity = "value"
    layout: ThermoelasticPlotLayout = "auto"
    colormap: str | None = None
    levels: int = 12
    isolines: bool = True
    isoline_labels: bool = True
    show_extrapolation: bool = True
    panel_columns: int = 2

    def __post_init__(self) -> None:
        """Validate pressure-temperature map controls."""
        if self.tensor_condition not in {"isothermal", "adiabatic"}:
            raise ValueError("invalid thermoelastic tensor condition")
        if self.quantity not in {
            "value",
            "uncertainty",
            "relative-uncertainty",
        }:
            raise ValueError("invalid thermoelastic P-T quantity")
        if self.layout not in {"auto", "facets", "separate"}:
            raise ValueError("P-T maps support auto, facets, or separate layouts")
        if self.levels < 2:
            raise ValueError("levels must be at least 2")
        if self.panel_columns < 1:
            raise ValueError("panel_columns must be positive")


@dataclass(slots=True)
class ThermoelasticProfilePlotOptions:
    """Options controlling elastic-component depth-profile plots.

    Parameters
    ----------
    style : ThermoelasticPlotStyleOptions, optional
        Common presentation controls.
    tensor_condition : {"isothermal", "adiabatic"}, optional
        Thermodynamic stiffness field represented along depth.
    mode : {"absolute", "relative"}, optional
        Absolute stiffness or percentage change relative to a reference depth.
    layout : {"auto", "overlay", "facets", "separate"}, optional
        Multi-component layout. ``"auto"`` overlays up to six components and
        uses facets above that threshold.
    uncertainty : {"auto", "none", "band", "bars"}, optional
        Representation of one-sigma component uncertainty. ``"auto"`` uses
        a band for one selected component and omits uncertainty when several
        components are requested together.
    color_by : {"component", "temperature", "none"}, optional
        Variable controlling line color. Temperature coloring uses a neutral
        scalar-colored path specification.
    background : {"none", "temperature"}, optional
        Optional temperature field painted behind the profile. When enabled,
        temperature is not also used to color the path.
    reference_depth : float or None, optional
        Reference depth in km for relative profiles. ``None`` uses the first
        profile point.
    show_pressure_axis : bool, optional
        Add a prepared pressure axis when pressure is monotonic with depth.
    show_extrapolation : bool, optional
        Mark contiguous QHA and elastic extrapolation intervals.
    max_overlay_components : int, optional
        Maximum number of components selected by ``layout="auto"`` for one
        shared axis.
    panel_columns : int, optional
        Preferred number of columns for faceted profiles.
    temperature_colormap : str, optional
        Portable colormap used for temperature encoding.
    """

    style: ThermoelasticPlotStyleOptions = field(
        default_factory=ThermoelasticPlotStyleOptions
    )
    tensor_condition: ThermoelasticTensorCondition = "isothermal"
    mode: ThermoelasticProfileMode = "absolute"
    layout: ThermoelasticPlotLayout = "auto"
    uncertainty: ThermoelasticUncertaintyMode = "auto"
    color_by: ThermoelasticProfileColor = "temperature"
    background: ThermoelasticProfileBackground = "none"
    reference_depth: float | None = None
    show_pressure_axis: bool = True
    show_extrapolation: bool = True
    max_overlay_components: int = 6
    panel_columns: int = 2
    temperature_colormap: str = "plasma"

    def __post_init__(self) -> None:
        """Validate depth-profile controls."""
        if self.tensor_condition not in {"isothermal", "adiabatic"}:
            raise ValueError("invalid thermoelastic tensor condition")
        if self.mode not in {"absolute", "relative"}:
            raise ValueError("invalid profile mode")
        if self.layout not in {"auto", "overlay", "facets", "separate"}:
            raise ValueError("invalid profile layout")
        if self.uncertainty not in {"auto", "none", "band", "bars"}:
            raise ValueError("invalid profile uncertainty mode")
        if self.color_by not in {"component", "temperature", "none"}:
            raise ValueError("invalid profile color mode")
        if self.background not in {"none", "temperature"}:
            raise ValueError("invalid profile background")
        if self.reference_depth is not None:
            self.reference_depth = float(self.reference_depth)
        if self.max_overlay_components < 1:
            raise ValueError("max_overlay_components must be positive")
        if self.panel_columns < 1:
            raise ValueError("panel_columns must be positive")


@dataclass(slots=True)
class ThermoelasticDomainPlotOptions:
    """Options controlling pressure-temperature domain diagnostic plots.

    Parameters
    ----------
    style : ThermoelasticPlotStyleOptions, optional
        Common presentation controls.
    colormap : str, optional
        Portable colormap used for equilibrium volume.
    levels : int, optional
        Number of principal volume contours.
    show_profiles : bool, optional
        Overlay all or selected archived depth profiles.
    color_profiles_by_depth : bool, optional
        Color profile paths by depth instead of using ordinary component-like
        line styles.
    profile_colormap : str, optional
        Colormap used for depth along overlaid geological paths. This is kept
        independent from the equilibrium-volume colormap.
    show_extrapolation : bool, optional
        Overlay QHA-coordinate and elastic-volume extrapolation masks.
    """

    style: ThermoelasticPlotStyleOptions = field(
        default_factory=ThermoelasticPlotStyleOptions
    )
    colormap: str = "viridis"
    levels: int = 12
    show_profiles: bool = True
    color_profiles_by_depth: bool = True
    profile_colormap: str = "plasma"
    show_extrapolation: bool = True

    def __post_init__(self) -> None:
        """Validate domain-map controls."""
        if self.levels < 2:
            raise ValueError("levels must be at least 2")


__all__ = [
    "ThermoelasticDomainPlotOptions",
    "ThermoelasticFitPlotOptions",
    "ThermoelasticPTPlotOptions",
    "ThermoelasticPTQuantity",
    "ThermoelasticPlotLayout",
    "ThermoelasticPlotPreset",
    "ThermoelasticPlotStyleOptions",
    "ThermoelasticProfileBackground",
    "ThermoelasticProfileColor",
    "ThermoelasticProfileMode",
    "ThermoelasticProfilePlotOptions",
    "ThermoelasticTensorCondition",
    "ThermoelasticUncertaintyMode",
]


@dataclass(slots=True)
class ThermoelasticComparePlotOptions:
    """Options for isothermal-versus-adiabatic line comparisons.

    Parameters
    ----------
    style : ThermoelasticPlotStyleOptions, optional
        Common frontend-neutral presentation controls.
    fixed_pressure : float or None, optional
        Pressure in GPa held fixed while temperature varies.
    fixed_temperature : float or None, optional
        Temperature in K held fixed while pressure varies.
    layout : {"auto", "facets", "separate"}, optional
        Multi-component arrangement.
    panel_columns : int, optional
        Preferred number of columns for faceted output.
    """

    style: ThermoelasticPlotStyleOptions = field(
        default_factory=ThermoelasticPlotStyleOptions
    )
    fixed_pressure: float | None = None
    fixed_temperature: float | None = None
    layout: ThermoelasticPlotLayout = "auto"
    panel_columns: int = 2

    def __post_init__(self) -> None:
        """Validate that exactly one thermodynamic coordinate is fixed."""
        if (self.fixed_pressure is None) == (self.fixed_temperature is None):
            raise ValueError("set exactly one of fixed_pressure or fixed_temperature")
        if self.fixed_pressure is not None:
            self.fixed_pressure = float(self.fixed_pressure)
        if self.fixed_temperature is not None:
            self.fixed_temperature = float(self.fixed_temperature)
        if self.layout not in {"auto", "facets", "separate"}:
            raise ValueError("compare plots support auto, facets, or separate layouts")
        if self.panel_columns < 1:
            raise ValueError("panel_columns must be positive")
