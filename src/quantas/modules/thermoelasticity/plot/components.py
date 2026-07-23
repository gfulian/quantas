# -*- coding: utf-8 -*-

"""Component selection, labels, and stable styles for thermoelastic plots."""

from __future__ import annotations

from dataclasses import dataclass

from quantas.models import LineStyle, PlotSeriesStyle
from quantas.modules.thermoelasticity.components import (
    ThermoelasticComponentGroup,
    VOIGT_COMPONENTS,
    component_indices,
    component_symbol,
    normalize_component_label,
    resolve_components,
)


_COMPONENT_MARKERS: tuple[str, ...] = (
    "o",
    "s",
    "^",
    "v",
    "D",
    "P",
    "X",
    "<",
    ">",
    "h",
    "p",
    "d",
    "8",
)
_COMPONENT_LINE_STYLES: tuple[LineStyle, ...] = (
    "solid",
    "dashed",
    "dotted",
    "dashdot",
)


@dataclass(frozen=True, slots=True)
class ThermoelasticComponentStyle:
    """Stable portable style assigned to one Voigt stiffness component.

    Parameters
    ----------
    label : str
        Canonical component label.
    marker : str or None
        Stable marker symbol. Components beyond the compact portable marker
        palette use line style and color without a marker.
    line_style : LineStyle
        Stable portable line style.
    color_index : int
        Deterministic index that renderers or builders may map to a palette.
    """

    label: str
    marker: str | None
    line_style: LineStyle
    color_index: int


def component_style(label: str) -> ThermoelasticComponentStyle:
    """Return the stable style assignment for one component.

    Parameters
    ----------
    label : str
        Component label.

    Returns
    -------
    ThermoelasticComponentStyle
        Deterministic marker, line style, and palette index.
    """
    canonical = normalize_component_label(label)
    index = VOIGT_COMPONENTS.index(canonical)
    return ThermoelasticComponentStyle(
        label=canonical,
        marker=(_COMPONENT_MARKERS[index] if index < len(_COMPONENT_MARKERS) else None),
        line_style=_COMPONENT_LINE_STYLES[index % len(_COMPONENT_LINE_STYLES)],
        color_index=index,
    )


def plot_series_style(
    label: str,
    *,
    preset: str,
    line_width: float,
    marker_size: float,
    color: str | None = None,
    marker_edge_color: str | None = "black",
    marker_edge_width: float = 0.7,
    line_only: bool = False,
) -> PlotSeriesStyle:
    """Build portable series style hints for one component.

    Parameters
    ----------
    label : str
        Component label.
    preset : str
        Thermoelastic style preset.
    line_width, marker_size : float
        Requested line and marker dimensions.
    color : str or None, optional
        Explicit portable color.
    marker_edge_color : str or None, optional
        Portable marker border color.
    marker_edge_width : float, optional
        Marker border width.
    line_only : bool, optional
        Omit markers while retaining the stable line style.

    Returns
    -------
    PlotSeriesStyle
        Frontend-neutral style object.
    """
    assigned = component_style(label)
    resolved_color = "black" if preset == "monochrome" else color
    return PlotSeriesStyle(
        color=resolved_color,
        line_style=assigned.line_style,
        line_width=float(line_width),
        marker=None if line_only else assigned.marker,
        marker_size=float(marker_size),
        marker_edge_color=marker_edge_color,
        marker_edge_width=float(marker_edge_width),
    )


__all__ = [
    "ThermoelasticComponentGroup",
    "ThermoelasticComponentStyle",
    "VOIGT_COMPONENTS",
    "component_indices",
    "component_style",
    "component_symbol",
    "normalize_component_label",
    "plot_series_style",
    "resolve_components",
]
