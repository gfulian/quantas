# -*- coding: utf-8 -*-

"""Public frontend-neutral plotting contracts.

This namespace exposes the passive plot specifications produced by Quantas
scientific builders.  The contracts contain prepared numerical data,
scientific axes, masks, overlays, uncertainty information, and portable style
hints without depending on Matplotlib, Plotly, Dash, Rich, or another concrete
renderer.

Applications should import these types from :mod:`quantas.api.plotting`
rather than from implementation packages.  The objects are the same
authoritative contracts used internally by Quantas; this module does not copy
or wrap plot data.
"""

from __future__ import annotations

from quantas.models.plot import (
    AxisFieldLayer,
    AxisLocation,
    AxisOrientation,
    BandOrientation,
    ColoredPathSeries,
    ColoredPathStyle,
    ContourPlotSpec,
    LinePlotSpec,
    LineStyle,
    PanelPlotSpec,
    PlotAxis,
    PlotBand,
    PlotBandStyle,
    PlotCollection,
    PlotMask,
    PlotSeries,
    PlotSeriesStyle,
    PlotSpan,
    PlotSpec,
    PolarPlotPanel,
    PolarPlotSpec,
    ScalarBackground,
    SecondaryAxis,
    SphericalMapSpec,
    SphericalMarker,
    SphericalProjection,
    SphericalSummarySpec,
    SurfaceLayer,
    SurfacePlotSpec,
    SurfaceStyle,
    VectorFieldLayer,
    VectorFieldStyle,
)
from quantas.models.plot_inventory import (
    PlotContextDescriptor,
    PlotContextValue,
    PlotInventory,
    PlotKind,
    PlotPropertyDescriptor,
    PlotRepresentationDescriptor,
)


def __dir__() -> list[str]:
    """Return only names in the supported public plotting contract."""
    return sorted(__all__)


__all__ = [
    "AxisFieldLayer",
    "AxisLocation",
    "AxisOrientation",
    "BandOrientation",
    "ColoredPathSeries",
    "ColoredPathStyle",
    "ContourPlotSpec",
    "LinePlotSpec",
    "LineStyle",
    "PanelPlotSpec",
    "PlotAxis",
    "PlotBand",
    "PlotBandStyle",
    "PlotCollection",
    "PlotContextDescriptor",
    "PlotContextValue",
    "PlotInventory",
    "PlotKind",
    "PlotPropertyDescriptor",
    "PlotRepresentationDescriptor",
    "PlotMask",
    "PlotSeries",
    "PlotSeriesStyle",
    "PlotSpan",
    "PlotSpec",
    "PolarPlotPanel",
    "PolarPlotSpec",
    "ScalarBackground",
    "SecondaryAxis",
    "SphericalMapSpec",
    "SphericalMarker",
    "SphericalProjection",
    "SphericalSummarySpec",
    "SurfaceLayer",
    "SurfacePlotSpec",
    "SurfaceStyle",
    "VectorFieldLayer",
    "VectorFieldStyle",
]
