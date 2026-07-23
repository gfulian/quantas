# -*- coding: utf-8 -*-

"""Neutral two- and three-dimensional plot builders for elasticity."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from quantas.core.physics.elasticity import ElasticitySurfaceCollection
from quantas.models import (
    PlotAxis,
    PlotCollection,
    PlotSeries,
    PlotSeriesStyle,
    PolarPlotPanel,
    PolarPlotSpec,
    SurfaceLayer,
    SurfacePlotSpec,
    SurfaceStyle,
)
from quantas.modules.elasticity.models import ElasticityResult


PLANES = ("xy", "xz", "yz")
ElasticityPlotProperty = Literal["young", "compressibility", "shear", "poisson"]
ElasticitySurfaceGeometry = Literal["physical", "unit_sphere"]


@dataclass(frozen=True, slots=True)
class _PolarProperty:
    """Internal description of one two-dimensional elasticity property."""

    stored_name: str
    title: str
    filename_stem: str
    labels: tuple[str, ...]
    styles: tuple[PlotSeriesStyle, ...]
    radial_column: int


_POLAR_PROPERTIES: dict[ElasticityPlotProperty, _PolarProperty] = {
    "young": _PolarProperty(
        stored_name="young_modulus",
        title="Young's modulus",
        filename_stem="2d_young",
        labels=("Young's modulus",),
        styles=(PlotSeriesStyle(color="#009010", line_width=2.0),),
        radial_column=0,
    ),
    "compressibility": _PolarProperty(
        stored_name="linear_compressibility",
        title="Linear compressibility",
        filename_stem="2d_compressibility",
        labels=("Positive", "Negative"),
        styles=(
            PlotSeriesStyle(color="#009010", line_width=2.0),
            PlotSeriesStyle(color="red", line_width=2.0),
        ),
        radial_column=0,
    ),
    "shear": _PolarProperty(
        stored_name="shear_modulus",
        title="Shear modulus",
        filename_stem="2d_shear",
        labels=("Minimum", "Maximum"),
        styles=(
            PlotSeriesStyle(color="#009010", line_style="dashed", line_width=2.0),
            PlotSeriesStyle(color="red", line_width=2.0),
        ),
        radial_column=1,
    ),
    "poisson": _PolarProperty(
        stored_name="poisson_ratio",
        title="Poisson's ratio",
        filename_stem="2d_poisson",
        labels=("Negative", "Minimum", "Maximum"),
        styles=(
            PlotSeriesStyle(color="blue", line_style="dashed", line_width=2.0),
            PlotSeriesStyle(color="#009010", line_width=2.0),
            PlotSeriesStyle(color="red", line_width=2.0),
        ),
        radial_column=2,
    ),
}


@dataclass(frozen=True, slots=True)
class _SurfaceProperty:
    """Internal display description of one three-dimensional surface branch."""

    title: str
    filename_stem: str
    label: str
    color: str


_SURFACE_PROPERTIES: dict[str, _SurfaceProperty] = {
    "young": _SurfaceProperty(
        "Young's modulus",
        "3d_young",
        "Young's modulus",
        "#009010",
    ),
    "compressibility_positive": _SurfaceProperty(
        "Positive linear compressibility",
        "3d_compressibility_positive",
        "Positive",
        "#009010",
    ),
    "compressibility_negative": _SurfaceProperty(
        "Negative linear compressibility",
        "3d_compressibility_negative",
        "Negative",
        "red",
    ),
    "shear_minimum": _SurfaceProperty(
        "Minimum shear modulus",
        "3d_shear_minimum",
        "Minimum",
        "#009010",
    ),
    "shear_maximum": _SurfaceProperty(
        "Maximum shear modulus",
        "3d_shear_maximum",
        "Maximum",
        "red",
    ),
    "poisson_negative": _SurfaceProperty(
        "Negative Poisson's ratio",
        "3d_poisson_negative",
        "Negative",
        "blue",
    ),
    "poisson_minimum": _SurfaceProperty(
        "Minimum positive Poisson's ratio",
        "3d_poisson_minimum",
        "Minimum",
        "#009010",
    ),
    "poisson_maximum": _SurfaceProperty(
        "Maximum Poisson's ratio",
        "3d_poisson_maximum",
        "Maximum",
        "red",
    ),
}


def build_elasticity_2d_plot_spec(
    result: ElasticityResult,
    property_name: ElasticityPlotProperty,
) -> PolarPlotSpec:
    """Build one neutral three-plane polar plot specification.

    Parameters
    ----------
    result : ElasticityResult
        Elasticity result containing principal-plane data.
    property_name : str
        Property family to plot.

    Returns
    -------
    PolarPlotSpec
        Frontend-neutral three-panel polar specification.

    Raises
    ------
    KeyError
        If the property is missing from a principal plane.
    ValueError
        If stored arrays have incompatible dimensions.
    """
    description = _POLAR_PROPERTIES[property_name]
    panels: list[PolarPlotPanel] = []
    for plane in PLANES:
        try:
            plane_data = result.properties_2d[plane]
            raw_values = plane_data[description.stored_name]
            raw_angle = plane_data["phi" if plane == "xy" else "theta"]
        except KeyError as exc:
            raise KeyError(
                f"Elasticity 2D property '{property_name}' is unavailable "
                f"on plane '{plane}'"
            ) from exc

        angle = np.degrees(np.asarray(raw_angle, dtype=np.float64))
        values = np.asarray(raw_values, dtype=np.float64)
        if values.ndim == 1:
            values = values[:, np.newaxis]
        if values.ndim != 2 or values.shape[0] != angle.size:
            raise ValueError(
                f"Elasticity 2D property '{property_name}' on plane '{plane}' "
                f"has shape {values.shape}, incompatible with {angle.size} angles"
            )
        if values.shape[1] != len(description.labels):
            raise ValueError(
                f"Elasticity 2D property '{property_name}' on plane '{plane}' has "
                f"{values.shape[1]} series, expected {len(description.labels)}"
            )

        series = [
            PlotSeries(
                key=f"{property_name}_{index}",
                label=label,
                x=angle.copy(),
                y=values[:, index].copy(),
                style=style,
                metadata={"plane": plane, "series_index": index},
            )
            for index, (label, style) in enumerate(
                zip(description.labels, description.styles, strict=True)
            )
        ]
        radial_data = np.abs(values[:, description.radial_column])
        radial_max = float(np.nanmax(radial_data)) * 1.1
        panels.append(
            PolarPlotPanel(
                key=plane,
                title=f"{description.title} ({plane})",
                series=series,
                angle_unit="degree",
                radial_limit=radial_max,
                theta_zero_location="N",
                theta_direction=1,
                metadata={"plane": plane},
            )
        )

    return PolarPlotSpec(
        key=f"2d_{property_name}",
        title=description.title,
        filename_stem=description.filename_stem,
        panels=panels,
        show_legend=len(description.labels) > 1,
        legend_columns=len(description.labels),
        metadata={
            "module": "elasticity",
            "dimension": 2,
            "property": property_name,
        },
    )


def build_elasticity_2d_plot_collection(
    result: ElasticityResult,
    properties: tuple[ElasticityPlotProperty, ...] | None = None,
) -> PlotCollection:
    """Build neutral polar specifications for available 2D properties."""
    collection = PlotCollection()
    if not result.properties_2d:
        return collection
    selected = tuple(_POLAR_PROPERTIES) if properties is None else properties
    for property_name in selected:
        description = _POLAR_PROPERTIES[property_name]
        if not _has_2d_property(result, description.stored_name):
            continue
        try:
            collection.plots.append(
                build_elasticity_2d_plot_spec(result, property_name)
            )
        except (KeyError, ValueError) as exc:
            collection.warnings.append(str(exc))
    return collection


def build_elasticity_surface_plot_collection(
    surfaces: ElasticitySurfaceCollection,
    *,
    geometry: str = "physical",
    color_mode: Literal["solid", "property"] = "property",
    colormap: str = "viridis",
    show_mesh: bool = False,
    mesh_color: str = "black",
    mesh_line_width: float = 0.5,
) -> PlotCollection:
    """Build neutral 3D surface specifications from numerical meshes.

    Parameters
    ----------
    surfaces : ElasticitySurfaceCollection
        Numerical elasticity surfaces.
    geometry : {"normalized", "physical", "unit_sphere"}, optional
        Radial geometry. Physical surfaces use the property magnitude directly
        and unit-sphere maps encode the property only through color. The legacy
        ``"normalized"`` geometry is accepted as an alias of ``"physical"``.
    color_mode : {"solid", "property"}, optional
        Use fixed semantic colors or map physical values through ``colormap``.
    colormap : str, optional
        Matplotlib-compatible colormap name used in property-color mode.

    Returns
    -------
    PlotCollection
        One independent plot specification for every available branch.

    Raises
    ------
    ValueError
        If a geometry or color mode is unsupported.
    """
    geometry = "physical" if geometry == "normalized" else geometry
    if geometry not in {"physical", "unit_sphere"}:
        raise ValueError("geometry must be 'physical' or 'unit_sphere'.")
    if color_mode not in {"solid", "property"}:
        raise ValueError("color_mode must be 'solid' or 'property'.")
    collection = PlotCollection(warnings=list(surfaces.warnings))
    limits = _shared_property_limits(surfaces)
    for key, surface in surfaces.surfaces.items():
        description = _SURFACE_PROPERTIES[key]
        x, y, z, values, theta, phi, radial_scale = _surface_plot_arrays(
            surface,
            geometry=geometry,
        )
        style = SurfaceStyle(
            color=description.color if color_mode == "solid" else None,
            colormap=colormap if color_mode == "property" else None,
            opacity=1.0,
            show_colorbar=color_mode == "property",
            value_limits=limits.get(surface.property_name),
            show_mesh=show_mesh,
            mesh_color=mesh_color,
            mesh_line_width=mesh_line_width,
        )
        layer = SurfaceLayer(
            key=key,
            label=description.label,
            x=x,
            y=y,
            z=z,
            values=values,
            theta=theta,
            phi=phi,
            style=style,
            metadata={
                "property": surface.property_name,
                "branch": surface.branch,
                "unit": surface.unit,
                "geometry": geometry,
                "radial_scale": radial_scale,
            },
        )
        filename_stem = (
            description.filename_stem
            if geometry == "physical"
            else f"{description.filename_stem}_{geometry}"
        )
        collection.plots.append(
            SurfacePlotSpec(
                key=key,
                title=description.title,
                filename_stem=filename_stem,
                value_axis=PlotAxis(
                    key=surface.property_name,
                    label=_quantity_label(description.title, surface.unit),
                    unit=surface.unit or None,
                ),
                layers=[layer],
                equal_aspect=True,
                show_axes=True,
                metadata={
                    "module": "elasticity",
                    "dimension": 3,
                    "property": surface.property_name,
                    "branch": surface.branch,
                    "geometry": geometry,
                    "radial_scale": radial_scale,
                    "axis_labels": ("X", "Y", "Z"),
                    "axis_draw_mode": "full_equal",
                    "axis_half_length_factor": 1.14,
                    "axis_label_factor": 1.07,
                },
            )
        )
    return collection


def build_elasticity_plot_collection(result: ElasticityResult) -> PlotCollection:
    """Build the default persisted elasticity plot collection.

    This function implements the uniform module contract and therefore uses
    only the 2D data already present in ``result``.  Transient 3D surfaces are
    built explicitly with :func:`build_elasticity_surface_plot_collection`.
    """
    return build_elasticity_2d_plot_collection(result)


def _has_2d_property(result: ElasticityResult, stored_name: str) -> bool:
    """Return whether a stored 2D property exists on all principal planes."""
    return all(
        plane in result.properties_2d and stored_name in result.properties_2d[plane]
        for plane in PLANES
    )


def _surface_plot_arrays(
    surface: object,
    *,
    geometry: str,
) -> tuple[np.ndarray, ...]:
    """Return closed Cartesian arrays for one selected geometry."""
    theta = np.asarray(getattr(surface, "theta"), dtype=float)
    phi = np.asarray(getattr(surface, "phi"), dtype=float)
    values = np.asarray(getattr(surface, "values"), dtype=float)
    physical_radius = np.asarray(getattr(surface, "radius"), dtype=float)
    if geometry == "physical":
        radius = physical_radius
        radial_scale = 1.0
    else:
        radius = np.where(np.isfinite(values), 1.0, np.nan)
        radial_scale = 1.0

    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)
    arrays = [x, y, z, values, theta, phi]
    closed = [np.concatenate((array, array[:, :1]), axis=1) for array in arrays]
    closed[-1][:, -1] = 2.0 * np.pi
    return (*closed, float(radial_scale))


def _format_unit(unit: str) -> str:
    """Return a Matplotlib-friendly unit string."""
    normalized = unit.strip()
    if not normalized or normalized == "dimensionless":
        return ""
    replacements = {
        "GPa": r"$\mathrm{GPa}$",
        "TPa^-1": r"$\mathrm{TPa^{-1}}$",
        "dimensionless": "",
    }
    return replacements.get(normalized, normalized)


def _quantity_label(title: str, unit: str) -> str:
    """Return a publication-style quantity label with unit."""
    rendered = _format_unit(unit)
    return title if not rendered else f"{title} ({rendered})"


def _shared_property_limits(
    surfaces: ElasticitySurfaceCollection,
) -> dict[str, tuple[float, float]]:
    """Return common finite color limits for branches of one property."""
    grouped: dict[str, list[np.ndarray]] = {}
    for surface in surfaces.surfaces.values():
        finite = np.asarray(surface.values, dtype=float)
        finite = finite[np.isfinite(finite)]
        if finite.size:
            grouped.setdefault(surface.property_name, []).append(finite)
    limits: dict[str, tuple[float, float]] = {}
    for property_name, chunks in grouped.items():
        values = np.concatenate(chunks)
        minimum = float(np.min(values))
        maximum = float(np.max(values))
        if np.isclose(minimum, maximum):
            delta = max(abs(minimum) * 1.0e-6, 1.0e-12)
            minimum -= delta
            maximum += delta
        limits[property_name] = (minimum, maximum)
    return limits
