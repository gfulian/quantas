# -*- coding: utf-8 -*-

"""Frontend-neutral plot specifications for seismic-wave results."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from quantas.core.geometry import close_periodic_seam
from quantas.core.physics.seismic import MODE_INDEX, MODE_SYMBOLS, WaveMode
from quantas.models import (
    AxisFieldLayer,
    PlotAxis,
    PlotCollection,
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
from quantas.modules.seismic.models import SeismicResult
from quantas.modules.seismic.plot.inventory import available_seismic_plot_properties
from quantas.modules.seismic.report import sampled_variation


SurfaceType = Literal["phase", "slowness", "group"]
SurfaceGeometry = Literal["physical", "unit_sphere"]


@dataclass(frozen=True, slots=True)
class SeismicPlotOptions:
    """Options controlling neutral spherical-map preparation.

    Parameters
    ----------
    properties : tuple of str or None, optional
        Stable property identifiers to include. ``None`` selects all
        quantities available in the result.
    projection : {"equal_area", "stereographic"}, optional
        Preferred spherical-map projection.
    colormap : str or None, optional
        Portable colormap name. ``None`` selects the default seismic map.
    levels : int, optional
        Number of principal contour levels.
    polarization_stride : int, optional
        Row and column stride used to decimate polarization axes.
    include_extrema_markers : bool, optional
        Whether sampled minima and maxima are attached to map specs.
    include_polarizations : bool, optional
        Whether tracked polarization axes associated with the represented
        acoustic mode are attached when available. Shear splitting and shear
        anisotropy use the :math:`V_{S1}` axis field.
    """

    properties: tuple[str, ...] | None = None
    projection: SphericalProjection = "equal_area"
    colormap: str | None = None
    levels: int = 16
    polarization_stride: int = 6
    include_extrema_markers: bool = True
    include_polarizations: bool = False

    def __post_init__(self) -> None:
        """Validate map-preparation options."""
        if self.polarization_stride < 1:
            raise ValueError("polarization_stride must be a positive integer.")
        if self.levels < 2:
            raise ValueError("levels must be greater than or equal to two.")


@dataclass(frozen=True, slots=True)
class SeismicSurfaceOptions:
    """Options controlling neutral three-dimensional surface preparation.

    Parameters
    ----------
    properties : tuple of str or None, optional
        Scalar properties represented in three dimensions. ``None`` selects
        every property available in the result when ``surface_types`` is
        empty.
    modes : tuple of WaveMode, optional
        Acoustic modes included in explicitly requested acoustic surfaces.
    surface_types : tuple of {"phase", "slowness", "group"}, optional
        Acoustic quantities used when a particular surface family is
        requested.
    geometry : {"unit_sphere", "physical"}, optional
        Geometrical carrier. Unit-sphere geometry displays a scalar field on
        the unit sphere. Physical geometry uses phase-velocity, slowness, or
        group-wavefront coordinates whenever the selected property has a
        natural physical surface.
    colormap : str or None, optional
        Portable colormap name. ``None`` selects a property-dependent default.
    complete_antipodal_surface : bool, optional
        Whether a hemispherical result is mirrored to display the complete
        antipodal surface.
    include_polarizations : bool, optional
        Whether decimated polarization axes are attached to each surface.
    polarization_stride : int, optional
        Row and column stride used for polarization-axis decimation.
    polarization_color : str, optional
        Portable color specification for polarization axes.
    polarization_line_width : float, optional
        Requested line width for polarization axes.
    polarization_scale : float, optional
        Polarization-axis length relative to the maximum surface extent.
    """

    properties: tuple[str, ...] | None = None
    modes: tuple[WaveMode, ...] = (
        WaveMode.V_P,
        WaveMode.V_S1,
        WaveMode.V_S2,
    )
    surface_types: tuple[SurfaceType, ...] = ()
    geometry: SurfaceGeometry = "unit_sphere"
    colormap: str | None = None
    complete_antipodal_surface: bool = True
    include_polarizations: bool = False
    polarization_stride: int = 8
    polarization_color: str = "black"
    polarization_line_width: float = 1.0
    polarization_scale: float = 0.08

    def __post_init__(self) -> None:
        """Validate surface-preparation options."""
        normalized_modes = tuple(WaveMode(mode) for mode in self.modes)
        normalized_types = tuple(self.surface_types)
        allowed = {"phase", "slowness", "group"}
        if not normalized_modes:
            raise ValueError("At least one acoustic mode must be selected.")
        if any(item not in allowed for item in normalized_types):
            raise ValueError("surface_types must contain phase, slowness, or group.")
        if self.geometry not in {"unit_sphere", "physical"}:
            raise ValueError("geometry must be 'unit_sphere' or 'physical'.")
        if self.polarization_stride < 1:
            raise ValueError("polarization_stride must be a positive integer.")
        if self.polarization_line_width <= 0.0:
            raise ValueError("polarization_line_width must be positive.")
        if self.polarization_scale <= 0.0:
            raise ValueError("polarization_scale must be positive.")
        object.__setattr__(self, "modes", normalized_modes)
        object.__setattr__(self, "surface_types", normalized_types)


def build_seismic_plot_collection(
    result: SeismicResult,
    options: SeismicPlotOptions | None = None,
) -> PlotCollection:
    """Build neutral spherical-map specifications from a seismic result.

    Parameters
    ----------
    result : SeismicResult
        Complete sampled seismic-wave result.
    options : SeismicPlotOptions or None, optional
        Map selection and preparation options.

    Returns
    -------
    PlotCollection
        Spherical scalar maps, optional polarization layers, and warnings.
    """
    selected = options or SeismicPlotOptions()
    available = _available_properties(result)
    requested = selected.properties or tuple(available)
    collection = PlotCollection()
    polarization_properties = {
        key for key in requested if _polarization_mode_for_property(key) is not None
    }
    if (
        selected.include_polarizations
        and result.field.tracking is None
        and any(key in polarization_properties for key in requested)
    ):
        collection.warnings.append(
            "Tracked polarization axes are unavailable; scalar maps were "
            "prepared without a V_S1 polarization overlay."
        )
    for key in requested:
        if key not in available:
            collection.warnings.append(
                f"Seismic property {key!r} is not available in this result."
            )
            continue
        collection.plots.append(_build_property_map(result, key, selected))
    return collection


def build_seismic_summary_spec(
    result: SeismicResult,
    options: SeismicPlotOptions | None = None,
) -> SphericalSummarySpec:
    """Build a six-panel seismic summary specification.

    Parameters
    ----------
    result : SeismicResult
        Complete sampled seismic-wave result.
    options : SeismicPlotOptions or None, optional
        Map projection, colormap, extrema, and polarization options.

    Returns
    -------
    SphericalSummarySpec
        Phase velocities and principal derived diagnostics arranged in a
        compact summary.
    """
    selected = options or SeismicPlotOptions()
    keys = (
        "phase_v_p",
        "phase_v_s1",
        "phase_v_s2",
        "shear_anisotropy",
        "phase_v_p_over_v_s1",
        "phase_v_p_over_v_s2",
    )
    maps = [_build_property_map(result, key, selected) for key in keys]
    has_polarization = any(map_spec.axis_layers for map_spec in maps)
    suffix = "_polarization" if has_polarization else ""
    return SphericalSummarySpec(
        key=f"seismic_summary{suffix}",
        title="Seismic summary",
        filename_stem=f"seismic_summary{suffix}",
        maps=maps,
        columns=3,
        metadata={
            "module": "seismic",
            "direction_frame": "elastic-tensor Cartesian frame",
        },
    )


def build_seismic_surface_collection(
    result: SeismicResult,
    options: SeismicSurfaceOptions | None = None,
) -> PlotCollection:
    """Build neutral three-dimensional seismic-property specifications.

    Parameters
    ----------
    result : SeismicResult
        Complete sampled seismic-wave result.
    options : SeismicSurfaceOptions or None, optional
        Property selection, geometry, coloring, antipodal completion, and
        polarization options.

    Returns
    -------
    PlotCollection
        Three-dimensional surface specifications and non-fatal warnings.
    """
    selected = options or SeismicSurfaceOptions()
    collection = PlotCollection()
    available = _available_properties(result)

    requested_properties = selected.properties
    if requested_properties is None and not selected.surface_types:
        requested_properties = tuple(available)
    if selected.include_polarizations and selected.geometry == "physical":
        collection.warnings.append(
            "3D polarization overlays are supported only on unit-sphere "
            "geometry. Physical-geometry surfaces were prepared without "
            "polarization overlays."
        )
    for key in requested_properties or ():
        if key not in available:
            collection.warnings.append(
                f"Seismic property {key!r} is not available in this result."
            )
            continue
        collection.plots.append(_build_property_surface_spec(result, key, selected))

    for surface_type in selected.surface_types:
        if surface_type == "group" and result.field.group is None:
            collection.warnings.append(
                "Group-velocity surfaces are unavailable in this result."
            )
            continue
        for mode in selected.modes:
            collection.plots.append(
                _build_surface_spec(result, mode, surface_type, selected)
            )
    if selected.include_polarizations and result.field.tracking is None:
        collection.warnings.append(
            "Tracked polarization axes are unavailable; surfaces were "
            "prepared without polarization overlays."
        )
    return collection


def _available_properties(result: SeismicResult) -> dict[str, tuple[str, str]]:
    """Return renderer labels derived from the shared scientific inventory."""
    return {
        descriptor.key: (f"${descriptor.symbol_math}$", descriptor.unit or "")
        for descriptor in available_seismic_plot_properties(result)
    }


def _mode_math(mode: WaveMode) -> str:
    """Return the publication-style phase-speed symbol for one mode."""
    return {
        WaveMode.V_P: r"$V_P$",
        WaveMode.V_S1: r"$V_{S1}$",
        WaveMode.V_S2: r"$V_{S2}$",
    }[mode]


def _group_math(mode: WaveMode) -> str:
    """Return the publication-style group-speed symbol for one mode."""
    return {
        WaveMode.V_P: r"$V_{g,P}$",
        WaveMode.V_S1: r"$V_{g,S1}$",
        WaveMode.V_S2: r"$V_{g,S2}$",
    }[mode]


def _power_flow_math(mode: WaveMode) -> str:
    """Return the publication-style power-flow symbol for one mode."""
    return {
        WaveMode.V_P: r"$\psi_P$",
        WaveMode.V_S1: r"$\psi_{S1}$",
        WaveMode.V_S2: r"$\psi_{S2}$",
    }[mode]


def _enhancement_math(mode: WaveMode) -> str:
    """Return the publication-style enhancement symbol for one mode."""
    return {
        WaveMode.V_P: r"$\log_{10}(A_P)$",
        WaveMode.V_S1: r"$\log_{10}(A_{S1})$",
        WaveMode.V_S2: r"$\log_{10}(A_{S2})$",
    }[mode]


def _quantity_label(symbol: str, unit: str) -> str:
    """Return a symbol followed by a publication-style unit."""
    units = {
        "km s^-1": r"$\mathrm{km\,s^{-1}}$",
        "s km^-1": r"$\mathrm{s\,km^{-1}}$",
        "%": r"$\%$",
        "degree": r"$^\circ$",
        "dimensionless": "",
        "": "",
    }
    rendered = units.get(unit, unit)
    return symbol if not rendered else f"{symbol} ({rendered})"


def _default_colormap(key: str) -> str:
    """Return the default colormap family for one seismic property."""
    if key.startswith("power_flow_"):
        return "quantas_powerflow"
    if key.startswith("log10_enhancement_"):
        return "quantas_enhancement"
    return "seismic"


def _build_property_map(
    result: SeismicResult,
    key: str,
    options: SeismicPlotOptions,
) -> SphericalMapSpec:
    """Build one spherical scalar-map specification."""
    label, unit = _available_properties(result)[key]
    values, mask = _property_values(result, key)
    mapped = np.where(mask & np.isfinite(values), values, np.nan)
    grid_values = result.grid.reshape_scalar_field(mapped)
    markers = (
        _extrema_markers(values, mask, result.field.phase.directions)
        if options.include_extrema_markers
        else []
    )
    axis_layers: list[AxisFieldLayer] = []
    polarization_mode = _polarization_mode_for_property(key)
    if options.include_polarizations and polarization_mode is not None:
        layer = _polarization_axis_layer(
            result,
            polarization_mode,
            options.polarization_stride,
        )
        if layer is not None:
            axis_layers.append(layer)
    metadata = {
        "direction_frame": "elastic-tensor Cartesian frame",
        "theta_definition": "polar angle from +z",
        "phi_definition": "azimuth from +x toward +y",
        "azimuth_endpoint_included": False,
        "color_center": 0.0 if key.startswith("log10_enhancement_") else None,
    }
    polarization_suffix = "_polarization" if axis_layers else ""
    return SphericalMapSpec(
        key=key,
        title=label,
        filename_stem=f"seismic_2d_{key}{polarization_suffix}",
        theta=result.grid.theta,
        phi=result.grid.phi,
        values=grid_values,
        value_axis=PlotAxis(key=key, label=_quantity_label(label, unit), unit=unit),
        hemisphere=result.grid.hemisphere.value,
        projection=options.projection,
        colormap=options.colormap or _default_colormap(key),
        levels=options.levels,
        isolines=True,
        markers=markers,
        axis_layers=axis_layers,
        metadata=metadata,
    )


def _polarization_mode_for_property(key: str) -> WaveMode | None:
    """Return the polarization mode naturally associated with a property.

    Parameters
    ----------
    key : str
        Stable seismic property identifier.

    Returns
    -------
    WaveMode or None
        Associated acoustic mode, or ``None`` when no single polarization
        field is physically distinguished.
    """
    mode = _property_mode(key)
    if mode is not None:
        return mode
    if key in {"shear_anisotropy", "shear_splitting"}:
        return WaveMode.V_S1
    return None


def _property_values(
    result: SeismicResult,
    key: str,
) -> tuple[NDArray[np.float64], NDArray[np.bool_]]:
    """Return flat values and an eligibility mask for one map property."""
    phase = result.field.phase
    if key.startswith("phase_v_") and "over" not in key:
        mode = WaveMode(key.removeprefix("phase_"))
        phase_mode = phase.for_mode(mode)
        return phase_mode.phase_speeds, phase_mode.valid_mask

    i_s2 = MODE_INDEX[WaveMode.V_S2]
    i_s1 = MODE_INDEX[WaveMode.V_S1]
    i_p = MODE_INDEX[WaveMode.V_P]
    speeds = phase.phase_speeds
    if key == "shear_anisotropy":
        mask = phase.valid_mask[:, i_s1] & phase.valid_mask[:, i_s2]
        values = _safe_divide(
            200.0 * (speeds[:, i_s1] - speeds[:, i_s2]),
            speeds[:, i_s1] + speeds[:, i_s2],
        )
        return values, mask
    if key == "shear_splitting":
        mask = phase.valid_mask[:, i_s1] & phase.valid_mask[:, i_s2]
        return speeds[:, i_s1] - speeds[:, i_s2], mask
    if key == "phase_v_p_over_v_s1":
        mask = phase.valid_mask[:, i_p] & phase.valid_mask[:, i_s1]
        return _safe_divide(speeds[:, i_p], speeds[:, i_s1]), mask
    if key == "phase_v_p_over_v_s2":
        mask = phase.valid_mask[:, i_p] & phase.valid_mask[:, i_s2]
        return _safe_divide(speeds[:, i_p], speeds[:, i_s2]), mask

    if key.startswith("group_"):
        assert result.field.group is not None
        mode = WaveMode(key.removeprefix("group_"))
        group_mode = result.field.group.for_mode(mode)
        return (
            group_mode.group_speeds,
            group_mode.valid_mask & group_mode.resolved_mask,
        )
    if key.startswith("power_flow_"):
        assert result.field.group is not None
        mode = WaveMode(key.removeprefix("power_flow_"))
        group_mode = result.field.group.for_mode(mode)
        return (
            np.degrees(group_mode.power_flow_angles),
            group_mode.valid_mask & group_mode.resolved_mask,
        )
    if key.startswith("log10_enhancement_"):
        assert result.field.enhancement is not None
        mode = WaveMode(key.removeprefix("log10_enhancement_"))
        enhancement_mode = result.field.enhancement.for_mode(mode)
        mask = (
            enhancement_mode.valid_mask
            & enhancement_mode.resolved_mask
            & enhancement_mode.finite_mask
        )
        return enhancement_mode.log10_enhancement, mask
    raise KeyError(f"Unknown seismic plot property: {key}")


def _extrema_markers(
    values: NDArray[np.float64],
    mask: NDArray[np.bool_],
    directions: NDArray[np.float64],
) -> list[SphericalMarker]:
    """Build minimum and maximum markers for a scalar map."""
    variation = sampled_variation(values, mask, directions)
    markers: list[SphericalMarker] = []
    for key, label, extremum, symbol in (
        ("minimum", "Minimum", variation.minimum, "circle"),
        ("maximum", "Maximum", variation.maximum, "square"),
    ):
        if extremum.direction is None:
            continue
        markers.append(
            SphericalMarker(
                key=key,
                label=label,
                directions=np.asarray(extremum.direction, dtype=float)[None, :],
                marker=symbol,
                metadata={
                    "value": extremum.value,
                    "multiplicity": extremum.multiplicity,
                    "antipodal": True,
                },
            )
        )
    return markers


def _polarization_axis_layer(
    result: SeismicResult,
    mode: WaveMode,
    stride: int,
) -> AxisFieldLayer | None:
    """Return a decimated sign-aligned polarization-axis layer."""
    selected, directions, axes = _tracked_mode_axes(result, mode, stride)
    if selected.size == 0:
        return None
    mask = np.ones(selected.size, dtype=bool)
    for array in (directions, axes, mask):
        array.setflags(write=False)
    return AxisFieldLayer(
        key=f"{mode.value}_polarization",
        label=f"{MODE_SYMBOLS[mode]} polarization axis",
        directions=directions,
        axes=axes,
        resolved_mask=mask,
        metadata={
            "axial": True,
            "antipodal": True,
            "stride": stride,
            "source": "tracked polarization field",
        },
    )


def _tracked_mode_axes(
    result: SeismicResult,
    mode: WaveMode,
    stride: int,
) -> tuple[NDArray[np.int64], NDArray[np.float64], NDArray[np.float64]]:
    """Return decimated directions and tracked axes for one local mode."""
    tracking = result.field.tracking
    if tracking is None:
        return (
            np.empty(0, dtype=np.int64),
            np.empty((0, 3), dtype=float),
            np.empty((0, 3), dtype=float),
        )
    target = MODE_INDEX[mode]
    matches = tracking.branch_mode_indices == target
    has_match = np.any(matches, axis=1)
    branch = np.argmax(matches, axis=1)
    point = np.arange(tracking.directions.shape[0])
    axes = tracking.polarizations[point, branch]
    resolved = tracking.resolved_mask[point, branch] & has_match

    index_grid = np.arange(result.grid.size).reshape(result.grid.shape)
    theta_stride = min(stride, max(1, result.grid.shape[0] // 4))
    phi_stride = min(stride, max(1, result.grid.shape[1] // 8))
    sampled = index_grid[::theta_stride, ::phi_stride].reshape(-1)
    selected: NDArray[np.int64] = np.asarray(
        sampled[resolved[sampled] & np.all(np.isfinite(axes[sampled]), axis=1)],
        dtype=np.int64,
    )
    if selected.size == 0:
        candidates = np.flatnonzero(resolved & np.all(np.isfinite(axes), axis=1))
        if candidates.size == 0:
            return (
                np.empty(0, dtype=np.int64),
                np.empty((0, 3), dtype=float),
                np.empty((0, 3), dtype=float),
            )
        fallback_stride = max(1, int(np.ceil(candidates.size / 200)))
        selected = np.asarray(
            candidates[::fallback_stride],
            dtype=np.int64,
        )
    directions = np.array(tracking.directions[selected], dtype=float, copy=True)
    axis_values = np.array(axes[selected], dtype=float, copy=True)
    return selected, directions, axis_values


def _build_property_surface_spec(
    result: SeismicResult,
    key: str,
    options: SeismicSurfaceOptions,
) -> SurfacePlotSpec:
    """Build one scalar property on a unit or natural physical surface."""
    label, unit = _available_properties(result)[key]
    values, mask = _property_values(result, key)
    mode = _property_mode(key)
    coordinates, geometry_name, geometry_unit = _property_surface_coordinates(
        result,
        key,
        mode,
        options.geometry,
    )
    valid = mask & np.isfinite(values) & np.all(np.isfinite(coordinates), axis=1)
    clean_values = np.where(valid, values, np.nan)
    clean_coordinates = np.where(valid[:, None], coordinates, np.nan)
    filename = f"seismic_3d_{key}_{options.geometry}"
    layers = _generic_surface_mesh_layers(
        result,
        clean_coordinates,
        clean_values,
        key,
        label,
        options,
    )
    vectors: list[VectorFieldLayer] = []
    if (
        options.include_polarizations
        and options.geometry == "unit_sphere"
        and mode is not None
    ):
        layer = _surface_polarization_layer(
            result,
            mode,
            coordinates,
            valid,
            options,
        )
        if layer is not None:
            vectors.append(layer)
    polarization_suffix = "_polarization" if vectors else ""
    return SurfacePlotSpec(
        key=f"property_{key}_{options.geometry}",
        title=label,
        filename_stem=f"{filename}{polarization_suffix}",
        value_axis=PlotAxis(
            key=key,
            label=_quantity_label(label, unit),
            unit=unit or None,
        ),
        layers=layers,
        vector_layers=vectors,
        equal_aspect=True,
        show_axes=True,
        metadata={
            "module": "seismic",
            "property": key,
            "geometry": geometry_name,
            "geometry_unit": geometry_unit,
            "mode": None if mode is None else mode.value,
            "direction_frame": "elastic-tensor Cartesian frame",
            "complete_antipodal_surface": options.complete_antipodal_surface,
            "axis_labels": ("x", "y", "z"),
            "color_center": 0.0 if key.startswith("log10_enhancement_") else None,
        },
    )


def _property_mode(key: str) -> WaveMode | None:
    """Return the acoustic mode encoded in a mode-resolved property key."""
    prefixes = ("phase_", "group_", "power_flow_", "log10_enhancement_")
    for prefix in prefixes:
        if not key.startswith(prefix):
            continue
        suffix = key.removeprefix(prefix)
        try:
            return WaveMode(suffix)
        except ValueError:
            return None
    return None


def _property_surface_coordinates(
    result: SeismicResult,
    key: str,
    mode: WaveMode | None,
    geometry: SurfaceGeometry,
) -> tuple[NDArray[np.float64], str, str]:
    """Return coordinates for a scalar property surface."""
    directions = result.field.phase.directions
    if geometry == "unit_sphere" or mode is None:
        return directions, "unit_sphere", "dimensionless"
    index = MODE_INDEX[mode]
    if key.startswith("phase_"):
        radius = result.field.phase.phase_speeds[:, index]
        return directions * radius[:, None], "phase_velocity", "km s^-1"
    if key.startswith(("group_", "power_flow_", "log10_enhancement_")):
        if result.field.group is None:
            return directions, "unit_sphere", "dimensionless"
        return (
            result.field.group.group_velocities[:, index, :],
            "group_wavefront",
            "km s^-1",
        )
    return directions, "unit_sphere", "dimensionless"


def _generic_surface_mesh_layers(
    result: SeismicResult,
    coordinates: NDArray[np.float64],
    values: NDArray[np.float64],
    key: str,
    label: str,
    options: SeismicSurfaceOptions,
) -> list[SurfaceLayer]:
    """Return closed meshes for a generic scalar property."""
    shape = result.grid.shape
    xyz = coordinates.reshape((*shape, 3))
    mapped = values.reshape(shape)
    colormap = options.colormap or _default_colormap(key)
    style = SurfaceStyle(colormap=colormap, opacity=0.96, show_colorbar=True)
    layers = [
        SurfaceLayer(
            key=f"{key}_sampled",
            label=label,
            x=close_periodic_seam(xyz[..., 0]),
            y=close_periodic_seam(xyz[..., 1]),
            z=close_periodic_seam(xyz[..., 2]),
            values=close_periodic_seam(mapped),
            theta=close_periodic_seam(result.grid.theta_grid),
            phi=close_periodic_seam(result.grid.phi_grid),
            style=style,
            metadata={"antipodal_copy": False},
        )
    ]
    if options.complete_antipodal_surface and result.grid.hemisphere.value != "full":
        layers.append(
            SurfaceLayer(
                key=f"{key}_antipodal",
                label=label,
                x=close_periodic_seam(-xyz[..., 0]),
                y=close_periodic_seam(-xyz[..., 1]),
                z=close_periodic_seam(-xyz[..., 2]),
                values=close_periodic_seam(mapped),
                style=SurfaceStyle(
                    colormap=colormap,
                    opacity=0.96,
                    show_colorbar=False,
                ),
                metadata={"antipodal_copy": True},
            )
        )
    return layers


def _build_surface_spec(
    result: SeismicResult,
    mode: WaveMode,
    surface_type: SurfaceType,
    options: SeismicSurfaceOptions,
) -> SurfacePlotSpec:
    """Build one acoustic surface with optional polarization axes."""
    index = MODE_INDEX[mode]
    phase = result.field.phase
    directions = phase.directions
    valid = phase.valid_mask[:, index]

    if surface_type == "phase":
        values = phase.phase_speeds[:, index]
        coordinates = (
            directions
            if options.geometry == "unit_sphere"
            else directions * values[:, None]
        )
        label = _mode_math(mode)
        unit = "km s^-1"
        filename = f"seismic_3d_phase_{mode.value}_{options.geometry}"
    elif surface_type == "slowness":
        values = _safe_divide(np.ones(phase.n_points), phase.phase_speeds[:, index])
        coordinates = (
            directions
            if options.geometry == "unit_sphere"
            else directions * values[:, None]
        )
        label = {
            WaveMode.V_P: r"$s_P$",
            WaveMode.V_S1: r"$s_{S1}$",
            WaveMode.V_S2: r"$s_{S2}$",
        }[mode]
        unit = "s km^-1"
        filename = f"seismic_3d_slowness_{mode.value}_{options.geometry}"
    else:
        assert result.field.group is not None
        group = result.field.group
        values = group.group_speeds[:, index]
        coordinates = (
            directions
            if options.geometry == "unit_sphere"
            else group.group_velocities[:, index, :]
        )
        valid = valid & group.valid_mask[:, index] & group.resolved_mask[:, index]
        label = _group_math(mode)
        unit = "km s^-1"
        filename = f"seismic_3d_group_{mode.value}_{options.geometry}"

    clean_values = np.where(valid & np.isfinite(values), values, np.nan)
    clean_coordinates = np.where(
        (valid & np.all(np.isfinite(coordinates), axis=1))[:, None],
        coordinates,
        np.nan,
    )
    layers = _surface_mesh_layers(
        result,
        clean_coordinates,
        clean_values,
        mode,
        surface_type,
        options,
    )
    vectors: list[VectorFieldLayer] = []
    if options.include_polarizations and options.geometry == "unit_sphere":
        layer = _surface_polarization_layer(
            result,
            mode,
            coordinates,
            valid,
            options,
        )
        if layer is not None:
            vectors.append(layer)
    polarization_suffix = "_polarization" if vectors else ""
    return SurfacePlotSpec(
        key=f"{surface_type}_{mode.value}",
        title=label,
        filename_stem=f"{filename}{polarization_suffix}",
        value_axis=PlotAxis(
            key=f"{surface_type}_{mode.value}",
            label=_quantity_label(label, unit),
            unit=unit,
        ),
        layers=layers,
        vector_layers=vectors,
        equal_aspect=True,
        show_axes=True,
        metadata={
            "module": "seismic",
            "surface_type": surface_type,
            "geometry": options.geometry,
            "mode": mode.value,
            "mode_symbol": MODE_SYMBOLS[mode],
            "direction_frame": "elastic-tensor Cartesian frame",
            "complete_antipodal_surface": options.complete_antipodal_surface,
            "axis_labels": ("x", "y", "z"),
        },
    )


def _surface_mesh_layers(
    result: SeismicResult,
    coordinates: NDArray[np.float64],
    values: NDArray[np.float64],
    mode: WaveMode,
    surface_type: SurfaceType,
    options: SeismicSurfaceOptions,
) -> list[SurfaceLayer]:
    """Return one or two closed surface meshes."""
    shape = result.grid.shape
    xyz = coordinates.reshape((*shape, 3))
    mapped = values.reshape(shape)
    style = SurfaceStyle(
        colormap=options.colormap or "seismic",
        opacity=0.96,
        show_colorbar=True,
    )
    layers = [
        SurfaceLayer(
            key=f"{surface_type}_{mode.value}_sampled",
            label=f"{MODE_SYMBOLS[mode]} sampled hemisphere",
            x=close_periodic_seam(xyz[..., 0]),
            y=close_periodic_seam(xyz[..., 1]),
            z=close_periodic_seam(xyz[..., 2]),
            values=close_periodic_seam(mapped),
            theta=close_periodic_seam(result.grid.theta_grid),
            phi=close_periodic_seam(result.grid.phi_grid),
            style=style,
            metadata={"antipodal_copy": False},
        )
    ]
    if options.complete_antipodal_surface and result.grid.hemisphere.value != "full":
        layers.append(
            SurfaceLayer(
                key=f"{surface_type}_{mode.value}_antipodal",
                label=f"{MODE_SYMBOLS[mode]} antipodal hemisphere",
                x=close_periodic_seam(-xyz[..., 0]),
                y=close_periodic_seam(-xyz[..., 1]),
                z=close_periodic_seam(-xyz[..., 2]),
                values=close_periodic_seam(mapped),
                style=SurfaceStyle(
                    colormap=options.colormap or "seismic",
                    opacity=0.96,
                    show_colorbar=False,
                ),
                metadata={"antipodal_copy": True},
            )
        )
    return layers


def _surface_polarization_layer(
    result: SeismicResult,
    mode: WaveMode,
    coordinates: NDArray[np.float64],
    valid: NDArray[np.bool_],
    options: SeismicSurfaceOptions,
) -> VectorFieldLayer | None:
    """Build one decimated Cartesian polarization overlay."""
    selected, _directions, axes = _tracked_mode_axes(
        result,
        mode,
        options.polarization_stride,
    )
    if selected.size == 0:
        return None
    eligible = valid[selected] & np.all(np.isfinite(coordinates[selected]), axis=1)
    selected = selected[eligible]
    axes = axes[eligible]
    if selected.size == 0:
        return None
    origins = np.array(coordinates[selected], dtype=float, copy=True)
    vectors = np.array(axes, dtype=float, copy=True)
    if options.complete_antipodal_surface and result.grid.hemisphere.value != "full":
        origins = np.concatenate((origins, -origins), axis=0)
        vectors = np.concatenate((vectors, vectors), axis=0)
    mask = np.ones(origins.shape[0], dtype=bool)
    for array in (origins, vectors, mask):
        array.setflags(write=False)
    return VectorFieldLayer(
        key=f"{mode.value}_polarization",
        label=f"{MODE_SYMBOLS[mode]} polarization axis",
        origins=origins,
        vectors=vectors,
        axial=True,
        resolved_mask=mask,
        style=VectorFieldStyle(
            color=options.polarization_color,
            line_width=options.polarization_line_width,
            scale=options.polarization_scale,
            opacity=0.9,
        ),
        metadata={
            "mode": mode.value,
            "antipodal": True,
            "source": "tracked polarization field",
        },
    )


def _safe_divide(
    numerator: NDArray[np.float64],
    denominator: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Return a finite-denominator division with NaN elsewhere."""
    result = np.full(np.broadcast_shapes(numerator.shape, denominator.shape), np.nan)
    np.divide(
        numerator,
        denominator,
        out=result,
        where=np.isfinite(denominator) & (denominator != 0.0),
    )
    return result
