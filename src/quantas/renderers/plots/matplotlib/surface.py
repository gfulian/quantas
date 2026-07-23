# -*- coding: utf-8 -*-

"""Three-dimensional surface renderer for Matplotlib."""

from __future__ import annotations

from typing import Any, Literal

import numpy as np

from quantas.models.plot import SurfacePlotSpec

from .colors import resolve_colormap
from .common import SURFACE_AZIMUTH, SURFACE_ELEVATION
from .options import MatplotlibOptions


def _render_surface_plot(
    plt: Any,
    spec: SurfacePlotSpec,
    options: MatplotlibOptions,
) -> Any:
    """Render a smooth three-dimensional surface with custom frame axes.

    Parameters
    ----------
    plt : module
        Imported :mod:`matplotlib.pyplot` module.
    spec : SurfacePlotSpec
        Neutral three-dimensional plot specification.
    options : MatplotlibOptions
        Renderer and file-output options.

    Returns
    -------
    matplotlib.figure.Figure
        Rendered figure.
    """
    from matplotlib import cm, colors

    figure = plt.figure(figsize=(9.2, 7.0))
    axis = figure.add_axes((0.01, 0.02, 0.80, 0.96), projection="3d")
    maximum_extent = 0.0
    positive_extents = np.zeros(3, dtype=float)
    colorbar_payload: tuple[Any, np.ndarray] | None = None

    for layer in spec.layers:
        x = np.asarray(layer.x, dtype=float)
        y = np.asarray(layer.y, dtype=float)
        z = np.asarray(layer.z, dtype=float)
        values = np.asarray(layer.values, dtype=float)
        coordinates = np.stack((x, y, z))
        finite_coordinates = coordinates[np.isfinite(coordinates)]
        if finite_coordinates.size == 0:
            continue
        maximum_extent = max(maximum_extent, float(np.nanmax(np.abs(coordinates))))
        for index, component in enumerate((x, y, z)):
            finite_component = component[np.isfinite(component)]
            if finite_component.size:
                positive_extents[index] = max(
                    positive_extents[index],
                    float(np.max(finite_component)),
                )

        kwargs: dict[str, Any] = {
            "linewidth": layer.style.mesh_line_width if layer.style.show_mesh else 0.0,
            "antialiased": bool(layer.style.show_mesh),
            "alpha": layer.style.opacity,
            "shade": False,
            "rcount": x.shape[0],
            "ccount": x.shape[1],
            "edgecolor": layer.style.mesh_color if layer.style.show_mesh else "none",
        }
        if layer.style.colormap is not None:
            finite = values[np.isfinite(values)]
            if finite.size == 0:
                continue
            limits = layer.style.value_limits
            vmin = float(np.min(finite)) if limits is None else limits[0]
            vmax = float(np.max(finite)) if limits is None else limits[1]
            norm: Any
            center = spec.metadata.get("color_center")
            if center is not None and vmin < float(center) < vmax:
                norm = colors.TwoSlopeNorm(
                    vmin=vmin,
                    vcenter=float(center),
                    vmax=vmax,
                )
            else:
                norm = colors.Normalize(vmin=vmin, vmax=vmax)
            cmap = resolve_colormap(plt, layer.style.colormap)
            mapped = np.where(np.isfinite(values), values, vmin)
            kwargs["facecolors"] = cmap(norm(mapped))
            axis.plot_surface(x, y, z, **kwargs)
            if layer.style.show_colorbar and colorbar_payload is None:
                scalar = cm.ScalarMappable(norm=norm, cmap=cmap)
                scalar.set_array(finite)
                colorbar_payload = scalar, finite
        else:
            if layer.style.color is not None:
                kwargs["color"] = layer.style.color
            axis.plot_surface(x, y, z, **kwargs)

    if spec.vector_layers and maximum_extent > 0.0:
        _render_cartesian_vector_layers(axis, spec, options, maximum_extent)

    if maximum_extent > 0.0:
        metadata = spec.metadata
        if spec.show_axes and metadata.get("axis_draw_mode") == "full_equal":
            factor = float(metadata.get("axis_half_length_factor", 1.14))
            label_factor = float(metadata.get("axis_label_factor", 1.07))
            half_length = factor * maximum_extent
            limit = max(
                half_length * (1.0 + 0.09 * label_factor), 1.09 * maximum_extent
            )
            lower = -limit
            upper = limit
        else:
            lower = -1.06 * maximum_extent
            upper = (1.48 if spec.show_axes else 1.06) * maximum_extent
        axis.set_xlim(lower, upper)
        axis.set_ylim(lower, upper)
        axis.set_zlim(lower, upper)
        if spec.equal_aspect:
            axis.set_box_aspect((1.0, 1.0, 1.0), zoom=1.48)
    axis.set_axis_off()
    if spec.show_axes:
        _decorate_cartesian_axes(
            axis,
            positive_extents,
            maximum_extent,
            options.axis_label_mode,
            metadata=spec.metadata,
        )
    axis.view_init(elev=SURFACE_ELEVATION, azim=SURFACE_AZIMUTH)
    axis.set_proj_type("persp", focal_length=0.92)
    if options.show_title and spec.title:
        figure.suptitle(spec.title, fontsize=15, y=0.965)

    if colorbar_payload is not None:
        scalar, _finite = colorbar_payload
        colorbar_axis = figure.add_axes((0.84, 0.14, 0.026, 0.72))
        colorbar = figure.colorbar(scalar, cax=colorbar_axis)
        colorbar.set_label(spec.value_axis.label, fontsize=13, labelpad=9)
        colorbar.ax.tick_params(labelsize=11)

    return figure


def _render_cartesian_vector_layers(
    axis: Any,
    spec: SurfacePlotSpec,
    options: MatplotlibOptions,
    maximum_extent: float,
) -> None:
    """Render vector and axial overlays on a Cartesian surface."""
    for layer in spec.vector_layers:
        origins = np.asarray(layer.origins, dtype=float)
        vectors = np.asarray(layer.vectors, dtype=float)
        mask = np.all(np.isfinite(origins), axis=1) & np.all(
            np.isfinite(vectors), axis=1
        )
        if layer.resolved_mask is not None:
            mask &= np.asarray(layer.resolved_mask, dtype=bool)
        origins = origins[mask]
        vectors = vectors[mask]
        if origins.size == 0:
            continue
        if spec.metadata.get("geometry") == "unit_sphere":
            visible = _unit_sphere_front_side_mask(origins)
            origins = origins[visible]
            vectors = vectors[visible]
            if origins.size == 0:
                continue
        norms = np.linalg.norm(vectors, axis=1)
        usable = norms > 0.0
        origins = origins[usable]
        vectors = vectors[usable] / norms[usable, None]
        if origins.size == 0:
            continue
        color = options.polarization_color or layer.style.color or "black"
        line_width = (
            options.polarization_line_width
            if options.polarization_line_width is not None
            else layer.style.line_width
        )
        scale_fraction = (
            options.polarization_scale
            if options.polarization_scale is not None
            else layer.style.scale
        )
        length = scale_fraction * maximum_extent
        if layer.axial:
            half = 0.5 * length * vectors
            from mpl_toolkits.mplot3d.art3d import Line3DCollection

            segments = np.stack((origins - half, origins + half), axis=1)
            collection = Line3DCollection(
                segments,
                colors=color,
                linewidths=line_width,
                alpha=layer.style.opacity,
            )
            axis.add_collection3d(collection)
        else:
            axis.quiver(
                origins[:, 0],
                origins[:, 1],
                origins[:, 2],
                vectors[:, 0],
                vectors[:, 1],
                vectors[:, 2],
                length=length,
                normalize=True,
                color=color,
                linewidth=line_width,
                alpha=layer.style.opacity,
            )


def _unit_sphere_front_side_mask(origins: np.ndarray) -> np.ndarray:
    """Return origins located on the visible side of the 3D unit sphere.

    Parameters
    ----------
    origins : ndarray
        Cartesian origins of the vector or axial segments.

    Returns
    -------
    ndarray
        Boolean mask selecting origins whose radial direction points toward
        the current fixed camera direction.
    """
    norms = np.linalg.norm(origins, axis=1)
    valid = np.isfinite(norms) & (norms > 0.0)
    unit = np.zeros_like(origins, dtype=float)
    unit[valid] = origins[valid] / norms[valid, None]
    view = _surface_view_direction()
    visibility = np.einsum("ij,j->i", unit, view)
    return valid & (visibility >= -1.0e-3)


def _surface_view_direction() -> np.ndarray:
    """Return the unit vector pointing from the origin to the camera."""
    elevation = np.radians(SURFACE_ELEVATION)
    azimuth = np.radians(SURFACE_AZIMUTH)
    direction = np.array(
        [
            np.cos(elevation) * np.cos(azimuth),
            np.cos(elevation) * np.sin(azimuth),
            np.sin(elevation),
        ],
        dtype=float,
    )
    return direction / np.linalg.norm(direction)


def _decorate_cartesian_axes(
    axis: Any,
    surface_extents: np.ndarray,
    maximum_extent: float,
    label_mode: Literal["cartesian", "crystal"],
    *,
    metadata: dict[str, Any] | None = None,
) -> None:
    """Draw clean Cartesian axes without Matplotlib's 3D box."""
    if maximum_extent <= 0.0:
        return
    settings = {} if metadata is None else metadata
    labels = _cartesian_axis_labels(label_mode, settings)
    if settings.get("axis_draw_mode") == "full_equal":
        _draw_full_cartesian_axes(axis, labels, maximum_extent, settings)
        return
    starts, ends, labels_at = _cartesian_axis_positions(
        surface_extents,
        maximum_extent,
        settings,
    )
    origins = np.diag(starts)
    vectors = np.diag(ends - starts)
    axis.quiver(
        origins[:, 0],
        origins[:, 1],
        origins[:, 2],
        vectors[:, 0],
        vectors[:, 1],
        vectors[:, 2],
        color=("black", "black", "black"),
        arrow_length_ratio=0.18,
        linewidth=1.8,
    )
    _draw_cartesian_axis_labels(axis, labels, labels_at)


def _draw_full_cartesian_axes(
    axis: Any,
    labels: tuple[str, str, str],
    maximum_extent: float,
    metadata: dict[str, Any],
) -> None:
    """Draw full X/Y/Z axes with equal length through the origin."""
    half_length = float(metadata.get("axis_half_length_factor", 1.14)) * maximum_extent
    label_factor = float(metadata.get("axis_label_factor", 1.07))
    line_width = float(metadata.get("axis_line_width", 1.8))
    coordinates = (
        ((-half_length, half_length), (0.0, 0.0), (0.0, 0.0)),
        ((0.0, 0.0), (-half_length, half_length), (0.0, 0.0)),
        ((0.0, 0.0), (0.0, 0.0), (-half_length, half_length)),
    )
    for xs, ys, zs in coordinates:
        axis.plot(
            xs, ys, zs, color="black", linewidth=line_width, solid_capstyle="round"
        )
    positions = np.full(3, half_length * label_factor, dtype=float)
    _draw_cartesian_axis_labels(axis, labels, positions)


def _cartesian_axis_labels(
    label_mode: Literal["cartesian", "crystal"],
    metadata: dict[str, Any],
) -> tuple[str, str, str]:
    """Resolve Cartesian or crystallographic labels for a surface."""
    default_labels = (
        (r"$x$", r"$y$", r"$z$")
        if label_mode == "cartesian"
        else (r"$[100]$", r"$[010]$", r"$[001]$")
    )
    if label_mode == "crystal":
        return default_labels
    labels = tuple(metadata.get("axis_labels", default_labels))
    if len(labels) != 3:
        raise ValueError("Surface axis_labels metadata must contain three labels.")
    return labels  # type: ignore[return-value]


def _cartesian_axis_positions(
    surface_extents: np.ndarray,
    maximum_extent: float,
    metadata: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return axis origins, arrow ends, and label positions."""
    starts = np.where(
        np.isfinite(surface_extents) & (surface_extents > 0.0),
        surface_extents,
        maximum_extent,
    )
    if metadata.get("axis_anchor_mode", "global") != "surface_relative":
        return (
            starts,
            np.full(3, 1.37 * maximum_extent, dtype=float),
            np.full(3, 1.45 * maximum_extent, dtype=float),
        )

    extension = float(metadata.get("axis_extension_fraction", 0.12))
    label_extension = float(metadata.get("axis_label_extension_fraction", 0.19))
    minimum_length = 0.04 * maximum_extent
    ends = np.maximum(starts * (1.0 + extension), starts + minimum_length)
    labels_at = np.maximum(
        starts * (1.0 + label_extension),
        ends + 0.03 * maximum_extent,
    )
    return starts, ends, labels_at


def _draw_cartesian_axis_labels(
    axis: Any,
    labels: tuple[str, str, str],
    positions: np.ndarray,
) -> None:
    """Draw X/Y/Z labels without restoring Matplotlib's 3D box."""
    label_box = {
        "facecolor": "white",
        "edgecolor": "none",
        "alpha": 0.72,
        "pad": 0.35,
    }
    coordinates = (
        (positions[0], 0.0, 0.0, "right", "center"),
        (0.0, positions[1], 0.0, "left", "center"),
        (0.0, 0.0, positions[2], "center", "bottom"),
    )
    for label, (x, y, z, horizontal, vertical) in zip(
        labels,
        coordinates,
        strict=True,
    ):
        axis.text(
            x,
            y,
            z,
            label,
            fontsize=16,
            weight="bold",
            ha=horizontal,
            va=vertical,
            bbox=label_box,
        )


__all__ = ["_render_surface_plot"]
