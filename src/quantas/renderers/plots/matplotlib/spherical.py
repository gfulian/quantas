# -*- coding: utf-8 -*-

"""Projected spherical-map renderers for Matplotlib."""

from __future__ import annotations

from typing import Any, Literal

import numpy as np

from quantas.models.plot import SphericalMapSpec, SphericalSummarySpec

from .colors import resolve_colormap
from .options import MatplotlibOptions


def _render_spherical_map(
    plt: Any,
    spec: SphericalMapSpec,
    options: MatplotlibOptions,
) -> Any:
    """Render a projected spherical scalar map."""
    panels = _spherical_panels(spec)
    figure, axes = plt.subplots(
        1,
        len(panels),
        figsize=(7.1 * len(panels), 6.7),
        squeeze=False,
    )
    for axis, (hemisphere, theta, values) in zip(axes.ravel(), panels, strict=True):
        _render_spherical_panel(
            figure,
            axis,
            spec,
            hemisphere,
            theta,
            values,
            options,
            compact=False,
        )
    figure.subplots_adjust(left=0.025, right=0.975, bottom=0.04, top=0.98, wspace=0.08)
    return figure


def _render_spherical_summary(
    plt: Any,
    spec: SphericalSummarySpec,
    options: MatplotlibOptions,
) -> Any:
    """Render a compact multi-panel spherical-map summary."""
    if not spec.maps:
        raise ValueError("a spherical summary requires at least one map")
    columns = max(1, int(spec.columns))
    rows = int(np.ceil(len(spec.maps) / columns))
    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=(5.45 * columns, 5.0 * rows),
        squeeze=False,
    )
    flat = axes.ravel()
    for axis, map_spec in zip(flat, spec.maps, strict=False):
        panels = _spherical_panels(map_spec)
        preferred = next(
            (panel for panel in panels if panel[0] == "upper"),
            panels[0],
        )
        hemisphere, theta, values = preferred
        _render_spherical_panel(
            figure,
            axis,
            map_spec,
            hemisphere,
            theta,
            values,
            options,
            compact=True,
        )
    for axis in flat[len(spec.maps) :]:
        axis.set_visible(False)
    figure.subplots_adjust(
        left=0.025, right=0.975, bottom=0.035, top=0.985, wspace=0.06, hspace=0.12
    )
    return figure


def _spherical_panels(
    spec: SphericalMapSpec,
) -> list[tuple[Literal["upper", "lower"], np.ndarray, np.ndarray]]:
    """Split a spherical specification into renderable hemispheres."""
    theta = np.asarray(spec.theta, dtype=float)
    values = np.asarray(spec.values, dtype=float)
    if spec.hemisphere == "upper":
        return [("upper", theta, values)]
    if spec.hemisphere == "lower":
        return [("lower", theta, values)]
    tolerance = 16.0 * np.finfo(float).eps
    upper = theta <= 0.5 * np.pi + tolerance
    lower = theta >= 0.5 * np.pi - tolerance
    return [
        ("upper", theta[upper], values[upper]),
        ("lower", theta[lower], values[lower]),
    ]


def _render_spherical_panel(
    figure: Any,
    axis: Any,
    spec: SphericalMapSpec,
    hemisphere: Literal["upper", "lower"],
    theta: np.ndarray,
    values: np.ndarray,
    options: MatplotlibOptions,
    *,
    compact: bool,
) -> None:
    """Render one projected hemisphere into an existing Matplotlib axis."""
    from matplotlib import colors
    import matplotlib.pyplot as plt

    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise ValueError(f"spherical map {spec.key!r} contains no finite values")
    phi = np.asarray(spec.phi, dtype=float)
    phi_closed = np.concatenate((phi, [2.0 * np.pi]))
    values_closed = np.concatenate((values, values[:, :1]), axis=1)
    theta_grid, phi_grid = np.meshgrid(theta, phi_closed, indexing="ij")
    x, y = _project_angles(
        theta_grid,
        phi_grid,
        hemisphere=hemisphere,
        projection=spec.projection,
    )
    levels = _mapped_levels(finite, spec.levels)
    center = spec.metadata.get("color_center")
    norm = None
    if center is not None:
        vmin = float(np.nanmin(finite))
        vmax = float(np.nanmax(finite))
        if vmin < float(center) < vmax:
            norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=float(center), vmax=vmax)
    filled = axis.contourf(
        x,
        y,
        values_closed,
        levels=max(64, spec.levels * 6),
        cmap=resolve_colormap(plt, spec.colormap),
        norm=norm,
        antialiased=True,
    )
    if spec.isolines and float(np.nanmax(finite) - np.nanmin(finite)) > 0.0:
        axis.contour(
            x,
            y,
            values_closed,
            levels=levels,
            colors="black",
            linewidths=0.45,
            alpha=0.55,
        )
    colorbar = figure.colorbar(
        filled,
        ax=axis,
        orientation="horizontal",
        fraction=0.040 if compact else 0.042,
        pad=0.038 if compact else 0.042,
        shrink=0.78 if compact else 0.74,
        aspect=30,
    )
    colorbar.set_label(
        spec.value_axis.label, fontsize=10 if compact else 12, labelpad=7
    )
    colorbar.ax.tick_params(labelsize=9 if compact else 10)

    _decorate_spherical_axis(
        axis,
        hemisphere,
        options.axis_label_mode,
        compact=compact,
    )
    _render_spherical_markers(
        axis,
        spec,
        hemisphere,
        options,
        compact=compact,
    )
    _render_axis_fields(
        axis,
        spec,
        hemisphere,
        options,
        compact=compact,
    )


def _project_angles(
    theta: np.ndarray,
    phi: np.ndarray,
    *,
    hemisphere: Literal["upper", "lower"],
    projection: Literal["equal_area", "stereographic"],
) -> tuple[np.ndarray, np.ndarray]:
    """Project polar coordinates to a unit-disk hemisphere map."""
    alpha = theta if hemisphere == "upper" else np.pi - theta
    if projection == "equal_area":
        radius = np.sqrt(2.0) * np.sin(0.5 * alpha)
    else:
        radius = np.tan(0.5 * alpha)
    return radius * np.cos(phi), radius * np.sin(phi)


def _project_directions(
    directions: np.ndarray,
    hemisphere: Literal["upper", "lower"],
    projection: Literal["equal_area", "stereographic"],
    *,
    antipodal: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project Cartesian directions into one selected hemisphere."""
    vectors = np.asarray(directions, dtype=float).copy()
    if antipodal:
        wrong = vectors[:, 2] < 0.0 if hemisphere == "upper" else vectors[:, 2] > 0.0
        vectors[wrong] *= -1.0
    eligible = (
        vectors[:, 2] >= -1.0e-12 if hemisphere == "upper" else vectors[:, 2] <= 1.0e-12
    )
    norms = np.linalg.norm(vectors, axis=1)
    eligible &= np.isfinite(norms) & (norms > 0.0)
    unit = vectors[eligible] / norms[eligible, None]
    theta = np.arccos(np.clip(unit[:, 2], -1.0, 1.0))
    phi = np.mod(np.arctan2(unit[:, 1], unit[:, 0]), 2.0 * np.pi)
    x, y = _project_angles(
        theta,
        phi,
        hemisphere=hemisphere,
        projection=projection,
    )
    return x, y, eligible


def _mapped_levels(values: np.ndarray, count: int) -> np.ndarray:
    """Return stable contour levels for variable or constant data."""
    minimum = float(np.nanmin(values))
    maximum = float(np.nanmax(values))
    if np.isclose(minimum, maximum, rtol=1.0e-12, atol=1.0e-14):
        width = max(abs(minimum), 1.0) * 1.0e-6
        return np.linspace(minimum - width, maximum + width, count)
    return np.linspace(minimum, maximum, count)


def _decorate_spherical_axis(
    axis: Any,
    hemisphere: Literal["upper", "lower"],
    label_mode: Literal["cartesian", "crystal"],
    *,
    compact: bool,
) -> None:
    """Draw the projection boundary and visible tensor-frame directions."""
    from matplotlib.patches import Circle

    axis.add_patch(Circle((0.0, 0.0), 1.0, fill=False, color="black", linewidth=1.1))
    axis.axhline(0.0, color="black", linewidth=0.45, alpha=0.55)
    axis.axvline(0.0, color="black", linewidth=0.45, alpha=0.55)
    axis.annotate(
        "",
        xy=(0.98, 0.0),
        xytext=(0.0, 0.0),
        arrowprops={"arrowstyle": "->", "linewidth": 0.9, "color": "black"},
    )
    axis.annotate(
        "",
        xy=(0.0, 0.98),
        xytext=(0.0, 0.0),
        arrowprops={"arrowstyle": "->", "linewidth": 0.9, "color": "black"},
    )
    if label_mode == "cartesian":
        right, left, top, bottom = "+x", "−x", "+y", "−y"
        center = "+z" if hemisphere == "upper" else "−z"
    else:
        right, left, top, bottom = "[100]", "[-100]", "[010]", "[0-10]"
        center = "[001]" if hemisphere == "upper" else "[00-1]"
    fontsize = 9 if compact else 11
    offset = 1.065
    axis.text(
        offset, 0.0, right, ha="left", va="center", fontsize=fontsize, weight="bold"
    )
    axis.text(-offset, 0.0, left, ha="right", va="center", fontsize=fontsize)
    axis.text(
        0.0, offset, top, ha="center", va="bottom", fontsize=fontsize, weight="bold"
    )
    axis.text(0.0, -offset, bottom, ha="center", va="top", fontsize=fontsize)
    axis.text(
        0.0,
        -0.055,
        center,
        ha="center",
        va="top",
        fontsize=fontsize,
        weight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.65, "pad": 0.5},
    )
    axis.set_xlim(-1.12, 1.12)
    axis.set_ylim(-1.12, 1.12)
    axis.set_aspect("equal")
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)


def _render_spherical_markers(
    axis: Any,
    spec: SphericalMapSpec,
    hemisphere: Literal["upper", "lower"],
    options: MatplotlibOptions,
    *,
    compact: bool,
) -> None:
    """Render extrema markers in a consistent minimum/maximum style."""
    marker_map = {"circle": "o", "square": "s", "triangle": "^", "diamond": "D"}
    for marker in spec.markers:
        x, y, eligible = _project_directions(
            marker.directions,
            hemisphere,
            spec.projection,
            antipodal=bool(marker.metadata.get("antipodal", False)),
        )
        if x.size == 0:
            continue
        is_minimum = marker.key == "minimum"
        axis.scatter(
            x,
            y,
            marker=marker_map.get(marker.marker, marker.marker),
            s=38 if compact else 58,
            facecolor="white" if is_minimum else "black",
            edgecolor="black" if is_minimum else "white",
            linewidth=1.1,
            zorder=8,
        )
        if options.annotate_extrema and not compact:
            value = marker.metadata.get("value")
            text = "min" if is_minimum else "max"
            if isinstance(value, (float, int)):
                text += f"\n{float(value):.4g}"
            axis.annotate(
                text,
                xy=(x[0], y[0]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=10,
                weight="bold",
                zorder=9,
            )


def _render_axis_fields(
    axis: Any,
    spec: SphericalMapSpec,
    hemisphere: Literal["upper", "lower"],
    options: MatplotlibOptions,
    *,
    compact: bool,
) -> None:
    """Render projected axial fields such as shear-wave polarizations."""
    from matplotlib.collections import LineCollection

    for layer in spec.axis_layers:
        directions = np.asarray(layer.directions, dtype=float)
        axes = np.asarray(layer.axes, dtype=float)
        mask = np.all(np.isfinite(directions), axis=1) & np.all(
            np.isfinite(axes), axis=1
        )
        if layer.resolved_mask is not None:
            mask &= np.asarray(layer.resolved_mask, dtype=bool)
        directions = directions[mask]
        axes = axes[mask]
        segments: list[np.ndarray] = []
        length = options.polarization_scale or (0.045 if compact else 0.065)
        for direction, polarization in zip(directions, axes, strict=True):
            prepared_direction = direction.copy()
            if hemisphere == "upper" and prepared_direction[2] < 0.0:
                prepared_direction *= -1.0
            elif hemisphere == "lower" and prepared_direction[2] > 0.0:
                prepared_direction *= -1.0
            if hemisphere == "upper" and prepared_direction[2] < -1.0e-12:
                continue
            if hemisphere == "lower" and prepared_direction[2] > 1.0e-12:
                continue
            tangent = (
                polarization
                - np.dot(polarization, prepared_direction) * prepared_direction
            )
            norm = float(np.linalg.norm(tangent))
            if not np.isfinite(norm) or norm <= 1.0e-10:
                continue
            tangent /= norm
            epsilon = 1.0e-4
            plus = prepared_direction + epsilon * tangent
            minus = prepared_direction - epsilon * tangent
            plus /= np.linalg.norm(plus)
            minus /= np.linalg.norm(minus)
            px, py, _ = _project_directions(
                np.stack((plus, minus)),
                hemisphere,
                spec.projection,
                antipodal=False,
            )
            cx, cy, _ = _project_directions(
                prepared_direction[None, :],
                hemisphere,
                spec.projection,
                antipodal=False,
            )
            if px.size != 2 or cx.size != 1:
                continue
            direction_2d = np.array([px[0] - px[1], py[0] - py[1]])
            direction_norm = float(np.linalg.norm(direction_2d))
            if direction_norm <= 1.0e-12:
                continue
            direction_2d /= direction_norm
            center = np.array([cx[0], cy[0]])
            half = 0.5 * length * direction_2d
            segments.append(np.stack((center - half, center + half)))
        if not segments:
            continue
        collection = LineCollection(
            segments,
            colors=options.polarization_color or "black",
            linewidths=options.polarization_line_width or (0.8 if compact else 1.2),
            alpha=0.9,
            zorder=7,
        )
        axis.add_collection(collection)


__all__ = ["_render_spherical_map", "_render_spherical_summary"]
