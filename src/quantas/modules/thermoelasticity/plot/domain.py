# -*- coding: utf-8 -*-

"""Pressure-temperature coverage and geological-path diagnostic plots."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from quantas.models import (
    LineStyle,
    ColoredPathSeries,
    ColoredPathStyle,
    ContourPlotSpec,
    PlotAxis,
    PlotCollection,
    PlotSeries,
    PlotSeriesStyle,
    ResultData,
)
from quantas.modules.thermoelasticity.models import ThermoelasticResult
from quantas.modules.thermoelasticity.plot.common import (
    extrapolation_masks,
    extract_thermoelastic_result,
    resolve_style,
)
from quantas.modules.thermoelasticity.plot.options import (
    ThermoelasticDomainPlotOptions,
)

_PROFILE_MARKERS = ("o", "s", "^", "v", "D", "P", "X", "h")
_PROFILE_LINE_STYLES: tuple[LineStyle, ...] = (
    "solid",
    "dashed",
    "dotted",
    "dashdot",
)


def build_thermoelastic_domain_plot(
    result: ResultData | ThermoelasticResult,
    *,
    profile_names: Sequence[str] | None = None,
    options: ThermoelasticDomainPlotOptions | None = None,
) -> PlotCollection:
    """Build a P-T domain map with volume coverage and optional profiles.

    Parameters
    ----------
    result : ResultData or ThermoelasticResult
        Thermoelastic result containing the QHA grid and optional profiles.
    profile_names : sequence of str or None, optional
        Explicit archived profiles to overlay. ``None`` selects all profiles
        when ``options.show_profiles`` is enabled.
    options : ThermoelasticDomainPlotOptions or None, optional
        Contour, profile, and extrapolation controls.

    Returns
    -------
    PlotCollection
        Collection containing one contour specification.

    Raises
    ------
    ValueError
        If the archive lacks a two-dimensional QHA grid or a requested profile.
    """
    payload = extract_thermoelastic_result(result)
    opts = options or ThermoelasticDomainPlotOptions()
    style = resolve_style(opts.style)
    if payload.temperature.size < 2 or payload.pressure.size < 2:
        raise ValueError(
            "thermoelastic domain plots require at least two temperatures and "
            "two pressures"
        )

    selected_profiles = _resolve_profiles(payload, profile_names, opts.show_profiles)
    ordinary: list[PlotSeries] = []
    colored: list[ColoredPathSeries] = []
    depth_limits = None
    if selected_profiles:
        depth_limits = (
            float(min(np.min(profile.depth) for profile in selected_profiles)),
            float(max(np.max(profile.depth) for profile in selected_profiles)),
        )
    for index, profile in enumerate(selected_profiles):
        if opts.color_profiles_by_depth:
            colored.append(
                ColoredPathSeries(
                    key=f"profile_{profile.name}",
                    label=profile.name,
                    x=profile.pressure,
                    y=profile.temperature,
                    values=profile.depth,
                    value_axis=PlotAxis(
                        key="depth",
                        label="Depth (km)",
                        unit="km",
                    ),
                    style=ColoredPathStyle(
                        colormap=opts.profile_colormap,
                        line_style=_PROFILE_LINE_STYLES[
                            index % len(_PROFILE_LINE_STYLES)
                        ],
                        line_width=style.line_width + 0.4,
                        marker=_PROFILE_MARKERS[index % len(_PROFILE_MARKERS)],
                        marker_size=style.marker_size,
                        marker_edge_color=style.marker_edge_color,
                        marker_edge_width=style.marker_edge_width,
                        show_colorbar=index == 0,
                        value_limits=depth_limits,
                    ),
                    metadata={"profile": profile.name},
                )
            )
        else:
            ordinary.append(
                PlotSeries(
                    key=f"profile_{profile.name}",
                    label=profile.name,
                    x=profile.pressure,
                    y=profile.temperature,
                    style=PlotSeriesStyle(
                        color="black" if style.preset == "monochrome" else None,
                        line_style=_PROFILE_LINE_STYLES[
                            index % len(_PROFILE_LINE_STYLES)
                        ],
                        line_width=style.line_width + 0.4,
                        marker=_PROFILE_MARKERS[index % len(_PROFILE_MARKERS)],
                        marker_size=style.marker_size,
                        marker_edge_color=style.marker_edge_color,
                        marker_edge_width=style.marker_edge_width,
                    ),
                    metadata={"profile": profile.name},
                )
            )

    masks = extrapolation_masks(payload) if opts.show_extrapolation else []
    spec = ContourPlotSpec(
        key="thermoelastic_domain",
        title="Thermoelastic P-T domain" if style.show_title else "",
        filename_stem="thermoelastic_PT_domain",
        x_axis=PlotAxis(
            key="pressure",
            label="Pressure (GPa)",
            unit="GPa",
        ),
        y_axis=PlotAxis(
            key="temperature",
            label="Temperature (K)",
            unit="K",
        ),
        value_axis=PlotAxis(
            key="equilibrium_volume",
            label=r"Equilibrium volume (Å$^3$)",
            unit="angstrom^3",
        ),
        x=payload.pressure.copy(),
        y=payload.temperature.copy(),
        z=np.asarray(payload.equilibrium_volume, dtype=np.float64).copy(),
        colormap=opts.colormap,
        mode="smooth",
        levels=opts.levels,
        isolines=True,
        isoline_labels=True,
        masks=masks,
        series=ordinary,
        colored_paths=colored,
        metadata={
            "module": "thermoelasticity",
            "plot_family": "domain",
            "profiles": [profile.name for profile in selected_profiles],
            "elastic_volume_bounds_A3": (
                float(payload.metadata.get("elastic_volume_min_A3", np.nan)),
                float(payload.metadata.get("elastic_volume_max_A3", np.nan)),
            ),
        },
    )
    return PlotCollection(plots=[spec])


def _resolve_profiles(
    result: ThermoelasticResult,
    names: Sequence[str] | None,
    show_profiles: bool,
):
    """Resolve archived profile objects in deterministic order."""
    if not show_profiles:
        return []
    if names is None:
        return [result.profiles[name] for name in sorted(result.profiles)]
    selected = []
    for name in names:
        try:
            selected.append(result.profiles[str(name)])
        except KeyError as exc:
            available = ", ".join(sorted(result.profiles)) or "none"
            raise ValueError(
                f"unknown thermoelastic profile '{name}'; available: {available}"
            ) from exc
    return selected


__all__ = ["build_thermoelastic_domain_plot"]
