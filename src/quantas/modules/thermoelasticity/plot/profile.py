# -*- coding: utf-8 -*-

"""Depth-profile plot builders for quasi-static thermoelastic stiffness."""

from __future__ import annotations

from typing import Any, Literal, Sequence

import numpy as np

from quantas.core.physics.elasticity import (
    cold_finite_strain_component_jacobian,
    stiffness_component_linear_coefficients,
)
from quantas.models import (
    ColoredPathSeries,
    ColoredPathStyle,
    ContourPlotSpec,
    LinePlotSpec,
    PanelPlotSpec,
    PlotAxis,
    PlotBand,
    PlotBandStyle,
    PlotCollection,
    PlotSeries,
    PlotSeriesStyle,
    ResultData,
    ScalarBackground,
    SecondaryAxis,
)
from quantas.modules.thermoelasticity.models import (
    ThermoelasticProfileResult,
    ThermoelasticResult,
)
from quantas.modules.thermoelasticity.plot.common import (
    extract_thermoelastic_result,
    profile_component,
    profile_extrapolation_spans,
    resolve_style,
)
from quantas.modules.thermoelasticity.plot.components import (
    ThermoelasticComponentGroup,
    component_indices,
    component_style,
    component_symbol,
    plot_series_style,
    resolve_components,
)
from quantas.modules.thermoelasticity.plot.options import (
    ThermoelasticProfilePlotOptions,
    ThermoelasticUncertaintyMode,
)


def build_thermoelastic_profile_plots(
    result: ResultData | ThermoelasticResult,
    *,
    profile_name: str | None = None,
    components: Sequence[str] | None = None,
    component_group: ThermoelasticComponentGroup = "independent",
    options: ThermoelasticProfilePlotOptions | None = None,
) -> PlotCollection:
    """Build absolute or relative stiffness plots along a depth profile.

    Parameters
    ----------
    result : ResultData or ThermoelasticResult
        Thermoelastic result containing one or more archived depth profiles.
    profile_name : str or None, optional
        Name of the profile to plot. When omitted, the sole archived profile is
        selected; an archive with several profiles requires an explicit name.
    components : sequence of str or None, optional
        Explicit stiffness components.
    component_group : ThermoelasticComponentGroup, optional
        Semantic group used when explicit components are absent.
    options : ThermoelasticProfilePlotOptions or None, optional
        Profile mode, layout, uncertainty, color, and annotation controls.

    Returns
    -------
    PlotCollection
        Overlay, faceted, or separate neutral profile specifications.

    Raises
    ------
    ValueError
        If the requested profile or component data are unavailable.
    """
    payload = extract_thermoelastic_result(result)
    opts = options or ThermoelasticProfilePlotOptions()
    profile = _resolve_profile(payload, profile_name)
    selected = resolve_components(
        payload,
        None if components is None else tuple(components),
        group=component_group,
    )
    collection = PlotCollection()
    uncertainty_mode = opts.uncertainty
    if uncertainty_mode == "auto":
        uncertainty_mode = "band" if len(selected) == 1 else "none"

    color_by = opts.color_by
    if opts.background == "temperature" and color_by == "temperature":
        color_by = "component"
        collection.warnings.append(
            "Temperature background and temperature-colored paths encode the "
            "same variable; paths were colored by component instead."
        )

    layout = opts.layout
    if layout == "auto":
        layout = "overlay" if len(selected) <= opts.max_overlay_components else "facets"

    if layout == "overlay":
        spec, warnings = _build_profile_spec(
            payload,
            profile,
            selected,
            opts,
            color_by=color_by,
            uncertainty_mode=uncertainty_mode,
        )
        collection.plots.append(spec)
        collection.warnings.extend(warnings)
        return collection

    panels: list[LinePlotSpec] = []
    for label in selected:
        spec, warnings = _build_profile_spec(
            payload,
            profile,
            (label,),
            opts,
            color_by=color_by,
            uncertainty_mode=uncertainty_mode,
        )
        panels.append(spec)
        collection.warnings.extend(warnings)

    if layout == "facets" and len(panels) > 1:
        style = resolve_style(opts.style)
        panel_specs: list[LinePlotSpec | ContourPlotSpec] = list(panels)
        collection.plots.append(
            PanelPlotSpec(
                key=f"profile_{profile.name}_{opts.mode}",
                title=(
                    f"{payload.jobname}: {profile.name}" if style.show_title else ""
                ),
                filename_stem=f"profile_{profile.name}_{opts.mode}",
                panels=panel_specs,
                columns=opts.panel_columns,
                share_y=True,
                metadata={
                    "module": "thermoelasticity",
                    "plot_family": "profile",
                    "profile": profile.name,
                    "mode": opts.mode,
                    "tensor_condition": opts.tensor_condition,
                    "components": list(selected),
                },
            )
        )
    else:
        collection.plots.extend(panels)
    return collection


def _build_profile_spec(
    result: ThermoelasticResult,
    profile: ThermoelasticProfileResult,
    labels: tuple[str, ...],
    options: ThermoelasticProfilePlotOptions,
    *,
    color_by: str,
    uncertainty_mode: ThermoelasticUncertaintyMode,
) -> tuple[LinePlotSpec, list[str]]:
    """Build one profile axis for one or several components."""
    style = resolve_style(options.style)
    reference_index = _reference_index(profile, options.reference_depth)
    series: list[PlotSeries] = []
    colored_paths: list[ColoredPathSeries] = []
    bands: list[PlotBand] = []
    x_axis: PlotAxis | None = None

    for label in labels:
        values, sigma, component_axis = _profile_component_values(
            result,
            profile,
            label,
            options,
            reference_index,
        )
        x_axis = component_axis
        _append_profile_component_layers(
            profile=profile,
            label=label,
            labels=labels,
            values=values,
            sigma=sigma,
            options=options,
            style=style,
            color_by=color_by,
            uncertainty_mode=uncertainty_mode,
            series=series,
            colored_paths=colored_paths,
            bands=bands,
        )

    if x_axis is None:
        raise ValueError("at least one thermoelastic component is required")
    secondary_axes, warnings = _profile_secondary_axes(profile, options)
    title = _profile_title(profile, options, style.show_title, reference_index)
    component_key = "_".join(labels)
    return (
        LinePlotSpec(
            key=f"profile_{profile.name}_{component_key}_{options.mode}",
            title=title,
            filename_stem=f"profile_{profile.name}_{component_key}_{options.mode}",
            x_axis=x_axis,
            y_axis=PlotAxis(key="depth", label="Depth (km)", unit="km"),
            series=series,
            bands=bands,
            colored_paths=colored_paths,
            secondary_axes=secondary_axes,
            spans=(
                profile_extrapolation_spans(profile)
                if options.show_extrapolation
                else []
            ),
            backgrounds=_profile_backgrounds(profile, options, color_by),
            legend_columns=1 if len(labels) < 4 else 2,
            grid=style.grid,
            invert_y_axis=True,
            metadata={
                "module": "thermoelasticity",
                "plot_family": "profile",
                "profile": profile.name,
                "mode": options.mode,
                "tensor_condition": options.tensor_condition,
                "components": list(labels),
                "reference_index": reference_index,
                "reference_depth_km": float(profile.depth[reference_index]),
            },
        ),
        warnings,
    )


def _profile_component_values(
    result: ThermoelasticResult,
    profile: ThermoelasticProfileResult,
    label: str,
    options: ThermoelasticProfilePlotOptions,
    reference_index: int,
) -> tuple[np.ndarray, np.ndarray, PlotAxis]:
    """Return one component series, uncertainty, and matching x axis."""
    absolute, absolute_sigma = profile_component(
        profile,
        label,
        tensor_condition=options.tensor_condition,
    )
    if options.mode == "relative":
        values, sigma = _relative_component_profile(
            result,
            profile,
            label,
            reference_index,
            tensor_condition=options.tensor_condition,
        )
        axis = PlotAxis(
            key="relative_stiffness",
            label=r"Relative change in $C_{IJ}$ (%)",
            unit="%",
        )
        return values, sigma, axis
    return (
        absolute,
        absolute_sigma,
        PlotAxis(
            key="stiffness",
            label=r"Elastic stiffness $C_{IJ}$ (GPa)",
            unit="GPa",
        ),
    )


def _append_profile_component_layers(
    *,
    profile: ThermoelasticProfileResult,
    label: str,
    labels: tuple[str, ...],
    values: np.ndarray,
    sigma: np.ndarray,
    options: ThermoelasticProfilePlotOptions,
    style: Any,
    color_by: str,
    uncertainty_mode: ThermoelasticUncertaintyMode,
    series: list[PlotSeries],
    colored_paths: list[ColoredPathSeries],
    bands: list[PlotBand],
) -> None:
    """Append line, scalar-color, and uncertainty layers for one component."""
    assigned = component_style(label)
    ordinary_style = plot_series_style(
        label,
        preset=style.preset,
        line_width=style.line_width,
        marker_size=style.marker_size,
        marker_edge_color=style.marker_edge_color,
        marker_edge_width=style.marker_edge_width,
    )
    if color_by == "temperature":
        colored_paths.append(
            ColoredPathSeries(
                key=f"{label}_{options.mode}",
                label=component_symbol(label),
                x=values,
                y=profile.depth,
                values=profile.temperature,
                value_axis=PlotAxis(
                    key="temperature",
                    label="Temperature (K)",
                    unit="K",
                ),
                style=ColoredPathStyle(
                    colormap=options.temperature_colormap,
                    line_style=assigned.line_style,
                    line_width=style.line_width,
                    marker=assigned.marker,
                    marker_size=style.marker_size,
                    marker_edge_color=style.marker_edge_color,
                    marker_edge_width=style.marker_edge_width,
                    show_colorbar=label == labels[0],
                    value_limits=(
                        float(np.nanmin(profile.temperature)),
                        float(np.nanmax(profile.temperature)),
                    ),
                ),
                metadata={"component": label, "profile": profile.name},
            )
        )
    else:
        if color_by == "none":
            ordinary_style.color = "black"
        series.append(
            PlotSeries(
                key=f"{label}_{options.mode}",
                label=component_symbol(label),
                x=values,
                y=profile.depth,
                style=ordinary_style,
                metadata={"component": label, "profile": profile.name},
            )
        )
    _append_profile_uncertainty(
        profile=profile,
        label=label,
        labels=labels,
        values=values,
        sigma=sigma,
        options=options,
        style=style,
        ordinary_style=ordinary_style,
        color_by=color_by,
        uncertainty_mode=uncertainty_mode,
        series=series,
        bands=bands,
    )


def _append_profile_uncertainty(
    *,
    profile: ThermoelasticProfileResult,
    label: str,
    labels: tuple[str, ...],
    values: np.ndarray,
    sigma: np.ndarray,
    options: ThermoelasticProfilePlotOptions,
    style: Any,
    ordinary_style: PlotSeriesStyle,
    color_by: str,
    uncertainty_mode: ThermoelasticUncertaintyMode,
    series: list[PlotSeries],
    bands: list[PlotBand],
) -> None:
    """Append the requested one-sigma presentation for one profile component."""
    if not np.any(np.isfinite(sigma)):
        return
    if uncertainty_mode == "band":
        bands.append(
            PlotBand(
                key=f"{label}_{options.mode}_uncertainty",
                label=(
                    rf"1σ {component_symbol(label)}"
                    if len(labels) == 1
                    else "_nolegend_"
                ),
                coordinates=profile.depth,
                lower=values - sigma,
                upper=values + sigma,
                orientation="horizontal",
                style=PlotBandStyle(
                    color=(
                        "0.75"
                        if style.preset == "monochrome"
                        else (
                            "0.65"
                            if color_by == "temperature"
                            else ordinary_style.color
                        )
                    ),
                    alpha=0.18,
                ),
                metadata={"component": label, "uncertainty": "one_sigma"},
            )
        )
    elif uncertainty_mode == "bars":
        series.append(
            PlotSeries(
                key=f"{label}_{options.mode}_errorbars",
                label="_nolegend_",
                x=values,
                y=profile.depth,
                x_error=sigma,
                style=PlotSeriesStyle(
                    color=(
                        "black"
                        if style.preset == "monochrome" or color_by == "temperature"
                        else ordinary_style.color
                    ),
                    line_style="none",
                    marker=None,
                    marker_size=0.1,
                    errorbar_line_width=style.errorbar_width,
                    errorbar_capsize=style.errorbar_capsize,
                ),
                metadata={"component": label, "uncertainty": "one_sigma"},
            )
        )


def _profile_backgrounds(
    profile: ThermoelasticProfileResult,
    options: ThermoelasticProfilePlotOptions,
    color_by: str,
) -> list[ScalarBackground]:
    """Return optional scalar backgrounds for one profile axis."""
    if options.background != "temperature":
        return []
    return [
        ScalarBackground(
            key="temperature_background",
            coordinates=profile.depth,
            values=profile.temperature,
            value_axis=PlotAxis(
                key="temperature",
                label="Temperature (K)",
                unit="K",
            ),
            axis="y",
            colormap=options.temperature_colormap,
            alpha=0.18,
            show_colorbar=color_by != "temperature",
        )
    ]


def _profile_secondary_axes(
    profile: ThermoelasticProfileResult,
    options: ThermoelasticProfilePlotOptions,
) -> tuple[list[SecondaryAxis], list[str]]:
    """Return prepared secondary axes and any omission warnings."""
    if not options.show_pressure_axis:
        return [], []
    pressure_axis = _pressure_secondary_axis(profile)
    if pressure_axis is None:
        return [], [
            f"Profile '{profile.name}' has non-monotonic pressure; "
            "the secondary pressure axis was omitted."
        ]
    return [pressure_axis], []


def _profile_title(
    profile: ThermoelasticProfileResult,
    options: ThermoelasticProfilePlotOptions,
    show_title: bool,
    reference_index: int,
) -> str:
    """Return the optional profile title."""
    if not show_title:
        return ""
    title = f"{profile.name}: {options.mode} thermoelastic profile"
    if options.mode == "relative":
        title += f" (reference {profile.depth[reference_index]:g} km)"
    return title


def _resolve_profile(
    result: ThermoelasticResult,
    name: str | None,
) -> ThermoelasticProfileResult:
    """Resolve one archived depth profile by explicit or unambiguous name."""
    if name is not None:
        try:
            return result.profiles[name]
        except KeyError as exc:
            available = ", ".join(sorted(result.profiles)) or "none"
            raise ValueError(
                f"unknown thermoelastic profile '{name}'; available: {available}"
            ) from exc
    if len(result.profiles) == 1:
        return next(iter(result.profiles.values()))
    if not result.profiles:
        raise ValueError("thermoelastic result does not contain depth profiles")
    raise ValueError(
        "profile_name is required when the archive contains several profiles"
    )


def _reference_index(
    profile: ThermoelasticProfileResult,
    reference_depth: float | None,
) -> int:
    """Return the closest profile index to the requested reference depth."""
    if reference_depth is None:
        return 0
    return int(np.argmin(np.abs(profile.depth - float(reference_depth))))


def _pressure_secondary_axis(
    profile: ThermoelasticProfileResult,
) -> SecondaryAxis | None:
    """Build a pressure axis when pressure varies monotonically with depth."""
    difference = np.diff(profile.pressure)
    if not (np.all(difference >= 0.0) or np.all(difference <= 0.0)):
        return None
    count = min(6, profile.depth.size)
    positions = np.linspace(
        float(profile.depth[0]),
        float(profile.depth[-1]),
        count,
        dtype=np.float64,
    )
    pressure = np.interp(positions, profile.depth, profile.pressure)
    return SecondaryAxis(
        key="pressure",
        label="Pressure (GPa)",
        orientation="y",
        location="right",
        positions=positions,
        labels=tuple(f"{value:.2f}" for value in pressure),
        metadata={"source": "profile_pressure"},
    )


def _relative_component_profile(
    result: ThermoelasticResult,
    profile: ThermoelasticProfileResult,
    label: str,
    reference_index: int,
    *,
    tensor_condition: str = "isothermal",
) -> tuple[np.ndarray, np.ndarray]:
    """Return relative change and correlated one-sigma uncertainty.

    Shared component-fit and reference-EOS covariance are retained between each
    profile state and the reference state. QHA volume uncertainties are already
    included in pointwise variances but are treated as uncorrelated between
    distinct profile states because cross-state QHA covariance is not archived.
    """
    values, _ = profile_component(profile, label, tensor_condition=tensor_condition)
    reference = float(values[reference_index])
    if not np.isfinite(reference) or abs(reference) <= 1.0e-14:
        raise ValueError(
            f"relative profile for {label} is undefined at a zero reference value"
        )
    relative = 100.0 * (values / reference - 1.0)

    coefficients = _component_coefficients(result, label)
    covariance = np.asarray(profile.independent_stiffness_covariance, dtype=np.float64)
    variance = np.einsum("i,nij,j->n", coefficients, covariance, coefficients)
    cross_independent = _cross_state_independent_covariance(
        result,
        profile,
        reference_index,
    )
    cross = np.einsum(
        "i,nij,j->n",
        coefficients,
        cross_independent,
        coefficients,
    )
    reference_variance = float(variance[reference_index])
    gradient_current = 100.0 / reference
    gradient_reference = -100.0 * values / (reference * reference)
    relative_variance = (
        np.square(gradient_current) * variance
        + np.square(gradient_reference) * reference_variance
        + 2.0 * gradient_current * gradient_reference * cross
    )
    relative_variance[reference_index] = 0.0
    sigma = np.sqrt(np.clip(relative_variance, 0.0, None))
    return np.asarray(relative, dtype=np.float64), np.asarray(sigma, dtype=np.float64)


def _component_coefficients(
    result: ThermoelasticResult,
    label: str,
) -> np.ndarray:
    """Return the linear independent-component coefficients for one matrix entry."""
    symmetry = str(result.metadata.get("elastic_symmetry", ""))
    coefficients = stiffness_component_linear_coefficients(
        symmetry,
        result.independent_labels,
    )
    first, second = component_indices(label)
    return np.asarray(coefficients[first, second, :], dtype=np.float64)


def _cross_state_independent_covariance(
    result: ThermoelasticResult,
    profile: ThermoelasticProfileResult,
    reference_index: int,
) -> np.ndarray:
    """Return covariance between every independent vector and one reference state."""
    labels = result.independent_labels
    count = len(labels)
    nstates = profile.depth.size
    beta_jacobian = np.zeros((nstates, count, 2), dtype=np.float64)
    eos_jacobian = np.zeros((nstates, count, 3), dtype=np.float64)
    order_value = int(result.metadata.get("finite_strain_order", 3))
    order: Literal[2, 3] = 2 if order_value == 2 else 3

    for component_index, independent_label in enumerate(labels):
        record = result.component_fits[independent_label]
        parameters = record.parameters
        if parameters is None or record.zero_by_tolerance:
            continue
        jacobian = cold_finite_strain_component_jacobian(
            profile.volume,
            reference_volume=result.reference_eos.reference_volume,
            bulk_modulus=result.reference_eos.bulk_modulus,
            bulk_modulus_derivative=result.reference_eos.bulk_modulus_derivative,
            reference_component=float(parameters[0]),
            component_pressure_derivative=float(parameters[1]),
            wallace_delta=record.wallace_delta,
            order=order,
        )
        local = np.asarray(jacobian[:, :2], dtype=np.float64)
        total_eos = np.asarray(jacobian[:, [2, 3, 4]], dtype=np.float64)
        sensitivity = record.metadata.get("reference_eos_parameter_sensitivity")
        if sensitivity is not None:
            sensitivity_array = np.asarray(sensitivity, dtype=np.float64)
            if sensitivity_array.shape == (2, 3):
                total_eos = total_eos + local @ sensitivity_array
        beta_jacobian[:, component_index, :] = local
        eos_jacobian[:, component_index, :] = total_eos

    cross = np.zeros((nstates, count, count), dtype=np.float64)
    for component_index, independent_label in enumerate(labels):
        fit_covariance = result.component_fits[independent_label].covariance
        if fit_covariance is None:
            continue
        cross[:, component_index, component_index] += np.einsum(
            "ni,ij,j->n",
            beta_jacobian[:, component_index, :],
            fit_covariance,
            beta_jacobian[reference_index, component_index, :],
        )
    if result.reference_eos.covariance is not None:
        cross += np.einsum(
            "nai,ij,bj->nab",
            eos_jacobian,
            result.reference_eos.covariance,
            eos_jacobian[reference_index],
        )
    cross[reference_index] = profile.independent_stiffness_covariance[reference_index]
    return np.asarray(cross, dtype=np.float64)


__all__ = ["build_thermoelastic_profile_plots"]
