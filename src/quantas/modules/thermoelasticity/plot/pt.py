# -*- coding: utf-8 -*-

"""Pressure-temperature contour builders for quasi-static stiffness data."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from quantas.models import (
    ContourPlotSpec,
    LinePlotSpec,
    PanelPlotSpec,
    PlotAxis,
    PlotCollection,
    ResultData,
)
from quantas.modules.thermoelasticity.models import ThermoelasticResult
from quantas.modules.thermoelasticity.plot.common import (
    component_grid,
    extrapolation_masks,
    extract_thermoelastic_result,
    resolve_style,
)
from quantas.modules.thermoelasticity.plot.components import (
    ThermoelasticComponentGroup,
    component_symbol,
    resolve_components,
)
from quantas.modules.thermoelasticity.plot.options import (
    ThermoelasticPTPlotOptions,
)


def build_thermoelastic_pt_plots(
    result: ResultData | ThermoelasticResult,
    components: Sequence[str] | None = None,
    *,
    component_group: ThermoelasticComponentGroup = "independent",
    options: ThermoelasticPTPlotOptions | None = None,
) -> PlotCollection:
    """Build pressure-temperature maps for selected stiffness components.

    Parameters
    ----------
    result : ResultData or ThermoelasticResult
        Reconstructed thermoelastic result.
    components : sequence of str or None, optional
        Explicit components. When omitted, ``component_group`` is resolved.
    component_group : ThermoelasticComponentGroup, optional
        Semantic group used when explicit components are absent.
    options : ThermoelasticPTPlotOptions or None, optional
        Quantity, contour, extrapolation, and layout controls.

    Returns
    -------
    PlotCollection
        Separate or faceted contour specifications.

    Raises
    ------
    ValueError
        If the archive does not contain a two-dimensional P-T grid.
    """
    payload = extract_thermoelastic_result(result)
    opts = options or ThermoelasticPTPlotOptions()
    selected = resolve_components(
        payload,
        None if components is None else tuple(components),
        group=component_group,
    )
    if payload.temperature.size < 2 or payload.pressure.size < 2:
        raise ValueError(
            "thermoelastic P-T contour maps require at least two temperatures "
            "and two pressures"
        )

    collection = PlotCollection()
    maps: list[ContourPlotSpec] = []
    for label in selected:
        try:
            spec, warnings = _build_component_map(payload, label, opts)
        except ValueError as exc:
            collection.warnings.append(f"{label}: {exc}")
            continue
        maps.append(spec)
        collection.warnings.extend(warnings)
    if not maps:
        details = "; ".join(collection.warnings) or "no components available"
        raise ValueError(f"No thermoelastic P-T plots could be prepared. {details}")

    layout = opts.layout
    if layout == "auto":
        layout = "facets" if 1 < len(maps) <= 4 else "separate"
    if layout == "facets" and len(maps) > 1:
        style = resolve_style(opts.style)
        panels: list[ContourPlotSpec | LinePlotSpec] = list(maps)
        collection.plots.append(
            PanelPlotSpec(
                key=f"pt_{opts.quantity}",
                title=(
                    f"{payload.jobname}: thermoelastic P-T maps"
                    if style.show_title
                    else ""
                ),
                filename_stem=f"thermoelastic_pt_{opts.quantity}",
                panels=panels,
                columns=opts.panel_columns,
                metadata={
                    "module": "thermoelasticity",
                    "plot_family": "pt",
                    "tensor_condition": opts.tensor_condition,
                    "quantity": opts.quantity,
                    "components": list(selected),
                },
            )
        )
    else:
        collection.plots.extend(maps)
    return collection


def _build_component_map(
    result: ThermoelasticResult,
    label: str,
    options: ThermoelasticPTPlotOptions,
) -> tuple[ContourPlotSpec, list[str]]:
    """Build one component map and return non-fatal warnings."""
    style = resolve_style(options.style)
    values, uncertainty = component_grid(
        result, label, tensor_condition=options.tensor_condition
    )
    warnings: list[str] = []

    if options.quantity == "value":
        mapped = values
        value_label = rf"{component_symbol(label)} (GPa)"
        value_unit = "GPa"
        suffix = "value"
        crosses_zero = bool(np.nanmin(mapped) < 0.0 < np.nanmax(mapped))
        center = 0.0 if crosses_zero else None
        colormap = options.colormap or ("coolwarm" if crosses_zero else "viridis")
    elif options.quantity == "uncertainty":
        if not np.any(np.isfinite(uncertainty)):
            raise ValueError("component uncertainty is unavailable")
        mapped = uncertainty
        value_label = rf"σ({component_symbol(label)}) (GPa)"
        value_unit = "GPa"
        suffix = "uncertainty"
        center = None
        colormap = options.colormap or "magma"
    else:
        if not np.any(np.isfinite(uncertainty)):
            raise ValueError("component uncertainty is unavailable")
        scale = np.abs(values)
        threshold = max(float(np.nanmax(scale)) * 1.0e-12, 1.0e-14)
        invalid = scale <= threshold
        mapped = np.full(values.shape, np.nan, dtype=np.float64)
        valid = ~invalid & np.isfinite(uncertainty)
        if not np.any(valid):
            raise ValueError(
                "relative uncertainty is undefined because no finite non-zero "
                "component values are available"
            )
        mapped[valid] = 100.0 * uncertainty[valid] / scale[valid]
        if np.any(invalid):
            warnings.append(
                f"{label}: relative uncertainty is undefined at "
                f"{int(np.count_nonzero(invalid))} near-zero grid point(s)."
            )
        value_label = rf"Relative σ({component_symbol(label)}) (%)"
        value_unit = "%"
        suffix = "relative_uncertainty"
        center = None
        colormap = options.colormap or "magma"

    masks = extrapolation_masks(result) if options.show_extrapolation else []
    return (
        ContourPlotSpec(
            key=f"{label}_{suffix}",
            title=(
                f"{component_symbol(label)}: {suffix.replace('_', ' ')}"
                if style.show_title
                else ""
            ),
            filename_stem=f"{label}_PT_{suffix}",
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
                key=f"{label}_{suffix}",
                label=value_label,
                unit=value_unit,
            ),
            x=result.pressure.copy(),
            y=result.temperature.copy(),
            z=np.asarray(mapped, dtype=np.float64).copy(),
            colormap=colormap,
            mode="smooth",
            levels=options.levels,
            isolines=options.isolines,
            isoline_labels=options.isoline_labels,
            center=center,
            masks=masks,
            metadata={
                "module": "thermoelasticity",
                "plot_family": "pt",
                "tensor_condition": options.tensor_condition,
                "component": label,
                "quantity": options.quantity,
            },
        ),
        warnings,
    )


__all__ = ["build_thermoelastic_pt_plots"]
