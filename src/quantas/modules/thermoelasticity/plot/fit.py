# -*- coding: utf-8 -*-

"""Frontend-neutral fit and residual plots for thermoelastic components."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from quantas.models import (
    LinePlotSpec,
    PanelPlotSpec,
    PlotAxis,
    PlotBand,
    PlotBandStyle,
    PlotCollection,
    PlotSeries,
    PlotSeriesStyle,
    ResultData,
)
from quantas.modules.thermoelasticity.fitting import (
    evaluate_component_predictions,
)
from quantas.modules.thermoelasticity.models import (
    ThermoelasticOptions,
    ThermoelasticResult,
)
from quantas.modules.thermoelasticity.plot.common import (
    confidence_multiplier,
    extract_thermoelastic_result,
    resolve_style,
)
from quantas.modules.thermoelasticity.plot.components import (
    ThermoelasticComponentGroup,
    component_symbol,
    plot_series_style,
    resolve_components,
)
from quantas.modules.thermoelasticity.plot.options import (
    ThermoelasticFitPlotOptions,
)


def build_thermoelastic_fit_plots(
    result: ResultData | ThermoelasticResult,
    components: Sequence[str] | None = None,
    *,
    component_group: ThermoelasticComponentGroup = "independent",
    options: ThermoelasticFitPlotOptions | None = None,
) -> PlotCollection:
    """Build observed, fitted, confidence-band, and residual specifications.

    Parameters
    ----------
    result : ResultData or ThermoelasticResult
        Thermoelastic fit result.
    components : sequence of str or None, optional
        Explicit independent components. When omitted, ``component_group`` is
        resolved against the available result.
    component_group : ThermoelasticComponentGroup, optional
        Semantic component group used when ``components`` is omitted. Fit plots
        are emitted only for components carrying independent fit records.
    options : ThermoelasticFitPlotOptions or None, optional
        Fit-plot controls.

    Returns
    -------
    PlotCollection
        One line or two-panel fit specification per selected component.

    Raises
    ------
    ValueError
        If no selected component has a fit record.
    """
    payload = extract_thermoelastic_result(result)
    opts = options or ThermoelasticFitPlotOptions()
    style = resolve_style(opts.style)
    selected = resolve_components(
        payload,
        None if components is None else tuple(components),
        group=component_group,
    )
    selected = tuple(label for label in selected if label in payload.component_fits)
    if not selected:
        raise ValueError("no selected component has an independent fit record")

    collection = PlotCollection()
    for label in selected:
        component = payload.component_fits[label]
        if component.fit is None or component.parameters is None:
            collection.warnings.append(
                f"{label}: component has no successful numerical fit; skipped."
            )
            continue
        fit_spec = _build_component_fit_spec(payload, label, opts)
        if opts.residuals:
            fit_spec.title = ""
            fit_spec.x_axis.label = ""
            residual_spec = _build_component_residual_spec(payload, label, opts)
            collection.plots.append(
                PanelPlotSpec(
                    key=f"{label}_fit_diagnostics",
                    title=(
                        f"{component_symbol(label)} cold finite-strain fit"
                        if style.show_title
                        else ""
                    ),
                    filename_stem=f"{label}_cold_finite_strain_fit",
                    panels=[fit_spec, residual_spec],
                    columns=1,
                    share_x=True,
                    metadata={
                        "module": "thermoelasticity",
                        "plot_family": "fit",
                        "component": label,
                    },
                )
            )
        else:
            collection.plots.append(fit_spec)
    if not collection.plots:
        details = "; ".join(collection.warnings) or "no successful component fits"
        raise ValueError(f"No thermoelastic fit plots could be prepared. {details}")
    return collection


def _build_component_fit_spec(
    result: ThermoelasticResult,
    label: str,
    options: ThermoelasticFitPlotOptions,
) -> LinePlotSpec:
    """Build the upper observed-versus-fitted panel for one component."""
    component = result.component_fits[label]
    style = resolve_style(options.style)
    curve_volume = np.linspace(
        float(np.min(component.volumes)),
        float(np.max(component.volumes)),
        options.curve_points,
        dtype=np.float64,
    )
    prediction_options = _prediction_options(result)
    predicted, sigma, _ = evaluate_component_predictions(
        result.component_fits,
        (label,),
        curve_volume,
        result.reference_eos,
        prediction_options,
    )
    curve = np.asarray(predicted[..., 0], dtype=np.float64)
    curve_sigma = np.asarray(sigma[..., 0], dtype=np.float64)

    observed_style = PlotSeriesStyle(
        color="black",
        line_style="none",
        marker="o",
        marker_size=style.marker_size,
        marker_edge_color=style.marker_edge_color,
        marker_edge_width=style.marker_edge_width,
        errorbar_line_width=style.errorbar_width,
        errorbar_capsize=style.errorbar_capsize,
    )
    fit_style = plot_series_style(
        label,
        preset=style.preset,
        line_width=style.line_width,
        marker_size=style.marker_size,
        marker_edge_color=style.marker_edge_color,
        marker_edge_width=style.marker_edge_width,
        line_only=True,
    )
    fit_style.color = "black" if style.preset == "monochrome" else None
    observed_error = (
        component.symmetry_spread.copy() if options.show_symmetry_spread else None
    )
    series = [
        PlotSeries(
            key=f"{label}_observed",
            label="Observed",
            x=component.volumes.copy(),
            y=component.observed.copy(),
            y_error=observed_error,
            style=observed_style,
            metadata={
                "component": label,
                "symmetry_spread_is_statistical_uncertainty": False,
            },
        ),
        PlotSeries(
            key=f"{label}_fitted",
            label="Finite-strain fit",
            x=curve_volume,
            y=curve,
            style=fit_style,
            metadata={"component": label, "model": "cold_finite_strain"},
        ),
    ]
    bands: list[PlotBand] = []
    if options.confidence is not None and np.any(np.isfinite(curve_sigma)):
        multiplier = confidence_multiplier(options.confidence)
        width = multiplier * curve_sigma
        bands.append(
            PlotBand(
                key=f"{label}_confidence",
                label=f"{100.0 * options.confidence:.1f}% propagated band",
                coordinates=curve_volume,
                lower=curve - width,
                upper=curve + width,
                style=PlotBandStyle(
                    color="0.75" if style.preset == "monochrome" else None,
                    alpha=0.22,
                ),
                metadata={
                    "probability": float(options.confidence),
                    "includes_component_parameter_covariance": True,
                    "includes_reference_eos_covariance": bool(
                        prediction_options.propagate_reference_eos_covariance
                    ),
                    "includes_volume_uncertainty": False,
                },
            )
        )

    return LinePlotSpec(
        key=f"{label}_fit",
        title=(f"{component_symbol(label)} versus volume" if style.show_title else ""),
        filename_stem=f"{label}_cold_finite_strain_fit",
        x_axis=PlotAxis(
            key="volume",
            label=r"Volume (Å$^3$)",
            unit="angstrom^3",
        ),
        y_axis=PlotAxis(
            key=label,
            label=rf"{component_symbol(label)} (GPa)",
            unit="GPa",
        ),
        series=series,
        bands=bands,
        grid=style.grid,
        metadata={
            "module": "thermoelasticity",
            "plot_family": "fit",
            "component": label,
        },
    )


def _build_component_residual_spec(
    result: ThermoelasticResult,
    label: str,
    options: ThermoelasticFitPlotOptions,
) -> LinePlotSpec:
    """Build the lower residual panel for one component."""
    component = result.component_fits[label]
    style = resolve_style(options.style)
    return LinePlotSpec(
        key=f"{label}_residuals",
        title="Residuals" if style.show_title else "",
        filename_stem=f"{label}_cold_finite_strain_residuals",
        x_axis=PlotAxis(
            key="volume",
            label=r"Volume (Å$^3$)",
            unit="angstrom^3",
        ),
        y_axis=PlotAxis(
            key="residual",
            label=r"Observed − fitted (GPa)",
            unit="GPa",
        ),
        series=[
            PlotSeries(
                key=f"{label}_residual",
                label="Residual",
                x=component.volumes.copy(),
                y=component.residuals.copy(),
                style=PlotSeriesStyle(
                    color="black" if style.preset == "monochrome" else None,
                    line_style="none",
                    marker=plot_series_style(
                        label,
                        preset=style.preset,
                        line_width=style.line_width,
                        marker_size=style.marker_size,
                        marker_edge_color=style.marker_edge_color,
                        marker_edge_width=style.marker_edge_width,
                    ).marker,
                    marker_size=style.marker_size,
                    marker_edge_color=style.marker_edge_color,
                    marker_edge_width=style.marker_edge_width,
                ),
            ),
            PlotSeries(
                key="zero",
                label="Zero residual",
                x=np.asarray(
                    [component.volumes.min(), component.volumes.max()],
                    dtype=np.float64,
                ),
                y=np.zeros(2, dtype=np.float64),
                style=PlotSeriesStyle(
                    color="0.35",
                    line_style="dashed",
                    line_width=1.0,
                ),
            ),
        ],
        show_legend=False,
        grid=style.grid,
        metadata={
            "module": "thermoelasticity",
            "plot_family": "fit_residuals",
            "component": label,
        },
    )


def _prediction_options(result: ThermoelasticResult) -> ThermoelasticOptions:
    """Reconstruct uncertainty-propagation options stored in result metadata."""
    uncertainty = result.metadata.get("uncertainty_propagation", {})
    if not isinstance(uncertainty, dict):
        uncertainty = {}
    order = int(result.metadata.get("finite_strain_order", 3))
    return ThermoelasticOptions(
        finite_strain_order=2 if order == 2 else 3,
        propagate_reference_eos_covariance=bool(
            uncertainty.get("shared_reference_eos_covariance", True)
        ),
        propagate_volume_uncertainty=False,
    )


__all__ = ["build_thermoelastic_fit_plots"]
