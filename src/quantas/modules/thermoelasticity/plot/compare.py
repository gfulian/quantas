# -*- coding: utf-8 -*-

"""Isothermal-versus-adiabatic thermoelastic comparison plots."""

from __future__ import annotations

import numpy as np

from quantas.models import (
    ContourPlotSpec,
    LinePlotSpec,
    PanelPlotSpec,
    PlotAxis,
    PlotCollection,
    PlotSeries,
    PlotSeriesStyle,
    ResultData,
)
from quantas.modules.thermoelasticity.models import ThermoelasticResult
from quantas.modules.thermoelasticity.plot.common import resolve_style
from quantas.modules.thermoelasticity.plot.components import (
    ThermoelasticComponentGroup,
    component_indices,
    component_symbol,
    resolve_components,
)
from quantas.modules.thermoelasticity.plot.options import (
    ThermoelasticComparePlotOptions,
)
from quantas.modules.thermoelasticity.analysis_engine import (
    ThermoelasticAnalysisEngine,
)


def build_thermoelastic_compare_plots(
    result_data: ResultData,
    *,
    components: tuple[str, ...] | list[str] | None = None,
    component_group: ThermoelasticComponentGroup = "independent",
    options: ThermoelasticComparePlotOptions | None = None,
) -> PlotCollection:
    """Build line plots comparing :math:`C^T_{IJ}` and :math:`C^S_{IJ}`.

    One pressure or one temperature is held fixed exactly, while the opposite
    archived coordinate is varied.  The model is re-evaluated at the requested
    coordinate rather than selecting the nearest stored grid plane.

    Parameters
    ----------
    result_data : ResultData
        Fit or analysis archive containing a thermoelastic payload.
    components : sequence of str or None, optional
        Explicit stiffness components.
    component_group : str, optional
        Semantic component group when explicit labels are absent.
    options : ThermoelasticComparePlotOptions or None, optional
        Fixed coordinate and presentation controls.

    Returns
    -------
    PlotCollection
        Frontend-neutral line or panel specifications.

    Raises
    ------
    ValueError
        If adiabatic stiffness is unavailable at any requested state.
    """
    source = result_data.results.get("thermoelasticity")
    if not isinstance(source, ThermoelasticResult):
        raise ValueError("result does not contain a thermoelasticity payload")
    resolved = options or ThermoelasticComparePlotOptions(fixed_pressure=0.0)
    labels = resolve_components(source, components, group=component_group)
    if resolved.fixed_pressure is not None:
        payload = ThermoelasticAnalysisEngine(result_data).constant_pressure_section(
            resolved.fixed_pressure
        )
        coordinate = payload.temperature
        fixed_text = f"P={resolved.fixed_pressure:g} GPa"
        x_axis = PlotAxis(key="temperature", label="Temperature (K)", unit="K")
        varying_temperature = True
        filename_suffix = f"P{resolved.fixed_pressure:g}GPa"
    else:
        assert resolved.fixed_temperature is not None
        payload = ThermoelasticAnalysisEngine(result_data).constant_temperature_section(
            resolved.fixed_temperature
        )
        coordinate = payload.pressure
        fixed_text = f"T={resolved.fixed_temperature:g} K"
        x_axis = PlotAxis(key="pressure", label="Pressure (GPa)", unit="GPa")
        varying_temperature = False
        filename_suffix = f"T{resolved.fixed_temperature:g}K"
    if payload.stiffness_isothermal is None or payload.stiffness_adiabatic is None:
        raise ValueError("both isothermal and adiabatic stiffness are required")
    if payload.adiabatic_valid_mask is not None and not np.all(
        payload.adiabatic_valid_mask
    ):
        raise ValueError("adiabatic stiffness is invalid along the requested path")
    style = resolve_style(resolved.style)
    panels: list[LinePlotSpec | ContourPlotSpec] = []
    for label in labels:
        i, j = component_indices(label)
        isothermal_values = (
            payload.stiffness_isothermal[:, 0, i, j]
            if varying_temperature
            else payload.stiffness_isothermal[0, :, i, j]
        )
        adiabatic_values = (
            payload.stiffness_adiabatic[:, 0, i, j]
            if varying_temperature
            else payload.stiffness_adiabatic[0, :, i, j]
        )
        panels.append(
            LinePlotSpec(
                key=f"{label}_isothermal_adiabatic_compare",
                title=(
                    f"{component_symbol(label)} comparison at {fixed_text}"
                    if style.show_title
                    else ""
                ),
                filename_stem=f"{label}_T_S_compare_{filename_suffix}",
                x_axis=x_axis,
                y_axis=PlotAxis(
                    key=label,
                    label=rf"{component_symbol(label)} (GPa)",
                    unit="GPa",
                ),
                series=[
                    PlotSeries(
                        key=f"{label}_isothermal",
                        label=r"Isothermal $C^T$",
                        x=coordinate.copy(),
                        y=np.asarray(isothermal_values, dtype=np.float64),
                        style=PlotSeriesStyle(
                            line_style="solid",
                            line_width=style.line_width,
                            marker="o",
                            marker_size=style.marker_size,
                            marker_edge_color=style.marker_edge_color,
                            marker_edge_width=style.marker_edge_width,
                        ),
                    ),
                    PlotSeries(
                        key=f"{label}_adiabatic",
                        label=r"Adiabatic $C^S$",
                        x=coordinate.copy(),
                        y=np.asarray(adiabatic_values, dtype=np.float64),
                        style=PlotSeriesStyle(
                            line_style="dashed",
                            line_width=style.line_width,
                            marker="s",
                            marker_size=style.marker_size,
                            marker_edge_color=style.marker_edge_color,
                            marker_edge_width=style.marker_edge_width,
                        ),
                    ),
                ],
                grid=style.grid,
                metadata={
                    "module": "thermoelasticity",
                    "plot_family": "compare",
                    "component": label,
                    "fixed_coordinate": fixed_text,
                },
            )
        )
    collection = PlotCollection()
    layout = resolved.layout
    if layout == "auto":
        layout = "facets" if len(panels) <= 6 else "separate"
    if layout == "facets" and len(panels) > 1:
        collection.plots.append(
            PanelPlotSpec(
                key="thermoelastic_isothermal_adiabatic_compare",
                title=(
                    f"Isothermal and adiabatic stiffness at {fixed_text}"
                    if style.show_title
                    else ""
                ),
                filename_stem=f"thermoelastic_T_S_compare_{filename_suffix}",
                panels=panels,
                columns=resolved.panel_columns,
                share_x=True,
                metadata={"module": "thermoelasticity", "plot_family": "compare"},
            )
        )
    else:
        collection.plots.extend(panels)
    return collection


def _payload(result_data: ResultData) -> ThermoelasticResult:
    value = result_data.results.get("thermoelasticity")
    if not isinstance(value, ThermoelasticResult):
        raise ValueError("result does not contain a thermoelasticity payload")
    return value


__all__ = ["build_thermoelastic_compare_plots"]
