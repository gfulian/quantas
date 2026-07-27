# -*- coding: utf-8 -*-

"""Neutral plot-specification builders for QHA results.

Unit conversions requested for plotting are applied only to the prepared plot
specifications. The scientific result object and stored HDF5 data remain
unchanged, and no concrete plotting backend is imported here.
"""

from __future__ import annotations

from typing import Literal, cast

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
from scipy import constants as cs

from quantas.core.physics.units import convert_pressure, convert_temperature
from quantas.models import (
    ContourPlotSpec,
    LinePlotSpec,
    PlotAxis,
    PlotCollection,
    PlotSeries,
    PlotSeriesStyle,
    ResultData,
)
from quantas.modules.qha.models import QHAResult
from quantas.modules.qha.plot.labels import (
    QHAPlotProperty,
    available_plot_properties,
    resolve_plot_property,
    unit_for_property,
    ylabel_for_property,
)


_ALLOWED_CMAPS = ("viridis", "plasma", "inferno", "magma", "cividis", "turbo")
_ALLOWED_CONTOUR_MODES = ("discrete", "smooth")

QHACurveAxis = Literal["temperature", "pressure"]


@dataclass(slots=True)
class QHAPlotOptions:
    """Options controlling preparation of neutral QHA plot specifications.

    Parameters
    ----------
    include_contours : bool, optional
        If ``True``, prepare contour specifications in addition to line plots
        when the pressure-temperature grid is sufficiently large.
    cmap : str, optional
        Portable colormap name used by contour renderers.
    contour_mode : {"discrete", "smooth"}, optional
        Preferred filled-contour rendering mode.
    levels : int, optional
        Number of principal contour levels.
    isolines : bool, optional
        Whether contour lines should be requested.
    isoline_labels : bool, optional
        Whether contour-line labels should be requested.
    temperature_unit : str or None, optional
        Display temperature unit. If ``None``, use the result unit.
    pressure_unit : str or None, optional
        Display pressure unit. If ``None``, use the result unit.
    energy_unit : str or None, optional
        Display energy unit for energy, entropy, and heat-capacity plots.
    include_dulong_petit : bool, optional
        Whether the Dulong-Petit reference should be included when reliable
        atom-count metadata are available.
    curve_axis : {"temperature", "pressure"}, optional
        Natural variable placed on the independent axis of line plots.
    selected_pressures : tuple of float or None, optional
        Exact native pressure-grid values included in temperature sections.
        ``None`` selects all stored pressures.
    selected_temperatures : tuple of float or None, optional
        Exact native temperature-grid values included in pressure sections.
        ``None`` selects all stored temperatures.

    Notes
    -----
    Slice selection is performed on the native stored grid before unit
    conversion. Quantas does not interpolate missing pressure or temperature
    coordinates implicitly.
    """

    include_contours: bool = False
    cmap: str = "viridis"
    contour_mode: str = "smooth"
    levels: int = 12
    isolines: bool = True
    isoline_labels: bool = True
    temperature_unit: str | None = None
    pressure_unit: str | None = None
    energy_unit: str | None = None
    include_dulong_petit: bool = True
    curve_axis: QHACurveAxis = "temperature"
    selected_pressures: tuple[float, ...] | None = None
    selected_temperatures: tuple[float, ...] | None = None


@dataclass(frozen=True, slots=True)
class _PlotUnitContext:
    """Resolved native and display units used by QHA plot builders.

    Parameters
    ----------
    native_units : dict
        Units stored in the QHA result.
    display_units : dict
        Units requested for the prepared plot data.
    formula_units_per_cell : float, optional
        Number of formula units used for cell-to-mole normalization.
    """

    native_units: dict[str, object]
    display_units: dict[str, object]
    formula_units_per_cell: float = 1.0


DEFAULT_PLOT_PROPERTIES: tuple[str, ...] = (
    "VT",
    "KT",
    "KS",
    "alphaV",
    "Cp",
    "Cv",
    "Cp-Cv",
    "F",
    "G",
)


def build_qha_plot_collection(
    result: ResultData | QHAResult,
    property_names: Sequence[str] | None = None,
    options: QHAPlotOptions | None = None,
) -> PlotCollection:
    """Build neutral QHA plot specifications for selected properties.

    Parameters
    ----------
    result : ResultData or QHAResult
        QHA result container.
    property_names : sequence of str or None, optional
        Properties to prepare. If omitted, available standard properties are
        selected.
    options : QHAPlotOptions or None, optional
        Plot-preparation options.

    Returns
    -------
    PlotCollection
        Neutral specifications and non-fatal preparation warnings.

    Raises
    ------
    KeyError
        If generic result data do not contain a QHA result.
    ValueError
        If required pressure or temperature grids are missing.
    """
    qha_result = _extract_qha_result(result)
    opts = _validated_options(options or QHAPlotOptions())
    names = list(property_names or []) or [
        key for key in DEFAULT_PLOT_PROPERTIES if _has_named_property(qha_result, key)
    ]

    collection = PlotCollection()
    for name in names:
        if name.lower() in {"heat_capacities", "heat_capacity", "cp_cv"}:
            collection.plots.append(build_heat_capacity_spec(qha_result, opts))
            if opts.include_contours:
                collection.warnings.append(
                    "Contour maps are generated for single properties only; "
                    "skipped heat-capacity combined contour."
                )
            continue

        property_info = resolve_plot_property(name)
        if getattr(qha_result, property_info.attribute, None) is None:
            collection.warnings.append(
                f"Property {property_info.key} is not available; skipped."
            )
            continue
        collection.plots.append(
            build_property_curve_spec(qha_result, property_info, opts)
        )
        if opts.include_contours:
            try:
                collection.plots.append(
                    build_property_contour_spec(qha_result, property_info, opts)
                )
            except ValueError as exc:
                collection.warnings.append(str(exc))
    return collection


def build_property_curve_spec(
    result: QHAResult,
    property_info: QHAPlotProperty,
    options: QHAPlotOptions | None = None,
) -> LinePlotSpec:
    """Build one exact-grid QHA line-section specification.

    The default representation places temperature on the independent axis and
    emits one series per selected pressure.  ``curve_axis="pressure"`` places
    pressure on the independent axis and emits one series per selected stored
    temperature. Selection is exact on the native result grid and never
    interpolates missing coordinates.
    """
    opts = _validated_options(options or QHAPlotOptions())
    context = _unit_context(result, opts)
    native_temperature = _required_grid(result.temperature, "temperature")
    native_pressure = _required_grid(result.pressure, "pressure")
    temperature, pressure, values = _property_grid(result, property_info, context)
    values = values * property_info.scale

    if opts.curve_axis == "temperature":
        indices = _selected_grid_indices(
            native_pressure,
            opts.selected_pressures,
            name="pressure",
        )
        series = [
            PlotSeries(
                key=f"pressure_{index}",
                label=_pressure_label(pressure[index], context.display_units),
                x=temperature.copy(),
                y=values[:, index].copy(),
                metadata={
                    "pressure": float(pressure[index]),
                    "pressure_native": float(native_pressure[index]),
                    "pressure_index": int(index),
                    "curve_axis": "temperature",
                },
            )
            for index in indices
        ]
        x_axis = PlotAxis(
            key="temperature",
            label=_temperature_label(context.display_units),
            unit=str(context.display_units.get("temperature", "K")),
        )
        legend_title = _pressure_axis_label(context.display_units)
        filename_stem = f"{property_info.key}_1D"
    else:
        if native_pressure.size < 2:
            raise ValueError(
                "QHA pressure sections require at least two stored pressures"
            )
        indices = _selected_grid_indices(
            native_temperature,
            opts.selected_temperatures,
            name="temperature",
        )
        series = [
            PlotSeries(
                key=f"temperature_{index}",
                label=_temperature_value_label(
                    temperature[index], context.display_units
                ),
                x=pressure.copy(),
                y=values[index, :].copy(),
                metadata={
                    "temperature": float(temperature[index]),
                    "temperature_native": float(native_temperature[index]),
                    "temperature_index": int(index),
                    "curve_axis": "pressure",
                },
            )
            for index in indices
        ]
        x_axis = PlotAxis(
            key="pressure",
            label=_pressure_axis_label(context.display_units),
            unit=str(context.display_units.get("pressure", "GPa")),
        )
        legend_title = _temperature_label(context.display_units)
        filename_stem = f"{property_info.key}_P"

    y_limits = None
    if (
        property_info.key == "Cv"
        and opts.include_dulong_petit
        and opts.curve_axis == "temperature"
    ):
        reference = _dulong_petit_series(result, context, temperature)
        if reference is not None:
            series.append(reference)
            y_limits = _upper_limit(values, reference.y[0])

    return LinePlotSpec(
        key=property_info.key,
        title=property_info.description,
        filename_stem=filename_stem,
        x_axis=x_axis,
        y_axis=PlotAxis(
            key=property_info.attribute,
            label=ylabel_for_property(property_info, context.display_units),
            unit=unit_for_property(property_info, context.display_units),
            limits=y_limits,
        ),
        series=series,
        legend_title=legend_title,
        metadata={
            "module": "qha",
            "property": property_info.attribute,
            "curve_axis": opts.curve_axis,
        },
    )


def build_heat_capacity_spec(
    result: QHAResult,
    options: QHAPlotOptions | None = None,
) -> LinePlotSpec:
    """Build combined exact-grid ``C_P`` and ``C_V`` line sections."""
    opts = _validated_options(options or QHAPlotOptions())
    context = _unit_context(result, opts)
    native_temperature = _required_grid(result.temperature, "temperature")
    native_pressure = _required_grid(result.pressure, "pressure")
    temperature = _convert_temperature_grid(native_temperature, context)
    pressure = _convert_pressure_grid(native_pressure, context)
    cp = _optional_property_grid(
        result, "isobaric_heat_capacity", context, "heat_capacity"
    )
    cv = _optional_property_grid(
        result, "isochoric_heat_capacity", context, "heat_capacity"
    )
    if cp is None and cv is None:
        raise ValueError("neither Cp nor Cv is available")

    series: list[PlotSeries] = []
    if opts.curve_axis == "temperature":
        indices = _selected_grid_indices(
            native_pressure,
            opts.selected_pressures,
            name="pressure",
        )
        for index in indices:
            suffix = _pressure_label(pressure[index], context.display_units)
            common_metadata = {
                "pressure": float(pressure[index]),
                "pressure_native": float(native_pressure[index]),
                "pressure_index": int(index),
                "curve_axis": "temperature",
            }
            if cp is not None:
                series.append(
                    PlotSeries(
                        key=f"Cp_pressure_{index}",
                        label=rf"$C_P$, {suffix}",
                        x=temperature.copy(),
                        y=cp[:, index].copy(),
                        metadata={**common_metadata, "property": "Cp"},
                    )
                )
            if cv is not None:
                series.append(
                    PlotSeries(
                        key=f"Cv_pressure_{index}",
                        label=rf"$C_V$, {suffix}",
                        x=temperature.copy(),
                        y=cv[:, index].copy(),
                        style=PlotSeriesStyle(line_style="dashed"),
                        metadata={**common_metadata, "property": "Cv"},
                    )
                )
        x_axis = PlotAxis(
            key="temperature",
            label=_temperature_label(context.display_units),
            unit=str(context.display_units.get("temperature", "K")),
        )
        legend_title = _pressure_axis_label(context.display_units)
        filename_stem = "heat_capacities_1D"
    else:
        if native_pressure.size < 2:
            raise ValueError(
                "QHA pressure sections require at least two stored pressures"
            )
        indices = _selected_grid_indices(
            native_temperature,
            opts.selected_temperatures,
            name="temperature",
        )
        for index in indices:
            suffix = _temperature_value_label(
                temperature[index], context.display_units
            )
            common_metadata = {
                "temperature": float(temperature[index]),
                "temperature_native": float(native_temperature[index]),
                "temperature_index": int(index),
                "curve_axis": "pressure",
            }
            if cp is not None:
                series.append(
                    PlotSeries(
                        key=f"Cp_temperature_{index}",
                        label=rf"$C_P$, {suffix}",
                        x=pressure.copy(),
                        y=cp[index, :].copy(),
                        metadata={**common_metadata, "property": "Cp"},
                    )
                )
            if cv is not None:
                series.append(
                    PlotSeries(
                        key=f"Cv_temperature_{index}",
                        label=rf"$C_V$, {suffix}",
                        x=pressure.copy(),
                        y=cv[index, :].copy(),
                        style=PlotSeriesStyle(line_style="dashed"),
                        metadata={**common_metadata, "property": "Cv"},
                    )
                )
        x_axis = PlotAxis(
            key="pressure",
            label=_pressure_axis_label(context.display_units),
            unit=str(context.display_units.get("pressure", "GPa")),
        )
        legend_title = _temperature_label(context.display_units)
        filename_stem = "heat_capacities_P"

    y_limits = None
    if (
        cv is not None
        and opts.include_dulong_petit
        and opts.curve_axis == "temperature"
    ):
        reference = _dulong_petit_series(result, context, temperature)
        if reference is not None:
            series.append(reference)
            y_limits = _upper_limit(cv, reference.y[0])

    heat_unit = str(
        context.display_units.get(
            "heat_capacity",
            str(context.display_units.get("energy", "energy unit")) + "/K",
        )
    )
    return LinePlotSpec(
        key="heat_capacities",
        title="Heat capacities",
        filename_stem=filename_stem,
        x_axis=x_axis,
        y_axis=PlotAxis(
            key="heat_capacity",
            label=rf"$C$ ({heat_unit})",
            unit=heat_unit,
            limits=y_limits,
        ),
        series=series,
        legend_title=legend_title,
        legend_columns=2 if len(series) > 2 else 1,
        metadata={
            "module": "qha",
            "property": "heat_capacities",
            "curve_axis": opts.curve_axis,
        },
    )


def build_property_contour_spec(
    result: QHAResult,
    property_info: QHAPlotProperty,
    options: QHAPlotOptions | None = None,
) -> ContourPlotSpec:
    """Build a neutral pressure-temperature contour specification.

    Parameters
    ----------
    result : QHAResult
        QHA result object.
    property_info : QHAPlotProperty
        Resolved property metadata.
    options : QHAPlotOptions or None, optional
        Plot-preparation options.

    Returns
    -------
    ContourPlotSpec
        Neutral contour specification.

    Raises
    ------
    ValueError
        If fewer than two temperature or pressure points are available.
    """
    opts = _validated_options(options or QHAPlotOptions())
    context = _unit_context(result, opts)
    temperature, pressure, values = _property_grid(result, property_info, context)
    if temperature.size < 2 or pressure.size < 2:
        raise ValueError(
            f"Contour plot for {property_info.key} requires at least two "
            "temperatures and two pressures."
        )
    values = values * property_info.scale
    value_label = ylabel_for_property(property_info, context.display_units)
    return ContourPlotSpec(
        key=property_info.key,
        title=property_info.description,
        filename_stem=f"{property_info.key}_2D",
        x_axis=PlotAxis(
            key="temperature",
            label=_temperature_label(context.display_units),
            unit=str(context.display_units.get("temperature", "K")),
        ),
        y_axis=PlotAxis(
            key="pressure",
            label=_pressure_axis_label(context.display_units),
            unit=str(context.display_units.get("pressure", "GPa")),
        ),
        value_axis=PlotAxis(
            key=property_info.attribute,
            label=value_label,
            unit=unit_for_property(property_info, context.display_units),
        ),
        x=temperature.copy(),
        y=pressure.copy(),
        z=values.T.copy(),
        colormap=opts.cmap,
        mode=cast(Literal["discrete", "smooth"], opts.contour_mode),
        levels=opts.levels,
        isolines=opts.isolines,
        isoline_labels=opts.isoline_labels,
        metadata={"module": "qha", "property": property_info.attribute},
    )


def list_available_plot_properties(
    result: ResultData | QHAResult,
) -> list[tuple[str, str, str]]:
    """List properties available for plotting in a result.

    Parameters
    ----------
    result : ResultData or QHAResult
        QHA result container.

    Returns
    -------
    list of tuple
        Tuples containing historical key, attribute name, and description.
    """
    qha_result = _extract_qha_result(result)
    rows: list[tuple[str, str, str]] = []
    for prop in available_plot_properties().values():
        if getattr(qha_result, prop.attribute, None) is not None:
            rows.append((prop.key, prop.attribute, prop.description))
    if (
        qha_result.isobaric_heat_capacity is not None
        or qha_result.isochoric_heat_capacity is not None
    ):
        rows.append(("heat_capacities", "Cp/Cv", "Heat capacities"))
    return rows


def _extract_qha_result(result: ResultData | QHAResult) -> QHAResult:
    """Extract a QHA result from a supported result container.

    Parameters
    ----------
    result : ResultData or QHAResult
        Generic or module-specific result container.

    Returns
    -------
    QHAResult
        Extracted module result.

    Raises
    ------
    KeyError
        If generic result data do not contain a QHA result.
    """
    if isinstance(result, QHAResult):
        return result
    qha = result.results.get("qha")
    if not isinstance(qha, QHAResult):
        raise KeyError("ResultData does not contain a QHAResult under the 'qha' key")
    return qha


def _has_named_property(result: QHAResult, name: str) -> bool:
    """Return whether a named property is available in a QHA result.

    Parameters
    ----------
    result : QHAResult
        QHA result object.
    name : str
        Property key or alias.

    Returns
    -------
    bool
        ``True`` when the property resolves and contains data.
    """
    try:
        prop = resolve_plot_property(name)
    except KeyError:
        return False
    return getattr(result, prop.attribute, None) is not None


def _property_grid(
    result: QHAResult,
    property_info: QHAPlotProperty,
    context: _PlotUnitContext,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return converted grids for one QHA plot property.

    Parameters
    ----------
    result : QHAResult
        QHA result object.
    property_info : QHAPlotProperty
        Resolved plot-property metadata.
    context : _PlotUnitContext
        Native and display units.

    Returns
    -------
    tuple of ndarray
        Temperature, pressure, and converted property arrays.

    Raises
    ------
    ValueError
        If the property or a required grid is unavailable.
    """
    values = getattr(result, property_info.attribute, None)
    if values is None:
        raise ValueError(f"QHA property {property_info.key} is not available")
    temperature = _convert_temperature_grid(
        _required_grid(result.temperature, "temperature"), context
    )
    pressure = _convert_pressure_grid(
        _required_grid(result.pressure, "pressure"), context
    )
    converted = _convert_property_values(
        _as_pressure_temperature_array(
            result, np.asarray(values, dtype=np.float64), property_info.key
        ),
        property_info.unit_key,
        context,
    )
    return temperature, pressure, converted


def _optional_property_grid(
    result: QHAResult,
    attribute: str,
    context: _PlotUnitContext,
    unit_key: str,
) -> np.ndarray | None:
    """Return an optional converted pressure-temperature property grid.

    Parameters
    ----------
    result : QHAResult
        QHA result object.
    attribute : str
        Result attribute name.
    context : _PlotUnitContext
        Native and display units.
    unit_key : str
        Unit category used for conversion.

    Returns
    -------
    ndarray or None
        Converted property grid, or ``None`` when unavailable.
    """
    values = getattr(result, attribute, None)
    if values is None:
        return None
    grid = _as_pressure_temperature_array(
        result, np.asarray(values, dtype=np.float64), attribute
    )
    return _convert_property_values(grid, unit_key, context)


def _required_grid(values: np.ndarray | None, name: str) -> np.ndarray:
    """Return a validated required one-dimensional grid.

    Parameters
    ----------
    values : ndarray or None
        Candidate grid values.
    name : str
        Grid name used in error messages.

    Returns
    -------
    ndarray
        Non-empty one-dimensional grid.

    Raises
    ------
    ValueError
        If the grid is missing, empty, or not one-dimensional.
    """
    if values is None:
        raise ValueError(f"QHA {name} grid is not available")
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"QHA {name} grid must be a non-empty one-dimensional array")
    if np.unique(array).size != array.size:
        raise ValueError(f"QHA {name} grid contains duplicate coordinates")
    return array


def _selected_grid_indices(
    grid: np.ndarray,
    selected: tuple[float, ...] | None,
    *,
    name: str,
) -> tuple[int, ...]:
    """Resolve exact native-grid selections without interpolation.

    Parameters
    ----------
    grid : ndarray
        Native one-dimensional pressure or temperature grid.
    selected : tuple of float or None
        Requested exact grid values. ``None`` selects all points.
    name : str
        Grid name used in validation messages.

    Returns
    -------
    tuple of int
        Grid indices in the requested order.

    Raises
    ------
    ValueError
        If the selection is empty, duplicated, or absent from the native grid.
    """
    if selected is None:
        return tuple(range(grid.size))
    if not selected:
        raise ValueError(f"selected {name} values cannot be empty")
    if len(set(selected)) != len(selected):
        raise ValueError(f"selected {name} values must be unique")

    indices: list[int] = []
    for requested in selected:
        matches = np.flatnonzero(grid == float(requested))
        if matches.size == 0:
            available = ", ".join(f"{value:.12g}" for value in grid)
            raise ValueError(
                f"QHA {name} {requested!r} is not present in the native grid; "
                f"available values: {available}"
            )
        indices.append(int(matches[0]))
    return tuple(indices)


def _as_pressure_temperature_array(
    result: QHAResult,
    values: np.ndarray,
    property_key: str,
) -> np.ndarray:
    """Normalize property values to a ``(nT, nP)`` array.

    Parameters
    ----------
    result : QHAResult
        QHA result defining temperature and pressure dimensions.
    values : ndarray
        Candidate property values.
    property_key : str
        Property key used in error messages.

    Returns
    -------
    ndarray
        Pressure-temperature property grid.

    Raises
    ------
    ValueError
        If the array shape is incompatible with the result grids.
    """
    temperature = _required_grid(result.temperature, "temperature")
    pressure = _required_grid(result.pressure, "pressure")
    if values.ndim == 1 and values.shape[0] == temperature.shape[0]:
        return np.repeat(values[:, np.newaxis], pressure.shape[0], axis=1)
    if values.ndim == 2 and values.shape == (temperature.shape[0], pressure.shape[0]):
        return values
    raise ValueError(
        f"QHA property {property_key} has shape {values.shape}, expected "
        f"({temperature.shape[0]}, {pressure.shape[0]})"
    )


def _validated_options(options: QHAPlotOptions) -> QHAPlotOptions:
    """Validate QHA plot-preparation options.

    Parameters
    ----------
    options : QHAPlotOptions
        Candidate options.

    Returns
    -------
    QHAPlotOptions
        Validated options.

    Raises
    ------
    ValueError
        If the colormap, contour mode, or level count is unsupported.
    """
    if options.curve_axis not in {"temperature", "pressure"}:
        raise ValueError("QHA curve_axis must be 'temperature' or 'pressure'")
    if options.curve_axis == "temperature" and options.selected_temperatures:
        raise ValueError(
            "selected_temperatures is valid only when curve_axis='pressure'"
        )
    if options.curve_axis == "pressure" and options.selected_pressures:
        raise ValueError(
            "selected_pressures is valid only when curve_axis='temperature'"
        )
    if options.cmap not in _ALLOWED_CMAPS:
        raise ValueError(f"unsupported colormap '{options.cmap}'")
    if options.contour_mode not in _ALLOWED_CONTOUR_MODES:
        raise ValueError(f"unsupported contour mode '{options.contour_mode}'")
    if options.levels < 2:
        raise ValueError("contour levels must be at least 2")
    return options


def _unit_context(result: QHAResult, options: QHAPlotOptions) -> _PlotUnitContext:
    """Build the native and display unit context for QHA plots.

    Parameters
    ----------
    result : QHAResult
        QHA result with unit and normalization metadata.
    options : QHAPlotOptions
        Requested display units.

    Returns
    -------
    _PlotUnitContext
        Resolved conversion context.
    """
    native = dict(result.metadata.get("units", {}))
    display = dict(native)
    if options.temperature_unit is not None:
        display["temperature"] = options.temperature_unit
    if options.pressure_unit is not None:
        display["pressure"] = options.pressure_unit
    if options.energy_unit is not None:
        display["energy"] = options.energy_unit
        display["heat_capacity"] = _heat_capacity_unit(
            options.energy_unit, display.get("temperature", "K")
        )
        display["entropy"] = _heat_capacity_unit(
            options.energy_unit, display.get("temperature", "K")
        )
    normalization = result.metadata.get("normalization", {})
    formula_units = 1.0
    if isinstance(normalization, Mapping):
        try:
            formula_units = float(normalization.get("formula_units_per_cell", 1.0))
        except (TypeError, ValueError):
            formula_units = 1.0
    if formula_units <= 0.0:
        formula_units = 1.0
    return _PlotUnitContext(
        native_units=native,
        display_units=display,
        formula_units_per_cell=formula_units,
    )


def _convert_temperature_grid(
    values: np.ndarray, context: _PlotUnitContext
) -> np.ndarray:
    """Convert a temperature grid to display units.

    Parameters
    ----------
    values : ndarray
        Native temperature values.
    context : _PlotUnitContext
        Native and display units.

    Returns
    -------
    ndarray
        Temperature values in display units.
    """
    native = str(context.native_units.get("temperature", "K"))
    target = str(context.display_units.get("temperature", native))
    if _same_unit(native, target):
        return values
    return np.asarray(convert_temperature(values, native, target), dtype=np.float64)


def _convert_pressure_grid(values: np.ndarray, context: _PlotUnitContext) -> np.ndarray:
    """Convert a pressure grid to display units.

    Parameters
    ----------
    values : ndarray
        Native pressure values.
    context : _PlotUnitContext
        Native and display units.

    Returns
    -------
    ndarray
        Pressure values in display units.
    """
    native = str(context.native_units.get("pressure", "GPa"))
    target = str(context.display_units.get("pressure", native))
    if _same_unit(native, target):
        return values
    return np.asarray(convert_pressure(values, native, target), dtype=np.float64)


def _convert_property_values(
    values: np.ndarray,
    unit_key: str,
    context: _PlotUnitContext,
) -> np.ndarray:
    """Convert QHA property values to display units.

    Parameters
    ----------
    values : ndarray
        Native property values.
    unit_key : str
        Property unit category.
    context : _PlotUnitContext
        Native and display units.

    Returns
    -------
    ndarray
        Converted property values.
    """
    if unit_key == "pressure":
        native = str(context.native_units.get("pressure", "GPa"))
        target = str(context.display_units.get("pressure", native))
        if _same_unit(native, target):
            return values
        return np.asarray(convert_pressure(values, native, target), dtype=np.float64)
    if unit_key == "energy":
        native = str(context.native_units.get("energy", "energy unit"))
        target = str(context.display_units.get("energy", native))
        return _convert_energy_like(
            values, native, target, context.formula_units_per_cell
        )
    if unit_key in {"heat_capacity", "entropy"}:
        native = str(
            context.native_units.get(
                "heat_capacity",
                str(context.native_units.get("energy", "energy unit")) + "/K",
            )
        )
        target = str(context.display_units.get("heat_capacity", native))
        return _convert_heat_capacity_like(
            values, native, target, context.formula_units_per_cell
        )
    return values


def _same_unit(left: str, right: str) -> bool:
    """Compare two unit labels after case and whitespace normalization.

    Parameters
    ----------
    left, right : str
        Unit labels to compare.

    Returns
    -------
    bool
        ``True`` when the normalized labels are equal.
    """
    return left.strip().lower() == right.strip().lower()


def _canonical_energy_unit(unit: str) -> str:
    """Return a canonical energy unit label.

    Parameters
    ----------
    unit : str
        Candidate energy unit label.

    Returns
    -------
    str
        Canonical label understood by the local converters.
    """
    normalized = unit.strip().lower().replace(" ", "")
    normalized = normalized.replace("joule", "j").replace("kilojoule", "kj")
    normalized = normalized.replace("per", "/")
    normalized = normalized.replace("mol^-1", "/mol").replace("mol-1", "/mol")
    for token in ("cell^-1", "cell-1", "/cell"):
        normalized = normalized.replace(token, "")
    normalized = normalized.rstrip("/")
    if normalized in {"j/mol", "jmol", "jmol^-1", "jmol-1"}:
        return "J/mol"
    if normalized in {"kj/mol", "kjmol", "kjmol^-1", "kjmol-1"}:
        return "kJ/mol"
    if normalized in {"ha", "hartree"}:
        return "Ha"
    if normalized in {"ev", "electronvolt"}:
        return "eV"
    if normalized in {"ry", "rydberg"}:
        return "Ry"
    return unit


def _energy_factor_to_j_per_mol(unit: str) -> float:
    """Return the conversion factor from one energy unit to J/mol.

    Parameters
    ----------
    unit : str
        Energy unit label.

    Returns
    -------
    float
        Multiplicative conversion factor.

    Raises
    ------
    ValueError
        If the energy unit is unsupported.
    """
    canonical = _canonical_energy_unit(unit)
    if canonical == "J/mol":
        return 1.0
    if canonical == "kJ/mol":
        return 1000.0
    if canonical == "Ha":
        return cs.physical_constants["Hartree energy"][0] * cs.Avogadro
    if canonical == "eV":
        return cs.electron_volt * cs.Avogadro
    if canonical == "Ry":
        return cs.physical_constants["Rydberg constant times hc in J"][0] * cs.Avogadro
    raise ValueError(f"unsupported energy unit for plotting conversion: {unit}")


def _split_heat_capacity_unit(unit: str) -> tuple[str, str]:
    """Split a heat-capacity unit into energy and temperature units.

    Parameters
    ----------
    unit : str
        Heat-capacity or entropy unit, for example ``"J mol^-1 K^-1"`` or
        ``"kJ/mol/Celsius"``.

    Returns
    -------
    tuple of str
        Canonical energy numerator and temperature denominator. The
        temperature denominator defaults to kelvin when it is omitted.
    """
    normalized = unit.strip().lower()
    normalized = normalized.replace("−", "-").replace("⁻¹", "^-1")
    normalized = normalized.replace("degrees", "degree")
    normalized = normalized.replace(" ", "")
    normalized = normalized.replace("per", "/")
    normalized = normalized.replace("mol^-1", "/mol").replace("mol-1", "/mol")

    suffixes = (
        (
            "Celsius",
            (
                "/degreecelsius",
                "/celsius",
                "/degc",
                "/°c",
                "degreecelsius^-1",
                "degreecelsius-1",
                "celsius^-1",
                "celsius-1",
                "°c^-1",
                "°c-1",
            ),
        ),
        (
            "Fahrenheit",
            (
                "/degreefahrenheit",
                "/fahrenheit",
                "/degf",
                "/°f",
                "degreefahrenheit^-1",
                "degreefahrenheit-1",
                "fahrenheit^-1",
                "fahrenheit-1",
                "°f^-1",
                "°f-1",
            ),
        ),
        (
            "Rankine",
            (
                "/degreerankine",
                "/rankine",
                "/degr",
                "/°r",
                "degreerankine^-1",
                "degreerankine-1",
                "rankine^-1",
                "rankine-1",
                "°r^-1",
                "°r-1",
            ),
        ),
        ("K", ("/kelvin", "/k", "kelvin^-1", "kelvin-1", "k^-1", "k-1")),
    )

    temperature_unit = "K"
    energy_unit = normalized
    for candidate, endings in suffixes:
        ending = next((item for item in endings if normalized.endswith(item)), None)
        if ending is not None:
            energy_unit = normalized[: -len(ending)]
            temperature_unit = candidate
            break

    energy_unit = energy_unit.rstrip("/")
    return _canonical_energy_unit(energy_unit), temperature_unit


def _temperature_interval_factor_to_kelvin(unit: str) -> float:
    """Return the kelvin width of one display-temperature interval.

    Parameters
    ----------
    unit : str
        Temperature unit label.

    Returns
    -------
    float
        Interval width expressed in kelvin.

    Raises
    ------
    ValueError
        If the temperature unit is unsupported.
    """
    normalized = unit.strip().lower().replace("°", "")
    if normalized in {"k", "kelvin", "c", "celsius", "degc", "degreecelsius"}:
        return 1.0
    if normalized in {
        "f",
        "fahrenheit",
        "degf",
        "degreefahrenheit",
        "r",
        "rankine",
        "degr",
        "degreerankine",
    }:
        return 5.0 / 9.0
    raise ValueError(
        f"unsupported temperature unit in heat-capacity conversion: {unit}"
    )


def _convert_energy_like(
    values: np.ndarray,
    native_unit: str,
    target_unit: str,
    formula_units_per_cell: float = 1.0,
) -> np.ndarray:
    """Convert energy-like values between display units.

    Parameters
    ----------
    values : ndarray
        Native energy-like values.
    native_unit, target_unit : str
        Source and destination energy units.
    formula_units_per_cell : float, optional
        Formula-unit normalization factor.

    Returns
    -------
    ndarray
        Converted values.
    """
    if _same_unit(native_unit, target_unit):
        return values
    factor = _energy_factor_to_j_per_mol(native_unit) / _energy_factor_to_j_per_mol(
        target_unit
    )
    factor *= _molar_normalization_factor(
        native_unit, target_unit, formula_units_per_cell
    )
    return values * factor


def _convert_heat_capacity_like(
    values: np.ndarray,
    native_unit: str,
    target_unit: str,
    formula_units_per_cell: float = 1.0,
) -> np.ndarray:
    """Convert entropy or heat-capacity values between display units.

    Parameters
    ----------
    values : ndarray
        Native values.
    native_unit, target_unit : str
        Source and destination compound units.
    formula_units_per_cell : float, optional
        Formula-unit normalization factor.

    Returns
    -------
    ndarray
        Converted values.
    """
    if _same_unit(native_unit, target_unit):
        return values
    native_energy, native_temperature = _split_heat_capacity_unit(native_unit)
    target_energy, target_temperature = _split_heat_capacity_unit(target_unit)
    energy_factor = _energy_factor_to_j_per_mol(
        native_energy
    ) / _energy_factor_to_j_per_mol(target_energy)
    temperature_factor = _temperature_interval_factor_to_kelvin(
        target_temperature
    ) / _temperature_interval_factor_to_kelvin(native_temperature)
    normalization_factor = _molar_normalization_factor(
        native_energy, target_energy, formula_units_per_cell
    )
    return values * energy_factor * temperature_factor * normalization_factor


def _molar_normalization_factor(
    native_unit: str,
    target_unit: str,
    formula_units_per_cell: float,
) -> float:
    """Return the cell-to-formula-unit normalization factor.

    Parameters
    ----------
    native_unit, target_unit : str
        Source and destination energy labels.
    formula_units_per_cell : float
        Formula units contained in the computational cell.

    Returns
    -------
    float
        Multiplicative normalization factor.
    """
    native_is_molar = "mol" in native_unit.strip().lower()
    target_is_molar = "mol" in target_unit.strip().lower()
    if native_is_molar == target_is_molar:
        return 1.0
    if formula_units_per_cell <= 0.0:
        formula_units_per_cell = 1.0
    if target_is_molar:
        return 1.0 / formula_units_per_cell
    return formula_units_per_cell


def _heat_capacity_unit(energy_unit: object, temperature_unit: object) -> str:
    """Build a heat-capacity unit label.

    Parameters
    ----------
    energy_unit, temperature_unit : object
        Energy and temperature unit labels.

    Returns
    -------
    str
        Compound heat-capacity unit label.
    """
    return (
        f"{energy_unit}/mol/{temperature_unit}"
        if str(energy_unit) in {"J", "kJ"}
        else f"{energy_unit}/{temperature_unit}"
    )


def _natoms_for_dulong_petit(result: QHAResult) -> float | None:
    """Return the atom count used for the Dulong-Petit limit.

    Parameters
    ----------
    result : QHAResult
        QHA result with normalization metadata.

    Returns
    -------
    float or None
        Positive atom count, or ``None`` when no reliable value is available.
    """
    metadata = result.metadata or {}
    normalization = metadata.get("normalization")
    candidates: tuple[object, ...] = (metadata.get("natoms"),)
    if isinstance(normalization, Mapping):
        candidates = candidates + (normalization.get("natoms_per_cell"),)
    candidates = candidates + (
        metadata.get("natoms_per_formula_unit"),
        metadata.get("formula_natoms"),
    )
    input_metadata = metadata.get("input")
    if isinstance(input_metadata, Mapping):
        candidates = candidates + (
            input_metadata.get("natoms"),
            input_metadata.get("natoms_per_formula_unit"),
            input_metadata.get("formula_natoms"),
        )
    for candidate in candidates:
        try:
            value = float(candidate)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            continue
        if value > 0:
            return value
    return None


def _dulong_petit_limit(result: QHAResult, context: _PlotUnitContext) -> float | None:
    """Return the Dulong-Petit limit in display heat-capacity units.

    Parameters
    ----------
    result : QHAResult
        QHA result with atom-count metadata.
    context : _PlotUnitContext
        Native and display unit context.

    Returns
    -------
    float or None
        Converted Dulong-Petit limit, or ``None`` when unavailable.
    """
    natoms = _natoms_for_dulong_petit(result)
    if natoms is None:
        return None
    native_unit = str(
        context.native_units.get(
            "heat_capacity",
            f"{context.native_units.get('energy', 'Ha')}/K",
        )
    )
    target_unit = str(context.display_units.get("heat_capacity", native_unit))
    value_j_mol_cell_k = 3.0 * natoms * cs.R
    native_energy, native_temperature = _split_heat_capacity_unit(native_unit)
    native_value = value_j_mol_cell_k / _energy_factor_to_j_per_mol(native_energy)
    native_value *= _temperature_interval_factor_to_kelvin(
        native_temperature
    ) / _temperature_interval_factor_to_kelvin("K")
    converted = _convert_heat_capacity_like(
        np.asarray(native_value, dtype=np.float64),
        native_unit,
        target_unit,
        context.formula_units_per_cell,
    )
    return float(converted)


def _dulong_petit_series(
    result: QHAResult,
    context: _PlotUnitContext,
    temperature: np.ndarray,
) -> PlotSeries | None:
    """Return a neutral Dulong-Petit reference series when available.

    Parameters
    ----------
    result : QHAResult
        QHA result object.
    context : _PlotUnitContext
        Resolved plotting units.
    temperature : ndarray
        Display-temperature coordinates.

    Returns
    -------
    PlotSeries or None
        Constant reference line, or ``None`` if no reliable limit is available.
    """
    limit = _dulong_petit_limit(result, context)
    if limit is None or not np.isfinite(limit):
        return None
    x = np.asarray([temperature[0], temperature[-1]], dtype=np.float64)
    y = np.asarray([limit, limit], dtype=np.float64)
    return PlotSeries(
        key="dulong_petit",
        label="Dulong-Petit limit",
        x=x,
        y=y,
        style=PlotSeriesStyle(color="black", line_style="dotted", line_width=1.0),
        metadata={"reference": True},
    )


def _upper_limit(values: np.ndarray, reference: float) -> tuple[None, float] | None:
    """Return an upper y-axis limit including a reference value.

    Parameters
    ----------
    values : ndarray
        Plotted property values.
    reference : float
        Reference value that must remain visible.

    Returns
    -------
    tuple or None
        Optional ``(None, upper)`` y-axis limits.
    """
    try:
        ymax = float(np.nanmax(values))
    except ValueError:
        return None
    if not np.isfinite(ymax):
        return None
    return None, max(ymax, float(reference)) * 1.08


def _temperature_label(units: Mapping[str, object]) -> str:
    """Return the temperature axis label.

    Parameters
    ----------
    units : mapping
        Display unit metadata.

    Returns
    -------
    str
        Temperature axis label.
    """
    return f"Temperature, T ({units.get('temperature', 'K')})"


def _pressure_axis_label(units: Mapping[str, object]) -> str:
    """Return the pressure axis label.

    Parameters
    ----------
    units : mapping
        Display unit metadata.

    Returns
    -------
    str
        Pressure axis label.
    """
    return f"Pressure, P ({units.get('pressure', 'GPa')})"


def _pressure_label(value: float, units: Mapping[str, object]) -> str:
    """Return a compact pressure legend label.

    Parameters
    ----------
    value : float
        Pressure value in display units.
    units : mapping
        Display unit metadata.

    Returns
    -------
    str
        Pressure legend label.
    """
    return f"P = {value:g} {units.get('pressure', 'GPa')}"


def _temperature_value_label(value: float, units: Mapping[str, object]) -> str:
    """Return a compact temperature legend label."""
    return f"T = {value:g} {units.get('temperature', 'K')}"
