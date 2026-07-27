# -*- coding: utf-8 -*-

"""Neutral plot-specification builders for harmonic results."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Literal, cast

import numpy as np

from quantas.models import (
    ContourPlotSpec,
    LinePlotSpec,
    PlotAxis,
    PlotCollection,
    PlotSeries,
)
from quantas.modules.ha.io.export import convert_property_values
from quantas.modules.ha.models import HAResult
from quantas.modules.ha.plot.labels import (
    available_plot_properties,
    format_property_axis_label,
    format_temperature_axis_label,
    format_volume_legend_label,
    resolve_plot_property,
)


_HA_PLOT_PROPERTIES = available_plot_properties()
PLOT_LABELS = {
    attribute: (property_info.key, property_info.name)
    for attribute, property_info in _HA_PLOT_PROPERTIES.items()
}
PLOT_PROPERTIES = tuple(
    property_info.key for property_info in _HA_PLOT_PROPERTIES.values()
)

DEFAULT_PLOT_PROPERTIES = ("F", "Utot", "S", "Cv")

HACurveAxis = Literal["temperature", "volume"]

_ALLOWED_CMAPS = ("viridis", "plasma", "inferno", "magma", "cividis", "turbo")
_ALLOWED_CONTOUR_MODES = ("discrete", "smooth")


@dataclass(slots=True)
class HAPlotOptions:
    """Options controlling neutral HA section and map construction.

    Parameters
    ----------
    curve_axis : {"temperature", "volume"}, optional
        Natural variable placed on the independent axis of line plots.
        Temperature curves contain one series per selected sampled volume;
        volume curves contain one series per selected stored temperature.
    selected_volumes : tuple of float or None, optional
        Exact sampled volumes to include in temperature curves, expressed in
        the native volume unit stored in the result. ``None`` selects all.
    selected_temperatures : tuple of float or None, optional
        Exact stored temperatures to include in volume curves, expressed in
        the native temperature unit stored in the result. ``None`` selects all.
    include_contours : bool, optional
        Build volume-temperature contour maps in addition to line sections
        when both native grids contain at least two points.
    energy_unit : str or None, optional
        Display unit for energy-like properties.
    cmap : str, optional
        Portable colormap hint for contour renderers.
    contour_mode : {"discrete", "smooth"}, optional
        Preferred filled-contour mode.
    levels : int, optional
        Number of principal contour levels.
    isolines : bool, optional
        Whether contour isolines should be requested.
    isoline_labels : bool, optional
        Whether contour-line labels should be requested.

    Notes
    -----
    Selection is exact and never interpolates the stored HA grid.
    """

    curve_axis: HACurveAxis = "temperature"
    selected_volumes: tuple[float, ...] | None = None
    selected_temperatures: tuple[float, ...] | None = None
    include_contours: bool = False
    energy_unit: str | None = None
    cmap: str = "viridis"
    contour_mode: str = "smooth"
    levels: int = 12
    isolines: bool = True
    isoline_labels: bool = True


def resolve_plot_properties(
    properties: str | list[str] | tuple[str, ...] | None,
) -> tuple[str, ...]:
    """Resolve a user-facing HA plot-property selection.

    Parameters
    ----------
    properties : str, list of str, tuple of str, or None
        Requested plot properties. ``"all"`` expands to all available HA
        properties. If ``None``, the compact default set is returned.

    Returns
    -------
    tuple of str
        Requested property keys.
    """
    if properties is None:
        return DEFAULT_PLOT_PROPERTIES
    if isinstance(properties, str):
        items = [item.strip() for item in properties.split(",") if item.strip()]
    else:
        items = [str(item).strip() for item in properties if str(item).strip()]
    if not items:
        return DEFAULT_PLOT_PROPERTIES
    if len(items) == 1 and items[0].lower() == "all":
        return PLOT_PROPERTIES
    return tuple(items)


def build_thermodynamic_plot_spec(
    result: HAResult,
    property_name: str,
    *,
    options: HAPlotOptions | None = None,
    show_legend: bool = True,
    unit: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
) -> LinePlotSpec:
    """Build a neutral HA thermodynamic line-plot specification.

    Parameters
    ----------
    result : HAResult
        Harmonic result containing temperature and property arrays.
    property_name : str
        Result attribute name or Quantas short property key.
    options : HAPlotOptions or None, optional
        Scientific section selection and portable contour/style hints.
    show_legend : bool, optional
        Whether a volume legend should be requested for multi-volume data.
    unit : str or None, optional
        Display energy unit for energy-like values.
    xlabel, ylabel, title : str or None, optional
        Optional presentation overrides.

    Returns
    -------
    LinePlotSpec
        Frontend-neutral plot specification.

    Raises
    ------
    KeyError
        If the property name is unknown.
    ValueError
        If required data are missing or have incompatible shapes.
    """
    opts = _validated_options(options or HAPlotOptions(), unit=unit)
    attr = _resolve_plot_property_name(property_name)
    _, description = PLOT_LABELS[attr]
    temperature = _temperature_array(result)
    values = _property_array(result, attr, temperature.size)
    units = result.metadata.get("units", {})
    values, property_unit = convert_property_values(
        value=values,
        attr=attr,
        units=units,
        target_unit=opts.energy_unit,
        normalization=result.metadata.get("normalization"),
    )

    volume_unit = str(units.get("volume", "unknown"))
    temperature_unit = str(units.get("temperature", "K"))
    volume = _matched_volume_array(result, values.shape[1])

    if opts.curve_axis == "temperature":
        indices = _selected_indices(
            volume,
            opts.selected_volumes,
            name="volume",
            allow_missing=opts.selected_volumes is None,
            size=values.shape[1],
        )
        labels = _volume_labels(result.volume, values.shape[1])
        series = []
        for index in indices:
            metadata: dict[str, float | int | str] = {
                "volume_index": int(index),
                "curve_axis": "temperature",
            }
            if volume is not None:
                metadata["volume"] = float(volume[index])
                metadata["volume_native"] = float(volume[index])
            series.append(
                PlotSeries(
                    key=f"volume_{index}",
                    label=labels[index],
                    x=temperature.copy(),
                    y=values[:, index].copy(),
                    metadata=metadata,
                )
            )
        x_axis = PlotAxis(
            key="temperature",
            label=xlabel or format_temperature_axis_label(temperature_unit),
            unit=temperature_unit,
        )
        legend_title = format_volume_legend_label(volume_unit)
        filename_stem = attr
    else:
        if volume is None:
            raise ValueError(
                f"HA property '{attr}' cannot be plotted against volume because "
                "its sampled-volume coordinates are unavailable or incompatible"
            )
        if volume.size < 2:
            raise ValueError(
                f"HA property '{attr}' requires at least two sampled volumes "
                "for a volume section"
            )
        indices = _selected_indices(
            temperature,
            opts.selected_temperatures,
            name="temperature",
            allow_missing=False,
            size=temperature.size,
        )
        series = [
            PlotSeries(
                key=f"temperature_{index}",
                label=f"{temperature[index]:.6g}",
                x=volume.copy(),
                y=values[index, :].copy(),
                metadata={
                    "temperature_index": int(index),
                    "temperature": float(temperature[index]),
                    "temperature_native": float(temperature[index]),
                    "curve_axis": "volume",
                },
            )
            for index in indices
        ]
        x_axis = PlotAxis(
            key="volume",
            label=xlabel or format_volume_legend_label(volume_unit),
            unit=volume_unit,
        )
        legend_title = format_temperature_axis_label(temperature_unit)
        filename_stem = f"{attr}_vs_volume"

    return LinePlotSpec(
        key=attr,
        title=title or description,
        filename_stem=filename_stem,
        x_axis=x_axis,
        y_axis=PlotAxis(
            key=attr,
            label=ylabel or format_property_axis_label(attr, property_unit),
            unit=property_unit,
        ),
        series=series,
        legend_title=legend_title,
        show_legend=show_legend and len(series) > 1,
        metadata={
            "module": "ha",
            "property": attr,
            "symbol": PLOT_LABELS[attr][0],
            "curve_axis": opts.curve_axis,
        },
    )


def build_thermodynamic_contour_spec(
    result: HAResult,
    property_name: str,
    *,
    options: HAPlotOptions | None = None,
    unit: str | None = None,
) -> ContourPlotSpec:
    """Build one exact-grid volume-temperature contour specification.

    No interpolation or resampling is performed. Temperature-independent
    properties are represented by their broadcast values on the stored grid.
    """
    opts = _validated_options(options or HAPlotOptions(), unit=unit)
    attr = _resolve_plot_property_name(property_name)
    _, description = PLOT_LABELS[attr]
    temperature = _temperature_array(result)
    values = _property_array(result, attr, temperature.size)
    volume = _matched_volume_array(result, values.shape[1])
    if volume is None:
        raise ValueError(
            f"HA contour for '{attr}' requires sampled-volume coordinates "
            "matching the property grid"
        )
    if temperature.size < 2 or volume.size < 2:
        raise ValueError(
            f"HA contour for '{attr}' requires at least two temperatures "
            "and two sampled volumes"
        )

    units = result.metadata.get("units", {})
    converted, property_unit = convert_property_values(
        value=values,
        attr=attr,
        units=units,
        target_unit=opts.energy_unit,
        normalization=result.metadata.get("normalization"),
    )
    temperature_unit = str(units.get("temperature", "K"))
    volume_unit = str(units.get("volume", "unknown"))
    return ContourPlotSpec(
        key=attr,
        title=description,
        filename_stem=f"{attr}_VT",
        x_axis=PlotAxis(
            key="temperature",
            label=format_temperature_axis_label(temperature_unit),
            unit=temperature_unit,
        ),
        y_axis=PlotAxis(
            key="volume",
            label=format_volume_legend_label(volume_unit),
            unit=volume_unit,
        ),
        value_axis=PlotAxis(
            key=attr,
            label=format_property_axis_label(attr, property_unit),
            unit=property_unit,
        ),
        x=temperature.copy(),
        y=volume.copy(),
        z=converted.T.copy(),
        colormap=opts.cmap,
        mode=cast(Literal["discrete", "smooth"], opts.contour_mode),
        levels=opts.levels,
        isolines=opts.isolines,
        isoline_labels=opts.isoline_labels,
        metadata={
            "module": "ha",
            "property": attr,
            "representation": "volume_temperature_contour",
        },
    )


def build_ha_plot_collection(
    result: HAResult,
    properties: str | list[str] | tuple[str, ...] | None = None,
    *,
    unit: str | None = None,
    options: HAPlotOptions | None = None,
) -> PlotCollection:
    """Build all requested neutral HA plot specifications.

    Parameters
    ----------
    result : HAResult
        Harmonic result object.
    properties : str, list of str, tuple of str, or None, optional
        Requested property keys.
    unit : str or None, optional
        Display energy unit for energy-like values.
    options : HAPlotOptions or None, optional
        Scientific section selection and contour-preparation options.

    Returns
    -------
    PlotCollection
        Prepared plot specifications and warnings for skipped properties.

    Raises
    ------
    ValueError
        If none of the requested properties can be prepared.
    """
    opts = _validated_options(options or HAPlotOptions(), unit=unit)
    collection = PlotCollection()
    for property_name in resolve_plot_properties(properties):
        try:
            collection.plots.append(
                build_thermodynamic_plot_spec(
                    result,
                    property_name,
                    options=opts,
                )
            )
        except (KeyError, ValueError) as exc:
            collection.warnings.append(f"{property_name}: {exc}")
            continue
        if opts.include_contours:
            try:
                collection.plots.append(
                    build_thermodynamic_contour_spec(
                        result,
                        property_name,
                        options=opts,
                    )
                )
            except (KeyError, ValueError) as exc:
                collection.warnings.append(f"{property_name}: {exc}")

    if not collection.plots:
        details = "; ".join(collection.warnings)
        raise ValueError(f"No HA plots could be prepared. {details}")
    return collection


def _resolve_plot_property_name(name: str) -> str:
    """Resolve a HA result attribute or stable short plotting key."""
    return resolve_plot_property(name).attribute


def _temperature_array(result: HAResult) -> np.ndarray:
    """Return the validated HA temperature array.

    Parameters
    ----------
    result : HAResult
        Harmonic result object.

    Returns
    -------
    ndarray
        One-dimensional temperature array.

    Raises
    ------
    ValueError
        If temperatures are missing or not one-dimensional.
    """
    if result.temperature is None:
        raise ValueError("HA result temperatures are not available")
    temperature = np.asarray(result.temperature, dtype=np.float64)
    if temperature.ndim != 1 or temperature.size == 0:
        raise ValueError(
            "HA result temperatures must be a non-empty one-dimensional array"
        )
    if np.unique(temperature).size != temperature.size:
        raise ValueError("HA result temperatures contain duplicate coordinates")
    return temperature


def _property_array(result: HAResult, attr: str, ntemp: int) -> np.ndarray:
    """Return a HA property array with shape ``(ntemp, nvol)``.

    Parameters
    ----------
    result : HAResult
        Harmonic result object.
    attr : str
        Property attribute name.
    ntemp : int
        Expected number of temperatures.

    Returns
    -------
    ndarray
        Two-dimensional property array.

    Raises
    ------
    ValueError
        If the property is unavailable or has an unsupported shape.
    """
    values = getattr(result, attr)
    if values is None:
        raise ValueError(f"HA result property '{attr}' is not available")
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 1:
        if attr in {"static_energy", "zero_point_energy"}:
            return np.tile(array[np.newaxis, :], (ntemp, 1))
        if array.shape[0] == ntemp:
            return array[:, np.newaxis]
        return np.tile(array[np.newaxis, :], (ntemp, 1))
    if array.ndim == 2:
        if array.shape[0] == ntemp:
            return array
        if attr in {"static_energy", "zero_point_energy"} and array.shape[0] == 1:
            return np.tile(array, (ntemp, 1))
        raise ValueError(
            f"HA result property '{attr}' has incompatible temperature dimension"
        )
    raise ValueError("HA thermodynamic properties must be one- or two-dimensional")


def _matched_volume_array(result: HAResult, nvol: int) -> np.ndarray | None:
    """Return sampled volumes only when they match one property grid."""
    if result.volume is None:
        return None
    volume = np.asarray(result.volume, dtype=np.float64)
    if volume.ndim != 1 or volume.size != nvol:
        return None
    if np.unique(volume).size != volume.size:
        return None
    return volume


def _selected_indices(
    grid: np.ndarray | None,
    selected: tuple[float, ...] | None,
    *,
    name: str,
    allow_missing: bool,
    size: int,
) -> tuple[int, ...]:
    """Resolve exact native-grid selections without interpolation."""
    if selected is None:
        return tuple(range(size))
    if not selected:
        raise ValueError(f"selected {name} values cannot be empty")
    if len(set(selected)) != len(selected):
        raise ValueError(f"selected {name} values must be unique")
    if grid is None:
        if allow_missing:
            return tuple(range(size))
        raise ValueError(f"HA {name} coordinates are not available")

    indices: list[int] = []
    for requested in selected:
        matches = np.flatnonzero(grid == float(requested))
        if matches.size == 0:
            available = ", ".join(f"{value:.12g}" for value in grid)
            raise ValueError(
                f"HA {name} {requested!r} is not present in the native grid; "
                f"available values: {available}"
            )
        indices.append(int(matches[0]))
    return tuple(indices)


def _validated_options(
    options: HAPlotOptions,
    *,
    unit: str | None = None,
) -> HAPlotOptions:
    """Validate scientific selections and portable contour hints."""
    if options.curve_axis not in {"temperature", "volume"}:
        raise ValueError("HA curve_axis must be 'temperature' or 'volume'")
    if options.curve_axis == "temperature" and options.selected_temperatures:
        raise ValueError(
            "selected_temperatures is valid only when curve_axis='volume'"
        )
    if options.curve_axis == "volume" and options.selected_volumes:
        raise ValueError("selected_volumes is valid only when curve_axis='temperature'")
    if options.cmap not in _ALLOWED_CMAPS:
        raise ValueError(f"unsupported colormap '{options.cmap}'")
    if options.contour_mode not in _ALLOWED_CONTOUR_MODES:
        raise ValueError(f"unsupported contour mode '{options.contour_mode}'")
    if options.levels < 2:
        raise ValueError("contour levels must be at least 2")
    if unit is not None:
        if options.energy_unit is not None and options.energy_unit != unit:
            raise ValueError(
                "unit and options.energy_unit must agree when both are set"
            )
        return replace(options, energy_unit=unit)
    return options


def _volume_labels(volume: np.ndarray | None, nvol: int) -> list[str]:
    """Return labels for volume-dependent HA curves.

    Parameters
    ----------
    volume : ndarray or None
        Unit-cell volumes.
    nvol : int
        Number of curves.

    Returns
    -------
    list of str
        Curve labels.
    """
    if volume is None:
        return [f"V{index}" for index in range(nvol)]
    volume_array = np.asarray(volume, dtype=np.float64)
    if volume_array.ndim != 1 or volume_array.size != nvol:
        return [f"V{index}" for index in range(nvol)]
    return [f"{value:.6g}" for value in volume_array]
