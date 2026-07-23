# -*- coding: utf-8 -*-

"""Neutral plot-specification builders for harmonic results."""

from __future__ import annotations

import numpy as np

from quantas.models import LinePlotSpec, PlotAxis, PlotCollection, PlotSeries
from quantas.modules.ha.io.export import convert_property_values
from quantas.modules.ha.models import HAResult
from quantas.modules.ha.plot.labels import (
    format_property_axis_label,
    format_temperature_axis_label,
    format_volume_legend_label,
)
from quantas.modules.ha.report import THERMODYNAMIC_LABELS


PLOT_LABELS = {
    "static_energy": ("U0", "Static energy"),
    **THERMODYNAMIC_LABELS,
}

PLOT_PROPERTIES = (
    "U0",
    "Uzp",
    "Uth",
    "Utot",
    "S",
    "Fvib",
    "F",
    "Cv",
)

DEFAULT_PLOT_PROPERTIES = ("F", "Utot", "S", "Cv")


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
    attr = _resolve_plot_property_name(property_name)
    _, description = PLOT_LABELS[attr]
    temperature = _temperature_array(result)
    values = _property_array(result, attr, temperature.size)
    units = result.metadata.get("units", {})
    values, property_unit = convert_property_values(
        value=values,
        attr=attr,
        units=units,
        target_unit=unit,
        normalization=result.metadata.get("normalization"),
    )

    volume_unit = str(units.get("volume", "unknown"))
    temperature_unit = str(units.get("temperature", "K"))
    labels = _volume_labels(result.volume, values.shape[1])
    series = [
        PlotSeries(
            key=f"volume_{index}",
            label=label,
            x=temperature.copy(),
            y=values[:, index].copy(),
            metadata={"volume_index": index},
        )
        for index, label in enumerate(labels)
    ]

    return LinePlotSpec(
        key=attr,
        title=title or description,
        filename_stem=attr,
        x_axis=PlotAxis(
            key="temperature",
            label=xlabel or format_temperature_axis_label(temperature_unit),
            unit=temperature_unit,
        ),
        y_axis=PlotAxis(
            key=attr,
            label=ylabel or format_property_axis_label(attr, property_unit),
            unit=property_unit,
        ),
        series=series,
        legend_title=format_volume_legend_label(volume_unit),
        show_legend=show_legend and values.shape[1] > 1,
        metadata={
            "module": "ha",
            "property": attr,
            "symbol": PLOT_LABELS[attr][0],
        },
    )


def build_ha_plot_collection(
    result: HAResult,
    properties: str | list[str] | tuple[str, ...] | None = None,
    *,
    unit: str | None = None,
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

    Returns
    -------
    PlotCollection
        Prepared plot specifications and warnings for skipped properties.

    Raises
    ------
    ValueError
        If none of the requested properties can be prepared.
    """
    collection = PlotCollection()
    for property_name in resolve_plot_properties(properties):
        try:
            collection.plots.append(
                build_thermodynamic_plot_spec(result, property_name, unit=unit)
            )
        except (KeyError, ValueError) as exc:
            collection.warnings.append(f"{property_name}: {exc}")

    if not collection.plots:
        details = "; ".join(collection.warnings)
        raise ValueError(f"No HA plots could be prepared. {details}")
    return collection


def _resolve_plot_property_name(name: str) -> str:
    """Resolve a HA plot key to a result attribute name.

    Parameters
    ----------
    name : str
        Result attribute name or Quantas short key.

    Returns
    -------
    str
        Result attribute name.

    Raises
    ------
    KeyError
        If the property is unknown.
    """
    if name in PLOT_LABELS:
        return name
    for attr, (symbol, _) in PLOT_LABELS.items():
        if name == symbol:
            return attr
    raise KeyError(f"unknown HA plot property: {name}")


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
    if temperature.ndim != 1:
        raise ValueError("HA result temperatures must be one-dimensional")
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
