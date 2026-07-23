# -*- coding: utf-8 -*-

"""Frontend-neutral EOS plot builders."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from quantas.models import (
    LinePlotSpec,
    PlotAxis,
    PlotCollection,
    PlotSeries,
    PlotSeriesStyle,
)
from quantas.modules.eos.calculator import EOSCalculator
from quantas.modules.eos.diagnostics import EOSDiagnosticResult, EOSDiagnostics
from quantas.modules.eos.history import EOSFitRecord, EOSResultSlot
from quantas.modules.eos.models import EOSDataset, EOSFitDomain

from .labels import format_unit, normalized_pressure_labels, property_label


EOS_PLOT_TYPES: tuple[str, ...] = (
    "fit",
    "residuals",
    "standardized-residuals",
    "normalized-pressure",
    "coverage",
    "isotherms",
    "isobars",
)


@dataclass(slots=True)
class EOSPlotOptions:
    """Control preparation of frontend-neutral EOS plot specifications.

    Parameters
    ----------
    show_uncertainties : bool, optional
        Include symmetric one-sigma error bars from the input dataset.
    show_excluded : bool, optional
        Display observations excluded from the fit.
    group_data : bool, optional
        Split included observations into one series per input group.
    point_size : float, optional
        Marker size in typographic points.
    curve_width : float, optional
        Width of fitted and calculated curves.
    errorbar_width : float, optional
        Width of error-bar strokes.
    errorbar_capsize : float, optional
        Error-bar cap size in typographic points.
    curve_points : int, optional
        Number of samples used for smooth fitted curves.
    curve_color : str or None, optional
        Portable color specification for fitted curves.
    excluded_color : str, optional
        Portable color specification for excluded observations.
    show_legend : bool, optional
        Request legends when multiple series are present.
    show_title : bool, optional
        Preserve human-readable figure titles in the neutral specifications.
        Disable this for publication-style figures that use captions instead.
    show_zero_pressure_point : bool, optional
        Display finite normalized-pressure observations whose measured
        pressure is statistically indistinguishable from zero. Such points are
        hidden by default because they carry little order-diagnostic
        information and can dominate the plotted strain range.
    grid : bool, optional
        Request Cartesian grids.
    isotherms, isobars : tuple of float, optional
        Explicit P--V--T curve coordinates. Empty tuples select representative
        values from the sampled data range.
    max_default_curves : int, optional
        Maximum number of automatically generated isotherms or isobars.
    """

    show_uncertainties: bool = True
    show_excluded: bool = True
    group_data: bool = True
    point_size: float = 5.0
    curve_width: float = 1.8
    errorbar_width: float = 1.0
    errorbar_capsize: float = 2.5
    curve_points: int = 300
    curve_color: str | None = "black"
    excluded_color: str = "0.55"
    show_legend: bool = True
    show_title: bool = True
    show_zero_pressure_point: bool = False
    grid: bool = True
    isotherms: tuple[float, ...] = ()
    isobars: tuple[float, ...] = ()
    max_default_curves: int = 5

    def __post_init__(self) -> None:
        if self.point_size <= 0.0:
            raise ValueError("point_size must be positive")
        if self.curve_width <= 0.0:
            raise ValueError("curve_width must be positive")
        if self.errorbar_width <= 0.0:
            raise ValueError("errorbar_width must be positive")
        if self.errorbar_capsize < 0.0:
            raise ValueError("errorbar_capsize must be non-negative")
        if self.curve_points < 20:
            raise ValueError("curve_points must be at least 20")
        if self.max_default_curves < 1:
            raise ValueError("max_default_curves must be positive")
        self.isotherms = tuple(float(item) for item in self.isotherms)
        self.isobars = tuple(float(item) for item in self.isobars)


class EOSPlotter:
    """Build EOS plots from one successful immutable fit record."""

    def __init__(self, record: EOSFitRecord, dataset: EOSDataset) -> None:
        if not record.result.fit.success or record.result.fit.parameters is None:
            raise ValueError("EOS plotting requires a successful fit record")
        self.record = record
        self.dataset = dataset
        self._diagnostics = EOSDiagnostics(record, dataset)
        self._calculator = EOSCalculator(record, dataset)

    @classmethod
    def from_archive(
        cls,
        path: str | Path,
        *,
        slot: str | EOSResultSlot | None = None,
        record_id: int | None = None,
    ) -> "EOSPlotter":
        """Construct a plotter from an explicit or accepted archive record."""
        calculator = EOSCalculator.from_archive(path, slot=slot, record_id=record_id)
        return cls(calculator.record, calculator.dataset)

    def available_plot_types(self) -> tuple[str, ...]:
        """Return plot types supported by the selected scientific domain."""
        domain = self.record.request.domain
        try:
            diagnostic = self._diagnostics.build(include_normalized_pressure=True)
        except (ValueError, NotImplementedError):
            diagnostic = None
        standardized_available = bool(
            diagnostic is not None
            and "standardized_residual" in diagnostic.columns
            and np.any(np.isfinite(diagnostic.columns["standardized_residual"]))
        )
        if domain is EOSFitDomain.PRESSURE_VOLUME:
            kinds = ["fit", "residuals"]
            if standardized_available:
                kinds.append("standardized-residuals")
            if diagnostic is not None and diagnostic.metadata[
                "normalized_pressure"
            ].get("available", False):
                kinds.append("normalized-pressure")
            return tuple(kinds)
        if domain is EOSFitDomain.VOLUME_TEMPERATURE:
            kinds = ["fit", "residuals"]
            if standardized_available:
                kinds.append("standardized-residuals")
            return tuple(kinds)
        if domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
            kinds = ["coverage", "isotherms", "isobars", "residuals"]
            if standardized_available:
                kinds.append("standardized-residuals")
            return tuple(kinds)
        return ()

    def build(
        self,
        plot_types: Sequence[str] | str | None = None,
        *,
        options: EOSPlotOptions | None = None,
    ) -> PlotCollection:
        """Build selected neutral EOS plot specifications.

        Parameters
        ----------
        plot_types : sequence of str, str, or None, optional
            Requested types. ``None`` and ``"all"`` select all types available
            for the current domain.
        options : EOSPlotOptions or None, optional
            Plot-preparation options.

        Returns
        -------
        PlotCollection
            Neutral specifications and non-fatal warnings.
        """
        opts = options or EOSPlotOptions()
        available = self.available_plot_types()
        selected = self._resolve_plot_types(plot_types, available)
        diagnostic = self._diagnostics.build(include_normalized_pressure=True)
        collection = PlotCollection(warnings=list(diagnostic.warnings))
        for kind in selected:
            try:
                prepared = self._build_kind(kind, diagnostic, opts)
                collection.plots.extend(prepared)
                for spec in prepared:
                    collection.warnings.extend(
                        str(item) for item in spec.metadata.get("warnings", [])
                    )
            except (ValueError, RuntimeError, NotImplementedError) as exc:
                collection.warnings.append(f"{kind}: {exc}")
        if not collection.plots:
            details = "; ".join(collection.warnings) or "no compatible plot types"
            raise ValueError(f"No EOS plots could be prepared. {details}")
        if not opts.show_title:
            for plot_spec in collection.plots:
                plot_spec.title = ""
        return collection

    def _build_kind(
        self,
        kind: str,
        diagnostic: EOSDiagnosticResult,
        options: EOSPlotOptions,
    ) -> list[LinePlotSpec]:
        if kind == "fit":
            return [self._build_fit_spec(diagnostic, options)]
        if kind == "residuals":
            return self._build_residual_specs(diagnostic, options, standardized=False)
        if kind == "standardized-residuals":
            return self._build_residual_specs(diagnostic, options, standardized=True)
        if kind == "normalized-pressure":
            return [self._build_normalized_pressure_spec(diagnostic, options)]
        if kind == "coverage":
            return [self._build_coverage_spec(diagnostic, options)]
        if kind == "isotherms":
            return [self._build_isotherm_spec(diagnostic, options)]
        if kind == "isobars":
            return [self._build_isobar_spec(diagnostic, options)]
        raise ValueError(f"unknown EOS plot type: {kind}")

    def _build_fit_spec(
        self,
        diagnostic: EOSDiagnosticResult,
        options: EOSPlotOptions,
    ) -> LinePlotSpec:
        request = self.record.request
        if request.domain is EOSFitDomain.PRESSURE_VOLUME:
            x_name = request.target
            y_name = "pressure"
            x = diagnostic.columns[x_name]
            y = diagnostic.columns["observed_pressure"]
            curve_x = self._curve_grid(
                x,
                diagnostic.columns["included"],
                points=options.curve_points,
            )
            calculated = self._calculator.calculate(
                volume=curve_x,
                propagate_uncertainty=False,
            )
            curve_y = calculated.columns["pressure"]
            x_error = self._input_error(x_name)
            y_error = self._input_error("pressure")
        elif request.domain is EOSFitDomain.VOLUME_TEMPERATURE:
            x_name = "temperature"
            y_name = request.target
            x = diagnostic.columns["temperature"]
            y = diagnostic.columns[f"observed_{request.target}"]
            curve_x = self._curve_grid(
                x,
                diagnostic.columns["included"],
                positive=True,
                points=options.curve_points,
            )
            calculated = self._calculator.calculate(
                temperature=curve_x,
                propagate_uncertainty=False,
            )
            curve_y = calculated.columns[request.target]
            x_error = self._input_error("temperature")
            y_error = self._input_error(request.target)
        else:
            raise ValueError("the generic fit plot is available only for P-V and V-T")

        series = [self._fit_curve(curve_x, curve_y, options)]
        series.extend(
            self._observation_series(
                x,
                y,
                diagnostic,
                options,
                x_error=x_error,
                y_error=y_error,
            )
        )
        return self._line_spec(
            key="fit",
            title=f"{self._title_name()}: fitted EOS",
            x_name=x_name,
            x_unit=diagnostic.units.get(x_name),
            y_name=y_name,
            y_unit=diagnostic.units.get(
                "observed_pressure" if y_name == "pressure" else f"observed_{y_name}"
            ),
            series=series,
            options=options,
        )

    def _build_residual_specs(
        self,
        diagnostic: EOSDiagnosticResult,
        options: EOSPlotOptions,
        *,
        standardized: bool,
    ) -> list[LinePlotSpec]:
        request = self.record.request
        y_name = "standardized_residual" if standardized else "residual"
        y = diagnostic.columns[y_name]
        if not np.any(np.isfinite(y)):
            raise ValueError(f"{y_name.replace('_', ' ')} is unavailable")
        if request.domain is EOSFitDomain.PRESSURE_VOLUME:
            x_names = ["observed_pressure"]
        elif request.domain is EOSFitDomain.VOLUME_TEMPERATURE:
            x_names = ["temperature"]
        elif request.domain is EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
            x_names = ["observed_pressure", "temperature"]
        else:
            raise ValueError("residual plots are unavailable for this domain")
        specs: list[LinePlotSpec] = []
        for x_column in x_names:
            axis_name = "pressure" if x_column == "observed_pressure" else x_column
            x = diagnostic.columns[x_column]
            series = self._observation_series(
                x,
                y,
                diagnostic,
                options,
                x_error=None,
                y_error=None,
            )
            finite = np.isfinite(x) & np.isfinite(y)
            if np.any(finite):
                xmin = float(np.min(x[finite]))
                xmax = float(np.max(x[finite]))
                series.insert(
                    0,
                    PlotSeries(
                        key="zero",
                        label="Zero",
                        x=np.array([xmin, xmax], dtype=np.float64),
                        y=np.zeros(2, dtype=np.float64),
                        style=PlotSeriesStyle(
                            color="black",
                            line_style="dotted",
                            line_width=max(0.8, 0.6 * options.curve_width),
                        ),
                    ),
                )
            specs.append(
                self._line_spec(
                    key=(
                        f"standardized_residuals_vs_{axis_name}"
                        if standardized
                        else f"residuals_vs_{axis_name}"
                    ),
                    title=(
                        f"{self._title_name()}: standardized residuals"
                        if standardized
                        else f"{self._title_name()}: fit residuals"
                    ),
                    x_name=axis_name,
                    x_unit=diagnostic.units.get(x_column),
                    y_name=y_name,
                    y_unit=diagnostic.units.get(y_name),
                    y_label=self._residual_axis_label(
                        standardized=standardized,
                        unit=diagnostic.units.get(y_name),
                    ),
                    series=series,
                    options=options,
                )
            )
        return specs

    def _build_normalized_pressure_spec(
        self,
        diagnostic: EOSDiagnosticResult,
        options: EOSPlotOptions,
    ) -> LinePlotSpec:
        metadata = diagnostic.metadata.get("normalized_pressure", {})
        if not metadata.get("available", False):
            raise ValueError("normalized pressure is unavailable for this EOS family")
        valid = (
            diagnostic.columns["normalized_pressure_valid"].astype(bool)
            & np.isfinite(diagnostic.columns["finite_strain"])
            & np.isfinite(diagnostic.columns["normalized_pressure"])
        )
        zero_pressure = np.zeros(valid.shape, dtype=np.bool_)
        if not options.show_zero_pressure_point:
            pressure = diagnostic.columns["observed_pressure"]
            sigma_pressure = self._input_error("pressure")
            scale = max(1.0, float(np.nanmax(np.abs(pressure))))
            numerical_tolerance = float(np.sqrt(np.finfo(np.float64).eps) * scale)
            if sigma_pressure is None:
                zero_pressure = np.abs(pressure) <= numerical_tolerance
            else:
                zero_pressure = np.abs(pressure) <= np.maximum(
                    np.asarray(sigma_pressure, dtype=np.float64),
                    numerical_tolerance,
                )
            valid &= ~zero_pressure
        if not np.any(valid):
            raise ValueError("no finite normalized-pressure observations are available")
        x = diagnostic.columns["finite_strain"]
        y = diagnostic.columns["normalized_pressure"]
        x_error = diagnostic.columns.get("sigma_finite_strain")
        y_error = diagnostic.columns.get("sigma_normalized_pressure")
        series = self._observation_series(
            x,
            y,
            diagnostic,
            options,
            x_error=x_error,
            y_error=y_error,
            extra_mask=valid,
        )
        included_valid = valid & diagnostic.columns["included"].astype(bool)
        order = np.argsort(x[included_valid])
        series.insert(
            0,
            PlotSeries(
                key="normalized_fit",
                label=f"Fit ({self._model_tag()})",
                x=x[included_valid][order],
                y=diagnostic.columns["calculated_normalized_pressure"][included_valid][
                    order
                ],
                style=PlotSeriesStyle(
                    color=options.curve_color,
                    line_width=options.curve_width,
                ),
            ),
        )
        x_label, y_label = normalized_pressure_labels(
            metadata,
            diagnostic.units.get("normalized_pressure"),
        )
        return LinePlotSpec(
            key="normalized_pressure",
            title=f"{self._title_name()}: normalized pressure",
            filename_stem=self._filename_stem("normalized_pressure"),
            x_axis=PlotAxis(key="finite_strain", label=x_label, unit="1"),
            y_axis=PlotAxis(
                key="normalized_pressure",
                label=y_label,
                unit=diagnostic.units.get("normalized_pressure"),
            ),
            series=series,
            show_legend=options.show_legend,
            grid=options.grid,
            metadata={
                **self._plot_metadata("normalized-pressure"),
                "zero_pressure_points_shown": options.show_zero_pressure_point,
                "zero_pressure_points_hidden": int(np.count_nonzero(zero_pressure)),
            },
        )

    def _build_coverage_spec(
        self,
        diagnostic: EOSDiagnosticResult,
        options: EOSPlotOptions,
    ) -> LinePlotSpec:
        if self.record.request.domain is not EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
            raise ValueError("P-T coverage is available only for P-V-T records")
        series = self._observation_series(
            diagnostic.columns["observed_pressure"],
            diagnostic.columns["temperature"],
            diagnostic,
            options,
            x_error=self._input_error("pressure"),
            y_error=self._input_error("temperature"),
        )
        return self._line_spec(
            key="coverage",
            title=f"{self._title_name()}: pressure-temperature coverage",
            x_name="pressure",
            x_unit=diagnostic.units.get("observed_pressure"),
            y_name="temperature",
            y_unit=diagnostic.units.get("temperature"),
            series=series,
            options=options,
        )

    def _build_isotherm_spec(
        self,
        diagnostic: EOSDiagnosticResult,
        options: EOSPlotOptions,
    ) -> LinePlotSpec:
        if self.record.request.domain is not EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
            raise ValueError("isotherms are available only for P-V-T records")
        volume = diagnostic.columns["volume"]
        selected = diagnostic.columns["included"].astype(bool)
        curve_volume = self._curve_grid(
            volume,
            selected,
            points=options.curve_points,
        )
        temperatures = options.isotherms or self._representative_values(
            diagnostic.columns["temperature"][selected], options.max_default_curves
        )
        series: list[PlotSeries] = []
        warnings: list[str] = []
        for index, temperature in enumerate(temperatures):
            try:
                result = self._calculator.calculate(
                    volume=curve_volume,
                    temperature=np.full(curve_volume.size, temperature),
                    propagate_uncertainty=False,
                )
            except (ValueError, RuntimeError) as exc:
                warnings.append(f"T={temperature:g} K: {exc}")
                continue
            series.append(
                PlotSeries(
                    key=f"isotherm_{index}",
                    label=f"{temperature:g} K",
                    x=curve_volume.copy(),
                    y=result.columns["pressure"].copy(),
                    style=PlotSeriesStyle(line_width=options.curve_width),
                    metadata={"temperature": float(temperature)},
                )
            )
        series.extend(
            self._observation_series(
                volume,
                diagnostic.columns["observed_pressure"],
                diagnostic,
                options,
                x_error=self._input_error("volume"),
                y_error=self._input_error("pressure"),
                generic_label="Observations",
            )
        )
        spec = self._line_spec(
            key="isotherms",
            title=f"{self._title_name()}: P-V-T isotherms",
            x_name="volume",
            x_unit=diagnostic.units.get("volume"),
            y_name="pressure",
            y_unit=diagnostic.units.get("observed_pressure"),
            series=series,
            options=options,
        )
        spec.metadata["warnings"] = warnings
        return spec

    def _build_isobar_spec(
        self,
        diagnostic: EOSDiagnosticResult,
        options: EOSPlotOptions,
    ) -> LinePlotSpec:
        if self.record.request.domain is not EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE:
            raise ValueError("isobars are available only for P-V-T records")
        temperature = diagnostic.columns["temperature"]
        selected = diagnostic.columns["included"].astype(bool)
        curve_temperature = self._curve_grid(
            temperature,
            selected,
            positive=False,
            points=options.curve_points,
        )
        pressures = options.isobars or self._representative_values(
            diagnostic.columns["observed_pressure"][selected],
            options.max_default_curves,
        )
        series: list[PlotSeries] = []
        warnings: list[str] = []
        for index, pressure in enumerate(pressures):
            try:
                result = self._calculator.calculate(
                    pressure=np.full(curve_temperature.size, pressure),
                    temperature=curve_temperature,
                    propagate_uncertainty=False,
                )
            except (ValueError, RuntimeError) as exc:
                warnings.append(f"P={pressure:g} GPa: {exc}")
                continue
            series.append(
                PlotSeries(
                    key=f"isobar_{index}",
                    label=f"{pressure:g} GPa",
                    x=curve_temperature.copy(),
                    y=result.columns["volume"].copy(),
                    style=PlotSeriesStyle(line_width=options.curve_width),
                    metadata={"pressure": float(pressure)},
                )
            )
        series.extend(
            self._observation_series(
                temperature,
                diagnostic.columns["volume"],
                diagnostic,
                options,
                x_error=self._input_error("temperature"),
                y_error=self._input_error("volume"),
                generic_label="Observations",
            )
        )
        spec = self._line_spec(
            key="isobars",
            title=f"{self._title_name()}: P-V-T isobars",
            x_name="temperature",
            x_unit=diagnostic.units.get("temperature"),
            y_name="volume",
            y_unit=diagnostic.units.get("volume"),
            series=series,
            options=options,
        )
        spec.metadata["warnings"] = warnings
        return spec

    def _observation_series(
        self,
        x: np.ndarray,
        y: np.ndarray,
        diagnostic: EOSDiagnosticResult,
        options: EOSPlotOptions,
        *,
        x_error: np.ndarray | None,
        y_error: np.ndarray | None,
        extra_mask: np.ndarray | None = None,
        generic_label: str = "Data",
    ) -> list[PlotSeries]:
        included = diagnostic.columns["included"].astype(bool)
        finite = np.isfinite(x) & np.isfinite(y)
        if extra_mask is not None:
            finite &= np.asarray(extra_mask, dtype=bool)
        groups = diagnostic.columns["group"].astype(int)
        series: list[PlotSeries] = []
        if options.group_data and np.any(groups[included & finite] > 0):
            identifiers: Iterable[int] = sorted(set(groups[included & finite].tolist()))
        else:
            identifiers = (0,)
        for identifier in identifiers:
            mask = included & finite
            if options.group_data and np.any(groups[included & finite] > 0):
                mask &= groups == identifier
                label = f"Group {identifier}" if identifier > 0 else "Ungrouped"
            else:
                label = generic_label
            if not np.any(mask):
                continue
            series.append(
                self._data_series(
                    f"data_group_{identifier}",
                    label,
                    x,
                    y,
                    mask,
                    options,
                    x_error=x_error,
                    y_error=y_error,
                )
            )
        if options.show_excluded:
            mask = (~included) & finite
            if np.any(mask):
                series.append(
                    self._data_series(
                        "excluded",
                        "Excluded",
                        x,
                        y,
                        mask,
                        options,
                        x_error=x_error,
                        y_error=y_error,
                        excluded=True,
                    )
                )
        return series

    def _data_series(
        self,
        key: str,
        label: str,
        x: np.ndarray,
        y: np.ndarray,
        mask: np.ndarray,
        options: EOSPlotOptions,
        *,
        x_error: np.ndarray | None,
        y_error: np.ndarray | None,
        excluded: bool = False,
    ) -> PlotSeries:
        style = PlotSeriesStyle(
            color=options.excluded_color if excluded else None,
            line_style="none",
            line_width=0.0,
            marker="x" if excluded else "o",
            marker_size=options.point_size,
            alpha=0.7 if excluded else 1.0,
            errorbar_line_width=options.errorbar_width,
            errorbar_capsize=options.errorbar_capsize,
        )
        return PlotSeries(
            key=key,
            label=label,
            x=np.asarray(x[mask], dtype=np.float64),
            y=np.asarray(y[mask], dtype=np.float64),
            x_error=(
                None
                if not options.show_uncertainties or x_error is None
                else np.asarray(x_error[mask], dtype=np.float64)
            ),
            y_error=(
                None
                if not options.show_uncertainties or y_error is None
                else np.asarray(y_error[mask], dtype=np.float64)
            ),
            style=style,
            metadata={"excluded": excluded},
        )

    def _fit_curve(
        self,
        x: np.ndarray,
        y: np.ndarray,
        options: EOSPlotOptions,
    ) -> PlotSeries:
        return PlotSeries(
            key="fit",
            label=f"Fit ({self._model_tag()})",
            x=np.asarray(x, dtype=np.float64),
            y=np.asarray(y, dtype=np.float64),
            style=PlotSeriesStyle(
                color=options.curve_color,
                line_width=options.curve_width,
            ),
        )

    def _line_spec(
        self,
        *,
        key: str,
        title: str,
        x_name: str,
        x_unit: str | None,
        y_name: str,
        y_unit: str | None,
        series: list[PlotSeries],
        x_label: str | None = None,
        y_label: str | None = None,
        options: EOSPlotOptions,
    ) -> LinePlotSpec:
        return LinePlotSpec(
            key=key,
            title=title,
            filename_stem=self._filename_stem(key),
            x_axis=PlotAxis(
                key=x_name,
                label=x_label or property_label(x_name, x_unit),
                unit=x_unit,
            ),
            y_axis=PlotAxis(
                key=y_name,
                label=y_label or property_label(y_name, y_unit),
                unit=y_unit,
            ),
            series=series,
            show_legend=options.show_legend and len(series) > 1,
            grid=options.grid,
            metadata=self._plot_metadata(key),
        )

    def _title_name(self) -> str:
        """Return a human-readable dataset title for figure headings."""
        return " ".join(self.dataset.jobname.replace("_", " ").split())

    def _residual_axis_label(
        self,
        *,
        standardized: bool,
        unit: str | None,
    ) -> str:
        if standardized:
            return property_label("standardized_residual", None)
        request = self.record.request
        if request.domain in {
            EOSFitDomain.PRESSURE_VOLUME,
            EOSFitDomain.PRESSURE_VOLUME_TEMPERATURE,
        }:
            base = "Pressure residual"
        elif request.target == "volume":
            base = "Volume residual"
        else:
            base = rf"Cell parameter ${request.target}$ residual"
        rendered_unit = format_unit(unit)
        return base if rendered_unit is None else f"{base} ({rendered_unit})"

    def _input_error(self, name: str) -> np.ndarray | None:
        sigma_name = f"sigma_{name}"
        if not self.dataset.has(sigma_name):
            return None
        values = self.dataset.column(sigma_name)
        return np.asarray(values, dtype=np.float64)

    def _curve_grid(
        self,
        values: np.ndarray,
        included: np.ndarray,
        *,
        positive: bool = True,
        points: int = 300,
    ) -> np.ndarray:
        mask = np.asarray(included, dtype=bool) & np.isfinite(values)
        if not np.any(mask):
            raise ValueError("no selected observations are available for a fit curve")
        lower = float(np.min(values[mask]))
        upper = float(np.max(values[mask]))
        span = upper - lower
        margin = 0.02 * span if span > 0.0 else 0.02 * max(abs(lower), 1.0)
        lower -= margin
        upper += margin
        if positive:
            lower = max(lower, float(np.finfo(np.float64).tiny))
        return np.linspace(
            lower,
            upper,
            int(points),
            dtype=np.float64,
        )

    @staticmethod
    def _representative_values(values: np.ndarray, count: int) -> tuple[float, ...]:
        array = np.asarray(values, dtype=np.float64)
        array = array[np.isfinite(array)]
        if array.size == 0:
            raise ValueError("no finite values are available for curve selection")
        unique = np.unique(array)
        if unique.size <= count:
            return tuple(float(item) for item in unique)
        return tuple(float(item) for item in np.linspace(unique[0], unique[-1], count))

    @staticmethod
    def _resolve_plot_types(
        plot_types: Sequence[str] | str | None,
        available: Sequence[str],
    ) -> tuple[str, ...]:
        if plot_types is None:
            return tuple(available)
        if isinstance(plot_types, str):
            requested = [item.strip() for item in plot_types.split(",") if item.strip()]
        else:
            requested = [str(item).strip() for item in plot_types if str(item).strip()]
        if not requested or any(item.lower() == "all" for item in requested):
            return tuple(available)
        normalized = tuple(item.lower().replace("_", "-") for item in requested)
        unsupported = [item for item in normalized if item not in available]
        if unsupported:
            raise ValueError(
                "unsupported EOS plot type(s): "
                + ", ".join(unsupported)
                + "; available: "
                + ", ".join(available)
            )
        return normalized

    def _model_tag(self) -> str:
        return str(getattr(self.record.request.model, "tag", self.record.request.model))

    def _filename_stem(self, key: str) -> str:
        return f"record_{self.record.record_id:06d}_{self.record.slot.key}_{key}"

    def _plot_metadata(self, kind: str) -> dict[str, object]:
        return {
            "module": "eos",
            "kind": kind,
            "record_id": self.record.record_id,
            "slot": self.record.slot.key,
            "model": self._model_tag(),
            "dataset_id": self.record.dataset_id,
        }


__all__ = ["EOS_PLOT_TYPES", "EOSPlotOptions", "EOSPlotter"]
