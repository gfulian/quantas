# -*- coding: utf-8 -*-

"""Neutral energy-volume plots for QHA input inspection."""

from __future__ import annotations

import numpy as np

from quantas.core.physics.eos import EnergyEOS
from quantas.models import (
    LinePlotSpec,
    PlotAxis,
    PlotCollection,
    PlotSeries,
    PlotSeriesStyle,
)
from quantas.modules.qha.formatting import (
    native_energy_unit_label,
    volume_unit_label,
)
from quantas.modules.qha.inspect import PressureEstimate, PressureVolumePreview


def _sampled_fit_series(
    estimate: PressureEstimate,
    volume: np.ndarray,
) -> PlotSeries | None:
    """Return a sampled energy-fit series when an estimate is usable."""
    if not estimate.success or estimate.fit.parameters is None:
        return None

    parameters = np.asarray(estimate.fit.parameters, dtype=np.float64)
    if estimate.method == "polynomial":
        energy = np.polynomial.polynomial.polyval(volume, parameters)
        key = "polynomial_fit"
        label = f"Polynomial fit (degree {estimate.metadata.get('degree', '?')})"
        style = PlotSeriesStyle(color="#d97706", line_width=2.0)
    elif estimate.method == "eos":
        eos = str(estimate.metadata.get("eos", ""))
        if not eos:
            return None
        energy = EnergyEOS().evaluate(eos, volume, parameters)
        key = "eos_fit"
        label = f"{eos} fit"
        style = PlotSeriesStyle(color="#2563eb", line_width=2.0)
    else:
        return None

    return PlotSeries(
        key=key,
        label=label,
        x=volume.copy(),
        y=np.asarray(energy, dtype=np.float64),
        style=style,
        metadata={
            "method": estimate.method,
            "fit_status": estimate.fit.status.value,
            "fit_quality": estimate.fit.quality.value,
        },
    )


def build_pressure_volume_preview_plots(
    preview: PressureVolumePreview,
    *,
    sample_points: int = 201,
) -> PlotCollection:
    """Build an energy-volume plot from a QHA inspection preview.

    The observed input values are always included. Successful polynomial and
    EOS fits are sampled only between the minimum and maximum input volumes.
    No extrapolation is performed.
    """
    if sample_points < 2:
        raise ValueError("sample_points must be at least 2")
    if sample_points > 10_001:
        raise ValueError("sample_points cannot exceed 10001")
    if preview.volume.size == 0 or preview.energy.size == 0:
        raise ValueError("inspection preview contains no energy-volume data")

    volume = np.linspace(
        float(np.min(preview.volume)),
        float(np.max(preview.volume)),
        int(sample_points),
        dtype=np.float64,
    )
    series = [
        PlotSeries(
            key="observed",
            label="Input data",
            x=preview.volume.copy(),
            y=preview.energy.copy(),
            style=PlotSeriesStyle(
                color="#111827",
                line_style="none",
                marker="circle",
                marker_size=6.0,
            ),
            metadata={"source": "input"},
        )
    ]
    warnings_ = list(preview.warnings)
    for estimate in (preview.polynomial, preview.eos):
        if estimate is None:
            continue
        try:
            fitted = _sampled_fit_series(estimate, volume)
        except (TypeError, ValueError, FloatingPointError) as exc:
            warnings_.append(
                f"could not sample the {estimate.method} energy fit: {exc}"
            )
            continue
        if fitted is not None:
            series.append(fitted)

    volume_unit = volume_unit_label(str(preview.metadata.get("volume_unit", "A")))
    energy_unit = native_energy_unit_label(
        str(preview.metadata.get("energy_unit", "Ha"))
    )
    spec = LinePlotSpec(
        key="qha_input_energy_volume",
        title="Static energy-volume inspection",
        filename_stem="qha_input_energy_volume",
        x_axis=PlotAxis(
            key="volume",
            label=f"V ({volume_unit})",
            unit=volume_unit,
        ),
        y_axis=PlotAxis(
            key="energy",
            label=f"E ({energy_unit})",
            unit=energy_unit,
        ),
        series=series,
        show_legend=len(series) > 1,
        metadata={
            "source": "qha_inspection",
            "sample_points": int(sample_points),
            "volume_interval": (
                float(np.min(preview.volume)),
                float(np.max(preview.volume)),
            ),
            "extrapolated": False,
        },
    )
    return PlotCollection(plots=[spec], warnings=warnings_)


__all__ = ["build_pressure_volume_preview_plots"]
