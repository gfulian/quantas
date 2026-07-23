# -*- coding: utf-8 -*-

"""Standard figure presets shared by all static Quantas plots."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, cast

FigurePreset = Literal["screen", "publication", "monochrome"]
FIGURE_PRESET_NAMES: tuple[FigurePreset, ...] = (
    "screen",
    "publication",
    "monochrome",
)


@dataclass(frozen=True, slots=True)
class FigurePresetSpec:
    """Resolved renderer defaults for one named figure preset.

    Parameters
    ----------
    name : {"screen", "publication", "monochrome"}
        Stable preset identifier exposed by the CLI.
    dpi : int
        Default raster resolution.
    figure_size : tuple of float or None
        Default figure size in inches. ``None`` preserves plot-specific sizes.
    tight_layout : bool
        Whether Matplotlib tight layout is requested after rendering.
    axis_label_font_size, legend_font_size, title_font_size,
    tick_label_font_size : float or None
        Renderer-level typography defaults.
    monochrome : bool
        Whether rendered colors are converted to a grayscale-safe style.
    savefig_kwargs : dict
        Additional deterministic save options.
    """

    name: FigurePreset
    dpi: int
    figure_size: tuple[float, float] | None
    tight_layout: bool
    axis_label_font_size: float | None
    legend_font_size: float | None
    title_font_size: float | None
    tick_label_font_size: float | None
    monochrome: bool = False
    savefig_kwargs: dict[str, Any] = field(default_factory=dict)


_PRESETS: dict[FigurePreset, FigurePresetSpec] = {
    "screen": FigurePresetSpec(
        name="screen",
        dpi=150,
        figure_size=None,
        tight_layout=False,
        axis_label_font_size=None,
        legend_font_size=None,
        title_font_size=None,
        tick_label_font_size=None,
        savefig_kwargs={},
    ),
    "publication": FigurePresetSpec(
        name="publication",
        dpi=300,
        figure_size=None,
        tight_layout=True,
        axis_label_font_size=13.0,
        legend_font_size=11.0,
        title_font_size=13.0,
        tick_label_font_size=11.0,
        savefig_kwargs={"bbox_inches": "tight"},
    ),
    "monochrome": FigurePresetSpec(
        name="monochrome",
        dpi=300,
        figure_size=None,
        tight_layout=True,
        axis_label_font_size=13.0,
        legend_font_size=11.0,
        title_font_size=13.0,
        tick_label_font_size=11.0,
        monochrome=True,
        savefig_kwargs={"bbox_inches": "tight"},
    ),
}


def figure_preset(name: str) -> FigurePresetSpec:
    """Return one validated standard figure preset.

    Parameters
    ----------
    name : str
        Case-insensitive preset identifier.

    Returns
    -------
    FigurePresetSpec
        Immutable resolved defaults.

    Raises
    ------
    ValueError
        If ``name`` is not a supported preset.
    """

    normalized = name.strip().lower()
    if normalized not in _PRESETS:
        choices = ", ".join(FIGURE_PRESET_NAMES)
        raise ValueError(f"unknown figure preset {name!r}; choose from {choices}")
    return _PRESETS[cast(FigurePreset, normalized)]


__all__ = [
    "FIGURE_PRESET_NAMES",
    "FigurePreset",
    "FigurePresetSpec",
    "figure_preset",
]
