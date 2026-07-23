# -*- coding: utf-8 -*-

"""File-output helpers for the Matplotlib plotting backend."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from quantas.models.plot import PlotSpec, SurfacePlotSpec

from .options import MatplotlibOptions


def save_figure(
    figure: Any,
    spec: PlotSpec,
    options: MatplotlibOptions,
) -> Path | None:
    """Save a rendered figure when file output is configured.

    Parameters
    ----------
    figure : matplotlib.figure.Figure
        Rendered figure.
    spec : PlotSpec
        Source neutral specification.
    options : MatplotlibOptions
        File-output options.

    Returns
    -------
    Path or None
        Written file path, or ``None`` when no output directory is configured.
    """
    if options.output_dir is None:
        return None
    output_dir = Path(options.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = options.image_format.lower().lstrip(".")
    stem = safe_stem(f"{options.filename_prefix}{spec.filename_stem}")
    path = output_dir / f"{stem}.{suffix}"
    savefig_kwargs = dict(options.savefig_kwargs)
    if isinstance(spec, SurfacePlotSpec):
        savefig_kwargs.pop("bbox_inches", None)
    figure.savefig(path, dpi=options.dpi, **savefig_kwargs)
    return path


def safe_stem(stem: str) -> str:
    """Return a filename-safe plot stem.

    Parameters
    ----------
    stem : str
        Raw filename stem.

    Returns
    -------
    str
        Sanitized filename stem.
    """
    return stem.replace("/", "_").replace(" ", "_").replace("$", "").replace("-", "_")


__all__ = ["safe_stem", "save_figure"]
