# -*- coding: utf-8 -*-

"""Path helpers shared by Quantas readers, exporters, and frontends."""

from __future__ import annotations

from pathlib import Path


def ensure_suffix(filename: str | Path, suffix: str) -> Path:
    """
    Return a path with the required suffix.

    If the filename already has the requested suffix, it is returned
    unchanged. Otherwise the suffix is appended.

    Parameters
    ----------
    filename
        Input filename.

    suffix
        Required suffix, including the leading dot.

    Returns
    -------
    Path
        Output filename with the requested suffix.
    """
    path = Path(filename)

    if path.suffix.lower() != suffix.lower():
        path = path.with_suffix(suffix)

    return path
