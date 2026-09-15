# -*- coding: utf-8 -*-

"""Shared normalization helpers for EOS passive contracts."""

from __future__ import annotations

import numpy as np

def as_float64_vector(values: np.ndarray, *, name: str) -> np.ndarray:
    """Return a finite one-dimensional ``float64`` copy.

    Parameters
    ----------
    values : np.ndarray
        Numerical values consumed by the operation.
    name : str
        Stable name or identifier.

    Returns
    -------
    np.ndarray
        A finite one-dimensional ``float64`` copy.

    Raises
    ------
    ValueError
        If the supplied data or workflow state violates the documented contract.
    """
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1:
        raise ValueError(f"EOS column '{name}' must be one-dimensional.")
    if array.size == 0:
        raise ValueError(f"EOS column '{name}' cannot be empty.")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"EOS column '{name}' must contain finite values.")
    return array.copy()


def optional_sigma(
    values: np.ndarray | None,
    shape: tuple[int, ...],
    name: str,
) -> np.ndarray | None:
    """Normalize an optional standard-uncertainty vector.

    Parameters
    ----------
    values : np.ndarray | None
        Numerical values consumed by the operation.
    shape : tuple[int, ...]
        Expected array shape.
    name : str
        Stable name or identifier.

    Returns
    -------
    np.ndarray | None
        Result described by the operation.

    Raises
    ------
    ValueError
        If the supplied data or workflow state violates the documented contract.
    """
    if values is None:
        return None
    array = as_float64_vector(values, name=name)
    if array.shape != shape:
        raise ValueError(f"EOS series {name} must match the data shape.")
    if np.any(array < 0.0):
        raise ValueError(f"EOS series {name} cannot contain negative values.")
    return array

__all__ = ["as_float64_vector", "optional_sigma"]
