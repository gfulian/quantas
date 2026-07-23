# -*- coding: utf-8 -*-

"""Fréchet derivatives of matrix functions used by Quantas workflows.

The primitives validate shapes, preserve ``float64`` numerical contracts, and
delegate the matrix-function evaluation to SciPy.
"""

from __future__ import annotations

from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import expm_frechet

FloatArray: TypeAlias = NDArray[np.float64]


def matrix_exponential_frechet(
    matrix: ArrayLike,
    direction: ArrayLike,
) -> FloatArray:
    r"""Return the Fréchet derivative of the matrix exponential.

    The returned matrix is

    .. math::

        L_{\exp}(A, E)
        = \left.\frac{d}{dt}\exp(A+tE)\right|_{t=0},

    where ``matrix`` is :math:`A` and ``direction`` is :math:`E`.

    Parameters
    ----------
    matrix : array_like
        Finite square matrix.
    direction : array_like
        Finite perturbation matrix with the same shape as ``matrix``.

    Returns
    -------
    ndarray
        Fréchet derivative in ``float64`` with the same shape as the inputs.

    Raises
    ------
    ValueError
        If the inputs are not finite square matrices of identical shape.
    """
    base = np.asarray(matrix, dtype=np.float64)
    perturbation = np.asarray(direction, dtype=np.float64)
    if base.ndim != 2 or base.shape[0] != base.shape[1]:
        raise ValueError("matrix must be a finite square matrix")
    if perturbation.shape != base.shape:
        raise ValueError("direction must have the same shape as matrix")
    if not np.all(np.isfinite(base)):
        raise ValueError("matrix must contain only finite values")
    if not np.all(np.isfinite(perturbation)):
        raise ValueError("direction must contain only finite values")
    return np.asarray(
        expm_frechet(base, perturbation, compute_expm=False),
        dtype=np.float64,
    )
