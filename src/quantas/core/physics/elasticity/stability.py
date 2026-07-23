# -*- coding: utf-8 -*-

"""Mechanical-stability checks for elastic stiffness matrices.

The generic criterion implemented here is positive definiteness of the
symmetric stress--strain stiffness representation.  For unstressed crystals
this is the Born elastic-stability criterion.  For a hydrostatically stressed
crystal it applies to the corresponding Wallace/Barron--Klein stress--strain
coefficients, not to an uncorrected energy Hessian.

References
----------
Canonical citation keys: ``mouhat_coudert_2014``, ``barron_klein_1965``,
and ``wallace_1972``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from quantas.references import method_citation_keys, render_citation_inline

from .tensor import ElasticTensor


@dataclass(frozen=True, slots=True)
class StabilityResult:
    """Result of a positive-definiteness stability check.

    Parameters
    ----------
    is_stable : bool
        Whether all stiffness eigenvalues exceed the selected tolerance.
    eigenvalues : ndarray
        Stiffness-matrix eigenvalues in GPa.
    tolerance : float
        Threshold used for the check.
    """

    is_stable: bool
    eigenvalues: np.ndarray
    tolerance: float = 0.0

    @property
    def failed_eigenvalues(self) -> np.ndarray:
        """Return eigenvalues that do not satisfy the stability criterion."""
        return self.eigenvalues[self.eigenvalues <= self.tolerance]


@dataclass(slots=True)
class StabilityFieldResult:
    """Mechanical-stability diagnostics for a grid of stiffness matrices.

    Parameters
    ----------
    eigenvalues : ndarray
        Ordered stiffness eigenvalues with shape ``field_shape + (6,)`` in GPa.
    minimum_eigenvalue : ndarray
        Minimum eigenvalue at every state in GPa.
    stable_mask : ndarray
        ``True`` where all eigenvalues exceed ``tolerance``.
    indeterminate_mask : ndarray
        ``True`` where the input matrix was non-finite or could not be assessed.
    tolerance : float
        Positive-definiteness threshold in GPa.
    criterion : str, optional
        Stable identifier for the applied criterion.
    metadata : dict, optional
        Tensor convention and reference provenance.
    """

    eigenvalues: NDArray[np.float64]
    minimum_eigenvalue: NDArray[np.float64]
    stable_mask: NDArray[np.bool_]
    indeterminate_mask: NDArray[np.bool_]
    tolerance: float = 0.0
    criterion: str = "positive_definite_wallace_stiffness"
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize arrays and validate aligned field shapes."""
        self.eigenvalues = np.asarray(self.eigenvalues, dtype=np.float64).copy()
        if self.eigenvalues.ndim < 1 or self.eigenvalues.shape[-1] != 6:
            raise ValueError("stability eigenvalues must have shape field + (6,)")
        shape = self.eigenvalues.shape[:-1]
        self.minimum_eigenvalue = np.asarray(
            self.minimum_eigenvalue, dtype=np.float64
        ).copy()
        self.stable_mask = np.asarray(self.stable_mask, dtype=np.bool_).copy()
        self.indeterminate_mask = np.asarray(
            self.indeterminate_mask, dtype=np.bool_
        ).copy()
        if (
            self.minimum_eigenvalue.shape != shape
            or self.stable_mask.shape != shape
            or self.indeterminate_mask.shape != shape
        ):
            raise ValueError("stability field arrays must share one field shape")
        if not np.isfinite(self.tolerance):
            raise ValueError("stability tolerance must be finite")
        self.tolerance = float(self.tolerance)
        self.criterion = str(self.criterion)
        self.metadata = dict(self.metadata)

    @property
    def unstable_mask(self) -> NDArray[np.bool_]:
        """Return assessed states that fail the stability criterion."""
        return np.asarray(~self.stable_mask & ~self.indeterminate_mask, dtype=np.bool_)

    @property
    def minimum_margin(self) -> NDArray[np.float64]:
        """Return ``lambda_min - tolerance`` in GPa."""
        return np.asarray(self.minimum_eigenvalue - self.tolerance, dtype=np.float64)


def check_positive_definiteness(
    tensor_or_matrix: ElasticTensor | Sequence[Sequence[float]],
    tolerance: float = 0.0,
) -> StabilityResult:
    """Check whether an elastic stiffness matrix is positive definite.

    Parameters
    ----------
    tensor_or_matrix : ElasticTensor or array_like
        Elastic tensor or ``6 x 6`` stiffness matrix.
    tolerance : float, optional
        Minimum accepted eigenvalue in GPa.

    Returns
    -------
    StabilityResult
        Stability flag and eigenvalues.

    Raises
    ------
    ValueError
        If an array input does not have shape ``(6, 6)``.
    """
    if isinstance(tensor_or_matrix, ElasticTensor):
        eigenvalues = tensor_or_matrix.eigenvalues
    else:
        matrix = np.asarray(tensor_or_matrix, dtype=float)
        if matrix.shape != (6, 6):
            raise ValueError("The elastic stiffness matrix must have shape (6, 6).")
        eigenvalues = np.linalg.eigvalsh(matrix)

    return StabilityResult(
        is_stable=bool(np.all(eigenvalues > tolerance)),
        eigenvalues=np.asarray(eigenvalues, dtype=float),
        tolerance=float(tolerance),
    )


def evaluate_stability_field(
    stiffness: ArrayLike,
    *,
    tolerance: float = 0.0,
    tensor_convention: str = "Wallace hydrostatic stress-strain stiffness",
) -> StabilityFieldResult:
    """Evaluate positive-definiteness over one or more stiffness matrices.

    Parameters
    ----------
    stiffness : array_like
        Matrix or field with shape ``(..., 6, 6)`` in GPa.
    tolerance : float, optional
        Minimum accepted eigenvalue in GPa.
    tensor_convention : str, optional
        Human-readable identification of the stiffness representation.

    Returns
    -------
    StabilityFieldResult
        Eigenvalues, minimum margins, and stable/indeterminate masks.

    Raises
    ------
    ValueError
        If the trailing matrix shape is not ``(6, 6)`` or the tolerance is
        non-finite.

    Notes
    -----
    Non-finite matrices are retained as indeterminate states.  They are never
    silently classified as stable or replaced by nearby values.
    """
    matrices = np.asarray(stiffness, dtype=np.float64)
    if matrices.shape[-2:] != (6, 6):
        raise ValueError("stiffness must have shape (..., 6, 6)")
    if not np.isfinite(tolerance):
        raise ValueError("stability tolerance must be finite")
    shape = matrices.shape[:-2]
    flattened = matrices.reshape((-1, 6, 6))
    nmatrices = int(flattened.size // 36)
    eigenvalues = np.full((nmatrices, 6), np.nan, dtype=np.float64)
    indeterminate = np.zeros(nmatrices, dtype=np.bool_)
    for index, matrix in enumerate(flattened):
        if np.any(~np.isfinite(matrix)):
            indeterminate[index] = True
            continue
        symmetric = 0.5 * (matrix + matrix.T)
        try:
            eigenvalues[index] = np.linalg.eigvalsh(symmetric)
        except np.linalg.LinAlgError:
            indeterminate[index] = True
    minimum = np.min(eigenvalues, axis=1, initial=np.nan)
    # np.min with initial does not ignore NaN; explicitly preserve indeterminate.
    for index in range(minimum.size):
        if indeterminate[index]:
            minimum[index] = np.nan
        else:
            minimum[index] = float(np.min(eigenvalues[index]))
    stable = np.asarray(
        (~indeterminate) & np.all(eigenvalues > float(tolerance), axis=1),
        dtype=np.bool_,
    )
    return StabilityFieldResult(
        eigenvalues=eigenvalues.reshape(shape + (6,)),
        minimum_eigenvalue=minimum.reshape(shape),
        stable_mask=stable.reshape(shape),
        indeterminate_mask=indeterminate.reshape(shape),
        tolerance=float(tolerance),
        metadata={
            "tensor_convention": str(tensor_convention),
            "criterion_reference": render_citation_inline("mouhat_coudert_2014"),
            "finite_pressure_reference": "; ".join(
                render_citation_inline(key)
                for key in ("barron_klein_1965", "wallace_1972")
            ),
            "citation_keys": list(method_citation_keys("mechanical_stability")),
        },
    )


__all__ = [
    "StabilityFieldResult",
    "StabilityResult",
    "check_positive_definiteness",
    "evaluate_stability_field",
]
