# -*- coding: utf-8 -*-

"""Passive helpers for solver failure analysis and optional debug traces.

The numerical fitting services use this module to collect backend-neutral
summaries without depending on terminal or GUI objects.  Detailed model-
evaluation traces are recorded only when ``solver_debug`` is true in the
solver metadata; compact input, uncertainty, bound, and termination summaries
are always retained.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np


SOLVER_DEBUG_METADATA_KEY = "solver_debug"
_DEFAULT_TRACE_LENGTH = 100


def solver_debug_enabled(metadata: Mapping[str, Any] | None) -> bool:
    """Return whether detailed numerical tracing was requested.

    Parameters
    ----------
    metadata : mapping or None
        Solver metadata supplied through :class:`~quantas.core.math.fitting.FitOptions`.

    Returns
    -------
    bool
        ``True`` when detailed traces should be retained.
    """
    return bool((metadata or {}).get(SOLVER_DEBUG_METADATA_KEY, False))


def array_summary(values: Any, *, positive_ratio: bool = False) -> dict[str, Any]:
    """Return a serialization-ready numerical summary.

    Parameters
    ----------
    values : array-like
        Numerical values to summarize.
    positive_ratio : bool, optional
        Include ``max/min`` for strictly positive finite entries.

    Returns
    -------
    dict
        Shape, size, finite count, extrema, median, Euclidean norm, and an
        optional positive dynamic-range ratio.
    """
    array = np.asarray(values, dtype=np.float64)
    flat = array.ravel()
    finite = flat[np.isfinite(flat)]
    summary: dict[str, Any] = {
        "shape": list(array.shape),
        "size": int(flat.size),
        "finite_count": int(finite.size),
    }
    if finite.size == 0:
        return summary
    summary.update(
        {
            "minimum": float(np.min(finite)),
            "maximum": float(np.max(finite)),
            "median": float(np.median(finite)),
            "norm": float(np.linalg.norm(finite)),
        }
    )
    if positive_ratio:
        positive = finite[finite > 0.0]
        if positive.size:
            minimum = float(np.min(positive))
            maximum = float(np.max(positive))
            summary["positive_minimum"] = minimum
            summary["positive_maximum"] = maximum
            summary["positive_dynamic_range"] = (
                float(maximum / minimum) if minimum > 0.0 else float("inf")
            )
    return summary


def problem_summary(
    x: Any,
    y: Any,
    initial: Sequence[float] | np.ndarray,
    bounds: tuple[Any, Any],
    sigma: Any | None,
) -> dict[str, Any]:
    """Return common diagnostics for one nonlinear fitting problem.

    Parameters
    ----------
    x, y : array-like
        Explanatory coordinates and observations.
    initial : sequence of float
        Initial free-parameter vector.
    bounds : tuple of array-like
        Lower and upper free-parameter bounds.
    sigma : array-like or None
        Dependent-variable standard uncertainty used by the objective.

    Returns
    -------
    dict
        Serialization-ready input and parameter summaries.
    """
    lower, upper = bounds
    payload: dict[str, Any] = {
        "x_summary": array_summary(x),
        "y_summary": array_summary(y),
        "initial_parameters": np.asarray(initial, dtype=np.float64).tolist(),
        "lower_bounds": np.asarray(lower, dtype=np.float64).tolist(),
        "upper_bounds": np.asarray(upper, dtype=np.float64).tolist(),
    }
    if sigma is not None:
        payload["sigma_summary"] = array_summary(sigma, positive_ratio=True)
    return payload


@dataclass(slots=True)
class ModelEvaluationRecorder:
    """Record compact nonlinear model-evaluation diagnostics.

    Parameters
    ----------
    observed : array-like
        Dependent observations used by the objective.
    sigma : array-like or None, optional
        Standard uncertainties used to standardize residuals.
    detailed : bool, optional
        Retain a bounded trace of individual evaluations.
    max_trace_length : int, optional
        Number of most-recent evaluations retained in addition to the first
        evaluation.  The total backend evaluation count is never truncated.
    """

    observed: np.ndarray
    sigma: np.ndarray | None = None
    detailed: bool = False
    max_trace_length: int = _DEFAULT_TRACE_LENGTH
    n_evaluations: int = 0
    first_record: dict[str, Any] | None = None
    last_record: dict[str, Any] | None = None
    _recent: deque[dict[str, Any]] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """Normalize arrays and initialize the bounded trace container."""
        self.observed = np.asarray(self.observed, dtype=np.float64).copy()
        if self.sigma is not None:
            self.sigma = np.asarray(self.sigma, dtype=np.float64).copy()
        self.max_trace_length = max(int(self.max_trace_length), 1)
        self._recent = deque(maxlen=self.max_trace_length)

    def evaluate(self, parameters: Sequence[float], fitted: Any) -> None:
        """Record one completed model evaluation.

        Parameters
        ----------
        parameters : sequence of float
            Current free-parameter vector.
        fitted : array-like
            Model response returned to the optimizer.
        """
        self.n_evaluations += 1
        predicted = np.asarray(fitted, dtype=np.float64)
        parameter_array = np.asarray(parameters, dtype=np.float64)
        record: dict[str, Any] = {
            "evaluation": self.n_evaluations,
            "parameters": parameter_array.tolist(),
        }
        if predicted.shape == self.observed.shape and np.all(np.isfinite(predicted)):
            residuals = self.observed - predicted
            scaled = residuals if self.sigma is None else residuals / self.sigma
            record.update(
                {
                    "objective": float(np.sum(np.square(scaled))),
                    "rmse": float(np.sqrt(np.mean(np.square(residuals)))),
                    "maximum_absolute_residual": float(np.max(np.abs(residuals))),
                }
            )
        else:
            record["objective"] = None
            record["non_finite_model_response"] = True
        if self.first_record is None:
            self.first_record = dict(record)
        self.last_record = dict(record)
        if self.detailed:
            self._recent.append(dict(record))

    def as_metadata(self) -> dict[str, Any]:
        """Return a serialization-ready recorder summary."""
        payload: dict[str, Any] = {
            "recorded_model_evaluations": self.n_evaluations,
            "first_evaluation": None
            if self.first_record is None
            else dict(self.first_record),
            "last_evaluation": None
            if self.last_record is None
            else dict(self.last_record),
            "detailed_trace_requested": self.detailed,
        }
        if self.detailed:
            recent = list(self._recent)
            if self.first_record is not None and (
                not recent or recent[0].get("evaluation") != 1
            ):
                recent.insert(0, dict(self.first_record))
            payload["evaluation_trace"] = recent
            payload["evaluation_trace_truncated"] = self.n_evaluations > len(recent)
        return payload


def termination_metadata(
    *,
    category: str,
    backend: str,
    exception: Exception | None = None,
    status_code: int | str | None = None,
    backend_version: str | None = None,
) -> dict[str, Any]:
    """Return normalized termination provenance.

    Parameters
    ----------
    category : str
        Stable backend-neutral category such as ``converged``,
        ``invalid_input``, ``iteration_limit``, or ``backend_exception``.
    backend : str
        Stable backend or algorithm name.
    exception : Exception or None, optional
        Exception retained by type and message only.
    status_code : int, str, or None, optional
        Native backend termination code.
    backend_version : str or None, optional
        Backend version provenance.

    Returns
    -------
    dict
        Serialization-ready termination details.
    """
    payload: dict[str, Any] = {
        "termination_category": str(category),
        "backend": str(backend),
    }
    if backend_version is not None:
        payload["backend_version"] = str(backend_version)
    if status_code is not None:
        payload["backend_status_code"] = status_code
    if exception is not None:
        payload["exception_type"] = type(exception).__name__
        payload["exception_message"] = str(exception)
    return payload


def detect_period_two_oscillation(history: Sequence[Mapping[str, Any]]) -> bool:
    """Return whether recent free-parameter vectors indicate period-two cycling.

    Parameters
    ----------
    history : sequence of mappings
        Iteration records containing a ``parameters`` entry.

    Returns
    -------
    bool
        ``True`` when the last vector is substantially closer to the vector two
        cycles earlier than to the immediately preceding vector.
    """
    if len(history) < 4:
        return False
    try:
        vectors = [
            np.asarray(item["parameters"], dtype=np.float64) for item in history[-4:]
        ]
    except (KeyError, TypeError, ValueError):
        return False
    if any(vector.shape != vectors[0].shape for vector in vectors):
        return False
    distance_one = float(np.linalg.norm(vectors[-1] - vectors[-2]))
    distance_two = float(np.linalg.norm(vectors[-1] - vectors[-3]))
    scale = max(float(np.linalg.norm(vectors[-1])), 1.0)
    return distance_one > 1.0e-8 * scale and distance_two < 0.1 * distance_one


__all__ = [
    "ModelEvaluationRecorder",
    "SOLVER_DEBUG_METADATA_KEY",
    "array_summary",
    "detect_period_two_oscillation",
    "problem_summary",
    "solver_debug_enabled",
    "termination_metadata",
]
