# -*- coding: utf-8 -*-

"""Scientific support diagnostics for quasi-static component fits."""

from __future__ import annotations

from typing import TypeAlias

import numpy as np
from numpy.typing import NDArray

from quantas.core.physics.elasticity import eulerian_finite_strain
from quantas.modules.thermoelasticity.models import (
    ThermoelasticFitQuality,
    ThermoelasticQualityLevel,
)


FloatArray: TypeAlias = NDArray[np.float64]


def assess_component_fit_quality(
    *,
    volumes: FloatArray,
    observed: FloatArray,
    symmetry_spread: FloatArray,
    reference_volume: float,
    design_matrix: FloatArray,
    fitted_parameters: FloatArray | None,
    leave_one_out_parameters: FloatArray | None,
    alternate_order_parameters: FloatArray | None,
    minimum_points: int,
    minimum_strain_span: float,
    maximum_design_condition: float,
    maximum_relative_symmetry_spread: float,
    maximum_leave_one_out_change: float,
    maximum_order_sensitivity: float,
    relative_floor: float,
) -> ThermoelasticFitQuality:
    """Classify whether one two-parameter QSA fit is scientifically supported.

    The assessment is deliberately diagnostic.  It does not modify observations,
    fitted parameters, or predicted tensors.  Thresholds are user-visible
    heuristics intended to distinguish numerical convergence from adequate data
    support.

    Parameters
    ----------
    volumes, observed, symmetry_spread : ndarray
        Aligned sampled data.
    reference_volume : float
        Static-EOS ``V0`` in the same unit as ``volumes``.
    design_matrix : ndarray
        Linear two-column design matrix for ``C0`` and ``Cprime``.
    fitted_parameters : ndarray or None
        Full-data parameter estimate.
    leave_one_out_parameters : ndarray or None
        Conditional estimates with each observation omitted in turn.
    alternate_order_parameters : ndarray or None
        Estimate obtained with the other supported finite-strain order.
    minimum_points, minimum_strain_span, maximum_design_condition : scalar
        User-configurable support thresholds.
    maximum_relative_symmetry_spread : float
        Maximum accepted symmetry-equivalent spread relative to component size.
    maximum_leave_one_out_change, maximum_order_sensitivity : float
        Maximum accepted relative parameter changes.
    relative_floor : float
        Positive GPa floor for relative comparisons.

    Returns
    -------
    ThermoelasticFitQuality
        Passive diagnostic record with level, issues, and point-level leverage.
    """
    volume = np.asarray(volumes, dtype=np.float64)
    values = np.asarray(observed, dtype=np.float64)
    spread = np.asarray(symmetry_spread, dtype=np.float64)
    design = np.asarray(design_matrix, dtype=np.float64)
    npoints = int(volume.size)
    nparameters = int(design.shape[1]) if design.ndim == 2 else 0
    rank = int(np.linalg.matrix_rank(design)) if design.ndim == 2 else 0
    condition = _normalized_design_condition(design)
    leverage = _leverage(design)
    strain = eulerian_finite_strain(volume, reference_volume)
    strain_min = float(np.min(strain))
    strain_max = float(np.max(strain))
    strain_span = strain_max - strain_min
    volume_min = float(np.min(volume))
    volume_max = float(np.max(volume))
    bracketed = volume_min <= reference_volume <= volume_max
    interval = float(max(volume_max - volume_min, float(np.finfo(np.float64).eps)))
    if reference_volume < volume_min:
        reference_distance = (volume_min - reference_volume) / interval
    elif reference_volume > volume_max:
        reference_distance = (reference_volume - volume_max) / interval
    else:
        reference_distance = 0.0
    scale = np.maximum(np.abs(values), relative_floor)
    max_relative_spread = float(np.max(np.abs(spread) / scale))
    loo_change = _maximum_parameter_change(
        fitted_parameters,
        leave_one_out_parameters,
        observed_scale=float(np.max(np.abs(values), initial=1.0)),
    )
    order_change = _maximum_parameter_change(
        fitted_parameters,
        alternate_order_parameters,
        observed_scale=float(np.max(np.abs(values), initial=1.0)),
    )

    unsupported: list[str] = []
    cautions: list[str] = []
    if npoints <= nparameters:
        unsupported.append("observations_do_not_exceed_parameter_count")
    if rank < nparameters:
        unsupported.append("rank_deficient_design_matrix")
    if not np.isfinite(condition):
        unsupported.append("singular_design_matrix")
    if npoints < minimum_points:
        cautions.append("few_observations")
    if strain_span < minimum_strain_span:
        cautions.append("narrow_eulerian_strain_span")
    if not bracketed:
        cautions.append("reference_volume_not_bracketed")
    if condition > maximum_design_condition:
        cautions.append("ill_conditioned_design_matrix")
    if max_relative_spread > maximum_relative_symmetry_spread:
        cautions.append("large_symmetry_equivalent_spread")
    if loo_change is not None and loo_change > maximum_leave_one_out_change:
        cautions.append("large_leave_one_out_parameter_sensitivity")
    if order_change is not None and order_change > maximum_order_sensitivity:
        cautions.append("large_finite_strain_order_sensitivity")

    if unsupported:
        level: ThermoelasticQualityLevel = "unsupported"
    elif cautions:
        level = "caution"
    else:
        level = "supported"
    return ThermoelasticFitQuality(
        level=level,
        issues=tuple(unsupported + cautions),
        n_observations=npoints,
        n_parameters=nparameters,
        degrees_of_freedom=max(npoints - nparameters, 0),
        eulerian_strain_min=strain_min,
        eulerian_strain_max=strain_max,
        eulerian_strain_span=strain_span,
        reference_volume_bracketed=bracketed,
        reference_volume_distance_fraction=float(reference_distance),
        design_rank=rank,
        design_condition_number=float(condition),
        leverage=leverage,
        maximum_leverage=(float(np.max(leverage)) if leverage.size else 0.0),
        maximum_relative_symmetry_spread=max_relative_spread,
        maximum_leave_one_out_parameter_change=loo_change,
        maximum_order_parameter_change=order_change,
        thresholds={
            "minimum_points": int(minimum_points),
            "minimum_eulerian_strain_span": float(minimum_strain_span),
            "maximum_design_condition_number": float(maximum_design_condition),
            "maximum_relative_symmetry_spread": float(maximum_relative_symmetry_spread),
            "maximum_leave_one_out_parameter_change": float(
                maximum_leave_one_out_change
            ),
            "maximum_order_parameter_change": float(maximum_order_sensitivity),
        },
    )


def _normalized_design_condition(design: FloatArray) -> float:
    """Return a column-scale-invariant two-norm condition number."""
    if design.ndim != 2 or design.shape[1] == 0:
        return float("inf")
    norms = np.linalg.norm(design, axis=0)
    if np.any(norms <= np.finfo(np.float64).eps):
        return float("inf")
    return float(np.linalg.cond(design / norms))


def _leverage(design: FloatArray) -> FloatArray:
    """Return diagonal entries of the ordinary least-squares hat matrix."""
    if design.ndim != 2 or design.shape[0] == 0:
        return np.empty(0, dtype=np.float64)
    hat = design @ np.linalg.pinv(design)
    return np.asarray(np.clip(np.diag(hat), 0.0, 1.0), dtype=np.float64)


def _maximum_parameter_change(
    reference: FloatArray | None,
    alternatives: FloatArray | None,
    *,
    observed_scale: float,
) -> float | None:
    """Return the largest scale-aware relative parameter change."""
    if reference is None or alternatives is None:
        return None
    base = np.asarray(reference, dtype=np.float64)
    other = np.asarray(alternatives, dtype=np.float64)
    if base.shape != (2,) or other.shape[-1:] != (2,) or np.any(~np.isfinite(other)):
        return None
    denominator = np.asarray(
        [max(abs(base[0]), observed_scale, 1.0), max(abs(base[1]), 1.0)],
        dtype=np.float64,
    )
    return float(np.max(np.abs(other - base) / denominator))


__all__ = ["assess_component_fit_quality"]
