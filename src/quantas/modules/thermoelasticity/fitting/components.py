# -*- coding: utf-8 -*-

"""Fitting services for cold quasi-static elastic components.

This module owns no file or frontend logic. It combines the general Quantas
fitting infrastructure with the Eulerian finite-strain relations implemented in
:mod:`quantas.core.physics.elasticity.quasistatic`.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from quantas.core.math.fitting import (
    LeastSquaresFitter,
    OLSOptions,
)
from quantas.core.physics.elasticity import (
    cold_finite_strain_component,
    cold_finite_strain_component_jacobian,
    elastic_component_definitions,
    extract_independent_stiffness_components,
    wallace_hydrostatic_delta_voigt,
)
from quantas.modules.thermoelasticity.models import (
    ElasticComponentFit,
    ElasticVolumeSeries,
    ReferenceEOSFit,
    ThermoelasticOptions,
    ThermoelasticFitQuality,
)
from quantas.modules.thermoelasticity.quality import assess_component_fit_quality

from .model import ColdFiniteStrainComponentModel

FloatArray: TypeAlias = NDArray[np.float64]


def collect_component_observations(
    series: ElasticVolumeSeries,
) -> tuple[dict[str, FloatArray], dict[str, FloatArray]]:
    """Return symmetry-averaged component series and equivalence spreads."""
    value_rows: list[dict[str, float]] = []
    spread_rows: list[dict[str, float]] = []
    for matrix in series.stiffness:
        values, spreads = extract_independent_stiffness_components(
            matrix,
            series.elastic_symmetry,
        )
        value_rows.append(values)
        spread_rows.append(spreads)
    labels = tuple(value_rows[0])
    component_values: dict[str, FloatArray] = {
        label: np.asarray([row[label] for row in value_rows], dtype=np.float64)
        for label in labels
    }
    component_spreads: dict[str, FloatArray] = {
        label: np.asarray([row[label] for row in spread_rows], dtype=np.float64)
        for label in labels
    }
    return component_values, component_spreads


def fit_elastic_components(
    series: ElasticVolumeSeries,
    reference_eos: ReferenceEOSFit,
    options: ThermoelasticOptions,
) -> dict[str, ElasticComponentFit]:
    """Fit every non-zero independent stiffness component.

    Symmetry-equivalent values are sign-corrected and averaged at each volume.
    Components below ``zero_tolerance`` remain exact zeros and receive a full
    point table but no numerical optimizer result.
    """
    observations, spreads = collect_component_observations(series)
    definitions = {
        definition.label: definition
        for definition in elastic_component_definitions(series.elastic_symmetry)
        if not definition.derived
    }
    wallace = wallace_hydrostatic_delta_voigt()
    results: dict[str, ElasticComponentFit] = {}
    for label, observed in observations.items():
        definition = definitions[label]
        representative = definition.entries[0]
        delta = float(wallace[representative[0], representative[1]])
        if float(np.max(np.abs(observed))) <= options.zero_tolerance:
            record = _zero_component_fit(
                label=label,
                entries=definition.entries,
                wallace_delta=delta,
                observed=observed,
                symmetry_spread=spreads[label],
                series=series,
                reference_eos=reference_eos,
                options=options,
            )
        else:
            record = _active_component_fit(
                label=label,
                entries=definition.entries,
                wallace_delta=delta,
                observed=observed,
                symmetry_spread=spreads[label],
                series=series,
                reference_eos=reference_eos,
                options=options,
            )
        _attach_reference_eos_sensitivity(record, reference_eos, options)
        results[label] = record
    return results


def _zero_component_fit(
    *,
    label: str,
    entries: tuple[tuple[int, int, float], ...],
    wallace_delta: float,
    observed: FloatArray,
    symmetry_spread: FloatArray,
    series: ElasticVolumeSeries,
    reference_eos: ReferenceEOSFit,
    options: ThermoelasticOptions,
) -> ElasticComponentFit:
    """Build one exact-zero component record without invoking a solver."""
    fitted = np.zeros_like(observed)
    residuals = observed - fitted
    relative_scale = np.maximum(np.abs(observed), options.relative_error_floor)
    bracketed = bool(
        series.volumes[0] <= reference_eos.reference_volume <= series.volumes[-1]
    )
    quality = ThermoelasticFitQuality(
        level="not_applicable",
        issues=("component_retained_as_exact_zero",),
        n_observations=int(series.npoints),
        n_parameters=0,
        degrees_of_freedom=int(series.npoints),
        eulerian_strain_min=0.0,
        eulerian_strain_max=0.0,
        eulerian_strain_span=0.0,
        reference_volume_bracketed=bracketed,
        reference_volume_distance_fraction=0.0,
        design_rank=0,
        design_condition_number=1.0,
        leverage=np.zeros(series.npoints, dtype=np.float64),
        maximum_leverage=0.0,
        maximum_relative_symmetry_spread=0.0,
        maximum_leave_one_out_parameter_change=None,
        maximum_order_parameter_change=None,
        thresholds={},
    )
    return ElasticComponentFit(
        label=label,
        entries=entries,
        wallace_delta=wallace_delta,
        volumes=series.volumes,
        pressures=series.pressures,
        observed=observed,
        fitted=fitted,
        residuals=residuals,
        relative_residuals=residuals / relative_scale,
        symmetry_spread=symmetry_spread,
        fit=None,
        active=False,
        zero_by_tolerance=True,
        metadata={
            "status": "numerically_zero",
            "zero_tolerance_GPa": options.zero_tolerance,
        },
        quality=quality,
    )


def _active_component_fit(
    *,
    label: str,
    entries: tuple[tuple[int, int, float], ...],
    wallace_delta: float,
    observed: FloatArray,
    symmetry_spread: FloatArray,
    series: ElasticVolumeSeries,
    reference_eos: ReferenceEOSFit,
    options: ThermoelasticOptions,
) -> ElasticComponentFit:
    """Fit one active component and build all scientific diagnostics."""
    model = ColdFiniteStrainComponentModel(
        reference_volume=reference_eos.reference_volume,
        bulk_modulus=reference_eos.bulk_modulus,
        bulk_modulus_derivative=reference_eos.bulk_modulus_derivative,
        wallace_delta=wallace_delta,
        order=options.finite_strain_order,
        label=label,
    )
    fit_options = OLSOptions(
        max_iterations=options.max_iterations,
        ftol=1.0e-15,
        xtol=1.0e-15,
        metadata={**model.metadata(), "solver_debug": options.solver_debug},
    )
    fit = LeastSquaresFitter().fit(
        model,
        series.volumes,
        observed,
        options=fit_options,
    )
    fitted, residuals, relative = _component_residual_arrays(
        fit.fitted,
        observed,
        relative_floor=options.relative_error_floor,
    )
    _record_component_fit_statistics(
        fit,
        residuals=residuals,
        relative_residuals=relative,
        symmetry_spread=symmetry_spread,
    )
    quality = _component_fit_quality(
        model=model,
        fitted_parameters=(
            None
            if fit.parameters is None
            else np.asarray(fit.parameters, dtype=np.float64)
        ),
        observed=observed,
        symmetry_spread=symmetry_spread,
        series=series,
        reference_eos=reference_eos,
        options=options,
    )
    return ElasticComponentFit(
        label=label,
        entries=entries,
        wallace_delta=wallace_delta,
        volumes=series.volumes,
        pressures=series.pressures,
        observed=observed,
        fitted=fitted,
        residuals=residuals,
        relative_residuals=relative,
        symmetry_spread=symmetry_spread,
        fit=fit,
        active=True,
        zero_by_tolerance=False,
        metadata={
            "model": model.name,
            "parameter_order": list(model.parameter_names),
        },
        quality=quality,
    )


def _component_residual_arrays(
    fitted_values: ArrayLike | None,
    observed: FloatArray,
    *,
    relative_floor: float,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Return fitted, residual, and relative-residual arrays."""
    if fitted_values is None:
        missing = np.full_like(observed, np.nan)
        return missing.copy(), missing.copy(), missing.copy()
    fitted = np.asarray(fitted_values, dtype=np.float64)
    residuals = observed - fitted
    scale = np.maximum(np.abs(observed), relative_floor)
    return fitted, residuals, residuals / scale


def _record_component_fit_statistics(
    fit,
    *,
    residuals: FloatArray,
    relative_residuals: FloatArray,
    symmetry_spread: FloatArray,
) -> None:
    """Attach deterministic residual diagnostics to one solver result."""
    if np.any(~np.isfinite(residuals)):
        return
    rss = float(np.sum(np.square(residuals)))
    reduced = float(rss / fit.dof) if fit.dof > 0 else np.nan
    fit.metadata.update(
        {
            "residual_sum_squares_GPa2": rss,
            "unweighted_chi_square_GPa2": rss,
            "reduced_unweighted_chi_square_GPa2": reduced,
            "maximum_relative_error": float(np.max(np.abs(relative_residuals))),
            "maximum_symmetry_spread_GPa": float(np.max(symmetry_spread)),
            "chi_square_interpretation": (
                "unweighted residual sum of squares; CRYSTAL SOEC "
                "uncertainties are unavailable"
            ),
        }
    )


def _component_fit_quality(
    *,
    model: ColdFiniteStrainComponentModel,
    fitted_parameters: FloatArray | None,
    observed: FloatArray,
    symmetry_spread: FloatArray,
    series: ElasticVolumeSeries,
    reference_eos: ReferenceEOSFit,
    options: ThermoelasticOptions,
) -> ThermoelasticFitQuality:
    """Evaluate scientific support for one active component fit."""
    return assess_component_fit_quality(
        volumes=series.volumes,
        observed=observed,
        symmetry_spread=symmetry_spread,
        reference_volume=reference_eos.reference_volume,
        design_matrix=component_design_matrix(model, series.volumes),
        fitted_parameters=fitted_parameters,
        leave_one_out_parameters=leave_one_out_component_parameters(
            model, series.volumes, observed
        ),
        alternate_order_parameters=alternate_order_component_parameters(
            model, series.volumes, observed
        ),
        minimum_points=options.minimum_fit_points,
        minimum_strain_span=options.minimum_eulerian_strain_span,
        maximum_design_condition=options.maximum_design_condition_number,
        maximum_relative_symmetry_spread=options.maximum_relative_symmetry_spread,
        maximum_leave_one_out_change=options.maximum_leave_one_out_parameter_change,
        maximum_order_sensitivity=options.maximum_order_parameter_change,
        relative_floor=options.relative_error_floor,
    )


def _attach_reference_eos_sensitivity(
    record: ElasticComponentFit,
    reference_eos: ReferenceEOSFit,
    options: ThermoelasticOptions,
) -> None:
    """Store the two-stage sensitivity to the shared static EOS fit."""
    sensitivity = component_parameter_sensitivity_to_reference_eos(
        record, reference_eos, options
    )
    if sensitivity is None:
        return
    record.metadata["reference_eos_parameter_sensitivity"] = sensitivity
    record.metadata["reference_eos_sensitivity_method"] = (
        "symmetric finite difference of exact conditional linear least-squares refits"
    )


def component_design_matrix(
    model: ColdFiniteStrainComponentModel,
    volume: ArrayLike,
) -> FloatArray:
    """Return the exact linear design matrix for ``C0`` and ``Cprime``.

    Parameters
    ----------
    model : ColdFiniteStrainComponentModel
        Fixed-EOS component model.
    volume : array_like
        Sampled volumes in angstrom cubed.

    Returns
    -------
    ndarray
        Matrix with shape ``(npoints, 2)``.
    """
    volumes = np.asarray(volume, dtype=np.float64)
    zero = model.evaluate(volumes, np.asarray([0.0, 0.0], dtype=np.float64))
    column_c0 = model.evaluate(volumes, np.asarray([1.0, 0.0], dtype=np.float64)) - zero
    column_cp = model.evaluate(volumes, np.asarray([0.0, 1.0], dtype=np.float64)) - zero
    return np.asarray(np.column_stack((column_c0, column_cp)), dtype=np.float64)


def exact_component_parameters(
    model: ColdFiniteStrainComponentModel,
    volume: ArrayLike,
    observed: ArrayLike,
) -> FloatArray | None:
    """Return the exact conditional OLS parameters or ``None`` if singular."""
    volumes = np.asarray(volume, dtype=np.float64)
    values = np.asarray(observed, dtype=np.float64)
    design = component_design_matrix(model, volumes)
    zero = model.evaluate(volumes, np.asarray([0.0, 0.0], dtype=np.float64))
    parameters, _, rank, _ = np.linalg.lstsq(design, values - zero, rcond=None)
    if rank < 2 or np.any(~np.isfinite(parameters)):
        return None
    return np.asarray(parameters, dtype=np.float64)


def leave_one_out_component_parameters(
    model: ColdFiniteStrainComponentModel,
    volume: ArrayLike,
    observed: ArrayLike,
) -> FloatArray | None:
    """Return exact leave-one-out estimates for support diagnostics."""
    volumes = np.asarray(volume, dtype=np.float64)
    values = np.asarray(observed, dtype=np.float64)
    if volumes.size <= 2:
        return None
    estimates = np.empty((volumes.size, 2), dtype=np.float64)
    for index in range(volumes.size):
        mask = np.ones(volumes.size, dtype=np.bool_)
        mask[index] = False
        estimate = exact_component_parameters(model, volumes[mask], values[mask])
        if estimate is None:
            return None
        estimates[index] = estimate
    return estimates


def alternate_order_component_parameters(
    model: ColdFiniteStrainComponentModel,
    volume: ArrayLike,
    observed: ArrayLike,
) -> FloatArray | None:
    """Return the exact estimate for the other supported strain order."""
    alternate = ColdFiniteStrainComponentModel(
        reference_volume=model.reference_volume,
        bulk_modulus=model.bulk_modulus,
        bulk_modulus_derivative=model.bulk_modulus_derivative,
        wallace_delta=model.wallace_delta,
        order=2 if model.order == 3 else 3,
        label=model.label,
    )
    return exact_component_parameters(alternate, volume, observed)


def component_parameter_sensitivity_to_reference_eos(
    record: ElasticComponentFit,
    reference_eos: ReferenceEOSFit,
    options: ThermoelasticOptions,
) -> FloatArray | None:
    """Return ``d(C0,Cprime)/d(V0,K0,Kprime)`` for a two-stage fit.

    The elastic parameters are conditionally re-fitted at symmetric
    perturbations of each fixed static-EOS parameter.  Because the component
    model is linear in ``C0`` and ``Cprime``, each conditional re-fit is an
    exact linear least-squares solve rather than a second nonlinear optimizer.
    This sensitivity is required for rigorous first-order propagation of the
    shared EOS covariance through the complete two-stage estimator.

    Parameters
    ----------
    record : ElasticComponentFit
        Successful active component fit.
    reference_eos : ReferenceEOSFit
        Fixed EOS parameter estimate.
    options : ThermoelasticOptions
        Finite-strain controls.

    Returns
    -------
    ndarray or None
        Matrix with shape ``(2, 3)`` and columns ``V0, K0, Kprime``.
        ``None`` is returned for exact zeros and failed fits.
    """
    if record.zero_by_tolerance:
        return np.zeros((2, 3), dtype=np.float64)
    if record.fit is None or not record.fit.success or record.fit.parameters is None:
        return None
    theta = np.asarray(
        [
            reference_eos.reference_volume,
            reference_eos.bulk_modulus,
            reference_eos.bulk_modulus_derivative,
        ],
        dtype=np.float64,
    )
    scale = np.maximum(np.abs(theta), np.asarray([1.0, 1.0, 1.0]))
    steps = np.cbrt(np.finfo(np.float64).eps) * scale
    sensitivity = np.empty((2, 3), dtype=np.float64)
    for index, step in enumerate(steps):
        plus = theta.copy()
        minus = theta.copy()
        plus[index] += step
        minus[index] -= step
        if index in (0, 1) and minus[index] <= 0.0:
            minus[index] = theta[index] - 0.5 * step
        beta_plus = _conditional_component_parameters(record, plus, options)
        beta_minus = _conditional_component_parameters(record, minus, options)
        sensitivity[:, index] = (beta_plus - beta_minus) / (plus[index] - minus[index])
    return sensitivity


def _conditional_component_parameters(
    record: ElasticComponentFit,
    theta: FloatArray,
    options: ThermoelasticOptions,
) -> FloatArray:
    """Return the exact conditional OLS component parameters for one EOS."""
    model = ColdFiniteStrainComponentModel(
        reference_volume=float(theta[0]),
        bulk_modulus=float(theta[1]),
        bulk_modulus_derivative=float(theta[2]),
        wallace_delta=record.wallace_delta,
        order=options.finite_strain_order,
        label=record.label,
    )
    return np.asarray(
        model.initial_guess(record.volumes, record.observed),
        dtype=np.float64,
    )


def evaluate_component_predictions(
    fits: Mapping[str, ElasticComponentFit],
    labels: tuple[str, ...],
    volume: ArrayLike,
    reference_eos: ReferenceEOSFit,
    options: ThermoelasticOptions,
    *,
    sigma_volume: ArrayLike | None = None,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Evaluate independent components and their full local covariance.

    Returns
    -------
    tuple of ndarray
        Values, one-sigma uncertainties, and covariance matrices with shapes
        ``volume.shape + (ncomponents,)`` and
        ``volume.shape + (ncomponents, ncomponents)``.
    """
    volumes = np.asarray(volume, dtype=np.float64)
    flat = volumes.ravel()
    ncomponents = len(labels)
    values = np.full((flat.size, ncomponents), np.nan, dtype=np.float64)
    covariance = np.zeros((flat.size, ncomponents, ncomponents), dtype=np.float64)
    eos_jacobian = np.zeros((flat.size, ncomponents, 3), dtype=np.float64)
    volume_jacobian = np.zeros((flat.size, ncomponents), dtype=np.float64)

    for component_index, label in enumerate(labels):
        record = fits[label]
        if record.zero_by_tolerance:
            values[:, component_index] = 0.0
            continue
        parameters = record.parameters
        if parameters is None:
            continue
        predicted = cold_finite_strain_component(
            flat,
            reference_volume=reference_eos.reference_volume,
            bulk_modulus=reference_eos.bulk_modulus,
            bulk_modulus_derivative=reference_eos.bulk_modulus_derivative,
            reference_component=float(parameters[0]),
            component_pressure_derivative=float(parameters[1]),
            wallace_delta=record.wallace_delta,
            order=options.finite_strain_order,
        )
        jacobian = cold_finite_strain_component_jacobian(
            flat,
            reference_volume=reference_eos.reference_volume,
            bulk_modulus=reference_eos.bulk_modulus,
            bulk_modulus_derivative=reference_eos.bulk_modulus_derivative,
            reference_component=float(parameters[0]),
            component_pressure_derivative=float(parameters[1]),
            wallace_delta=record.wallace_delta,
            order=options.finite_strain_order,
        )
        values[:, component_index] = predicted
        fit_covariance = record.covariance
        if fit_covariance is not None:
            local = np.asarray(jacobian[:, :2], dtype=np.float64)
            covariance[:, component_index, component_index] += np.einsum(
                "ni,ij,nj->n", local, fit_covariance, local
            )
        conditional_sensitivity = record.metadata.get(
            "reference_eos_parameter_sensitivity"
        )
        if conditional_sensitivity is None:
            conditional_sensitivity = component_parameter_sensitivity_to_reference_eos(
                record, reference_eos, options
            )
        total_eos_jacobian = np.asarray(jacobian[:, [2, 3, 4]], dtype=np.float64)
        if conditional_sensitivity is not None:
            sensitivity_array = np.asarray(conditional_sensitivity, dtype=np.float64)
            if sensitivity_array.shape != (2, 3):
                raise ValueError(
                    "reference EOS parameter sensitivity must have shape (2, 3)"
                )
            total_eos_jacobian += np.asarray(jacobian[:, :2]) @ sensitivity_array
        eos_jacobian[:, component_index, :] = total_eos_jacobian
        volume_jacobian[:, component_index] = jacobian[:, 5]

    if (
        options.propagate_reference_eos_covariance
        and reference_eos.covariance is not None
    ):
        covariance += np.einsum(
            "nai,ij,nbj->nab",
            eos_jacobian,
            reference_eos.covariance,
            eos_jacobian,
        )
    if options.propagate_volume_uncertainty and sigma_volume is not None:
        sigma = np.asarray(sigma_volume, dtype=np.float64)
        if sigma.shape != volumes.shape:
            raise ValueError("sigma_volume must match the evaluation-volume shape")
        flat_sigma = sigma.ravel()
        valid = np.isfinite(flat_sigma) & (flat_sigma >= 0.0)
        if np.any(valid):
            covariance[valid] += (
                volume_jacobian[valid, :, np.newaxis]
                * volume_jacobian[valid, np.newaxis, :]
                * np.square(flat_sigma[valid])[:, np.newaxis, np.newaxis]
            )
    covariance[...] = 0.5 * (
        covariance + np.swapaxes(covariance, -1, -2)
    )
    variances = np.clip(
        np.diagonal(covariance, axis1=-2, axis2=-1),
        0.0,
        None,
    )
    errors = np.sqrt(variances)
    value_shape = volumes.shape + (ncomponents,)
    covariance_shape = volumes.shape + (ncomponents, ncomponents)
    return (
        values.reshape(value_shape),
        errors.reshape(value_shape),
        covariance.reshape(covariance_shape),
    )


__all__ = [
    "alternate_order_component_parameters",
    "collect_component_observations",
    "component_design_matrix",
    "component_parameter_sensitivity_to_reference_eos",
    "evaluate_component_predictions",
    "exact_component_parameters",
    "fit_elastic_components",
    "leave_one_out_component_parameters",
]
