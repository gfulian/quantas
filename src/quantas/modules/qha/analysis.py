# -*- coding: utf-8 -*-

"""Orchestrate the numerical stages of quasi-harmonic analysis.

The functions combine normalized inputs, fitting, minimization, diagnostics,
validity masks, and partial-result policies into structured QHA results.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np

from quantas.core.physics.units import energy_to_pressure
from quantas.modules.qha.core.minimization import (
    VolumeMinimumResult,
    minimize_eos,
    minimize_polynomial,
    target_pressure_energy_density,
)
from quantas.modules.qha.models import (
    QHAFailedPoint,
    QHAInput,
    QHAOptions,
    QHAResult,
    result_metadata_from_options,
)

ArrayLike = np.ndarray | Sequence[float]


@dataclass(slots=True)
class QHAWorkflowState:
    """Mutable state used while evaluating a QHA pressure-temperature grid.

    Parameters
    ----------
    consecutive_failures : int, optional
        Number of consecutive failed pressure-temperature points.
    stopped : bool, optional
        Whether the workflow requested an early stop.
    stop_message : str, optional
        Explanation associated with an early stop.
    metadata : dict, optional
        Additional workflow details.
    """

    consecutive_failures: int = 0
    stopped: bool = False
    stop_message: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)


def initialize_result(input_data: QHAInput, options: QHAOptions) -> QHAResult:
    """Create an empty QHA result for the requested calculation grid.

    Parameters
    ----------
    input_data : QHAInput
        Normalized QHA input data.
    options : QHAOptions
        Workflow options.

    Returns
    -------
    QHAResult
        Result container initialized with input arrays, grids, metadata and an
        empty validity mask.

    Raises
    ------
    ValueError
        If input data or options are inconsistent.
    """
    input_data.validate_shapes()
    options.validate()

    temperature = options.temperature_grid()
    pressure = options.pressure_grid()
    shape = (temperature.size, pressure.size)
    metadata = result_metadata_from_options(options)
    metadata["input"] = {
        "natoms": int(input_data.natoms),
        "formula_units": int(input_data.formula_units),
        "natoms_per_formula_unit": float(input_data.natoms_per_formula_unit),
        "qpoints": int(input_data.qpoints),
        "nvol": int(input_data.nvol),
        "nmodes": int(input_data.nmodes),
        "mode_continuity": input_data.mode_continuity_status(),
    }
    metadata["thermodynamic_unit_convention"] = "native_energy_per_cell_per_kelvin"
    metadata["normalization"] = {
        "native_basis": "cell",
        "formula_units_per_cell": int(input_data.formula_units),
        "natoms_per_cell": int(input_data.natoms),
        "natoms_per_formula_unit": float(input_data.natoms_per_formula_unit),
        "molar_basis": "formula_unit",
    }

    return QHAResult(
        jobname=input_data.jobname,
        temperature=temperature,
        pressure=pressure,
        volume=np.asarray(input_data.volume, dtype=np.float64).copy(),
        static_energy=np.asarray(input_data.energy, dtype=np.float64).copy(),
        equilibrium_volume=np.full(shape, np.nan, dtype=np.float64),
        isothermal_bulk_modulus=np.full(shape, np.nan, dtype=np.float64),
        bulk_modulus_derivative=np.full(shape, np.nan, dtype=np.float64),
        valid_mask=np.zeros(shape, dtype=np.bool_),
        completed=True,
        metadata=metadata,
    )


def prepare_free_energy_grid(
    free_energy: ArrayLike,
    *,
    ntemperatures: int,
    nvolumes: int,
) -> np.ndarray:
    """Return free energies with shape ``(ntemperatures, nvolumes)``.

    Parameters
    ----------
    free_energy : array-like
        Free energies sampled on the volume grid.  A one-dimensional array is
        interpreted as temperature-independent and is broadcast over the
        requested temperature grid.
    ntemperatures : int
        Number of temperature points.
    nvolumes : int
        Number of volume points.

    Returns
    -------
    ndarray
        Two-dimensional free-energy array.

    Raises
    ------
    ValueError
        If the input shape is not compatible with the requested grid.
    """
    array = np.asarray(free_energy, dtype=np.float64)
    if array.ndim == 1:
        if array.size != nvolumes:
            raise ValueError("free_energy must match the number of volumes")
        return np.broadcast_to(array, (ntemperatures, nvolumes)).copy()
    if array.ndim == 2:
        if array.shape != (ntemperatures, nvolumes):
            raise ValueError("free_energy must have shape (ntemperatures, nvolumes)")
        return array.copy()
    raise ValueError("free_energy must be one- or two-dimensional")


def pressure_energy_density(pressure: float, options: QHAOptions) -> float:
    """Convert pressure to the energy-density scale used in minimization.

    Parameters
    ----------
    pressure : float
        Pressure value in ``options.pressure_unit``.
    options : QHAOptions
        Workflow options defining energy, volume and pressure units.

    Returns
    -------
    float
        Pressure expressed as energy per volume.
    """
    return target_pressure_energy_density(
        pressure,
        energy_unit=options.energy_unit,
        volume_unit=options.volume_unit,
        pressure_unit=options.pressure_unit,
    )


def minimize_volume_at_point(
    volume: ArrayLike,
    free_energy: ArrayLike,
    *,
    pressure: float,
    options: QHAOptions,
) -> VolumeMinimumResult:
    """Minimize a free-energy curve at a single pressure.

    Parameters
    ----------
    volume : array-like
        Sampled volume grid.
    free_energy : array-like
        Free-energy values sampled on ``volume``.
    pressure : float
        Target pressure.
    options : QHAOptions
        Workflow options selecting polynomial or EOS minimization.

    Returns
    -------
    VolumeMinimumResult
        Local minimization result and diagnostics.
    """
    p_edens = pressure_energy_density(pressure, options)
    if options.minimization == "poly":
        minimum = minimize_polynomial(
            volume,
            free_energy,
            pressure_energy_density=p_edens,
            degree=options.free_energy_degree,
            derivative_method=options.polynomial_derivative_method,
            local_grid_points=options.polynomial_grid_points,
            local_grid_separation=options.polynomial_grid_separation,
            local_degree=options.energy_degree,
        )
        if minimum.bulk_modulus is not None:
            minimum.bulk_modulus = float(
                energy_to_pressure(
                    minimum.bulk_modulus,
                    options.energy_unit,
                    options.volume_unit,
                    options.pressure_unit,
                )
            )
        return minimum
    if options.minimization == "eos":
        return minimize_eos(
            volume,
            free_energy,
            pressure_energy_density=p_edens,
            eos=options.eos,
            energy_unit=options.energy_unit,
            volume_unit=options.volume_unit,
            pressure_unit=options.pressure_unit,
        )
    return VolumeMinimumResult.failed(
        f"unsupported minimization method: {options.minimization}",
        method=str(options.minimization),
    )


def should_store_fit_record(
    minimum: VolumeMinimumResult,
    *,
    options: QHAOptions,
) -> bool:
    """Return whether a local fit diagnostic should be stored.

    Parameters
    ----------
    minimum : VolumeMinimumResult
        Local minimization result.
    options : QHAOptions
        Workflow options controlling diagnostics.

    Returns
    -------
    bool
        ``True`` if the fit record should be retained.
    """
    if not options.store_fit_diagnostics:
        return False
    return options.debug or (not minimum.success) or bool(minimum.warnings)


def apply_minimum_to_result(
    result: QHAResult,
    minimum: VolumeMinimumResult,
    *,
    temperature_index: int,
    pressure_index: int,
) -> None:
    """Store a successful local minimum in a result object.

    Parameters
    ----------
    result : QHAResult
        Result container to update.
    minimum : VolumeMinimumResult
        Successful local minimization result.
    temperature_index : int
        Temperature-grid index.
    pressure_index : int
        Pressure-grid index.
    """
    if result.equilibrium_volume is not None:
        result.equilibrium_volume[temperature_index, pressure_index] = (
            np.nan if minimum.volume is None else minimum.volume
        )
    if result.isothermal_bulk_modulus is not None:
        value = minimum.bulk_modulus
        result.isothermal_bulk_modulus[temperature_index, pressure_index] = (
            np.nan if value is None else value
        )
    if result.bulk_modulus_derivative is not None:
        value = minimum.bulk_modulus_derivative
        result.bulk_modulus_derivative[temperature_index, pressure_index] = (
            np.nan if value is None else value
        )
    if result.valid_mask is not None:
        result.valid_mask[temperature_index, pressure_index] = True

    uncertainty_values = {
        "sigma_VT": minimum.sigma_volume,
        "sigma_KT": minimum.sigma_bulk_modulus,
        "sigma_Kp": minimum.sigma_bulk_modulus_derivative,
    }
    shape = None if result.valid_mask is None else result.valid_mask.shape
    for key, value in uncertainty_values.items():
        if value is None or shape is None:
            continue
        if key not in result.uncertainties:
            result.uncertainties[key] = np.full(shape, np.nan, dtype=np.float64)
        result.uncertainties[key][temperature_index, pressure_index] = float(value)


def register_failure(
    result: QHAResult,
    state: QHAWorkflowState,
    *,
    temperature: float,
    pressure: float,
    stage: str,
    message: str,
    diagnostics: dict[str, Any] | None = None,
    options: QHAOptions,
) -> None:
    """Record a failed pressure-temperature point and update workflow state.

    Parameters
    ----------
    result : QHAResult
        Result container to update.
    state : QHAWorkflowState
        Mutable workflow state.
    temperature : float
        Temperature associated with the failure.
    pressure : float
        Pressure associated with the failure.
    stage : str
        Workflow stage where the failure occurred.
    message : str
        Failure message.
    diagnostics : dict, optional
        Additional failure details.
    options : QHAOptions
        Workflow options controlling failure policy.

    Raises
    ------
    RuntimeError
        If ``options.fit_failure_policy`` is ``"raise"``.
    """
    result.add_failed_point(
        QHAFailedPoint(
            temperature=float(temperature),
            pressure=float(pressure),
            stage=stage,
            message=message,
            diagnostics=dict(diagnostics or {}),
        )
    )
    state.consecutive_failures += 1
    if options.fit_failure_policy == "raise":
        raise RuntimeError(message)
    if (
        options.fit_failure_policy == "stop"
        and state.consecutive_failures >= options.max_consecutive_failures
    ):
        state.stopped = True
        state.stop_message = (
            f"stopped after {state.consecutive_failures} consecutive failed fits"
        )
        result.completed = False
        result.metadata.setdefault("stop", {})
        result.metadata["stop"].update(
            {
                "message": state.stop_message,
                "temperature": float(temperature),
                "pressure": float(pressure),
                "max_consecutive_failures": int(options.max_consecutive_failures),
            }
        )


def analyze_volume_minimization(
    input_data: QHAInput,
    options: QHAOptions,
    *,
    free_energy: ArrayLike | None = None,
) -> QHAResult:
    """Evaluate equilibrium volumes on the requested pressure-temperature grid.

    This convenience entry point delegates to the shared QHA workflow
    controller so that direct numerical use and calculator-driven execution
    follow the same fitting, minimization, diagnostic, and failure policies.

    Parameters
    ----------
    input_data : QHAInput
        Normalized input data containing sampled volumes and static energies.
    options : QHAOptions
        Workflow options controlling pressure and temperature grids,
        minimization method, diagnostics and failure policy.
    free_energy : array-like or None, optional
        Free-energy data sampled on the input volume grid. If omitted, the
        static energies stored in ``input_data`` are used at all temperatures.

    Returns
    -------
    QHAResult
        Result containing equilibrium volumes, thermoelastic quantities,
        validity masks, fit diagnostics and failed-point records.

    Raises
    ------
    ValueError
        If input arrays or options are inconsistent.
    RuntimeError
        If the selected failure policy requests exceptions on local failures.
    """
    from quantas.modules.qha.workflow import run_volume_minimization_workflow

    return run_volume_minimization_workflow(
        input_data,
        options,
        free_energy=(
            None if free_energy is None else np.asarray(free_energy, dtype=np.float64)
        ),
    )
