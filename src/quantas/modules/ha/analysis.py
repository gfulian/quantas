# -*- coding: utf-8 -*-

"""
High-level numerical workflow for harmonic-approximation calculations.

The functions defined here form the library-level HA workflow. They validate a
normalized :class:`~quantas.modules.ha.models.HAInput`, prepare unit-converted
arrays, call the low-level thermodynamic routines, and collect the numerical
results in a passive :class:`~quantas.modules.ha.models.HAResult` container.

No function in this module writes to the terminal, opens progress bars, or
performs frontend-specific rendering. Long or multi-step workflows may report
progress through an optional callback supplied by calculators or frontends.
"""

from __future__ import annotations

from collections.abc import Callable
from time import perf_counter
from typing import Any

import numpy as np

from quantas.core.physics.units import (
    convert_energy,
    convert_energy_per_temperature,
    convert_frequency,
    convert_temperature,
)
from quantas.core.physics.thermodynamics import (
    entropy,
    free_energy,
    internal_energy,
    isochoric_heat_capacity,
    thermal_energy,
    vibrational_free_energy,
    zero_point_energy,
)
from quantas.modules.ha.models import HAInput, HAOptions, HAResult


ProgressCallback = Callable[[str, int, int], None]
StepCallback = Callable[[str, str], None]
ResultCallback = Callable[[str, dict[str, Any]], None]
TimingCallback = Callable[[str, float, str], None]


THERMODYNAMIC_STEPS = (
    "zero_point_energy",
    "thermal_energy",
    "entropy",
    "isochoric_heat_capacity",
    "vibrational_free_energy",
    "total_energies",
)


def validate_input(input_data: HAInput) -> None:
    """
    Validate harmonic-approximation input data.

    Parameters
    ----------
    input_data : HAInput
        Normalized harmonic-approximation input data.

    Raises
    ------
    ValueError
        If required arrays are missing or have incompatible shapes.
    """
    if input_data.volume is None:
        raise ValueError("HA input volumes are missing.")
    if input_data.energy is None:
        raise ValueError("HA input static energies are missing.")
    if input_data.frequencies is None:
        raise ValueError("HA input phonon frequencies are missing.")
    if input_data.weights is None:
        raise ValueError("HA input q-point weights are missing.")

    volume = np.asarray(input_data.volume, dtype=np.float64)
    energy = np.asarray(input_data.energy, dtype=np.float64)
    frequencies = np.asarray(input_data.frequencies, dtype=np.float64)
    weights = np.asarray(input_data.weights, dtype=np.float64)

    if volume.ndim != 1:
        raise ValueError("HA volumes must be a one-dimensional array.")
    if energy.ndim != 1:
        raise ValueError("HA static energies must be a one-dimensional array.")
    if volume.shape[0] == 0:
        raise ValueError("HA input must contain at least one volume.")
    if energy.shape[0] != volume.shape[0]:
        raise ValueError("HA static energies and volumes must have the same length.")

    if input_data.natoms <= 0:
        raise ValueError("The number of atoms must be positive.")
    if input_data.formula_units <= 0:
        raise ValueError("The number of formula units must be positive.")
    expected_modes = int(input_data.natoms) * 3

    if frequencies.ndim != 3:
        raise ValueError(
            "HA phonon frequencies must have shape (qpoints, modes, volumes)."
        )
    if frequencies.shape[1] != expected_modes:
        raise ValueError(
            "The number of phonon modes must be equal to three times natoms."
        )
    if frequencies.shape[2] != volume.shape[0]:
        raise ValueError(
            "The phonon frequency volume dimension is incompatible with volumes."
        )

    if weights.ndim != 1:
        raise ValueError("HA q-point weights must be a one-dimensional array.")
    if weights.shape[0] != frequencies.shape[0]:
        raise ValueError(
            "The number of q-point weights must match the phonon q-points."
        )
    if np.sum(weights) <= 0.0:
        raise ValueError("The sum of q-point weights must be positive.")

    if input_data.qpoints and int(input_data.qpoints) != frequencies.shape[0]:
        raise ValueError(
            "The declared number of q-points is incompatible with phonon data."
        )


def prepare_temperature_grid(options: HAOptions) -> np.ndarray:
    """
    Build the HA temperature grid in Kelvin.

    Parameters
    ----------
    options : HAOptions
        Workflow options defining the temperature range and temperature unit.

    Returns
    -------
    ndarray
        One-dimensional temperature grid converted to Kelvin.

    Raises
    ------
    ValueError
        If the temperature range is invalid.
    NotImplementedError
        If the requested temperature unit is unsupported.
    """
    grid = options.temperature_grid()
    return np.asarray(
        convert_temperature(grid, options.temperature_unit, "K"),
        dtype=np.float64,
    )


def prepare_phonon_data(
    input_data: HAInput,
    options: HAOptions,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Return unit-converted arrays required by HA thermodynamics.

    Parameters
    ----------
    input_data : HAInput
        Normalized harmonic-approximation input data.
    options : HAOptions
        Workflow options defining frequency and temperature units.

    Returns
    -------
    tuple of ndarray
        ``temperature``, ``frequencies``, ``weights``, and ``static_energy``.
        Temperatures are expressed in Kelvin, frequencies in Hertz, weights are
        normalized by their sum, and static energies are kept in the requested
        output energy unit, matching the HA workflow output convention.

    Raises
    ------
    ValueError
        If the input data are incomplete or inconsistent.
    NotImplementedError
        If a requested unit conversion is unsupported.
    """
    validate_input(input_data)

    temperature = prepare_temperature_grid(options)
    frequencies = np.asarray(
        convert_frequency(input_data.frequencies, options.frequency_unit, "Hz"),
        dtype=np.float64,
    )
    weights = input_data.normalized_weights()
    static_energy = np.asarray(input_data.energy, dtype=np.float64)

    return temperature, frequencies, weights, static_energy


def calculate_thermodynamic_properties(
    input_data: HAInput,
    options: HAOptions | None = None,
    progress_callback: ProgressCallback | None = None,
    step_callback: StepCallback | None = None,
    result_callback: ResultCallback | None = None,
    timing_callback: TimingCallback | None = None,
) -> HAResult:
    """
    Calculate harmonic thermodynamic properties.

    Parameters
    ----------
    input_data : HAInput
        Normalized harmonic-approximation input data.
    options : HAOptions or None, optional
        Options controlling units and temperature range. If ``None``, default
        options are used.
    progress_callback : callable or None, optional
        Optional callback receiving ``label``, ``current`` and ``total`` after
        each major calculation step.
    step_callback : callable or None, optional
        Optional callback receiving ``label`` and ``status`` at the beginning
        and end of each major calculation step.
    result_callback : callable or None, optional
        Optional callback receiving a result label and a structured payload each
        time a thermodynamic quantity has been calculated.
    timing_callback : callable or None, optional
        Optional callback receiving step label, elapsed seconds, and backend name
        after each timed numerical operation.

    Returns
    -------
    HAResult
        Harmonic thermodynamic results.

    Raises
    ------
    ValueError
        If the input data are incomplete or inconsistent.
    NotImplementedError
        If a requested unit conversion is unsupported.
    """
    options = options or HAOptions()
    temperature, frequencies, weights, static_energy = prepare_phonon_data(
        input_data,
        options,
    )
    backend_name = "numpy"

    timing_records: list[dict[str, Any]] = []

    def record_timing(label: str, elapsed: float, backend: str) -> None:
        """Store one timing record and forward it to the caller callback.

        Parameters
        ----------
        label : str
            Name of the timed numerical operation.
        elapsed : float
            Elapsed wall-clock time in seconds.
        backend : str
            Numerical backend identifier.
        """
        timing_records.append(
            {
                "label": label,
                "elapsed_seconds": float(elapsed),
                "backend": backend,
            }
        )
        if timing_callback is not None:
            timing_callback(label, elapsed, backend)

    result = HAResult(
        jobname=input_data.jobname,
        temperature=np.asarray(
            options.temperature_grid(),
            dtype=np.float64,
        ),
        volume=np.asarray(input_data.volume, dtype=np.float64),
        static_energy=static_energy.copy(),
        metadata=_build_result_metadata(
            input_data,
            options,
            backend_name=backend_name,
        ),
    )

    _notify_result(
        "backend",
        {
            "kind": "thermodynamic_backend",
            "backend": backend_name,
        },
        result_callback,
    )

    if not options.calculate_thermodynamics:
        return result

    total = len(THERMODYNAMIC_STEPS)
    current = 0

    _notify_step("zero_point_energy", "started", step_callback)
    uzp_raw = _time_call(
        "zero_point_energy",
        backend_name,
        record_timing,
        zero_point_energy,
        np.zeros(1, dtype=np.float64),
        frequencies,
        weights,
    )
    result.zero_point_energy = np.asarray(
        convert_energy(uzp_raw, "kjmol", options.energy_unit),
        dtype=np.float64,
    )
    _notify_result(
        "zero_point_energy",
        {
            "kind": "thermodynamic_property",
            "property": "zero_point_energy",
            "symbol": "Uzp",
            "unit": options.energy_unit,
            "value": result.zero_point_energy.copy(),
        },
        result_callback,
    )
    current = _notify_progress(
        "zero_point_energy", current, total, progress_callback, step_callback
    )

    _notify_step("thermal_energy", "started", step_callback)
    uth_raw = _time_call(
        "thermal_energy",
        backend_name,
        record_timing,
        thermal_energy,
        temperature,
        frequencies,
        weights,
    )
    result.thermal_energy = np.asarray(
        convert_energy(uth_raw, "kjmol", options.energy_unit),
        dtype=np.float64,
    )
    _notify_result(
        "thermal_energy",
        {
            "kind": "thermodynamic_property",
            "property": "thermal_energy",
            "symbol": "Uth",
            "unit": options.energy_unit,
            "value": result.thermal_energy.copy(),
        },
        result_callback,
    )
    current = _notify_progress(
        "thermal_energy", current, total, progress_callback, step_callback
    )

    _notify_step("entropy", "started", step_callback)
    entropy_raw = _time_call(
        "entropy",
        backend_name,
        record_timing,
        entropy,
        temperature,
        frequencies,
        weights,
    )
    result.entropy = np.asarray(
        convert_energy_per_temperature(
            entropy_raw,
            "J mol^-1 K^-1",
            f"{options.energy_unit} cell^-1 K^-1",
        ),
        dtype=np.float64,
    )
    _notify_result(
        "entropy",
        {
            "kind": "thermodynamic_property",
            "property": "entropy",
            "symbol": "S",
            "unit": f"{options.energy_unit} cell^-1 K^-1",
            "value": result.entropy.copy(),
        },
        result_callback,
    )
    current = _notify_progress(
        "entropy", current, total, progress_callback, step_callback
    )

    _notify_step("isochoric_heat_capacity", "started", step_callback)
    cv_raw = _time_call(
        "isochoric_heat_capacity",
        backend_name,
        record_timing,
        isochoric_heat_capacity,
        temperature,
        frequencies,
        weights,
    )
    result.isochoric_heat_capacity = np.asarray(
        convert_energy_per_temperature(
            cv_raw,
            "J mol^-1 K^-1",
            f"{options.energy_unit} cell^-1 K^-1",
        ),
        dtype=np.float64,
    )
    _notify_result(
        "isochoric_heat_capacity",
        {
            "kind": "thermodynamic_property",
            "property": "isochoric_heat_capacity",
            "symbol": "Cv",
            "unit": f"{options.energy_unit} cell^-1 K^-1",
            "value": result.isochoric_heat_capacity.copy(),
        },
        result_callback,
    )
    current = _notify_progress(
        "isochoric_heat_capacity",
        current,
        total,
        progress_callback,
        step_callback,
    )

    _notify_step("vibrational_free_energy", "started", step_callback)
    fvib_raw = _time_call(
        "vibrational_free_energy",
        backend_name,
        record_timing,
        vibrational_free_energy,
        temperature,
        frequencies,
        weights,
    )
    result.vibrational_free_energy = np.asarray(
        convert_energy(fvib_raw, "kjmol", options.energy_unit),
        dtype=np.float64,
    )
    _notify_result(
        "vibrational_free_energy",
        {
            "kind": "thermodynamic_property",
            "property": "vibrational_free_energy",
            "symbol": "Fvib",
            "unit": options.energy_unit,
            "value": result.vibrational_free_energy.copy(),
        },
        result_callback,
    )
    current = _notify_progress(
        "vibrational_free_energy",
        current,
        total,
        progress_callback,
        step_callback,
    )

    _notify_step("total_energies", "started", step_callback)
    start = perf_counter()
    result.internal_energy = internal_energy(
        static_energy,
        result.zero_point_energy[0],
        result.thermal_energy,
    )
    result.free_energy = free_energy(
        static_energy,
        result.vibrational_free_energy,
    )
    record_timing("total_energies", perf_counter() - start, backend_name)
    _notify_result(
        "total_energies",
        {
            "kind": "thermodynamic_property",
            "property": "total_energies",
            "unit": options.energy_unit,
            "values": {
                "Utot": result.internal_energy.copy(),
                "F": result.free_energy.copy(),
            },
        },
        result_callback,
    )
    _notify_progress("total_energies", current, total, progress_callback, step_callback)
    result.metadata["timing"] = _timing_records_to_dict(timing_records)

    return result


def _build_result_metadata(
    input_data: HAInput,
    options: HAOptions,
    *,
    backend_name: str = "numpy",
) -> dict[str, Any]:
    """
    Build metadata describing a HA result.

    Parameters
    ----------
    input_data : HAInput
        Input data used in the calculation.
    options : HAOptions
        Calculation options.
    backend_name : str, optional
        Name of the thermodynamics backend selected for the workflow.

    Returns
    -------
    dict
        Metadata dictionary suitable for later HDF5 export.
    """
    return {
        "module": "ha",
        "method": "harmonic",
        "jobname": input_data.jobname,
        "source": str(input_data.source) if input_data.source is not None else None,
        "units": {
            "energy": options.energy_unit,
            "entropy": f"{options.energy_unit} cell^-1 K^-1",
            "heat_capacity": f"{options.energy_unit} cell^-1 K^-1",
            "volume": f"{options.volume_unit}^3",
            "frequency": options.frequency_unit,
            "temperature": options.temperature_unit,
        },
        "input": {
            "natoms": int(input_data.natoms),
            "formula_units": int(input_data.formula_units),
            "natoms_per_formula_unit": float(input_data.natoms_per_formula_unit),
            "qpoints": int(input_data.qpoints),
            "nvol": int(input_data.nvol),
            "nmodes": int(input_data.nmodes),
            "total_q_points": float(input_data.total_q_points),
        },
        "thermodynamic_unit_convention": "native_energy_per_cell_per_kelvin",
        "normalization": {
            "native_basis": "cell",
            "formula_units_per_cell": int(input_data.formula_units),
            "natoms_per_cell": int(input_data.natoms),
            "natoms_per_formula_unit": float(input_data.natoms_per_formula_unit),
            "molar_basis": "formula_unit",
        },
        "options": {
            "temperature_min": float(options.temperature_min),
            "temperature_max": float(options.temperature_max),
            "temperature_step": float(options.temperature_step),
            "calculate_thermodynamics": bool(options.calculate_thermodynamics),
        },
        "backend": {
            "thermodynamics": backend_name,
        },
    }


def _timing_records_to_dict(records: list[dict[str, Any]]) -> dict[str, Any]:
    """
    Convert timing records to an HDF5-friendly dictionary.

    Parameters
    ----------
    records : list of dict
        Timing records collected during the HA workflow.

    Returns
    -------
    dict
        Dictionary containing labels, elapsed times, and backend names.
    """
    return {
        "label": [str(record["label"]) for record in records],
        "elapsed_seconds": np.asarray(
            [float(record["elapsed_seconds"]) for record in records],
            dtype=np.float64,
        ),
        "backend": [str(record["backend"]) for record in records],
    }


def _time_call(
    label: str,
    backend_name: str,
    timing_callback: TimingCallback | None,
    function: Callable[..., np.ndarray],
    *args: Any,
) -> np.ndarray:
    """
    Execute a numerical function and report its elapsed time.

    Parameters
    ----------
    label : str
        Name of the numerical operation.
    backend_name : str
        Name of the backend used for the operation.
    timing_callback : callable or None
        Optional callback receiving label, elapsed seconds, and backend name.
    function : callable
        Numerical function to execute.
    *args : tuple
        Positional arguments passed to ``function``.

    Returns
    -------
    ndarray
        Numerical result returned by ``function``.
    """
    start = perf_counter()
    value = function(*args)
    elapsed = perf_counter() - start
    if timing_callback is not None:
        timing_callback(label, elapsed, backend_name)
    return value


def _notify_result(
    label: str,
    payload: dict[str, Any],
    result_callback: ResultCallback | None,
) -> None:
    """
    Notify availability of a structured intermediate result.

    Parameters
    ----------
    label : str
        Result label.
    payload : dict
        Structured numerical payload associated with the result.
    result_callback : callable or None
        Optional callback receiving ``label`` and ``payload``.
    """
    if result_callback is not None:
        result_callback(label, payload)


def _notify_step(
    label: str,
    status: str,
    step_callback: StepCallback | None,
) -> None:
    """
    Notify the beginning or completion of a workflow step.

    Parameters
    ----------
    label : str
        Step label.
    status : str
        Step status, for example ``"started"`` or ``"finished"``.
    step_callback : callable or None
        Optional callback receiving ``label`` and ``status``.
    """
    if step_callback is not None:
        step_callback(label, status)


def _notify_progress(
    label: str,
    current: int,
    total: int,
    progress_callback: ProgressCallback | None,
    step_callback: StepCallback | None,
) -> int:
    """
    Notify completion of a workflow step and return the updated counter.

    Parameters
    ----------
    label : str
        Completed step label.
    current : int
        Current completed step count before notification.
    total : int
        Total number of workflow steps.
    progress_callback : callable or None
        Optional callback receiving ``label``, ``current`` and ``total``.
    step_callback : callable or None
        Optional callback receiving ``label`` and ``status``.

    Returns
    -------
    int
        Updated completed step count.
    """
    updated = current + 1
    _notify_step(label, "finished", step_callback)
    if progress_callback is not None:
        progress_callback(label, updated, total)
    return updated
