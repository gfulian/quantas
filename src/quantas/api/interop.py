# -*- coding: utf-8 -*-

"""Explicit transformations between Quantas scientific workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np

from quantas.core.physics.units import convert_volume
from quantas.models import ResultData
from quantas.models.phonons import PhononInputData
from quantas.modules.qha.models import QHAInput, QHAOptions, QHAResult
from quantas.modules.seismic.models import SeismicInput
from quantas.modules.thermoelasticity.models import (
    ThermoelasticContext,
    ThermoelasticInput,
    ThermoelasticResult,
    select_stiffness_tensor,
)


from .common import _public_dir

def load_qha_result(
    source: ResultData | QHAInput | PhononInputData | str | Path,
    *,
    options: QHAOptions | None = None,
) -> ResultData:
    """Load or calculate a validated QHA result envelope.

    Parameters
    ----------
    source : ResultData, QHAInput, PhononInputData, str, or Path
        Existing result, QHA/phonon input, native QHA HDF5 path, or YAML input
        path.
    options : QHAOptions or None, optional
        QHA options used only when a calculation is required.

    Returns
    -------
    ResultData
        Validated QHA result envelope.

    Raises
    ------
    ValueError
        If a file or envelope does not contain a valid QHA result.
    """
    from . import qha

    if isinstance(source, ResultData):
        result = source
    elif isinstance(source, (str, Path)) and Path(source).suffix.lower() in {
        ".h5",
        ".hdf5",
    }:
        result = qha.read_result(source)
    else:
        result = qha.run(source, options=options)
    qha.get_result(result)
    return result


def qha_to_thermoelastic_context(
    input_data: ThermoelasticInput | str | Path,
    qha_source: ResultData | QHAInput | PhononInputData | str | Path,
    *,
    qha_options: QHAOptions | None = None,
) -> ThermoelasticContext:
    """Pair elastic-volume input with QHA data and validate coverage.

    Parameters
    ----------
    input_data : ThermoelasticInput, str, or Path
        Elastic-volume series or thermoelastic YAML path.
    qha_source : ResultData, QHAInput, PhononInputData, str, or Path
        Existing QHA result, native HDF5 result, or input requiring a QHA run.
    qha_options : QHAOptions or None, optional
        Options used when ``qha_source`` must be calculated.

    Returns
    -------
    ThermoelasticContext
        Validated coupling context including extrapolation and completeness
        diagnostics.

    Raises
    ------
    ValueError
        If QHA equilibrium volumes are missing, array shapes are inconsistent,
        or primitive atomic normalization differs.
    """
    from . import thermoelasticity

    normalized = thermoelasticity.normalize_input(input_data)
    result_data = load_qha_result(qha_source, options=qha_options)
    qha_payload = result_data.results.get("qha")
    if not isinstance(qha_payload, QHAResult):
        raise ValueError("the supplied source does not contain a valid QHA result")
    if qha_payload.equilibrium_volume is None:
        raise ValueError(
            "QHA equilibrium volumes are required for quasi-static coupling"
        )
    volume_unit = str(result_data.options.get("volume_unit", "A"))
    volumes = np.asarray(
        convert_volume(qha_payload.equilibrium_volume, volume_unit, "A"),
        dtype=np.float64,
    )
    if volumes.ndim != 2:
        raise ValueError("QHA equilibrium_volume must have shape (nT, nP)")
    lower, upper = normalized.elastic_series.volume_bounds
    extrapolation = np.asarray(
        (~np.isfinite(volumes)) | (volumes < lower) | (volumes > upper),
        dtype=np.bool_,
    )
    missing = tuple(
        name
        for name in (
            "temperature",
            "pressure",
            "volume",
            "static_energy",
            "equilibrium_volume",
            "isochoric_heat_capacity",
            "isothermal_bulk_modulus",
            "bulk_modulus_derivative",
            "thermal_expansion",
            "axial_thermal_expansion",
            "thermal_expansion_tensor",
            "equilibrium_lattice",
            "lattice_parameters",
        )
        if getattr(qha_payload, name) is None
    )
    _validate_atomic_normalization(normalized, result_data)
    return ThermoelasticContext(
        input_data=normalized,
        qha_result_data=result_data,
        qha=qha_payload,
        extrapolation_mask=extrapolation,
        missing_qha_fields=missing,
        metadata={
            "elastic_volume_min": lower,
            "elastic_volume_max": upper,
            "qha_volume_min": float(np.nanmin(volumes)),
            "qha_volume_max": float(np.nanmax(volumes)),
            "extrapolated_points": int(np.count_nonzero(extrapolation)),
            "total_points": int(extrapolation.size),
            "qha_source_mode": (
                "hdf5"
                if isinstance(qha_source, (str, Path))
                and Path(qha_source).suffix.lower() in {".h5", ".hdf5"}
                else "calculated-or-object"
            ),
        },
    )


def thermoelastic_to_seismic(
    source: ResultData | str | Path,
    *,
    pressure: float,
    temperature: float,
    tensor_condition: Literal["isothermal", "adiabatic"] = "adiabatic",
    extrapolation_policy: Literal["fail", "warn", "allow"] = "fail",
) -> SeismicInput:
    """Build a seismic input from one thermoelastic state.

    Parameters
    ----------
    source : ResultData, str, or Path
        Thermoelastic result envelope or native HDF5 path.
    pressure : float
        Requested pressure in GPa.
    temperature : float
        Requested temperature in K.
    tensor_condition : {"isothermal", "adiabatic"}, optional
        Stiffness tensor condition exported to SEISMIC.
    extrapolation_policy : {"fail", "warn", "allow"}, optional
        Behavior outside the calibrated volume/thermodynamic domain.

    Returns
    -------
    SeismicInput
        Density and stiffness tensor for the requested state.

    Raises
    ------
    ValueError
        If the state cannot be reconstructed or density/stiffness is invalid.
    TypeError
        If ``source`` is not a result envelope or path.
    """
    from . import thermoelasticity

    result_data = (
        thermoelasticity.read_result(source)
        if isinstance(source, (str, Path))
        else source
    )
    if not isinstance(result_data, ResultData):
        raise TypeError("source must be a ResultData object or HDF5 path")
    analyzed = thermoelasticity.analyze_grid(
        result_data,
        pressure=np.asarray([float(pressure)], dtype=np.float64),
        temperature=np.asarray([float(temperature)], dtype=np.float64),
        extrapolation_policy=extrapolation_policy,
    )
    payload = analyzed.results.get("thermoelasticity")
    if not isinstance(payload, ThermoelasticResult):
        raise ValueError("source does not contain a thermoelasticity result")
    stiffness, _ = select_stiffness_tensor(payload, tensor_condition)
    density = np.asarray(payload.density, dtype=np.float64)
    if density.shape != (1, 1) or not np.isfinite(density[0, 0]):
        raise ValueError("thermoelastic density is unavailable at the requested state")
    return SeismicInput(
        jobname=(
            f"{payload.jobname} at P={float(pressure):g} GPa, "
            f"T={float(temperature):g} K ({tensor_condition})"
        ),
        stiffness=np.asarray(stiffness[0, 0], dtype=np.float64),
        density=float(density[0, 0]),
        source=source if isinstance(source, (str, Path)) else None,
    )


def _validate_atomic_normalization(
    thermoelastic: ThermoelasticInput,
    qha_result: ResultData,
) -> None:
    """Validate primitive atomic-number ordering when QHA stores it."""
    if qha_result.input_data is None:
        return
    structure = qha_result.input_data.data.get("structure")
    if not isinstance(structure, dict):
        return
    numbers = structure.get("atomic_numbers")
    reference = structure.get("reference")
    if numbers is None and isinstance(reference, dict):
        numbers = reference.get("atomic_numbers")
    if numbers is None:
        return
    qha_numbers = np.asarray(numbers, dtype=np.int64)
    elastic_numbers = thermoelastic.elastic_series.reference_structure.atomic_numbers
    if not np.array_equal(qha_numbers, elastic_numbers):
        raise ValueError(
            "QHA and elastic inputs use different primitive atomic species/order"
        )


def __dir__() -> list[str]:
    """Return only names in the supported public contract."""
    return _public_dir(__all__)


__all__ = [
    "load_qha_result",
    "qha_to_thermoelastic_context",
    "thermoelastic_to_seismic",
]
