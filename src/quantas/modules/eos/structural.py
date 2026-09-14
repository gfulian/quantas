# -*- coding: utf-8 -*-

"""Shared structural response derived from theoretical energy-volume fits.

The module connects the public EnergyEOS workflow to the frontend-neutral
:class:`quantas.core.geometry.StructuralPathModel`.  No alternative structural
interpolator is implemented here: EOS and QHA deliberately share the same
volume-constrained lattice path.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
from numpy.typing import NDArray

from quantas.core.geometry import StructuralPathModel, lattice_from_parameters
from quantas.core.math.fitting import FitResult
from quantas.models.structures import LatticeVolumeSeries, SymmetryMetadata

from .dataset_models import (
    CrystalReference,
    EOSDataset,
    parse_crystal_reference,
    parse_eos_crystal_system,
)
from .domains.ev import EnergyEOSFitModel

FloatArray = NDArray[np.float64]

_CELL_PARAMETER_COLUMNS = ("a", "b", "c", "alpha", "beta", "gamma")
_AXIS_NAMES = ("a", "b", "c")
_PRIMARY_ORDER = (
    "a0",
    "b0",
    "c0",
    "eta_a",
    "eta_b",
    "eta_c",
    "M_a",
    "M_b",
    "M_c",
)


@dataclass(slots=True)
class EnergyStructuralResponse:
    """Structural properties derived from one successful EnergyEOS fit.

    Parameters
    ----------
    predictions : dict
        Structural quantities evaluated on the complete observed volume grid.
    derived : dict
        Equilibrium cell lengths, logarithmic volume responses, axial moduli,
        and one-sigma uncertainties where available.
    metadata : dict
        Structural-path, covariance, crystal-reference, and symmetry metadata.
    selected_pressure : ndarray
        Pressures derived from the EnergyEOS at selected observations.
    pressure_covariance : ndarray or None
        Full covariance among ``selected_pressure`` values due to EnergyEOS
        parameter covariance.
    independent_axes : tuple of str
        Axes that are independent in the declared crystal reference.
    """

    predictions: dict[str, FloatArray]
    derived: dict[str, float]
    metadata: dict[str, Any]
    selected_pressure: FloatArray
    pressure_covariance: FloatArray | None
    independent_axes: tuple[str, ...]


def build_lattice_volume_series(
    dataset: EOSDataset,
    *,
    mask: np.ndarray | None = None,
) -> LatticeVolumeSeries:
    """Build the shared lattice-only structural path from one EOS dataset.

    Parameters
    ----------
    dataset : EOSDataset
        EOS dataset containing volume and all six lattice parameters.
    mask : ndarray or None, optional
        Observation selection used to construct the structural path.

    Returns
    -------
    LatticeVolumeSeries
        Volume-aligned lattice matrices suitable for
        :class:`~quantas.core.geometry.StructuralPathModel`.

    Raises
    ------
    ValueError
        If the required structural columns are missing or inconsistent with
        the declared volumes.
    """
    missing = [name for name in ("volume", *_CELL_PARAMETER_COLUMNS) if not dataset.has(name)]
    if missing:
        raise ValueError(
            "structural EnergyEOS analysis requires columns: " + ", ".join(missing)
        )
    selection = dataset.selection_mask(mask)
    indices = np.flatnonzero(selection)
    if indices.size < 2:
        raise ValueError("structural EnergyEOS analysis requires at least two states")
    volumes = np.asarray(dataset.column("volume")[indices], dtype=np.float64)
    parameters = np.column_stack(
        [
            np.asarray(dataset.column(name)[indices], dtype=np.float64)
            for name in _CELL_PARAMETER_COLUMNS
        ]
    )
    lattices = np.asarray(
        [lattice_from_parameters(*row) for row in parameters],
        dtype=np.float64,
    )
    lattice_volumes = np.abs(np.linalg.det(lattices))
    if not np.allclose(
        lattice_volumes,
        volumes,
        rtol=5.0e-6,
        atol=5.0e-5,
    ):
        maximum = float(np.max(np.abs(lattice_volumes - volumes)))
        raise ValueError(
            "EOS lattice parameters are inconsistent with the declared cell "
            f"volumes; maximum |det(A)-V|={maximum:.6g} angstrom^3"
        )

    symmetry = None
    number_value = dataset.metadata.get("space_group_number")
    if number_value is not None:
        symmetry = SymmetryMetadata(
            space_group_number=int(number_value),
            international_symbol=str(dataset.metadata.get("space_group_symbol", "")),
        )
    if dataset.has("energy"):
        energies = np.asarray(dataset.column("energy")[indices], dtype=np.float64)
        reference_index = int(np.argmin(energies))
    else:
        reference_index = int(indices.size // 2)
    return LatticeVolumeSeries(
        lattices=lattices,
        volumes=volumes,
        symmetry=symmetry,
        orientation="eos_input_reference",
        reference_index=reference_index,
        metadata={
            "source": "eos_dataset",
            "crystal_reference": dataset.metadata.get("crystal_reference"),
        },
    )


def analyze_energy_structural_response(
    dataset: EOSDataset,
    adapter: EnergyEOSFitModel,
    fit: FitResult,
    *,
    mask: np.ndarray | None = None,
    structural_degree: int = 3,
    relative_step: float = 1.0e-5,
) -> EnergyStructuralResponse:
    r"""Derive equilibrium axes and axial moduli from an EnergyEOS fit.

    The primary axial modulus follows

    .. math::

        M_i(V) = K(V) / \eta_i(V),\qquad
        \eta_i = \partial\ln l_i / \partial\ln V.

    EnergyEOS parameter covariance and structural-path covariance are combined
    by first-order delta propagation at the equilibrium volume.

    Parameters
    ----------
    dataset : EOSDataset
        Dataset used for the E-V fit.
    adapter : EnergyEOSFitModel
        Fitted EnergyEOS model adapter.
    fit : FitResult
        Successful EnergyEOS fit result.
    mask : ndarray or None, optional
        Observation selection used by the fit.
    structural_degree : int, optional
        Polynomial degree for the shared structural-path model.
    relative_step : float, optional
        Relative finite-difference step used for EnergyEOS covariance
        propagation.

    Returns
    -------
    EnergyStructuralResponse
        Primary structural response and covariance information.

    Raises
    ------
    ValueError
        If the fit is unsuccessful, lattice data are incomplete, or an axial
        logarithmic response vanishes at equilibrium.
    """
    if not fit.success or fit.parameters is None:
        raise ValueError("structural response requires a successful EnergyEOS fit")
    series = build_lattice_volume_series(dataset, mask=mask)
    path = StructuralPathModel(series, degree=structural_degree, basis="sampled")
    parameters = np.asarray(fit.parameters, dtype=np.float64)
    parameter_names = tuple(fit.parameter_names)
    try:
        v0_index = parameter_names.index("V0")
    except ValueError as exc:
        raise ValueError("EnergyEOS fit does not report V0") from exc
    v0 = float(parameters[v0_index])

    equilibrium = path.log_volume_response(v0, include_fit_uncertainty=True)
    cell = np.asarray(equilibrium.lattice_parameters, dtype=np.float64).reshape(6)
    eta = np.asarray(equilibrium.logarithmic_length_response, dtype=np.float64).reshape(3)
    if np.any(np.isclose(eta, 0.0, rtol=0.0, atol=1.0e-12)):
        raise ValueError("axial logarithmic volume response is zero at equilibrium")
    k0 = float(adapter.bulk_modulus(np.asarray([v0]), parameters)[0])
    modulus = k0 / eta
    primary = np.concatenate((cell[:3], eta, modulus))

    energy_covariance = _propagate_covariance(
        lambda values: _equilibrium_primary_vector(path, adapter, values, parameter_names),
        parameters,
        fit.covariance,
        relative_step=relative_step,
    )
    structural_covariance = _equilibrium_structural_covariance(
        equilibrium.covariance,
        k0=k0,
        eta=eta,
    )
    total_covariance = _sum_covariances(energy_covariance, structural_covariance)

    derived: dict[str, float] = {
        name: float(value) for name, value in zip(_PRIMARY_ORDER, primary, strict=True)
    }
    if total_covariance is not None:
        errors = np.sqrt(np.clip(np.diag(total_covariance), 0.0, None))
        for name, error in zip(_PRIMARY_ORDER, errors, strict=True):
            derived[f"sigma_{name}"] = float(error)

    all_volumes = np.asarray(dataset.column("volume"), dtype=np.float64)
    sampled = path.log_volume_response(all_volumes, include_fit_uncertainty=False)
    sampled_cell = np.asarray(sampled.lattice_parameters, dtype=np.float64)
    sampled_eta = np.asarray(sampled.logarithmic_length_response, dtype=np.float64)
    sampled_bulk = np.asarray(adapter.bulk_modulus(all_volumes, parameters), dtype=np.float64)
    sampled_modulus = sampled_bulk[:, None] / sampled_eta
    predictions: dict[str, FloatArray] = {}
    for index, axis in enumerate(_AXIS_NAMES):
        predictions[f"structural_{axis}"] = np.asarray(sampled_cell[:, index], dtype=np.float64)
        predictions[f"eta_{axis}"] = np.asarray(sampled_eta[:, index], dtype=np.float64)
        predictions[f"M_{axis}"] = np.asarray(sampled_modulus[:, index], dtype=np.float64)

    selection = dataset.selection_mask(mask)
    selected_volumes = np.asarray(dataset.column("volume")[selection], dtype=np.float64)
    selected_pressure = np.asarray(adapter.pressure(selected_volumes, parameters), dtype=np.float64)
    pressure_covariance = _propagate_covariance(
        lambda values: np.asarray(adapter.pressure(selected_volumes, values), dtype=np.float64),
        parameters,
        fit.covariance,
        relative_step=relative_step,
    )

    independent_axes = _independent_axes(dataset)
    metadata: dict[str, Any] = {
        "available": True,
        "method": "energy_eos_plus_shared_structural_path",
        "relationship": "M_i=K/(dln(l_i)/dln(V))",
        "structural_path": equilibrium.metadata,
        "primary_order": list(_PRIMARY_ORDER),
        "independent_axes": list(independent_axes),
        "crystal_reference": dataset.metadata.get("crystal_reference", "unspecified"),
        "crystal_system": dataset.metadata.get("crystal_system"),
        "space_group_number": dataset.metadata.get("space_group_number"),
        "space_group_symbol": dataset.metadata.get("space_group_symbol"),
        "structural_degree": int(path.degree),
        "uncertainty_method": "first_order_delta_method",
        "uncertainty_sources": {
            "energy_eos_parameter_covariance": energy_covariance is not None,
            "structural_path_fit_covariance": structural_covariance is not None,
            "cross_covariance_energy_structural": False,
        },
        "uncertainty_assumption": (
            "EnergyEOS and structural-path fit covariances are treated as independent."
        ),
    }
    if total_covariance is not None:
        metadata["primary_covariance"] = total_covariance
    if pressure_covariance is not None:
        metadata["derived_pressure_covariance"] = pressure_covariance
        metadata["derived_pressure_uncertainty"] = np.sqrt(
            np.clip(np.diag(pressure_covariance), 0.0, None)
        )
    return EnergyStructuralResponse(
        predictions=predictions,
        derived=derived,
        metadata=metadata,
        selected_pressure=selected_pressure,
        pressure_covariance=pressure_covariance,
        independent_axes=independent_axes,
    )


def _equilibrium_primary_vector(
    path: StructuralPathModel,
    adapter: EnergyEOSFitModel,
    parameters: FloatArray,
    parameter_names: tuple[str, ...],
) -> FloatArray:
    """Return ``a,b,c,eta_a,eta_b,eta_c,M_a,M_b,M_c`` for parameters."""
    v0 = float(parameters[parameter_names.index("V0")])
    state = path.log_volume_response(v0, include_fit_uncertainty=False)
    cell = np.asarray(state.lattice_parameters, dtype=np.float64).reshape(6)
    eta = np.asarray(state.logarithmic_length_response, dtype=np.float64).reshape(3)
    bulk = float(adapter.bulk_modulus(np.asarray([v0]), parameters)[0])
    return np.concatenate((cell[:3], eta, bulk / eta))


def _equilibrium_structural_covariance(
    covariance: FloatArray | None,
    *,
    k0: float,
    eta: FloatArray,
) -> FloatArray | None:
    """Map structural covariance to the nine primary EOS response quantities."""
    if covariance is None:
        return None
    source = np.asarray(covariance, dtype=np.float64).reshape(6, 6)
    if not np.all(np.isfinite(source)):
        return None
    jacobian = np.zeros((9, 6), dtype=np.float64)
    jacobian[:6, :6] = np.eye(6, dtype=np.float64)
    for index in range(3):
        jacobian[6 + index, 3 + index] = -k0 / float(eta[index]) ** 2
    return jacobian @ source @ jacobian.T


def _propagate_covariance(
    evaluator: Callable[[FloatArray], FloatArray],
    parameters: FloatArray,
    covariance: FloatArray | None,
    *,
    relative_step: float,
) -> FloatArray | None:
    """Propagate one parameter covariance through a vector-valued evaluator."""
    if covariance is None:
        return None
    if not np.isfinite(relative_step) or relative_step <= 0.0:
        raise ValueError("relative_step must be finite and positive")
    values = np.asarray(parameters, dtype=np.float64)
    source = np.asarray(covariance, dtype=np.float64)
    if source.shape != (values.size, values.size):
        raise ValueError("parameter covariance has incompatible shape")
    base = np.asarray(evaluator(values), dtype=np.float64).reshape(-1)
    jacobian = np.zeros((base.size, values.size), dtype=np.float64)
    for index in range(values.size):
        if source[index, index] == 0.0:
            continue
        scale = max(abs(float(values[index])), 1.0)
        step = max(
            relative_step * scale,
            np.sqrt(np.finfo(np.float64).eps) * scale,
        )
        plus = values.copy()
        minus = values.copy()
        plus[index] += step
        minus[index] -= step
        plus_value: FloatArray | None = None
        minus_value: FloatArray | None = None
        try:
            plus_value = np.asarray(evaluator(plus), dtype=np.float64).reshape(-1)
        except (ValueError, FloatingPointError, OverflowError):
            pass
        try:
            minus_value = np.asarray(evaluator(minus), dtype=np.float64).reshape(-1)
        except (ValueError, FloatingPointError, OverflowError):
            pass
        if plus_value is not None and minus_value is not None:
            jacobian[:, index] = (plus_value - minus_value) / (2.0 * step)
        elif plus_value is not None:
            jacobian[:, index] = (plus_value - base) / step
        elif minus_value is not None:
            jacobian[:, index] = (base - minus_value) / step
    propagated = jacobian @ source @ jacobian.T
    return 0.5 * (propagated + propagated.T)


def _sum_covariances(
    first: FloatArray | None,
    second: FloatArray | None,
) -> FloatArray | None:
    """Return the sum of available covariance contributions."""
    if first is None:
        return None if second is None else np.asarray(second, dtype=np.float64)
    if second is None:
        return np.asarray(first, dtype=np.float64)
    return np.asarray(first, dtype=np.float64) + np.asarray(second, dtype=np.float64)


def _independent_axes(dataset: EOSDataset) -> tuple[str, ...]:
    """Return independent axes for the dataset's declared cell reference."""
    reference_value = dataset.metadata.get("crystal_reference")
    if reference_value is not None:
        reference = parse_crystal_reference(str(reference_value))
        if reference is CrystalReference.PRIMITIVE:
            return _AXIS_NAMES
    system_value = dataset.metadata.get("crystal_system")
    if system_value is None:
        return _AXIS_NAMES
    return parse_eos_crystal_system(str(system_value)).independent_axes


__all__ = [
    "EnergyStructuralResponse",
    "analyze_energy_structural_response",
    "build_lattice_volume_series",
]
