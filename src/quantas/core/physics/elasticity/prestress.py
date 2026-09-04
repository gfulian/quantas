# -*- coding: utf-8 -*-

"""Hydrostatic finite-stress correction of elastic-state tensors.

The functions in this module convert explicitly raw energy--strain stiffness
coefficients into the Wallace incremental coefficients required by acoustic
wave propagation.  Positive pressure denotes compression throughout.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from quantas.models.elastic_states import (
    ElasticState,
    ElasticStateSeries,
    ElasticTensorKind,
    PressureSource,
    PrestressProvenance,
)

from .quasistatic import wallace_hydrostatic_delta_voigt


def hydrostatic_wallace_stiffness(
    raw_stiffness: ArrayLike,
    pressure_gpa: float,
) -> NDArray[np.float64]:
    r"""Return hydrostatic Wallace coefficients from a raw stiffness tensor.

    Quantas uses the Eulerian finite-strain convention

    .. math::

       B_{ijkl}=C^{\mathrm{raw}}_{ijkl}-P\,\Delta_{ijkl},

    where ``P`` is positive in compression and ``Delta`` is returned by
    :func:`wallace_hydrostatic_delta_voigt`.

    Parameters
    ----------
    raw_stiffness : array_like
        Finite symmetric ``(6, 6)`` energy--strain stiffness matrix in GPa.
    pressure_gpa : float
        Hydrostatic pressure in GPa, positive in compression.

    Returns
    -------
    ndarray
        Corrected symmetric stiffness matrix in GPa.

    Raises
    ------
    ValueError
        If the matrix or pressure is invalid.
    """
    stiffness = np.asarray(raw_stiffness, dtype=np.float64)
    pressure = float(pressure_gpa)
    if stiffness.shape != (6, 6) or not np.all(np.isfinite(stiffness)):
        raise ValueError("raw_stiffness must be finite with shape (6, 6)")
    if not np.allclose(stiffness, stiffness.T, rtol=0.0, atol=1.0e-10):
        raise ValueError("raw_stiffness must be symmetric")
    if not np.isfinite(pressure):
        raise ValueError("pressure_gpa must be finite")
    corrected = stiffness - pressure * wallace_hydrostatic_delta_voigt()
    return np.asarray(0.5 * (corrected + corrected.T), dtype=np.float64)


def correct_hydrostatic_elastic_state(
    state: ElasticState,
    *,
    correction_applied_by: str = "quantas",
) -> ElasticState:
    """Return one raw elastic state corrected at its recorded pressure.

    Parameters
    ----------
    state : ElasticState
        State whose tensor is explicitly marked ``raw_energy_strain`` and
        whose pressure value and source are available.
    correction_applied_by : str, optional
        Provenance label for the software or caller applying the correction.

    Returns
    -------
    ElasticState
        Independent state containing Wallace hydrostatic coefficients.

    Raises
    ------
    TypeError
        If ``state`` is not an :class:`ElasticState`.
    ValueError
        If the source tensor is not raw, pressure is unavailable, or the
        provenance label is empty.
    """
    if not isinstance(state, ElasticState):
        raise TypeError("state must be an ElasticState")
    source_kind = ElasticTensorKind(state.prestress.tensor_kind)
    if source_kind is not ElasticTensorKind.RAW_ENERGY_STRAIN:
        raise ValueError(
            "hydrostatic correction requires an explicitly raw energy-strain tensor"
        )
    pressure = state.prestress.pressure_gpa
    pressure_source = PressureSource(state.prestress.pressure_source)
    if pressure is None or pressure_source is PressureSource.UNAVAILABLE:
        raise ValueError("hydrostatic correction requires pressure provenance")
    applied_by = str(correction_applied_by).strip()
    if not applied_by:
        raise ValueError("correction_applied_by must be non-empty")

    method = "barron-klein-wallace-hydrostatic"
    metadata = dict(state.metadata)
    metadata["prestress_correction"] = {
        "method": method,
        "pressure_gpa": pressure,
        "pressure_source": pressure_source.value,
        "applied_by": applied_by,
        "source_tensor_kind": source_kind.value,
        "target_tensor_kind": ElasticTensorKind.WALLACE_HYDROSTATIC.value,
    }
    return ElasticState(
        volume=state.volume,
        density=state.density,
        stiffness=hydrostatic_wallace_stiffness(state.stiffness, pressure),
        prestress=PrestressProvenance(
            tensor_kind=ElasticTensorKind.WALLACE_HYDROSTATIC,
            pressure_gpa=pressure,
            pressure_source=pressure_source,
            correction_method=method,
            correction_applied_by=applied_by,
            source_tensor_kind=source_kind,
        ),
        energy=state.energy,
        energy_unit=state.energy_unit,
        lattice=state.lattice,
        source=state.source,
        metadata=metadata,
    )


def correct_hydrostatic_elastic_series(
    series: ElasticStateSeries,
    *,
    correction_applied_by: str = "quantas",
) -> ElasticStateSeries:
    """Correct every raw state in a hydrostatic elastic volume series.

    Parameters
    ----------
    series : ElasticStateSeries
        Increasing raw elastic-state series with pressure provenance at every
        volume.
    correction_applied_by : str, optional
        Provenance label recorded in every corrected state.

    Returns
    -------
    ElasticStateSeries
        New series containing only Wallace hydrostatic tensors.

    Raises
    ------
    TypeError
        If ``series`` has an unsupported type.
    ValueError
        If any state cannot be corrected exactly once.
    """
    if not isinstance(series, ElasticStateSeries):
        raise TypeError("series must be an ElasticStateSeries")
    corrected_states: list[ElasticState] = []
    for index, state in enumerate(series.states):
        try:
            corrected_states.append(
                correct_hydrostatic_elastic_state(
                    state,
                    correction_applied_by=correction_applied_by,
                )
            )
        except ValueError as exc:
            raise ValueError(f"elastic state {index}: {exc}") from exc
    states = tuple(corrected_states)
    metadata = dict(series.metadata)
    metadata["prestress_correction"] = {
        "method": "barron-klein-wallace-hydrostatic",
        "applied_by": str(correction_applied_by).strip(),
        "state_count": len(states),
    }
    return ElasticStateSeries(
        states=states,
        reference_index=series.reference_index,
        orientation=series.orientation,
        metadata=metadata,
    )


__all__ = [
    "correct_hydrostatic_elastic_series",
    "correct_hydrostatic_elastic_state",
    "hydrostatic_wallace_stiffness",
]
