# -*- coding: utf-8 -*-

"""Hydrostatic finite-stress operators for backend-neutral elastic states.

The correction implemented here follows Quantas' Eulerian finite-strain
convention and is retained for the internal thermodynamic formulation.  It is
not a universal external-code ingestion rule: interfaces must apply the
finite-pressure transformation appropriate to the tensor definition emitted
by their backend.  In particular, raw CRYSTAL energy--strain tensors are
converted in :mod:`quantas.interfaces.crystal`.  Positive pressure denotes
compression throughout.
"""

from __future__ import annotations

from typing import Any, Mapping

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


def assign_hydrostatic_pressures(
    series: ElasticStateSeries,
    pressures_gpa: ArrayLike,
    *,
    pressure_source: PressureSource | str,
    assignment_method: str,
    metadata: Mapping[str, Any] | None = None,
) -> ElasticStateSeries:
    """Attach externally derived pressures to an explicitly raw series.

    This operation does not alter any stiffness coefficient.  It prepares raw
    energy--strain tensors for a subsequent, separately auditable finite-
    prestress conversion.  Existing pressure or correction provenance cannot
    be replaced.

    Parameters
    ----------
    series : ElasticStateSeries
        Increasing series of raw tensors without assigned pressures.
    pressures_gpa : array_like
        One finite hydrostatic pressure per state, positive in compression.
    pressure_source : PressureSource or str
        Origin of the supplied pressure values.
    assignment_method : str
        Stable description of the calculation that produced the values.
    metadata : mapping, optional
        Additional series-level provenance, for example EOS fit diagnostics.

    Returns
    -------
    ElasticStateSeries
        Independent raw series with complete pressure provenance.

    Raises
    ------
    TypeError
        If ``series`` has an unsupported type.
    ValueError
        If values, sources, or existing tensor provenance are incompatible.
    """
    if not isinstance(series, ElasticStateSeries):
        raise TypeError("series must be an ElasticStateSeries")
    pressures = np.asarray(pressures_gpa, dtype=np.float64)
    if pressures.shape != (series.nstates,) or not np.all(np.isfinite(pressures)):
        raise ValueError("pressures_gpa must contain one finite value per state")
    source = PressureSource(pressure_source)
    if source in {PressureSource.UNAVAILABLE, PressureSource.APPLIED_PRESTRESS}:
        raise ValueError("pressure_source must identify an external pressure value")
    method = str(assignment_method).strip()
    if not method:
        raise ValueError("assignment_method must be non-empty")

    states: list[ElasticState] = []
    for index, (state, pressure) in enumerate(
        zip(series.states, pressures, strict=True)
    ):
        tensor_kind = ElasticTensorKind(state.prestress.tensor_kind)
        state_source = PressureSource(state.prestress.pressure_source)
        if tensor_kind is not ElasticTensorKind.RAW_ENERGY_STRAIN:
            raise ValueError(
                f"elastic state {index}: pressure assignment requires an "
                "explicitly raw energy-strain tensor"
            )
        if state.prestress.pressure_gpa is not None or (
            state_source is not PressureSource.UNAVAILABLE
        ):
            raise ValueError(
                f"elastic state {index}: existing pressure provenance cannot "
                "be replaced"
            )
        state_metadata = dict(state.metadata)
        state_metadata["pressure_assignment"] = {
            "method": method,
            "pressure_gpa": float(pressure),
            "pressure_source": source.value,
        }
        states.append(
            ElasticState(
                volume=state.volume,
                density=state.density,
                stiffness=state.stiffness,
                prestress=PrestressProvenance(
                    tensor_kind=tensor_kind,
                    pressure_gpa=float(pressure),
                    pressure_source=source,
                ),
                energy=state.energy,
                energy_unit=state.energy_unit,
                lattice=state.lattice,
                source=state.source,
                metadata=state_metadata,
            )
        )

    series_metadata = dict(series.metadata)
    assignment: dict[str, Any] = {
        "method": method,
        "pressure_source": source.value,
        "pressure_unit": "GPa",
        "state_count": series.nstates,
    }
    assignment.update(dict(metadata or {}))
    series_metadata["pressure_assignment"] = assignment
    return ElasticStateSeries(
        states=tuple(states),
        reference_index=series.reference_index,
        orientation=series.orientation,
        metadata=series_metadata,
    )


EULERIAN_HYDROSTATIC_PRESTRESS_METHOD = "quantas-eulerian-hydrostatic-incremental"


def eulerian_hydrostatic_incremental_stiffness(
    raw_stiffness: ArrayLike,
    pressure_gpa: float,
) -> NDArray[np.float64]:
    r"""Apply Quantas' Eulerian hydrostatic incremental-stiffness operator.

    The internal finite-strain convention used by Quantas is

    .. math::

       B_{ijkl}=C^{\mathrm{raw}}_{ijkl}-P\,\Delta_{ijkl},

    where ``P`` is positive in compression and ``Delta`` is returned by
    :func:`wallace_hydrostatic_delta_voigt`.

    Parameters
    ----------
    raw_stiffness : array_like
        Finite symmetric ``(6, 6)`` stiffness matrix in GPa defined in the
        Eulerian convention expected by this operator.
    pressure_gpa : float
        Hydrostatic pressure in GPa, positive in compression.

    Returns
    -------
    ndarray
        Symmetric incremental stiffness matrix in GPa.

    Raises
    ------
    ValueError
        If the matrix or pressure is invalid.

    Notes
    -----
    This is an internal Quantas finite-strain operator, not a universal
    external-code ingestion rule.  A backend may define its raw elastic
    derivatives with a different strain measure and therefore require a
    backend-specific finite-pressure transformation.  In particular, raw
    CRYSTAL energy--strain coefficients are converted by
    :func:`quantas.interfaces.crystal.crystal_hydrostatic_stiffness` using the
    Erba/Barron--Klein relation implemented by CRYSTAL itself.
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


def convert_eulerian_hydrostatic_elastic_state(
    state: ElasticState,
    *,
    correction_applied_by: str = "quantas",
) -> ElasticState:
    """Convert one compatible raw state to Eulerian incremental coefficients.

    Parameters
    ----------
    state : ElasticState
        State whose tensor is explicitly marked ``raw_energy_strain`` and
        whose pressure value and source are available.  The tensor definition
        must be compatible with Quantas' Eulerian hydrostatic operator.
    correction_applied_by : str, optional
        Provenance label for the software or caller applying the conversion.

    Returns
    -------
    ElasticState
        Independent state containing hydrostatic incremental coefficients.

    Raises
    ------
    TypeError
        If ``state`` is not an :class:`ElasticState`.
    ValueError
        If the source tensor is not raw, pressure is unavailable, or the
        provenance label is empty.

    Notes
    -----
    External-code adapters must not call this conversion merely because a
    tensor is labelled ``raw_energy_strain``.  They must first establish that
    the backend's raw derivative convention matches this Eulerian operator.
    CRYSTAL raw tensors deliberately use the CRYSTAL-specific conversion in
    :mod:`quantas.interfaces.crystal` instead.
    """
    if not isinstance(state, ElasticState):
        raise TypeError("state must be an ElasticState")
    source_kind = ElasticTensorKind(state.prestress.tensor_kind)
    if source_kind is not ElasticTensorKind.RAW_ENERGY_STRAIN:
        raise ValueError(
            "hydrostatic conversion requires an explicitly raw energy-strain tensor"
        )
    pressure = state.prestress.pressure_gpa
    pressure_source = PressureSource(state.prestress.pressure_source)
    if pressure is None or pressure_source is PressureSource.UNAVAILABLE:
        raise ValueError("hydrostatic conversion requires pressure provenance")
    applied_by = str(correction_applied_by).strip()
    if not applied_by:
        raise ValueError("correction_applied_by must be non-empty")

    method = EULERIAN_HYDROSTATIC_PRESTRESS_METHOD
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
        stiffness=eulerian_hydrostatic_incremental_stiffness(
            state.stiffness,
            pressure,
        ),
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


def convert_eulerian_hydrostatic_elastic_series(
    series: ElasticStateSeries,
    *,
    correction_applied_by: str = "quantas",
) -> ElasticStateSeries:
    """Convert a compatible raw series to Eulerian incremental coefficients.

    Parameters
    ----------
    series : ElasticStateSeries
        Increasing raw elastic-state series with pressure provenance at every
        volume and a tensor convention compatible with Quantas' Eulerian
        hydrostatic operator.
    correction_applied_by : str, optional
        Provenance label recorded in every converted state.

    Returns
    -------
    ElasticStateSeries
        New series containing hydrostatic incremental tensors.

    Raises
    ------
    TypeError
        If ``series`` has an unsupported type.
    ValueError
        If any state cannot be converted exactly once.
    """
    if not isinstance(series, ElasticStateSeries):
        raise TypeError("series must be an ElasticStateSeries")
    corrected_states: list[ElasticState] = []
    for index, state in enumerate(series.states):
        try:
            corrected_states.append(
                convert_eulerian_hydrostatic_elastic_state(
                    state,
                    correction_applied_by=correction_applied_by,
                )
            )
        except ValueError as exc:
            raise ValueError(f"elastic state {index}: {exc}") from exc
    states = tuple(corrected_states)
    metadata = dict(series.metadata)
    metadata["prestress_correction"] = {
        "method": EULERIAN_HYDROSTATIC_PRESTRESS_METHOD,
        "applied_by": str(correction_applied_by).strip(),
        "state_count": len(states),
    }
    return ElasticStateSeries(
        states=states,
        reference_index=series.reference_index,
        orientation=series.orientation,
        metadata=metadata,
    )


def hydrostatic_wallace_stiffness(
    raw_stiffness: ArrayLike,
    pressure_gpa: float,
) -> NDArray[np.float64]:
    """Compatibility alias for the Eulerian hydrostatic stiffness operator.

    Parameters
    ----------
    raw_stiffness : array_like
        Finite symmetric ``(6, 6)`` stiffness matrix in GPa.
    pressure_gpa : float
        Hydrostatic pressure in GPa, positive in compression.

    Returns
    -------
    ndarray
        Symmetric incremental stiffness matrix in GPa.

    Notes
    -----
    New code should use :func:`eulerian_hydrostatic_incremental_stiffness`.
    The historical name is retained for source compatibility only; it must not
    be interpreted as the CRYSTAL/Erba finite-pressure transformation.
    """
    return eulerian_hydrostatic_incremental_stiffness(raw_stiffness, pressure_gpa)


def correct_hydrostatic_elastic_state(
    state: ElasticState,
    *,
    correction_applied_by: str = "quantas",
) -> ElasticState:
    """Compatibility alias for the Eulerian hydrostatic state conversion.

    Parameters
    ----------
    state : ElasticState
        Compatible raw elastic state with pressure provenance.
    correction_applied_by : str, optional
        Provenance label recorded on the converted state.

    Returns
    -------
    ElasticState
        Independent hydrostatic incremental state.

    Notes
    -----
    New code should use :func:`convert_eulerian_hydrostatic_elastic_state`.
    """
    return convert_eulerian_hydrostatic_elastic_state(
        state,
        correction_applied_by=correction_applied_by,
    )


def correct_hydrostatic_elastic_series(
    series: ElasticStateSeries,
    *,
    correction_applied_by: str = "quantas",
) -> ElasticStateSeries:
    """Compatibility alias for the Eulerian hydrostatic series conversion.

    Parameters
    ----------
    series : ElasticStateSeries
        Compatible raw elastic-state series with pressure provenance.
    correction_applied_by : str, optional
        Provenance label recorded on each converted state.

    Returns
    -------
    ElasticStateSeries
        Independent series of hydrostatic incremental tensors.

    Notes
    -----
    New code should use :func:`convert_eulerian_hydrostatic_elastic_series`.
    """
    return convert_eulerian_hydrostatic_elastic_series(
        series,
        correction_applied_by=correction_applied_by,
    )


__all__ = [
    "EULERIAN_HYDROSTATIC_PRESTRESS_METHOD",
    "assign_hydrostatic_pressures",
    "convert_eulerian_hydrostatic_elastic_series",
    "convert_eulerian_hydrostatic_elastic_state",
    "eulerian_hydrostatic_incremental_stiffness",
    "correct_hydrostatic_elastic_series",
    "correct_hydrostatic_elastic_state",
    "hydrostatic_wallace_stiffness",
]
