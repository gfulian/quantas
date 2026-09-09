# -*- coding: utf-8 -*-

"""Composition of elastic volume states into Kieffer acoustic cutoffs."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import (
    ChristoffelSolver,
    ElasticMedium,
    average_acoustic_phase_velocities,
    kieffer_cutoff_frequencies,
)
from quantas.models.elastic_states import (
    ElasticStateSeries,
    ElasticTensorKind,
    PressureSource,
)
from quantas.models.kieffer import KiefferCutoffState, KiefferVolumeSeries


ProgressCallback = Callable[[int, int], None]


def build_kieffer_volume_series(
    elastic_series: ElasticStateSeries,
    *,
    mu_order: int = 12,
    phi_order: int = 24,
    refinement_factor: int = 2,
    batch_size: int = 512,
    progress_callback: ProgressCallback | None = None,
) -> KiefferVolumeSeries:
    """Build direct Kieffer cutoffs from incremental elastic volume states.

    Parameters
    ----------
    elastic_series : ElasticStateSeries
        Increasing primitive-volume series. Every stiffness tensor must be
        explicitly marked as hydrostatic Wallace or full-stress incremental.
    mu_order, phi_order : int, optional
        Coarse spherical quadrature orders.
    refinement_factor : int, optional
        Integer multiplier applied to both final quadrature orders.
    batch_size : int, optional
        Maximum directions solved in one Christoffel batch.
    progress_callback : callable or None, optional
        Callback receiving the number of completed volume states and total
        number of states.

    Returns
    -------
    KiefferVolumeSeries
        Direct, volume-resolved effective velocities and cutoff frequencies.

    Raises
    ------
    TypeError
        If ``elastic_series`` has an unsupported type.
    ValueError
        If a tensor is not incremental or an acoustic calculation fails.
    """
    if not isinstance(elastic_series, ElasticStateSeries):
        raise TypeError("elastic_series must be an ElasticStateSeries")
    elastic_series.require_incremental()
    states: list[KiefferCutoffState] = []
    for index, elastic_state in enumerate(elastic_series.states):
        medium = ElasticMedium(
            ElasticTensor(np.asarray(elastic_state.stiffness, dtype=np.float64)),
            elastic_state.density,
        )
        average = average_acoustic_phase_velocities(
            ChristoffelSolver(medium),
            mu_order=mu_order,
            phi_order=phi_order,
            refinement_factor=refinement_factor,
            batch_size=batch_size,
        )
        cutoff = kieffer_cutoff_frequencies(
            average.effective_velocities,
            elastic_state.volume,
        )
        metadata: dict[str, Any] = {
            "tensor_kind": ElasticTensorKind(elastic_state.prestress.tensor_kind).value,
            "pressure_gpa": elastic_state.prestress.pressure_gpa,
            "pressure_source": PressureSource(
                elastic_state.prestress.pressure_source
            ).value,
            "quadrature": {
                "mu_order": average.mu_order,
                "phi_order": average.phi_order,
                "direction_count": average.direction_count,
                "relative_errors": average.relative_errors.copy(),
                "degenerate_direction_count": average.degenerate_direction_count,
                "clamped_mode_count": average.clamped_mode_count,
            },
            "angular_frequencies_rad_s": cutoff.angular_frequencies_rad_s.copy(),
            "wavenumbers_cm1": cutoff.wavenumbers_cm1.copy(),
        }
        states.append(
            KiefferCutoffState(
                volume=elastic_state.volume,
                frequencies_hz=cutoff.frequencies_hz,
                effective_velocities_km_s=average.effective_velocities,
                source_elastic_indices=(index,),
                metadata=metadata,
            )
        )
        if progress_callback is not None:
            progress_callback(index + 1, elastic_series.nstates)
    return KiefferVolumeSeries(
        states=tuple(states),
        metadata={
            "source": "elastic_state_series",
            "orientation": elastic_series.orientation,
            "reference_index": elastic_series.reference_index,
        },
    )


__all__ = ["build_kieffer_volume_series"]
