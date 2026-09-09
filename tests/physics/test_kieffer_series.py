"""End-to-end elastic-state to Kieffer-cutoff composition."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.kieffer import build_kieffer_volume_series
from quantas.models import (
    ElasticState,
    ElasticStateSeries,
    ElasticTensorKind,
    PressureSource,
    PrestressProvenance,
)


def _isotropic_stiffness(lame: float, shear: float) -> np.ndarray:
    """Return an isotropic stiffness matrix in GPa."""
    matrix = np.zeros((6, 6), dtype=np.float64)
    matrix[:3, :3] = lame
    np.fill_diagonal(matrix[:3, :3], lame + 2.0 * shear)
    np.fill_diagonal(matrix[3:, 3:], shear)
    return matrix


def _series(kind: ElasticTensorKind) -> ElasticStateSeries:
    """Return two isotropic elastic states with explicit provenance."""
    states = tuple(
        ElasticState(
            volume=volume,
            density=density,
            stiffness=_isotropic_stiffness(70.0, 50.0),
            prestress=PrestressProvenance(
                tensor_kind=kind,
                pressure_gpa=pressure,
                pressure_source=PressureSource.APPLIED_PRESTRESS,
            ),
            source=f"state-{index}.out",
        )
        for index, (volume, density, pressure) in enumerate(
            ((100.0, 3000.0, 2.0), (110.0, 2900.0, 1.0))
        )
    )
    return ElasticStateSeries(states=states, reference_index=0)


@pytest.mark.physics
@pytest.mark.seismic
def test_incremental_series_builds_direct_cutoffs() -> None:
    progress: list[tuple[int, int]] = []
    result = build_kieffer_volume_series(
        _series(ElasticTensorKind.WALLACE_HYDROSTATIC),
        mu_order=4,
        phi_order=8,
        refinement_factor=2,
        progress_callback=lambda current, total: progress.append((current, total)),
    )

    assert len(result.states) == 2
    assert result.frequencies_hz.shape == (3, 2)
    assert result.effective_velocities_km_s.shape == (3, 2)
    assert np.all(result.frequencies_hz > 0.0)
    assert result.states[0].source_elastic_indices == (0,)
    assert result.states[0].metadata["tensor_kind"] == "wallace_hydrostatic"
    assert result.states[0].metadata["quadrature"]["direction_count"] == 128
    assert progress == [(1, 2), (2, 2)]


@pytest.mark.physics
@pytest.mark.seismic
def test_density_and_volume_change_cutoffs_consistently() -> None:
    result = build_kieffer_volume_series(
        _series(ElasticTensorKind.WALLACE_HYDROSTATIC),
        mu_order=4,
        phi_order=8,
        refinement_factor=1,
    )
    velocity_ratio = np.sqrt(3000.0 / 2900.0)
    np.testing.assert_allclose(
        result.effective_velocities_km_s[:, 1] / result.effective_velocities_km_s[:, 0],
        velocity_ratio,
    )
    expected_cutoff_ratio = velocity_ratio * (100.0 / 110.0) ** (1.0 / 3.0)
    np.testing.assert_allclose(
        result.frequencies_hz[:, 1] / result.frequencies_hz[:, 0],
        expected_cutoff_ratio,
    )


@pytest.mark.physics
@pytest.mark.seismic
def test_raw_series_is_rejected_before_acoustic_work() -> None:
    calls: list[tuple[int, int]] = []
    with pytest.raises(ValueError, match="incremental"):
        build_kieffer_volume_series(
            _series(ElasticTensorKind.RAW_ENERGY_STRAIN),
            progress_callback=lambda current, total: calls.append((current, total)),
        )
    assert calls == []
