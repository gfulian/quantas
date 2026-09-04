"""Tests for explicit hydrostatic correction of raw elastic states."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import (
    correct_hydrostatic_elastic_series,
    correct_hydrostatic_elastic_state,
    hydrostatic_wallace_stiffness,
)
from quantas.core.physics.kieffer import build_kieffer_volume_series
from quantas.models import (
    ElasticState,
    ElasticStateSeries,
    ElasticTensorKind,
    PressureSource,
    PrestressProvenance,
)


def _state(
    volume: float,
    pressure: float | None,
    *,
    kind: ElasticTensorKind = ElasticTensorKind.RAW_ENERGY_STRAIN,
) -> ElasticState:
    """Return a simple elastic state with explicit pressure provenance."""
    edge = volume ** (1.0 / 3.0)
    source = (
        PressureSource.UNAVAILABLE if pressure is None else PressureSource.OUTPUT_STRESS
    )
    stiffness = np.diag([200.0, 210.0, 220.0, 70.0, 75.0, 80.0])
    return ElasticState(
        volume=volume,
        density=3200.0,
        stiffness=stiffness,
        prestress=PrestressProvenance(
            tensor_kind=kind,
            pressure_gpa=pressure,
            pressure_source=source,
        ),
        energy=-100.0,
        energy_unit="Ha",
        lattice=np.eye(3) * edge,
        source="elastic.out",
        metadata={"backend": "crystal"},
    )


def test_wallace_correction_uses_positive_compression_convention() -> None:
    """Known normal, coupling, and shear terms follow the documented tensor."""
    raw = np.diag([200.0, 210.0, 220.0, 70.0, 75.0, 80.0])
    corrected = hydrostatic_wallace_stiffness(raw, 2.0)
    assert corrected.dtype == np.float64
    assert corrected[0, 0] == pytest.approx(206.0)
    assert corrected[0, 1] == pytest.approx(2.0)
    assert corrected[3, 3] == pytest.approx(72.0)
    np.testing.assert_allclose(corrected, corrected.T)
    np.testing.assert_allclose(hydrostatic_wallace_stiffness(raw, 0.0), raw)


def test_state_correction_preserves_data_and_records_provenance() -> None:
    """Correction returns an independent state with a complete audit trail."""
    raw = _state(100.0, 2.0)
    corrected = correct_hydrostatic_elastic_state(raw)
    assert corrected is not raw
    assert corrected.prestress.tensor_kind is ElasticTensorKind.WALLACE_HYDROSTATIC
    assert corrected.prestress.source_tensor_kind is ElasticTensorKind.RAW_ENERGY_STRAIN
    assert corrected.prestress.pressure_source is PressureSource.OUTPUT_STRESS
    assert corrected.prestress.correction_method == ("barron-klein-wallace-hydrostatic")
    assert corrected.volume == raw.volume
    assert corrected.energy == raw.energy
    assert corrected.source == raw.source
    assert corrected.metadata["backend"] == "crystal"
    assert corrected.metadata["prestress_correction"]["pressure_gpa"] == 2.0
    np.testing.assert_allclose(raw.stiffness[0, 0], 200.0)


def test_correction_rejects_missing_pressure_and_second_application() -> None:
    """Neither unknown pressure nor an incremental source can be corrected."""
    with pytest.raises(ValueError, match="pressure provenance"):
        correct_hydrostatic_elastic_state(_state(100.0, None))
    with pytest.raises(ValueError, match="explicitly raw"):
        correct_hydrostatic_elastic_state(
            _state(100.0, 2.0, kind=ElasticTensorKind.WALLACE_HYDROSTATIC)
        )

    incomplete = ElasticStateSeries(
        states=(_state(100.0, 2.0), _state(110.0, None)),
        reference_index=0,
    )
    with pytest.raises(ValueError, match="elastic state 1:.*pressure provenance"):
        correct_hydrostatic_elastic_series(incomplete)


def test_series_correction_produces_acoustic_ready_states() -> None:
    """A complete raw volume series becomes uniformly incremental."""
    raw = ElasticStateSeries(
        states=(_state(100.0, 2.0), _state(110.0, -1.0)),
        reference_index=0,
        orientation="crystal-cartesian",
        metadata={"dataset": "synthetic"},
    )
    corrected = correct_hydrostatic_elastic_series(raw)
    corrected.require_incremental()
    assert corrected.reference_index == raw.reference_index
    assert corrected.orientation == raw.orientation
    assert corrected.metadata["dataset"] == "synthetic"
    assert corrected.metadata["prestress_correction"]["state_count"] == 2
    assert all(
        state.prestress.tensor_kind is ElasticTensorKind.WALLACE_HYDROSTATIC
        for state in corrected.states
    )
    cutoffs = build_kieffer_volume_series(
        corrected,
        mu_order=4,
        phi_order=8,
        refinement_factor=1,
    )
    assert len(cutoffs.states) == 2
    assert np.all(cutoffs.frequencies_hz > 0.0)
