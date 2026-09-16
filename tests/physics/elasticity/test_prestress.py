"""Tests for explicit hydrostatic correction of raw elastic states."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.elasticity import (
    EULERIAN_HYDROSTATIC_PRESTRESS_METHOD,
    assign_hydrostatic_pressures,
    convert_eulerian_hydrostatic_elastic_series,
    convert_eulerian_hydrostatic_elastic_state,
    eulerian_hydrostatic_incremental_stiffness,
    correct_hydrostatic_elastic_series,
    correct_hydrostatic_elastic_state,
    hydrostatic_wallace_stiffness,
    resolve_energy_derived_pressures,
)
from quantas.core.physics.kieffer import build_kieffer_volume_series
from quantas.core.physics.units import energy_to_pressure
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


def test_eulerian_operator_uses_positive_compression_convention() -> None:
    """Known normal, coupling, and shear terms follow the Eulerian operator."""
    raw = np.diag([200.0, 210.0, 220.0, 70.0, 75.0, 80.0])
    corrected = eulerian_hydrostatic_incremental_stiffness(raw, 2.0)
    assert corrected.dtype == np.float64
    assert corrected[0, 0] == pytest.approx(206.0)
    assert corrected[0, 1] == pytest.approx(2.0)
    assert corrected[3, 3] == pytest.approx(72.0)
    np.testing.assert_allclose(corrected, corrected.T)
    np.testing.assert_allclose(
        eulerian_hydrostatic_incremental_stiffness(raw, 0.0), raw
    )


def test_state_correction_preserves_data_and_records_provenance() -> None:
    """Correction returns an independent state with a complete audit trail."""
    raw = _state(100.0, 2.0)
    corrected = convert_eulerian_hydrostatic_elastic_state(raw)
    assert corrected is not raw
    assert corrected.prestress.tensor_kind is ElasticTensorKind.WALLACE_HYDROSTATIC
    assert corrected.prestress.source_tensor_kind is ElasticTensorKind.RAW_ENERGY_STRAIN
    assert corrected.prestress.pressure_source is PressureSource.OUTPUT_STRESS
    assert (
        corrected.prestress.correction_method
        == EULERIAN_HYDROSTATIC_PRESTRESS_METHOD
    )
    assert corrected.volume == raw.volume
    assert corrected.energy == raw.energy
    assert corrected.source == raw.source
    assert corrected.metadata["backend"] == "crystal"
    assert corrected.metadata["prestress_correction"]["pressure_gpa"] == 2.0
    np.testing.assert_allclose(raw.stiffness[0, 0], 200.0)


def test_correction_rejects_missing_pressure_and_second_application() -> None:
    """Neither unknown pressure nor an incremental source can be corrected."""
    with pytest.raises(ValueError, match="pressure provenance"):
        convert_eulerian_hydrostatic_elastic_state(_state(100.0, None))
    with pytest.raises(ValueError, match="explicitly raw"):
        convert_eulerian_hydrostatic_elastic_state(
            _state(100.0, 2.0, kind=ElasticTensorKind.WALLACE_HYDROSTATIC)
        )

    incomplete = ElasticStateSeries(
        states=(_state(100.0, 2.0), _state(110.0, None)),
        reference_index=0,
    )
    with pytest.raises(ValueError, match="elastic state 1:.*pressure provenance"):
        convert_eulerian_hydrostatic_elastic_series(incomplete)


def test_series_correction_produces_acoustic_ready_states() -> None:
    """A complete raw volume series becomes uniformly incremental."""
    raw = ElasticStateSeries(
        states=(_state(100.0, 2.0), _state(110.0, -1.0)),
        reference_index=0,
        orientation="crystal-cartesian",
        metadata={"dataset": "synthetic"},
    )
    corrected = convert_eulerian_hydrostatic_elastic_series(raw)
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


def test_external_pressure_assignment_precedes_wallace_correction() -> None:
    """Fitted E(V) pressures remain traceable through the separate correction."""
    raw = ElasticStateSeries(
        states=(_state(100.0, None), _state(110.0, None)),
        reference_index=0,
        metadata={"dataset": "synthetic"},
    )

    assigned = assign_hydrostatic_pressures(
        raw,
        [2.5, -1.0],
        pressure_source=PressureSource.ENERGY_EOS,
        assignment_method="energy_eos",
        metadata={"settings": {"eos": "BM3"}},
    )
    corrected = convert_eulerian_hydrostatic_elastic_series(assigned)

    assert assigned.states[0].stiffness[0, 0] == pytest.approx(200.0)
    assert assigned.states[0].prestress.pressure_gpa == pytest.approx(2.5)
    assert assigned.states[0].prestress.pressure_source is PressureSource.ENERGY_EOS
    assert assigned.metadata["pressure_assignment"]["settings"] == {"eos": "BM3"}
    assert corrected.states[0].stiffness[0, 0] == pytest.approx(207.5)
    assert corrected.states[0].prestress.pressure_source is PressureSource.ENERGY_EOS


def test_external_pressure_assignment_cannot_replace_provenance() -> None:
    """Pressure assignment accepts only raw states with no existing value."""
    series = ElasticStateSeries(states=(_state(100.0, 1.0),), reference_index=0)

    with pytest.raises(ValueError, match="cannot be replaced"):
        assign_hydrostatic_pressures(
            series,
            [2.0],
            pressure_source=PressureSource.ENERGY_POLYNOMIAL,
            assignment_method="energy_polynomial",
        )

def test_compatibility_names_match_eulerian_api() -> None:
    """Historical public names preserve numerics while new code uses Eulerian names."""
    raw_matrix = np.diag([200.0, 210.0, 220.0, 70.0, 75.0, 80.0])
    np.testing.assert_allclose(
        hydrostatic_wallace_stiffness(raw_matrix, 2.0),
        eulerian_hydrostatic_incremental_stiffness(raw_matrix, 2.0),
    )

    raw_state = _state(100.0, 2.0)
    legacy_state = correct_hydrostatic_elastic_state(raw_state)
    canonical_state = convert_eulerian_hydrostatic_elastic_state(raw_state)
    np.testing.assert_allclose(legacy_state.stiffness, canonical_state.stiffness)
    assert legacy_state.prestress == canonical_state.prestress

    raw_series = ElasticStateSeries(states=(raw_state,), reference_index=0)
    legacy_series = correct_hydrostatic_elastic_series(raw_series)
    canonical_series = convert_eulerian_hydrostatic_elastic_series(raw_series)
    np.testing.assert_allclose(
        legacy_series.states[0].stiffness,
        canonical_series.states[0].stiffness,
    )
    assert legacy_series.metadata == canonical_series.metadata


def test_energy_pressure_resolution_reuses_fit_matching_and_assignment() -> None:
    """One backend-neutral service owns E(V) fitting and volume assignment."""
    raw = ElasticStateSeries(
        states=(
            _state(90.0, None),
            _state(100.0, None),
            _state(110.0, None),
        ),
        reference_index=1,
        metadata={"dataset": "synthetic"},
    )
    source_volumes = np.asarray([110.0, 95.0, 90.0, 105.0, 100.0])
    source_energies = -100.0 + 1.0e-4 * (source_volumes - 100.0) ** 2

    resolution = resolve_energy_derived_pressures(
        raw,
        source_volumes,
        source_energies,
        pressure_source=PressureSource.ENERGY_POLYNOMIAL,
        source_dataset="synthetic-e-v",
        energy_unit="Ha",
        volume_length_unit="angstrom",
        volume_unit="angstrom^3",
        polynomial_degree=2,
    )

    expected = energy_to_pressure(
        -2.0e-4 * (raw.volumes - 100.0),
        "Ha",
        "angstrom",
        "GPa",
    )
    np.testing.assert_allclose(
        resolution.pressures_gpa, expected, rtol=2.0e-11, atol=1.0e-10
    )
    assert [match.source_index for match in resolution.matches] == [2, 4, 0]
    assert resolution.provenance["method"] == "energy_polynomial"
    assert resolution.provenance["relation"] == "P(V) = -dE/dV"
    assert resolution.provenance["source_dataset"] == "synthetic-e-v"
    assert resolution.provenance["settings"] == {"degree": 2}
    assert resolution.provenance["volume_unit"] == "angstrom^3"
    assert resolution.provenance["volume_matches"][0]["source_index"] == 2
    assignment = resolution.series.metadata["pressure_assignment"]
    assert assignment["fit"]["success"] is True
    assert assignment["volume_matches"][1]["source_index"] == 4
    np.testing.assert_allclose(resolution.series.stiffness, raw.stiffness)
    assert all(
        state.prestress.pressure_source is PressureSource.ENERGY_POLYNOMIAL
        for state in resolution.series.states
    )
