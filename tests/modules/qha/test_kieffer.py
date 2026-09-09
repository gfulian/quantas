"""Tests for multi-volume Kieffer enrichment of thermodynamic QHA."""

from __future__ import annotations

from copy import deepcopy

import numpy as np
import pytest

from quantas.api import qha
from quantas.core.physics.thermodynamics import kieffer_vibrational_free_energy
from quantas.core.physics.units import convert_energy
from quantas.models.kieffer import KiefferCutoffState, KiefferVolumeSeries
from quantas.modules.qha.kieffer import validate_kieffer_qha_applicability
from quantas.modules.qha.models import QHAInput, QHAOptions
from quantas.modules.qha.thermodynamics import (
    FrequencyThermodynamicEvaluator,
    calculate_sampled_thermodynamics,
)


def _input() -> QHAInput:
    """Return primitive Gamma phonons on a deliberately reordered volume grid."""
    volume = np.array([10.0, 9.0, 11.0, 9.5, 10.5], dtype=np.float64)
    return QHAInput(
        jobname="Kieffer QHA",
        natoms=1,
        supercell=np.eye(3),
        qpoints=1,
        qcoords=np.zeros((1, 3)),
        volume=volume,
        energy=0.02 * (volume - 10.0) ** 2 - 10.0,
        frequencies=np.broadcast_to(
            np.array([[[100.0], [200.0], [300.0]]]),
            (1, 3, volume.size),
        ).copy(),
        weights=np.array([1.0]),
    )


def _cutoffs() -> KiefferVolumeSeries:
    """Return direct cutoffs ordered independently by increasing volume."""
    states = []
    for index, volume in enumerate((9.0, 9.5, 10.0, 10.5, 11.0)):
        scale = (10.0 / volume) ** (1.0 / 3.0)
        states.append(
            KiefferCutoffState(
                volume=volume,
                frequencies_hz=scale * np.array([3.0e12, 4.0e12, 5.0e12]),
                effective_velocities_km_s=np.array([3.0, 4.0, 7.0]),
                source_elastic_indices=(index,),
            )
        )
    return KiefferVolumeSeries(states=tuple(states))


def _options(**overrides) -> QHAOptions:
    values = {
        "temperature_min": 0.0,
        "temperature_max": 300.0,
        "temperature_step": 300.0,
        "pressure_min": 0.0,
        "pressure_max": 0.0,
        "scheme": "td",
        "minimization": "poly",
        "energy_degree": 2,
        "free_energy_degree": 2,
        "frequency_unit": "cm^-1",
        "estimate_uncertainties": False,
        "calculate_mode_gruneisen": False,
    }
    values.update(overrides)
    return QHAOptions(**values)


def test_sampled_qha_adds_kieffer_in_explicit_volume_order() -> None:
    """Acoustic properties are matched, reordered, retained, and added."""
    data = _input()
    original = data.frequencies.copy()
    harmonic = calculate_sampled_thermodynamics(data, _options())
    enriched = calculate_sampled_thermodynamics(
        data,
        _options(),
        kieffer_cutoffs=_cutoffs(),
    )

    contribution = enriched.kieffer_contribution
    assert contribution is not None
    expected_cutoffs = _cutoffs().frequencies_hz[:, [2, 0, 4, 1, 3]]
    np.testing.assert_allclose(contribution.cutoff_frequencies_hz, expected_cutoffs)
    expected_fvib = convert_energy(
        kieffer_vibrational_free_energy(
            np.array([0.0, 300.0]),
            expected_cutoffs,
        ),
        "kjmol",
        "Ha",
    )
    np.testing.assert_allclose(
        enriched.vibrational_free_energy - harmonic.vibrational_free_energy,
        expected_fvib,
    )
    np.testing.assert_array_equal(data.frequencies, original)
    assert [
        item["cutoff_index"] for item in contribution.metadata["volume_matches"]
    ] == [2, 0, 4, 1, 3]


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda data: setattr(data, "qcoords", None), "explicit Gamma"),
        (
            lambda data: setattr(data, "qcoords", np.array([[0.1, 0.0, 0.0]])),
            "Gamma q-point",
        ),
        (
            lambda data: setattr(data, "supercell", np.diag([2.0, 1.0, 1.0])),
            "identity supercell",
        ),
        (
            lambda data: setattr(
                data,
                "volume",
                np.array([10.0, 9.0, 11.1, 9.5, 10.5]),
            ),
            "no source match",
        ),
    ],
)
def test_kieffer_qha_rejects_inapplicable_inputs(change, message: str) -> None:
    """QHA applicability failures are explicit and occur before fitting."""
    data = deepcopy(_input())
    change(data)
    with pytest.raises(ValueError, match=message):
        validate_kieffer_qha_applicability(data, _cutoffs())


def test_qha_api_td_round_trip_retains_sampled_kieffer_component(tmp_path) -> None:
    """The public TD workflow and native archive retain acoustic provenance."""
    envelope = qha.run(
        _input(),
        options=_options(),
        kieffer_cutoffs=_cutoffs(),
    )
    result = qha.get_result(envelope)
    assert result.completed is True
    assert result.kieffer_sampled_contribution is not None
    assert result.metadata["kieffer"]["composition"] == ("additional-acoustic-branches")

    archive = qha.write_result(envelope, tmp_path / "qha-kieffer")
    loaded = qha.get_result(qha.read_result(archive))
    assert loaded.kieffer_sampled_contribution is not None
    np.testing.assert_allclose(
        loaded.kieffer_sampled_contribution.cutoff_frequencies_hz,
        result.kieffer_sampled_contribution.cutoff_frequencies_hz,
    )


def test_frequency_evaluator_adds_fitted_kieffer_properties() -> None:
    """The frequency route evaluates acoustic cutoffs at arbitrary volumes."""
    options = _options(scheme="freq")
    harmonic = FrequencyThermodynamicEvaluator(_input(), options)
    enriched = FrequencyThermodynamicEvaluator(
        _input(),
        options,
        kieffer_cutoffs=_cutoffs(),
    )
    volumes = np.array([9.25, 10.25], dtype=np.float64)
    temperature = 300.0
    expected = convert_energy(
        kieffer_vibrational_free_energy(
            np.array([temperature]),
            enriched.kieffer_cutoffs_at(volumes),
        )[0],
        "kjmol",
        "Ha",
    )
    np.testing.assert_allclose(
        enriched.properties_at(volumes, temperature)["vibrational_free_energy"]
        - harmonic.properties_at(volumes, temperature)["vibrational_free_energy"],
        expected,
    )


def test_frequency_evaluator_rejects_nonpositive_extrapolated_cutoff() -> None:
    """A fitted acoustic branch cannot disappear during extrapolation."""
    evaluator = FrequencyThermodynamicEvaluator(
        _input(),
        _options(scheme="freq"),
        kieffer_cutoffs=_cutoffs(),
    )
    with pytest.raises(ValueError, match="non-positive or non-finite"):
        evaluator.kieffer_cutoffs_at(100.0)


def test_qha_frequency_scheme_includes_kieffer_at_equilibrium() -> None:
    """Local minimization and final frequency thermodynamics include Kieffer."""
    options = _options(scheme="freq")
    harmonic = qha.get_result(qha.run(_input(), options=options))
    enriched = qha.get_result(
        qha.run(
            _input(),
            options=options,
            kieffer_cutoffs=_cutoffs(),
        )
    )
    assert enriched.completed is True
    assert enriched.equilibrium_volume is not None
    assert harmonic.equilibrium_volume is not None
    assert not np.allclose(enriched.equilibrium_volume, harmonic.equilibrium_volume)
    assert enriched.metadata["thermodynamics"]["kieffer_cutoff_fit_status"] == (
        "success"
    )


def test_kieffer_frequency_scheme_rejects_mode_gruneisen_analysis() -> None:
    """An incomplete optical-only mode-Gruneisen average is never reported."""
    with pytest.raises(ValueError, match="does not yet support mode-Gruneisen"):
        qha.run(
            _input(),
            options=_options(scheme="freq", calculate_mode_gruneisen=True),
            kieffer_cutoffs=_cutoffs(),
        )
