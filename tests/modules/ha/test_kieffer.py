"""Tests for additive Kieffer enrichment of single-volume HA data."""

from __future__ import annotations

from copy import deepcopy

import h5py
import numpy as np
import pytest

from quantas.api import ha
from quantas.core.physics.thermodynamics import kieffer_vibrational_free_energy
from quantas.core.physics.units import convert_energy
from quantas.models.kieffer import KiefferCutoffState, KiefferVolumeSeries
from quantas.modules.ha.analysis import calculate_thermodynamic_properties
from quantas.modules.ha.io.hdf5_payload import (
    read_current_ha_payload,
    write_ha_payload,
)
from quantas.modules.ha.kieffer import validate_kieffer_ha_applicability
from quantas.modules.ha.models import HAInput, HAOptions


def _input() -> HAInput:
    """Return an applicable primitive Gamma-only HA dataset."""
    return HAInput(
        jobname="Kieffer HA",
        natoms=2,
        supercell=np.eye(3),
        qpoints=1,
        qcoords=np.zeros((1, 3)),
        volume=np.array([80.0]),
        energy=np.array([-100.0]),
        frequencies=np.array([[[1.0], [2.0], [3.0], [400.0], [500.0], [600.0]]]),
        weights=np.array([1.0]),
    )


def _cutoffs() -> KiefferVolumeSeries:
    """Return one direct cutoff state at the HA volume."""
    return KiefferVolumeSeries(
        states=(
            KiefferCutoffState(
                volume=80.0,
                frequencies_hz=np.array([3.0e12, 4.0e12, 5.0e12]),
                effective_velocities_km_s=np.array([3.0, 4.0, 7.0]),
                source_elastic_indices=(0,),
            ),
        )
    )


def test_kieffer_ha_is_additive_and_retains_components() -> None:
    """Kieffer adds acoustic branches without changing Gamma frequencies."""
    data = _input()
    original_frequencies = data.frequencies.copy()
    options = HAOptions(
        temperature_min=0.0,
        temperature_max=300.0,
        temperature_step=300.0,
        energy_unit="Ha",
        frequency_unit="cm^-1",
    )
    harmonic = calculate_thermodynamic_properties(data, options)
    enriched = calculate_thermodynamic_properties(
        data,
        options,
        kieffer_cutoffs=_cutoffs(),
    )

    expected = convert_energy(
        kieffer_vibrational_free_energy(
            np.array([0.0, 300.0]),
            _cutoffs().frequencies_hz,
        ),
        "kjmol",
        "Ha",
    )
    np.testing.assert_allclose(
        enriched.vibrational_free_energy - harmonic.vibrational_free_energy,
        expected,
    )
    np.testing.assert_array_equal(data.frequencies, original_frequencies)
    contribution = enriched.kieffer_contribution
    assert contribution is not None
    np.testing.assert_allclose(contribution.vibrational_free_energy, expected)
    assert contribution.metadata["composition"] == "additional-acoustic-branches"
    assert enriched.metadata["kieffer"]["volume_match"]["relative_difference"] == 0.0


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda data: setattr(data, "qpoints", 2), "qpoints == 1"),
        (
            lambda data: setattr(data, "qcoords", np.array([[0.25, 0.0, 0.0]])),
            "Gamma q-point",
        ),
        (lambda data: setattr(data, "qcoords", None), "explicit Gamma"),
        (
            lambda data: setattr(data, "supercell", np.diag([2.0, 1.0, 1.0])),
            "identity supercell",
        ),
        (lambda data: setattr(data, "volume", np.array([81.0])), "no source match"),
    ],
)
def test_kieffer_ha_rejects_inapplicable_inputs(change, message: str) -> None:
    """Every physical applicability condition is checked explicitly."""
    data = deepcopy(_input())
    change(data)

    with pytest.raises(ValueError, match=message):
        validate_kieffer_ha_applicability(data, _cutoffs())


def test_kieffer_contribution_round_trips_through_ha_payload(tmp_path) -> None:
    """The acoustic component remains separate in a native HA archive."""
    result = calculate_thermodynamic_properties(
        _input(),
        HAOptions(
            temperature_min=0.0,
            temperature_max=300.0,
            temperature_step=300.0,
            energy_unit="Ha",
            frequency_unit="cm^-1",
        ),
        kieffer_cutoffs=_cutoffs(),
    )
    filename = tmp_path / "ha-kieffer.hdf5"
    with h5py.File(filename, "w") as archive:
        write_ha_payload(archive, result)
    with h5py.File(filename, "r") as archive:
        loaded = read_current_ha_payload(archive["results"])

    assert loaded.kieffer_contribution is not None
    np.testing.assert_allclose(
        loaded.kieffer_contribution.cutoff_frequencies_hz,
        _cutoffs().frequencies_hz,
    )
    np.testing.assert_allclose(
        loaded.kieffer_contribution.vibrational_free_energy,
        result.kieffer_contribution.vibrational_free_energy,
    )


def test_public_ha_api_accepts_explicit_kieffer_cutoffs(tmp_path) -> None:
    """The public API exposes enrichment without introducing frontend logic."""
    envelope = ha.run(
        _input(),
        options=ha.Options(
            temperature_min=300.0,
            temperature_max=300.0,
            frequency_unit="cm^-1",
        ),
        kieffer_cutoffs=_cutoffs(),
    )

    result = ha.get_result(envelope)
    assert result.kieffer_contribution is not None
    assert result.metadata["kieffer"]["composition"] == ("additional-acoustic-branches")
    archive = ha.write_result(envelope, tmp_path / "public-ha-kieffer")
    loaded = ha.get_result(ha.read_result(archive))
    assert loaded.kieffer_contribution is not None
    np.testing.assert_allclose(
        loaded.kieffer_contribution.cutoff_frequencies_hz,
        _cutoffs().frequencies_hz,
    )
