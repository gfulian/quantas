"""OHAp regression for the complete Kieffer frequency-QHA composition."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from quantas.api import qha
from quantas.core.physics.thermodynamics import (
    entropy,
    isochoric_heat_capacity,
    kieffer_entropy,
    kieffer_isochoric_heat_capacity,
    kieffer_thermal_energy,
    thermal_energy,
)
from quantas.core.physics.units import convert_frequency
from quantas.models.kieffer import KiefferCutoffState, KiefferVolumeSeries
from quantas.modules.qha.models import QHAInput, QHAResult


_REFERENCE_FILE = Path(__file__).with_name("data") / "ohap_kieffer_reference.json"


def _reference() -> dict[str, Any]:
    """Load the compact values extracted from the complete QM dataset."""
    return json.loads(_REFERENCE_FILE.read_text(encoding="utf-8"))


def _input_and_cutoffs() -> tuple[QHAInput, KiefferVolumeSeries]:
    """Reconstruct the portable numerical subset used by the regression."""
    reference = _reference()
    raw_input = reference["input"]
    raw_kieffer = reference["kieffer"]
    volumes = np.asarray(raw_input["volumes_angstrom3"], dtype=np.float64)
    frequencies = np.asarray(raw_input["frequencies_cm1"], dtype=np.float64)
    cutoffs = np.asarray(raw_kieffer["cutoffs_hz"], dtype=np.float64)
    velocities = np.asarray(raw_kieffer["velocities_km_s"], dtype=np.float64)
    input_data = qha.Input(
        jobname="OHAp Kieffer regression",
        natoms=int(raw_input["natoms"]),
        formula_units=int(raw_input["formula_units"]),
        supercell=np.eye(3, dtype=np.float64),
        qpoints=1,
        qcoords=np.zeros((1, 3), dtype=np.float64),
        volume=volumes,
        energy=np.asarray(raw_input["energies_hartree"], dtype=np.float64),
        frequencies=frequencies[np.newaxis, :, :],
        weights=np.ones(1, dtype=np.float64),
        units={
            "energy": "Ha",
            "volume": "angstrom^3",
            "frequency": "cm^-1",
            "length": "angstrom",
        },
        mode_continuity="verified",
    )
    states = tuple(
        KiefferCutoffState(
            volume=float(volume),
            frequencies_hz=cutoffs[:, index],
            effective_velocities_km_s=velocities[:, index],
            source_elastic_indices=(index,),
        )
        for index, volume in enumerate(volumes)
    )
    return input_data, KiefferVolumeSeries(states=states)


def _acoustic_percent(optical: np.ndarray, acoustic: np.ndarray) -> np.ndarray:
    """Return the acoustic share of a positive additive property."""
    return 100.0 * acoustic[:, 0] / (optical[:, 0] + acoustic[:, 0])


def test_ohap_acoustic_share_follows_temperature_limits() -> None:
    """The real spectrum has the expected low- and high-temperature limits."""
    input_data, series = _input_and_cutoffs()
    reference = _reference()["direct_reference"]
    temperatures = np.asarray(reference["temperature_k"], dtype=np.float64)
    reference_index = int(np.argmin(input_data.energy))
    assert input_data.volume[reference_index] == reference["reference_volume_angstrom3"]

    frequency_hz = np.asarray(
        convert_frequency(
            input_data.frequencies[:, :, reference_index : reference_index + 1],
            "cm^-1",
            "Hz",
        ),
        dtype=np.float64,
    )
    cutoff_hz = series.frequencies_hz[:, reference_index : reference_index + 1]
    weights = input_data.normalized_weights()
    percentages = {
        "thermal_energy": _acoustic_percent(
            thermal_energy(temperatures, frequency_hz, weights),
            kieffer_thermal_energy(temperatures, cutoff_hz),
        ),
        "entropy": _acoustic_percent(
            entropy(temperatures, frequency_hz, weights),
            kieffer_entropy(temperatures, cutoff_hz),
        ),
        "isochoric_heat_capacity": _acoustic_percent(
            isochoric_heat_capacity(temperatures, frequency_hz, weights),
            kieffer_isochoric_heat_capacity(temperatures, cutoff_hz),
        ),
    }
    for name, actual in percentages.items():
        np.testing.assert_allclose(
            actual,
            reference["acoustic_percent"][name],
            rtol=5.0e-8,
            atol=5.0e-8,
        )

    cv_share = percentages["isochoric_heat_capacity"]
    classical_share = 100.0 * 3.0 / (3.0 * input_data.natoms)
    assert cv_share[0] > cv_share[4] > cv_share[5] > classical_share
    optical_share = 100.0 - cv_share
    assert optical_share[0] < optical_share[4] < optical_share[5]
    np.testing.assert_allclose(cv_share[-1], classical_share, rtol=1.0e-8)
    # The 10--20 K crossover is not monotonic because low optical modes enter
    # the same energy window.  The broad low-to-high trend remains physical.
    assert cv_share[3] > cv_share[2]


def _run_reference_pair() -> tuple[QHAInput, QHAResult, QHAResult]:
    """Run the controlled zero-pressure frequency-polynomial comparison."""
    input_data, series = _input_and_cutoffs()
    options = qha.Options(
        temperature_min=0.0,
        temperature_max=1500.0,
        temperature_step=100.0,
        pressure_min=0.0,
        pressure_max=0.0,
        pressure_step=1.0,
        scheme="freq",
        minimization="poly",
        eos="BM3",
        energy_degree=3,
        free_energy_degree=3,
        frequency_degree=3,
        thermal_expansion_method="mixed_derivative",
        calculate_mode_gruneisen=False,
        estimate_uncertainties=False,
        fit_failure_policy="raise",
        extrapolation_policy="fail",
    )
    normal = qha.get_result(qha.run(input_data, options=options))
    enriched = qha.get_result(
        qha.run(input_data, options=options, kieffer_cutoffs=series)
    )
    return input_data, normal, enriched


def test_ohap_frequency_qha_matches_frozen_reference() -> None:
    """The paired QHA calculation reproduces selected OHAp observables."""
    input_data, normal, enriched = _run_reference_pair()
    reference = _reference()["qha_reference"]
    selected = np.asarray(reference["temperature_k"], dtype=np.float64)
    indices = [
        int(np.flatnonzero(normal.temperature == value)[0]) for value in selected
    ]

    for label, result in (("phonons_only", normal), ("with_kieffer", enriched)):
        expected = reference[label]
        for name, values in expected.items():
            actual = np.asarray(getattr(result, name), dtype=np.float64)[indices, 0]
            # K_T is a second derivative of a local polynomial free-energy fit.
            # Supported NumPy/LAPACK combinations differ at about 1e-6 relative
            # for this derived quantity, while primary thermodynamic observables
            # reproduce the frozen reference much more tightly.
            rtol = 2.0e-6 if name == "isothermal_bulk_modulus" else 2.0e-7
            np.testing.assert_allclose(
                actual,
                np.asarray(values, dtype=np.float64),
                rtol=rtol,
                atol=2.0e-10,
            )
        validation = qha.validate_result(result, input_data)
        assert validation.completed is True
        assert validation.finite_properties is True
        assert validation.valid_points == validation.total_points
        assert validation.volumes_below_sampled_range == 0
        assert validation.volumes_above_sampled_range == 0

    assert normal.kieffer_sampled_contribution is None
    assert enriched.kieffer_sampled_contribution is not None
    assert np.all(enriched.equilibrium_volume - normal.equilibrium_volume > 0.0)
    assert np.all(
        enriched.isothermal_bulk_modulus - normal.isothermal_bulk_modulus < 0.0
    )
    finite = normal.temperature > 0.0
    cv_gain = (
        enriched.isochoric_heat_capacity[finite, 0]
        / normal.isochoric_heat_capacity[finite, 0]
        - 1.0
    )
    finite_temperatures = normal.temperature[finite]
    room_index = int(np.flatnonzero(finite_temperatures == 300.0)[0])
    assert cv_gain[room_index] > cv_gain[-1]
