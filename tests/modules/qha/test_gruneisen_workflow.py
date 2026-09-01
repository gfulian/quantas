"""Tests for mode-continuity and Grüneisen QHA workflows."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from quantas.modules.qha.calculator import QHACalculator
from quantas.modules.qha.core.gruneisen import (
    ModeGruneisenEvaluator,
    gruneisen_from_power_law,
    thermal_gruneisen_from_modes,
)
from quantas.modules.qha.io.export import QHAHDF5Export
from quantas.modules.qha.io.hdf5 import read_qha_hdf5
from quantas.modules.qha.io.reader import read_qha_input
from quantas.modules.qha.models import QHAInput, QHAOptions
from quantas.modules.qha.thermodynamics import FrequencyThermodynamicEvaluator


def _power_law_input(*, continuity: str = "verified") -> QHAInput:
    volumes = np.linspace(9.4, 10.8, 8)
    energy = 2.0e-3 * (volumes - 10.0) ** 2
    gamma = np.array([1.2, 1.8, 2.4])
    frequencies = np.stack(
        [
            gruneisen_from_power_law(
                volumes,
                value,
                400.0 + 100.0 * index,
                10.0,
            )
            for index, value in enumerate(gamma)
        ],
        axis=0,
    )[None, :, :]
    return QHAInput(
        jobname="gruneisen-workflow",
        natoms=1,
        formula_units=1,
        qpoints=1,
        mode_continuity=continuity,  # type: ignore[arg-type]
        volume=volumes,
        energy=energy,
        frequencies=frequencies,
        weights=np.array([1.0]),
    )


def _options() -> QHAOptions:
    return QHAOptions(
        temperature_min=0.0,
        temperature_max=600.0,
        temperature_step=300.0,
        pressure_min=0.0,
        pressure_max=1.0,
        pressure_step=1.0,
        scheme="freq",
        minimization="poly",
        energy_degree=2,
        free_energy_degree=3,
        frequency_degree=2,
        polynomial_derivative_method="analytic",
        calculate_gruneisen=True,
        calculate_mode_gruneisen=True,
    )


def test_mode_gruneisen_evaluator_interpolates_power_law() -> None:
    volumes = np.linspace(8.0, 12.0, 7)
    gamma_values = np.array([0.7, 1.4, 2.1])
    frequencies = np.vstack(
        [
            gruneisen_from_power_law(volumes, gamma, 300.0, 10.0)
            for gamma in gamma_values
        ]
    )
    evaluator = ModeGruneisenEvaluator(frequencies, volumes, degree=2)

    calculated = evaluator.gamma_at(np.array([9.25, 10.75]))

    np.testing.assert_allclose(
        calculated,
        np.repeat(gamma_values[:, None], 2, axis=1),
        rtol=1.0e-11,
        atol=1.0e-11,
    )


def test_mode_gruneisen_matches_numeric_derivative() -> None:
    input_data = _power_law_input()
    evaluator = FrequencyThermodynamicEvaluator(input_data, _options())
    volumes = np.array([9.7, 10.2, 10.6])
    step = 1.0e-5
    plus = evaluator.frequencies_at(volumes + step)
    minus = evaluator.frequencies_at(volumes - step)
    derivative = (plus - minus) / (2.0 * step)
    frequencies = evaluator.frequencies_at(volumes)
    expected = -(volumes * derivative) / frequencies

    np.testing.assert_allclose(
        evaluator.mode_gruneisen_at(volumes),
        expected,
        rtol=2.0e-9,
        atol=2.0e-9,
    )


def test_mode_weighted_gruneisen_has_correct_high_temperature_limit() -> None:
    gamma = np.array(
        [
            [[1.0], [2.0]],
            [[3.0], [5.0]],
        ]
    )
    frequencies = np.full_like(gamma, 1.0e12)
    result = thermal_gruneisen_from_modes(
        gamma,
        frequencies,
        1.0e8,
        np.array([1.0, 3.0]),
    )

    assert result.success
    assert result.gamma is not None
    expected = (1.0 + 2.0 + 3.0 * 3.0 + 3.0 * 5.0) / 8.0
    assert result.gamma[0] == pytest.approx(expected, rel=1.0e-8)

    zero = thermal_gruneisen_from_modes(
        gamma,
        frequencies,
        0.0,
        np.array([1.0, 3.0]),
    )
    assert zero.success
    np.testing.assert_allclose(zero.gamma, 0.0)


def test_frequency_workflow_calculates_both_gruneisen_definitions() -> None:
    result = QHACalculator(_power_law_input(), _options()).execute().results["qha"]

    assert result.gruneisen is not None
    assert result.mode_weighted_gruneisen is not None
    assert result.mode_gruneisen is not None
    assert result.gruneisen.shape == (3, 2)
    assert result.mode_weighted_gruneisen.shape == (3, 2)
    assert result.mode_gruneisen.shape == (1, 3, 8)
    evaluator = FrequencyThermodynamicEvaluator(_power_law_input(), _options())
    expected = evaluator.mode_gruneisen_at(_power_law_input().volume)
    np.testing.assert_allclose(
        result.mode_gruneisen,
        expected,
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(result.gruneisen[0], 0.0)
    np.testing.assert_allclose(result.mode_weighted_gruneisen[0], 0.0)
    assert result.metadata["gruneisen"]["low_heat_capacity_policy"] == "nan"
    assert result.metadata["gruneisen"]["mode_continuity"] == "verified"
    assert result.metadata["gruneisen"]["n_usable_modes"] == 3
    assert result.metadata["gruneisen"]["mode_source"] == (
        "derivative of fitted frequency-volume polynomials"
    )
    assert "consistency" in result.metadata["gruneisen"]


def test_unknown_mode_continuity_is_rejected_for_frequency_qha() -> None:
    calculator = QHACalculator(
        _power_law_input(continuity="unknown"),
        _options(),
    )

    with pytest.raises(ValueError, match="verified or assumed"):
        calculator.execute()


def test_qha_reader_preserves_explicit_mode_continuity(tmp_path: Path) -> None:
    source = Path(__file__).parent / "data" / "mgo_b3lyp_qha.yaml"
    payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    payload["mode_continuity"] = "verified"
    payload["mode_continuity_metadata"] = {
        "method": "eigenvector_overlap",
        "minimum_overlap": 0.97,
    }
    filename = tmp_path / "verified.yaml"
    filename.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")

    data = read_qha_input(filename)

    assert data.mode_continuity_status() == "verified"
    assert data.has_verified_mode_continuity()
    assert data.metadata["mode_continuity_metadata"]["method"] == (
        "eigenvector_overlap"
    )


def test_gruneisen_arrays_survive_hdf5_round_trip(tmp_path: Path) -> None:
    generic = QHACalculator(_power_law_input(), _options()).execute()
    filename = tmp_path / "gruneisen.hdf5"
    QHAHDF5Export().export(generic, filename)

    loaded = read_qha_hdf5(filename).results["qha"]

    np.testing.assert_allclose(loaded.gruneisen, generic.results["qha"].gruneisen)
    np.testing.assert_allclose(
        loaded.mode_weighted_gruneisen,
        generic.results["qha"].mode_weighted_gruneisen,
    )
    np.testing.assert_allclose(
        loaded.mode_gruneisen,
        generic.results["qha"].mode_gruneisen,
    )


def test_assumed_mode_continuity_emits_warning() -> None:
    generic = QHACalculator(
        _power_law_input(continuity="assumed"),
        _options(),
    ).execute()

    assert any("continuity is assumed" in item for item in generic.warnings)


def test_td_workflow_calculates_only_thermodynamic_gruneisen() -> None:
    options = _options()
    options.scheme = "td"
    input_data = _power_law_input(continuity="unknown")
    result = QHACalculator(input_data, options).execute().results["qha"]

    assert result.gruneisen is not None
    assert result.mode_weighted_gruneisen is None
    assert result.mode_gruneisen is None
