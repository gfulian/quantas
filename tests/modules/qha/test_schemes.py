"""Tests for distinct frequency and thermodynamic QHA schemes."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from quantas.modules.qha.calculator import QHACalculator
from quantas.modules.qha.models import QHAInput, QHAOptions
from quantas.modules.qha.report import phonon_frequency_fit_tables


def _input() -> QHAInput:
    volume = np.linspace(9.2, 10.8, 7)
    energy = 2.0e-3 * (volume - 10.0) ** 2
    base = 450.0 * (10.0 / volume) ** 2.7
    frequencies = np.stack(
        [
            base,
            1.1 * base + 0.8 * (volume - 10.0) ** 3,
            1.4 * base - 0.5 * (volume - 10.0) ** 3,
        ],
        axis=0,
    )[None, :, :]
    return QHAInput(
        jobname="scheme-separation",
        natoms=1,
        qpoints=1,
        volume=volume,
        energy=energy,
        frequencies=frequencies,
        weights=np.array([1.0]),
    )


def _options(scheme: str) -> QHAOptions:
    return QHAOptions(
        temperature_min=300.0,
        temperature_max=1000.0,
        temperature_step=700.0,
        pressure_min=0.0,
        pressure_max=1.0,
        pressure_step=1.0,
        scheme=scheme,
        minimization="poly",
        energy_degree=2,
        free_energy_degree=3,
        frequency_degree=3,
    )


def test_frequency_and_thermodynamic_schemes_use_distinct_evaluators() -> None:
    freq = QHACalculator(_input(), _options("freq")).execute().results["qha"]
    td = QHACalculator(_input(), _options("td")).execute().results["qha"]

    np.testing.assert_allclose(freq.equilibrium_volume, td.equilibrium_volume)
    assert not np.array_equal(
        freq.isochoric_heat_capacity,
        td.isochoric_heat_capacity,
    )
    assert freq.metadata["thermodynamics"]["scheme"] == "freq"
    assert "mode-resolved frequencies" in freq.metadata["thermodynamics"]["source"]
    assert td.metadata["thermodynamics"]["scheme"] == "td"
    assert "sampled volume grid" in td.metadata["thermodynamics"]["source"]


def test_frequency_degree_controls_frequency_fit_diagnostics() -> None:
    options = _options("freq")
    options.frequency_degree = 4
    _, summary = phonon_frequency_fit_tables(_input(), options)

    assert ["Polynomial degree", 4] in summary.rows


def test_frequency_scheme_uses_single_eos_workflow_with_uncertainties() -> None:
    from quantas.modules.qha.io.reader import read_qha_input

    input_data = read_qha_input(Path(__file__).parent / "data" / "mgo_b3lyp_qha.yaml")
    options = QHAOptions(
        temperature_min=300.0,
        temperature_max=300.0,
        temperature_step=1.0,
        pressure_min=0.0,
        pressure_max=1.0,
        pressure_step=1.0,
        scheme="freq",
        minimization="eos",
        eos="BM",
        energy_degree=3,
        free_energy_degree=3,
        frequency_degree=3,
        estimate_uncertainties=True,
        uncertainty_method="covariance",
    )

    result = QHACalculator(input_data, options).execute().results["qha"]

    assert result.completed is True
    assert result.metadata["eos_workflow"]["fit_count"] == 1
    assert result.metadata["eos_workflow"]["state_count"] == 2
    assert result.metadata["thermodynamics"]["scheme"] == "freq"
    assert set(result.uncertainties) >= {"sigma_VT", "sigma_KT", "sigma_Kp"}
    assert np.all(np.isfinite(result.isochoric_heat_capacity))
