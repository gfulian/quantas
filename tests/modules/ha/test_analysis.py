"""Tests for HA high-level analysis functions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.modules.ha.analysis import (
    calculate_thermodynamic_properties,
    prepare_phonon_data,
    validate_input,
)
from quantas.modules.ha.io.reader import read_ha_input
from quantas.modules.ha.models import HAInput, HAOptions


DATA = Path(__file__).parent / "data/mgo_b3lyp_qha.yaml"


def minimal_input() -> HAInput:
    """Return a minimal valid HA input object."""
    return HAInput(
        jobname="minimal",
        natoms=1,
        qpoints=1,
        volume=np.array([10.0, 11.0]),
        energy=np.array([-1.0, -0.9]),
        frequencies=np.array([[[100.0, 110.0], [200.0, 210.0], [300.0, 310.0]]]),
        weights=np.array([2.0]),
    )


def test_validate_input_accepts_minimal_dataset() -> None:
    """Valid HA input data should pass validation."""
    validate_input(minimal_input())


def test_validate_input_rejects_inconsistent_modes() -> None:
    """The number of modes must be three times the number of atoms."""
    data = minimal_input()
    data.frequencies = np.ones((1, 2, 2), dtype=float)

    with pytest.raises(ValueError, match="three times natoms"):
        validate_input(data)


def test_prepare_phonon_data_normalizes_weights_and_converts_frequency() -> None:
    """Prepared data should use normalized q weights and Hz frequencies."""
    data = minimal_input()
    options = HAOptions(
        temperature_min=0.0,
        temperature_max=100.0,
        temperature_step=100.0,
        frequency_unit="cm^-1",
    )

    temperature, frequencies, weights, static_energy = prepare_phonon_data(
        data,
        options,
    )

    assert np.allclose(temperature, [0.0, 100.0])
    assert np.allclose(weights, [1.0])
    assert frequencies[0, 0, 0] > data.frequencies[0, 0, 0]
    assert np.allclose(static_energy, data.energy)


def test_calculate_thermodynamic_properties_shapes_and_progress() -> None:
    """HA analysis should fill all harmonic thermodynamic result arrays."""
    data = minimal_input()
    options = HAOptions(
        temperature_min=0.0,
        temperature_max=200.0,
        temperature_step=100.0,
        energy_unit="Ha",
        frequency_unit="cm^-1",
    )
    progress: list[tuple[str, int, int]] = []

    result = calculate_thermodynamic_properties(
        data,
        options,
        progress_callback=lambda label, current, total: progress.append(
            (label, current, total)
        ),
    )

    assert result.temperature.shape == (3,)
    assert result.volume.shape == (2,)
    assert result.zero_point_energy.shape == (1, 2)
    assert result.thermal_energy.shape == (3, 2)
    assert result.internal_energy.shape == (3, 2)
    assert result.entropy.shape == (3, 2)
    assert result.vibrational_free_energy.shape == (3, 2)
    assert result.free_energy.shape == (3, 2)
    assert result.isochoric_heat_capacity.shape == (3, 2)
    assert progress[-1] == ("total_energies", 6, 6)
    assert [item[1] for item in progress] == [1, 2, 3, 4, 5, 6]


def test_calculate_thermodynamic_properties_supports_metadata() -> None:
    """HA results should carry library-level metadata for later export."""
    result = calculate_thermodynamic_properties(
        minimal_input(),
        HAOptions(temperature_min=298.15, temperature_max=298.15),
    )

    assert result.metadata["module"] == "ha"
    assert result.metadata["method"] == "harmonic"
    assert result.metadata["input"]["natoms"] == 1
    assert result.has_thermodynamic_data()


def test_reader_and_analysis_process_mgo_yaml_example() -> None:
    """The uploaded MgO HA/QHA YAML example should be usable directly."""
    data = read_ha_input(DATA)
    options = HAOptions(
        temperature_min=0.0,
        temperature_max=300.0,
        temperature_step=300.0,
        energy_unit="Ha",
        frequency_unit="cm^-1",
    )

    result = calculate_thermodynamic_properties(data, options)

    assert data.jobname.startswith("Quasi-Harmonic analysis")
    assert data.natoms == 2
    assert data.qpoints == 32
    assert data.nvol == 11
    assert data.frequencies.shape == (32, 6, 11)
    assert result.free_energy.shape == (2, 11)
    assert np.all(np.isfinite(result.free_energy))


def test_heat_capacity_is_stored_in_true_hartree_per_cell_kelvin() -> None:
    """The high-temperature MgO heat capacity approaches 3NkB in native units."""
    from scipy import constants as cs

    data = read_ha_input(DATA)
    result = calculate_thermodynamic_properties(
        data,
        HAOptions(
            temperature_min=3000.0,
            temperature_max=3000.0,
            temperature_step=1.0,
            energy_unit="Ha",
            frequency_unit="cm^-1",
        ),
    )

    cv_native = float(np.max(result.isochoric_heat_capacity))
    cv_j_mol_cell_k = (
        cv_native * cs.physical_constants["Hartree energy"][0] * cs.Avogadro
    )
    dulong_petit = 3.0 * data.natoms * cs.R

    assert cv_j_mol_cell_k / dulong_petit == pytest.approx(1.0, rel=0.03)
    assert result.metadata["units"]["heat_capacity"] == "Ha cell^-1 K^-1"
