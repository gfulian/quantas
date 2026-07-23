from __future__ import annotations

import numpy as np
import pytest

from quantas.modules.ha.models import HAInput, HAOptions, HAResult


def test_ha_input_exposes_shared_derived_quantities():
    ha_input = HAInput(
        jobname="HA model test",
        natoms=2,
        supercell=np.diag([2, 2, 1]),
        qpoints=2,
        volume=np.array([10.0, 11.0]),
        energy=np.array([-20.0, -19.5]),
        frequencies=np.zeros((2, 6, 2)),
        weights=np.array([1.0, 3.0]),
    )

    assert ha_input.jobname == "HA model test"
    assert ha_input.nvol == 2
    assert ha_input.nmodes == 6
    assert ha_input.kpoints == 4
    assert ha_input.total_q_points == 4.0
    assert ha_input.has_phonons() is True
    np.testing.assert_allclose(ha_input.normalized_weights(), [0.25, 0.75])


def test_ha_input_rejects_invalid_weight_normalization():
    ha_input = HAInput(weights=np.array([0.0, 0.0]))

    with pytest.raises(ValueError, match="sum of q-point weights"):
        ha_input.normalized_weights()


def test_ha_options_temperature_grid_includes_the_upper_bound():
    options = HAOptions(
        temperature_min=100.0,
        temperature_max=300.0,
        temperature_step=100.0,
    )

    np.testing.assert_allclose(options.temperature_grid(), [100.0, 200.0, 300.0])


def test_ha_options_rejects_invalid_temperature_range():
    with pytest.raises(ValueError, match="temperature_step"):
        HAOptions(temperature_step=0.0).temperature_grid()

    with pytest.raises(ValueError, match="temperature_max"):
        HAOptions(temperature_min=300.0, temperature_max=100.0).temperature_grid()


def test_ha_result_returns_compact_property_names():
    result = HAResult(
        zero_point_energy=np.array([[1.0, 2.0]]),
        thermal_energy=np.array([[0.1, 0.2]]),
        entropy=np.array([[0.01, 0.02]]),
    )

    properties = result.as_property_dict()

    assert result.has_thermodynamic_data() is True
    assert set(properties) == {"Uzp", "Uth", "S"}
    np.testing.assert_allclose(properties["Uzp"], [[1.0, 2.0]])
