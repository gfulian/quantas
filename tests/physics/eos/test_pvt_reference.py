"""Independent formula-reference tests for compositional P--V--T EOS."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from quantas.core.physics.eos import PVTEOS, PVTModel


_REFERENCE = (
    Path(__file__).resolve().parents[3]
    / "tests"
    / "reference"
    / "eos"
    / "eos_pvt_formula_reference.json"
)


def _reference() -> dict[str, object]:
    return json.loads(_REFERENCE.read_text(encoding="utf-8"))


def _common_parameters() -> tuple[dict[str, float], dict[str, float]]:
    data = _reference()["reference_parameters"]
    assert isinstance(data, dict)
    pressure = {
        "K0": float(data["K0_GPa"]),
        "KP": float(data["KP"]),
        "KPP": -0.02,
        "V0": float(data["V0_A3"]),
    }
    thermal = {
        "V0": float(data["V0_A3"]),
        "temperature_ref": float(data["temperature_K"]),
        "alpha0": float(data["alpha0_K-1"]),
        "alpha1": float(data["alpha1_K-2"]),
    }
    return pressure, thermal


@pytest.mark.parametrize(
    ("coupling_name", "temperature_model", "coupling_parameters"),
    [
        (
            "linear-bulk-modulus",
            "berman:quadratic",
            {"dK0_dT": -0.02},
        ),
        (
            "anderson-gruneisen",
            "berman:quadratic",
            {"delta": 4.4},
        ),
        (
            "thermal-pressure",
            None,
            {
                "temperature_ref": 300.0,
                "alpha_ref": 3.0e-5,
                "theta_e": 500.0,
            },
        ),
    ],
)
def test_pvt_core_matches_independent_reference(
    coupling_name: str,
    temperature_model: str | None,
    coupling_parameters: dict[str, float],
) -> None:
    """Compare pressure and coupling properties with independent formulae."""
    pressure_parameters, thermal_parameters = _common_parameters()
    model = PVTModel("BM3", coupling_name, temperature_model)
    core = PVTEOS()
    entry = _reference()["couplings"]
    assert isinstance(entry, dict)
    points_entry = entry[coupling_name]
    assert isinstance(points_entry, dict)
    points = points_entry["points"]
    assert isinstance(points, list)

    volume = np.asarray([point["volume_A3"] for point in points], dtype=np.float64)
    temperature = np.asarray(
        [point["temperature_K"] for point in points], dtype=np.float64
    )
    expected_pressure = np.asarray(
        [point["pressure_GPa"] for point in points], dtype=np.float64
    )
    calculated = core.pressure(
        model,
        pressure_parameters,
        None if temperature_model is None else thermal_parameters,
        coupling_parameters,
        volume,
        temperature,
    )
    np.testing.assert_allclose(
        calculated, expected_pressure, rtol=2.0e-13, atol=2.0e-13
    )

    unique_temperature = np.unique(temperature)
    expected_k0 = np.asarray(
        [
            next(
                point["zero_pressure_bulk_modulus_GPa"]
                for point in points
                if point["temperature_K"] == item
            )
            for item in unique_temperature
        ],
        dtype=np.float64,
    )
    expected_derivative = np.asarray(
        [
            next(
                point["dK0_dT_GPa_K-1"]
                for point in points
                if point["temperature_K"] == item
            )
            for item in unique_temperature
        ],
        dtype=np.float64,
    )
    if coupling_name != "thermal-pressure":
        np.testing.assert_allclose(
            core.zero_pressure_bulk_modulus(
                model,
                pressure_parameters,
                thermal_parameters,
                coupling_parameters,
                unique_temperature,
            ),
            expected_k0,
            rtol=2.0e-13,
            atol=2.0e-13,
        )
        np.testing.assert_allclose(
            core.zero_pressure_bulk_modulus_temperature_derivative(
                model,
                pressure_parameters,
                thermal_parameters,
                coupling_parameters,
                unique_temperature,
            ),
            expected_derivative,
            rtol=2.0e-12,
            atol=2.0e-13,
        )
    else:
        expected_thermal_pressure = np.asarray(
            [
                next(
                    point["thermal_pressure_GPa"]
                    for point in points
                    if point["temperature_K"] == item
                )
                for item in unique_temperature
            ],
            dtype=np.float64,
        )
        np.testing.assert_allclose(
            core.thermal_pressure(
                pressure_parameters,
                coupling_parameters,
                unique_temperature,
            ),
            expected_thermal_pressure,
            rtol=2.0e-13,
            atol=2.0e-13,
        )
        np.testing.assert_allclose(
            core.thermal_pressure_temperature_derivative(
                pressure_parameters,
                coupling_parameters,
                unique_temperature,
            ),
            expected_derivative,
            rtol=2.0e-13,
            atol=2.0e-13,
        )
