"""Tests for analytical volume-temperature equations of state."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from quantas.core.physics.eos import (
    TemperatureEOS,
    TemperatureEOSFamily,
    TemperatureEOSParameters,
    TemperatureEOSVariant,
    available_temperature_eos_models,
    parse_temperature_eos_model,
)

ROOT = Path(__file__).resolve().parents[3]
GOLDEN = ROOT / "tests/reference/eos/eos_vt_formula_reference.json"


def _golden() -> dict[str, object]:
    return json.loads(GOLDEN.read_text(encoding="utf-8"))


@pytest.mark.parametrize(
    ("name", "model", "parameters"),
    [
        (
            "berman",
            "berman",
            {
                "V0": 100.0,
                "temperature_ref": 298.15,
                "alpha0": 3.0e-5,
                "alpha1": 1.0e-8,
            },
        ),
        (
            "fei",
            "fei",
            {
                "V0": 100.0,
                "temperature_ref": 298.15,
                "alpha0": 2.5e-5,
                "alpha1": 1.0e-8,
                "alpha2": 0.1,
            },
        ),
        (
            "modified_holland_powell",
            "mhp",
            {
                "V0": 100.0,
                "temperature_ref": 298.15,
                "alpha0": 3.0e-5,
                "alpha1": 5.0e-5,
            },
        ),
        (
            "salje",
            "salje",
            {"V0": 100.0, "temperature_ref": 0.0, "p1": 1.0e-4, "theta_sat": 300.0},
        ),
        (
            "kroll_holland_powell",
            "khp",
            {
                "V0": 100.0,
                "temperature_ref": 298.15,
                "alpha_ref": 3.0e-5,
                "theta_e": 500.0,
                "kp": 4.0,
            },
        ),
    ],
)
def test_temperature_models_match_formula_reference(
    name: str, model: str, parameters: dict[str, float]
) -> None:
    reference = _golden()["models"][name]
    temperatures = np.asarray([point["temperature_K"] for point in reference["points"]])
    expected_value = np.asarray([point["value"] for point in reference["points"]])
    expected_alpha = np.asarray([point["alpha_K-1"] for point in reference["points"]])
    eos = TemperatureEOS()
    np.testing.assert_allclose(
        eos.value(model, parameters, temperatures),
        expected_value,
        rtol=2e-13,
        atol=2e-13,
    )
    np.testing.assert_allclose(
        eos.expansion_coefficient(model, parameters, temperatures),
        expected_alpha,
        rtol=3e-13,
        atol=3e-13,
    )
    np.testing.assert_allclose(
        eos.derivative(model, parameters, temperatures),
        expected_value * expected_alpha,
        rtol=3e-13,
        atol=3e-13,
    )


def test_reference_identities() -> None:
    eos = TemperatureEOS()
    cases = [
        (
            "berman",
            {"V0": 10.0, "temperature_ref": 300.0, "alpha0": 2e-5, "alpha1": 1e-8},
            300.0,
        ),
        (
            "fei",
            {
                "V0": 10.0,
                "temperature_ref": 300.0,
                "alpha0": 2e-5,
                "alpha1": 1e-8,
                "alpha2": 0.1,
            },
            300.0,
        ),
        (
            "mhp",
            {"V0": 10.0, "temperature_ref": 300.0, "alpha0": 2e-5, "alpha1": 3e-5},
            300.0,
        ),
        (
            "salje",
            {"V0": 10.0, "temperature_ref": 0.0, "p1": 1e-4, "theta_sat": 250.0},
            0.0,
        ),
        (
            "khp",
            {
                "V0": 10.0,
                "temperature_ref": 300.0,
                "alpha_ref": 2e-5,
                "theta_e": 500.0,
                "kp": 4.0,
            },
            300.0,
        ),
    ]
    for model, parameters, tref in cases:
        assert eos.value(model, parameters, tref)[0] == pytest.approx(10.0)
    assert eos.expansion_coefficient("salje", cases[3][1], 0.0)[0] == 0.0
    assert eos.expansion_coefficient("khp", cases[4][1], 300.0)[0] == pytest.approx(
        2e-5
    )


def test_variants_set_implicit_coefficients() -> None:
    eos = TemperatureEOS()
    t = np.array([300.0, 500.0])
    blinear = {"V0": 10.0, "temperature_ref": 300.0, "alpha0": 2e-5}
    bfull = {**blinear, "alpha1": 0.0}
    np.testing.assert_allclose(
        eos.value("berman", blinear, t, variant="linear"),
        eos.value("berman", bfull, t, variant="quadratic"),
    )
    fei_linear = {"V0": 10.0, "temperature_ref": 300.0, "alpha0": 2e-5, "alpha1": 1e-8}
    fei_full = {**fei_linear, "alpha2": 0.0}
    np.testing.assert_allclose(
        eos.value("fei", fei_linear, t, variant="linear"),
        eos.value("fei", fei_full, t, variant="inverse-square"),
    )


def test_model_registry_and_aliases() -> None:
    assert len(available_temperature_eos_models()) == 8
    assert (
        parse_temperature_eos_model("MHP").family
        is TemperatureEOSFamily.MODIFIED_HOLLAND_POWELL
    )
    assert (
        parse_temperature_eos_model("berman", "linear").variant
        is TemperatureEOSVariant.LINEAR
    )


def test_linear_expansion_conversion() -> None:
    np.testing.assert_allclose(
        TemperatureEOS.linear_expansion_coefficient([3e-5, 6e-5]), [1e-5, 2e-5]
    )


@pytest.mark.parametrize(
    ("model", "parameters", "temperature", "match"),
    [
        (
            "fei",
            {
                "V0": 10.0,
                "temperature_ref": 300.0,
                "alpha0": 2e-5,
                "alpha1": 1e-8,
                "alpha2": 0.1,
            },
            0.0,
            "T > 0",
        ),
        (
            "mhp",
            {"V0": 10.0, "temperature_ref": 300.0, "alpha0": 2e-5, "alpha1": 0.0},
            0.0,
            "T > 0",
        ),
        (
            "salje",
            {"V0": 10.0, "temperature_ref": 1.0, "p1": 1e-4, "theta_sat": 250.0},
            100.0,
            "T_ref = 0",
        ),
        (
            "khp",
            {
                "V0": 10.0,
                "temperature_ref": 300.0,
                "alpha_ref": 2e-5,
                "theta_e": 0.0,
                "kp": 4.0,
            },
            300.0,
            "theta_e",
        ),
    ],
)
def test_model_domains_are_enforced(
    model: str, parameters: dict[str, float], temperature: float, match: str
) -> None:
    with pytest.raises(ValueError, match=match):
        TemperatureEOS().value(model, parameters, temperature)


def test_parameter_mapping_is_strict() -> None:
    with pytest.raises(ValueError, match="missing"):
        TemperatureEOS().value("berman", {"V0": 10.0, "temperature_ref": 300.0}, 300.0)
    with pytest.raises(ValueError, match="unknown"):
        TemperatureEOS().value(
            "berman",
            {"V0": 10.0, "temperature_ref": 300.0, "alpha0": 1e-5, "banana": 1.0},
            300.0,
        )


def test_dataclass_parameters_are_supported() -> None:
    parameters = TemperatureEOSParameters(
        V0=10.0, temperature_ref=300.0, alpha0=2e-5, alpha1=0.0
    )
    assert TemperatureEOS().value("berman", parameters, 300.0)[0] == pytest.approx(10.0)
