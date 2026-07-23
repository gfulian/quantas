from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.eos import PressureEOSDiagnostics


def test_birch_normalized_pressure_has_expected_bm3_line():
    diagnostics = PressureEOSDiagnostics()
    k0 = 37.0
    kp = 6.0
    v0 = 113.0
    volume = np.array([110.0, 105.0, 100.0])
    f = 0.5 * ((v0 / volume) ** (2.0 / 3.0) - 1.0)
    pressure = 3.0 * k0 * f * (1.0 + 2.0 * f) ** 2.5 * (1.0 + 1.5 * (kp - 4.0) * f)

    result = diagnostics.transform("BM3", pressure, volume, v0)

    np.testing.assert_allclose(result.strain, f)
    np.testing.assert_allclose(
        result.normalized_pressure,
        k0 * (1.0 + 1.5 * (kp - 4.0) * f),
    )
    assert np.all(result.valid)


def test_natural_normalized_pressure_matches_definition():
    diagnostics = PressureEOSDiagnostics()
    v0 = 100.0
    volume = np.array([98.0, 95.0])
    pressure = np.array([2.0, 5.0])
    strain = np.log(v0 / volume) / 3.0

    result = diagnostics.transform("NS3", pressure, volume, v0)

    np.testing.assert_allclose(result.strain, strain)
    np.testing.assert_allclose(
        result.normalized_pressure,
        pressure / (3.0 * (v0 / volume) * strain),
    )


def test_vinet_normalized_pressure_matches_definition():
    diagnostics = PressureEOSDiagnostics()
    v0 = 100.0
    volume = np.array([98.0, 92.0])
    pressure = np.array([2.0, 8.0])
    strain = 1.0 - (volume / v0) ** (1.0 / 3.0)

    result = diagnostics.transform("V3", pressure, volume, v0)

    np.testing.assert_allclose(result.strain, strain)
    np.testing.assert_allclose(
        result.normalized_pressure,
        pressure * (1.0 - strain) ** 2 / (3.0 * strain),
    )


def test_reference_state_is_marked_singular_not_fabricated():
    result = PressureEOSDiagnostics().transform("BM3", [0.0], [100.0], 100.0)

    assert not result.valid[0]
    assert np.isnan(result.normalized_pressure[0])


def test_transform_derivatives_match_finite_differences():
    diagnostics = PressureEOSDiagnostics()
    p = np.array([4.0])
    x = np.array([96.0])
    x0 = 100.0
    step = 1.0e-5

    for model in ("BM3", "NS3", "V3"):
        result = diagnostics.transform(model, p, x, x0)
        plus_x = diagnostics.transform(model, p, x + step, x0)
        minus_x = diagnostics.transform(model, p, x - step, x0)
        numerical_x = (plus_x.normalized_pressure - minus_x.normalized_pressure) / (
            2.0 * step
        )
        plus_x0 = diagnostics.transform(model, p, x, x0 + step)
        minus_x0 = diagnostics.transform(model, p, x, x0 - step)
        numerical_x0 = (plus_x0.normalized_pressure - minus_x0.normalized_pressure) / (
            2.0 * step
        )
        np.testing.assert_allclose(
            result.dnormalized_dcoordinate, numerical_x, rtol=2e-5, atol=1e-8
        )
        np.testing.assert_allclose(
            result.dnormalized_dreference, numerical_x0, rtol=2e-5, atol=1e-8
        )


def test_unsupported_eos_rejected():
    with pytest.raises(NotImplementedError, match="normalized-pressure"):
        PressureEOSDiagnostics().transform("M", [1.0], [99.0], 100.0)
