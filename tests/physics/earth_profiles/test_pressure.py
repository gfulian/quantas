"""Scientific regression tests for terrestrial pressure-depth models."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.earth_profiles import (
    LayeredLithostaticPressureModel,
    LithostaticLayer,
    PremPressureModel,
)

pytestmark = [pytest.mark.physics, pytest.mark.scientific, pytest.mark.fast]


def test_prem_reproduces_global_mass_gravity_and_mantle_pressures() -> None:
    """Integrated PREM density reproduces standard global and mantle anchors."""
    model = PremPressureModel(integration_step_km=0.25)
    assert model.earth_mass_kg == pytest.approx(5.9722e24, rel=3.0e-4)
    assert model.surface_gravity_m_s2 == pytest.approx(9.82, rel=3.0e-3)

    pressure = model.pressure(np.asarray([0.0, 410.0, 660.0, 2891.0]))
    np.testing.assert_allclose(
        pressure,
        np.asarray([0.0, 13.73625, 23.45266, 135.85045]),
        rtol=2.0e-4,
        atol=2.0e-4,
    )


def test_prem_pressure_is_float64_monotonic_and_continuous() -> None:
    """PREM pressure grows continuously through density discontinuities."""
    model = PremPressureModel(integration_step_km=0.5)
    depth = np.linspace(0.0, 2891.0, 3000, dtype=np.float64)
    pressure = model.pressure(depth)
    assert pressure.dtype == np.float64
    assert np.all(np.diff(pressure) > 0.0)
    for boundary in (3.0, 15.0, 24.4, 80.0, 220.0, 400.0, 600.0, 670.0):
        values = model.pressure(np.asarray([boundary - 1.0e-5, boundary + 1.0e-5]))
        assert abs(values[1] - values[0]) < 1.0e-4


def test_prem_provenance_is_complete() -> None:
    """PREM metadata records complete article provenance and numerical method."""
    metadata = PremPressureModel().metadata()
    assert metadata["doi"] == "10.1016/0031-9201(81)90046-7"
    assert "Dziewonski" in metadata["citation"]
    assert "Physics of the Earth and Planetary Interiors" in metadata["citation"]
    assert metadata["integration_step_km"] == 0.25


def test_layered_lithostatic_model_matches_analytic_pressure() -> None:
    """Constant-density layers reproduce the direct rho-g-z expression."""
    model = LayeredLithostaticPressureModel(
        (
            LithostaticLayer(10.0, 2700.0, "upper crust"),
            LithostaticLayer(20.0, 3000.0, "lower crust"),
        ),
        gravity_m_s2=10.0,
    )
    depth = np.asarray([0.0, 5.0, 10.0, 20.0, 30.0])
    expected = np.asarray([0.0, 0.135, 0.270, 0.570, 0.870])
    np.testing.assert_allclose(model.pressure(depth), expected, atol=1.0e-14)
