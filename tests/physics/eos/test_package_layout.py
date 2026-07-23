"""Package-layout tests for shared equation-of-state physics."""

from __future__ import annotations

from quantas.core.physics import eos
from quantas.core.physics.eos import pvt, temperature, thermal_pressure


def test_shared_eos_public_api_is_available():
    assert eos.EnergyEOS is not None
    assert eos.PressureEOS is not None
    assert eos.FittedEnergyEOS is not None
    assert eos.EOSState is not None
    assert eos.EOSStateUncertainty is not None


def test_temperature_and_pvt_cores_are_public():
    assert "TemperatureEOS" in temperature.__all__
    assert "TemperatureEOSModel" in temperature.__all__
    assert "PVTEOS" in pvt.__all__
    assert "PVTModel" in pvt.__all__
    assert "PVTCouplingFamily" in pvt.__all__


def test_thermal_pressure_contracts_are_public():
    assert "ThermalPressureModel" in thermal_pressure.__all__
    assert "MGDNormalization" in thermal_pressure.__all__
    assert "MGDParameters" in thermal_pressure.__all__
