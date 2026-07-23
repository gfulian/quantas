"""Tests for composed, tabulated, preset, and YAML Earth profiles."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.physics.earth_profiles import (
    EarthProfileModel,
    LinearThermalBoundaryLayer,
    PiecewiseTemperatureModel,
    PremPressureModel,
    TabulatedTemperatureModel,
    TemperatureSegment,
    build_earth_profile_preset,
    earth_profile_presets,
    read_earth_profile_spec,
)

pytestmark = [pytest.mark.physics, pytest.mark.scientific, pytest.mark.fast]


def test_composed_profile_preserves_depths_and_provenance() -> None:
    """Generated grids retain source knots and recursively record both models."""
    temperature = TabulatedTemperatureModel(
        np.asarray([0.0, 35.0, 100.0]),
        np.asarray([288.15, 800.0, 1500.0]),
        citation="User reference.",
    )
    model = EarthProfileModel(PremPressureModel(), temperature, name="custom")
    profile = model.regular_profile(depth_step_km=20.0)
    assert 35.0 in profile.depth
    assert profile.depth.dtype == np.float64
    assert profile.metadata["pressure"]["model"] == "prem"
    assert profile.metadata["temperature"]["citation"] == "User reference."


def test_piecewise_continuous_offset_records_applied_transform() -> None:
    """A declared offset removes a join jump and remains visible in metadata."""
    first = LinearThermalBoundaryLayer(
        depth_top_km=0.0,
        depth_bottom_km=50.0,
        temperature_top_K=300.0,
        temperature_bottom_K=1000.0,
    )
    second = LinearThermalBoundaryLayer(
        depth_top_km=50.0,
        depth_bottom_km=100.0,
        temperature_top_K=1200.0,
        temperature_bottom_K=1700.0,
    )
    model = PiecewiseTemperatureModel(
        (
            TemperatureSegment(0.0, 50.0, first),
            TemperatureSegment(
                50.0,
                100.0,
                second,
                join_mode="continuous-offset",
            ),
        )
    )
    values = model.temperature(np.asarray([49.999999, 50.0, 50.000001]))
    assert np.max(np.abs(np.diff(values))) < 1.0e-3
    assert model.metadata()["segments"][1]["applied_offset_K"] == pytest.approx(-200.0)


def test_yaml_spec_combines_prem_with_user_temperature_table(tmp_path: Path) -> None:
    """A profile specification independently composes PREM and a table."""
    table = tmp_path / "temperature.dat"
    table.write_text(
        "depth_km T_K\n0 288.15\n35 800\n100 1500\n",
        encoding="utf-8",
    )
    specification = tmp_path / "profile.yaml"
    specification.write_text(
        """schema_version: 1
name: prem-custom-temperature
depth:
  min_km: 0
  max_km: 100
  step_km: 10
pressure:
  model: prem
temperature:
  source: table
  file: temperature.dat
  interpolation: pchip
  citation: User temperature profile.
""",
        encoding="utf-8",
    )
    profile = read_earth_profile_spec(specification)
    assert profile.name == "prem-custom-temperature"
    assert 35.0 in profile.depth
    assert profile.pressure[-1] > profile.pressure[1] > 0.0
    assert profile.temperature[-1] == pytest.approx(1500.0)
    assert profile.metadata["profile_schema_version"] == 1


def test_presets_are_bounded_and_referenced() -> None:
    """Every scientific preset has explicit parameters and full citations."""
    presets = earth_profile_presets()
    names = {preset.name for preset in presets}
    assert {
        "continental-cratonic",
        "continental-reference",
        "continental-active",
        "oceanic-10ma",
        "oceanic-50ma",
        "oceanic-100ma",
        "mantle-katsura-2022",
    } <= names
    for preset in presets:
        assert preset.depth_max_km > preset.depth_min_km
        assert preset.references
        assert all("https://doi.org/" in reference for reference in preset.references)

    mantle = build_earth_profile_preset("mantle-katsura-2022")
    assert mantle.depth[0] == 50.0
    assert mantle.depth[-1] == 2800.0
    assert mantle.temperature[0] == pytest.approx(1646.0)
    assert mantle.temperature[-1] == pytest.approx(2587.0)
    assert mantle.metadata["pressure"]["model"] == "prem"
