"""Scientific regression tests for terrestrial temperature-depth models."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.physics.earth_profiles import (
    ConductiveLayer,
    ContinentalConductiveGeotherm,
    Katsura2022MantleAdiabat,
    OceanicHalfSpaceGeotherm,
    OceanicPlateGeotherm,
)

pytestmark = [pytest.mark.physics, pytest.mark.scientific, pytest.mark.fast]


def test_single_layer_conduction_matches_analytic_solution() -> None:
    """The numerical evaluator follows the closed steady conductive solution."""
    model = ContinentalConductiveGeotherm(
        (ConductiveLayer(40.0, 2.5, 1.0, "crust"),),
        surface_temperature_K=280.0,
        surface_heat_flow_mW_m2=60.0,
    )
    depth = np.asarray([0.0, 10.0, 25.0, 40.0])
    z_m = depth * 1000.0
    expected = 280.0 + 0.060 * z_m / 2.5 - 1.0e-6 * z_m**2 / 5.0
    np.testing.assert_allclose(model.temperature(depth), expected)


def test_layered_conduction_is_continuous_at_interfaces() -> None:
    """Temperature remains continuous when conductivity and heat production change."""
    model = ContinentalConductiveGeotherm(
        (
            ConductiveLayer(20.0, 2.5, 1.0, "upper crust"),
            ConductiveLayer(30.0, 3.0, 0.2, "lower crust"),
        ),
        surface_heat_flow_mW_m2=55.0,
    )
    values = model.temperature(np.asarray([19.999999, 20.0, 20.000001]))
    assert np.max(np.abs(np.diff(values))) < 1.0e-3


def test_half_space_cooling_boundary_and_age_dependence() -> None:
    """The surface is fixed and older oceanic lithosphere is cooler at depth."""
    young = OceanicHalfSpaceGeotherm(age_Ma=10.0)
    old = OceanicHalfSpaceGeotherm(age_Ma=100.0)
    assert young.temperature(np.asarray([0.0]))[0] == pytest.approx(273.15)
    assert (
        old.temperature(np.asarray([40.0]))[0]
        < young.temperature(np.asarray([40.0]))[0]
    )


def test_plate_cooling_satisfies_surface_and_basal_boundaries() -> None:
    """Finite plate cooling preserves both imposed thermal boundaries."""
    model = OceanicPlateGeotherm(age_Ma=50.0, plate_thickness_km=125.0)
    values = model.temperature(np.asarray([0.0, 125.0]))
    np.testing.assert_allclose(values, np.asarray([273.15, 1623.15]))


def test_katsura_matches_published_constraints() -> None:
    """The deterministic reconstruction reproduces every published node."""
    model = Katsura2022MantleAdiabat()
    shallow = np.asarray([50.0, 409.999999, 519.999999, 659.999999, 2400.0, 2800.0])
    expected_shallow = np.asarray([1646.0, 1799.0, 1899.0, 1994.0, 2490.0, 2587.0])
    np.testing.assert_allclose(
        model.temperature(shallow), expected_shallow, atol=5.0e-5
    )
    deep = model.temperature(np.asarray([410.0, 520.0, 660.0]))
    np.testing.assert_allclose(deep, np.asarray([1860.0, 1942.0, 1960.0]))


def test_katsura_metadata_records_provenance() -> None:
    """The mantle profile carries both primary and software citations."""
    metadata = Katsura2022MantleAdiabat().metadata()
    assert metadata["doi"] == "10.1029/2021JB023562"
    assert metadata["software_doi"] == "10.5281/zenodo.5903286"
    assert "Journal of Geophysical Research" in metadata["citation"]
    assert "Version 1.1.0" in metadata["software_citation"]
