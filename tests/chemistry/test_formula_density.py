"""Tests for chemical formulas, molar masses and crystal densities."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.chemistry import (
    convert_density,
    density_from_formula,
    molar_mass,
    parse_formula,
)


def test_parse_formula_handles_nested_groups_and_counts():
    composition = parse_formula("CaMg(CO3)2")

    assert composition == {"Ca": 1.0, "Mg": 1.0, "C": 2.0, "O": 6.0}


def test_parse_formula_handles_hydrate_dot_notation():
    composition = parse_formula("CuSO4·5H2O")

    assert composition["Cu"] == 1.0
    assert composition["S"] == 1.0
    assert composition["O"] == 9.0
    assert composition["H"] == 10.0


def test_molar_mass_returns_expected_value_for_mgo():
    np.testing.assert_allclose(molar_mass("MgO"), 40.304, atol=1.0e-6)


def test_density_from_formula_uses_z_formula_units():
    density = density_from_formula(
        "MgO",
        74.7,
        z=4,
        volume_unit="angstrom",
        density_unit="g cm^-3",
    )

    np.testing.assert_allclose(density, 3.584, rtol=5.0e-4)


def test_density_from_full_cell_formula_matches_z_formula_units():
    formula_unit = density_from_formula(
        "MgO",
        74.7,
        z=4,
        volume_unit="angstrom",
        density_unit="kg m^-3",
    )
    full_cell = density_from_formula(
        "Mg4O4",
        74.7,
        volume_unit="angstrom",
        density_unit="kg m^-3",
    )

    np.testing.assert_allclose(full_cell, formula_unit)


def test_density_conversion_between_supported_units():
    assert convert_density(1.0, "g cm^-3", "kg m^-3") == 1000.0
    assert convert_density(1000.0, "kg/m^3", "g/cm^3") == 1.0


def test_formula_and_density_reject_invalid_input():
    with pytest.raises(ValueError):
        parse_formula("Mg(OH")
    with pytest.raises(ValueError):
        parse_formula("Xx2")
    with pytest.raises(ValueError):
        density_from_formula("MgO", 0.0)
