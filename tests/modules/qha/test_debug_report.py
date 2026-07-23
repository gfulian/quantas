"""Tests for QHA debug and final text reporting."""

from __future__ import annotations

import numpy as np

from quantas.modules.qha.core.minimization import minimize_polynomial
from quantas.modules.qha.models import QHAResult
from quantas.modules.qha.report import (
    debug_thermodynamic_property_tables,
    selected_property_tables,
)


def test_polynomial_minimization_reports_bulk_modulus_and_derivative() -> None:
    """Polynomial minimization exposes KT and Kp estimates."""
    volume = np.linspace(9.0, 11.0, 7)
    energy = 0.5 * (volume - 10.0) ** 2 + 0.01 * (volume - 10.0) ** 3

    minimum = minimize_polynomial(volume, energy, degree=3)

    assert minimum.success
    assert minimum.bulk_modulus is not None
    assert minimum.bulk_modulus > 0.0
    assert minimum.bulk_modulus_derivative is not None
    assert np.isfinite(minimum.bulk_modulus_derivative)


def test_selected_property_tables_use_standard_qha_columns() -> None:
    """Final QHA tables use the compact pressure-grouped layout."""
    result = _fake_result()

    tables = selected_property_tables(result)

    assert len(tables) == 2
    assert tables[0].columns == [
        "T",
        "V",
        "KT",
        "KS",
        "alphaV x 10^5",
        "Cp-Cv",
    ]
    assert tables[0].rows[0][0] == "0.00"
    assert tables[0].rows[0][1] == "10.000000"
    assert tables[0].metadata["column_units"][-1] == "J mol^-1 K^-1"


def test_debug_thermodynamic_tables_include_requested_properties() -> None:
    """Debug thermodynamic tables expose detailed values in scientific format."""
    result = _fake_result()

    tables = debug_thermodynamic_property_tables(result)

    assert len(tables) == 2
    assert "U_zp" in tables[0].columns
    assert "C_P-C_V" in tables[0].columns
    assert "alpha_V" in tables[0].columns
    assert "E" in tables[0].rows[1][1]


def _fake_result() -> QHAResult:
    """Return a small QHA result with all final-report properties."""
    shape = (2, 2)
    return QHAResult(
        temperature=np.array([0.0, 300.0]),
        pressure=np.array([0.0, 1.0]),
        equilibrium_volume=np.array([[10.0, 9.8], [10.2, 10.0]]),
        zero_point_energy=np.full(shape, 1.0),
        thermal_energy=np.full(shape, 2.0),
        internal_energy=np.full(shape, 3.0),
        entropy=np.full(shape, 4.0),
        isochoric_heat_capacity=np.full(shape, 5.0),
        isobaric_heat_capacity=np.full(shape, 6.0),
        heat_capacity_difference=np.full(shape, 1.0e-4),
        vibrational_free_energy=np.full(shape, 7.0),
        free_energy=np.full(shape, 8.0),
        enthalpy=np.full(shape, 9.0),
        gibbs_free_energy=np.full(shape, 10.0),
        isothermal_bulk_modulus=np.full(shape, 100.0),
        adiabatic_bulk_modulus=np.full(shape, 101.0),
        thermal_expansion=np.full(shape, 1.0e-5),
        metadata={"units": {"energy": "Ha"}},
    )
