"""Systematic validation of QHA workflow combinations and table formatting."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.modules.qha.calculator import QHACalculator
from quantas.cli.qha_render import render_table
from quantas.modules.qha.io.reader import read_qha_input
from quantas.modules.qha.models import QHAOptions, QHAResult
from quantas.modules.qha.report import selected_property_tables
from quantas.core.physics.units import convert_energy_per_temperature
from quantas.modules.qha.validation import compare_qha_results, validate_qha_result


@pytest.fixture(scope="module")
def workflow_results() -> tuple[object, dict[tuple[str, str], QHAResult]]:
    """Return MgO results for the four QHA workflow combinations."""
    input_data = read_qha_input(Path(__file__).parent / "data" / "mgo_b3lyp_qha.yaml")
    results: dict[tuple[str, str], QHAResult] = {}
    for scheme in ("freq", "td"):
        for minimization in ("poly", "eos"):
            options = QHAOptions(
                temperature_min=0.0,
                temperature_max=2000.0,
                temperature_step=1000.0,
                pressure_min=0.0,
                pressure_max=30.0,
                pressure_step=15.0,
                scheme=scheme,
                minimization=minimization,
                eos="BM",
                energy_degree=3,
                free_energy_degree=3,
                frequency_degree=3,
                polynomial_derivative_method="local_grid",
                estimate_uncertainties=minimization == "eos",
                uncertainty_method="covariance",
            )
            results[(scheme, minimization)] = (
                QHACalculator(
                    input_data,
                    options,
                )
                .execute()
                .results["qha"]
            )
    return input_data, results


def test_all_four_workflows_pass_basic_physical_validation(workflow_results) -> None:
    input_data, results = workflow_results

    for result in results.values():
        summary = validate_qha_result(result, input_data)
        assert summary.completed
        assert summary.valid_points == summary.total_points
        assert summary.finite_properties
        assert summary.volume_decreases_with_pressure
        assert summary.positive_bulk_moduli
        assert summary.zero_kelvin_consistency
        assert summary.cp_not_below_cv
        assert summary.dulong_petit_ratio is not None
        assert 0.90 <= summary.dulong_petit_ratio <= 1.05


def test_frequency_and_td_schemes_share_equilibrium_volumes(workflow_results) -> None:
    _, results = workflow_results

    for minimization in ("poly", "eos"):
        np.testing.assert_allclose(
            results[("freq", minimization)].equilibrium_volume,
            results[("td", minimization)].equilibrium_volume,
            rtol=0.0,
            atol=0.0,
        )


def test_eos_thermoelastic_properties_are_scheme_independent(workflow_results) -> None:
    _, results = workflow_results
    differences = {
        item.property_name: item
        for item in compare_qha_results(
            results[("freq", "eos")], results[("td", "eos")]
        )
    }

    assert differences["isothermal_bulk_modulus"].maximum_absolute == 0.0
    assert differences["bulk_modulus_derivative"].maximum_absolute == 0.0
    assert differences["free_energy"].maximum_absolute > 0.0


def test_polynomial_extrapolation_is_reported_by_validation(workflow_results) -> None:
    input_data, results = workflow_results
    summary = validate_qha_result(results[("freq", "poly")], input_data)

    assert summary.volumes_above_sampled_range > 0
    assert any("outside the sampled interval" in item for item in summary.warnings)


def test_compact_report_formats_and_molar_cp() -> None:
    result = QHAResult(
        temperature=np.array([300.0]),
        pressure=np.array([0.0]),
        equilibrium_volume=np.array([[19.123456789]]),
        isothermal_bulk_modulus=np.array([[160.1234567]]),
        bulk_modulus_derivative=np.array([[4.1234567]]),
        adiabatic_bulk_modulus=np.array([[161.1234567]]),
        thermal_expansion=np.array([[3.1234567e-5]]),
        heat_capacity_difference=np.array([[1.0e-6]]),
        uncertainties={
            "sigma_VT": np.array([[0.0001234567]]),
            "sigma_KT": np.array([[0.01234567]]),
            "sigma_Kp": np.array([[0.001234567]]),
        },
        metadata={
            "units": {
                "temperature": "K",
                "pressure": "GPa",
                "volume": "A",
                "energy": "Ha",
                "heat_capacity": "Ha cell^-1 K^-1",
            },
            "normalization": {"formula_units_per_cell": 1},
        },
    )

    table = selected_property_tables(result)[0]
    text = render_table(table)

    assert table.rows[0][0] == "300.00"
    assert table.rows[0][1] == "19.123457"
    assert table.rows[0][2] == "0.000123"
    assert table.rows[0][3] == "160.12346"
    assert table.rows[0][5] == "4.12346"
    assert table.rows[0][8] == "3.123457"
    expected = float(
        convert_energy_per_temperature(
            1.0e-6,
            "Ha cell^-1 K^-1",
            "J mol^-1 K^-1",
        )
    )
    assert table.rows[0][9] == f"{expected:.12E}"
    assert "(A^3)" in text
    assert "(J mol^-1 K^-1)" in text
