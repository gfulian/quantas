"""Tests for QHA data containers."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.math.fitting import FitQuality, FitResult, FitStatus
from quantas.modules.qha.models import (
    QHAFailedPoint,
    QHAFitRecord,
    QHAInput,
    QHAOptions,
    QHAResult,
    diagnostics_to_dict,
    merge_uncertainties,
    result_metadata_from_options,
)


def test_qha_input_reports_dimensions_and_normalized_weights() -> None:
    """QHAInput exposes dimensions and normalized q-point weights."""
    data = QHAInput(
        jobname="MgO",
        natoms=2,
        supercell=np.eye(3) * 2.0,
        qpoints=2,
        volume=np.array([70.0, 75.0, 80.0]),
        energy=np.array([-100.0, -101.0, -100.5]),
        frequencies=np.ones((2, 6, 3)),
        weights=np.array([1.0, 3.0]),
    )

    assert data.nvol == 3
    assert data.nmodes == 6
    assert data.kpoints == 8
    assert data.total_q_points == pytest.approx(4.0)
    np.testing.assert_allclose(data.normalized_weights(), [0.25, 0.75])
    assert data.has_phonons()
    assert data.has_mode_continuity()


def test_qha_input_validate_shapes_rejects_inconsistent_frequency_axis() -> None:
    """QHAInput validates the frequency volume axis."""
    data = QHAInput(
        qpoints=1,
        volume=np.array([70.0, 75.0, 80.0]),
        energy=np.array([-100.0, -101.0, -100.5]),
        frequencies=np.ones((1, 6, 2)),
        weights=np.array([1.0]),
    )

    with pytest.raises(ValueError, match="frequency volume axis"):
        data.validate_shapes()


def test_qha_input_validate_shapes_accepts_consistent_arrays() -> None:
    """QHAInput accepts consistent volume, energy, phonon, and q-point arrays."""
    data = QHAInput(
        qpoints=2,
        volume=np.array([70.0, 75.0, 80.0]),
        energy=np.array([-100.0, -101.0, -100.5]),
        frequencies=np.ones((2, 6, 3)),
        weights=np.array([1.0, 1.0]),
        qcoords=np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]),
    )

    data.validate_shapes()


def test_qha_options_builds_temperature_and_pressure_grids() -> None:
    """QHAOptions builds inclusive temperature and pressure grids."""
    options = QHAOptions(
        temperature_min=100.0,
        temperature_max=300.0,
        temperature_step=100.0,
        pressure_min=0.0,
        pressure_max=2.0,
        pressure_step=1.0,
    )

    np.testing.assert_allclose(options.temperature_grid(), [100.0, 200.0, 300.0])
    np.testing.assert_allclose(options.pressure_grid(), [0.0, 1.0, 2.0])
    assert options.requires_mode_continuity()
    options.validate()


def test_qha_options_validation_rejects_invalid_settings() -> None:
    """QHAOptions rejects unsupported schemes and invalid failure limits."""
    with pytest.raises(ValueError, match="scheme"):
        QHAOptions(scheme="bad").validate()  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="max_consecutive_failures"):
        QHAOptions(max_consecutive_failures=0).validate()


def test_fit_record_serializes_fit_result() -> None:
    """QHAFitRecord serializes nested FitResult diagnostics."""
    fit = FitResult(
        success=True,
        status=FitStatus.SUCCESS,
        quality=FitQuality.GOOD,
        message="ok",
        parameters=np.array([1.0, 2.0]),
        n_points=5,
        n_parameters=2,
        dof=3,
    )
    record = QHAFitRecord(
        quantity="F",
        method="poly",
        temperature=300.0,
        pressure=1.0,
        fit=fit,
    )

    payload = record.as_dict()

    assert payload["quantity"] == "F"
    assert payload["fit"]["status"] == "success"
    assert payload["fit"]["parameters"] == [1.0, 2.0]


def test_failed_point_serializes_diagnostics() -> None:
    """QHAFailedPoint keeps pressure-temperature failure metadata."""
    point = QHAFailedPoint(
        temperature=500.0,
        pressure=10.0,
        stage="minimization",
        message="minimum outside volume range",
        diagnostics={"volume": 100.0},
    )

    payload = point.as_dict()

    assert payload["stage"] == "minimization"
    assert payload["diagnostics"]["volume"] == 100.0


def test_qha_result_exposes_property_symbols_and_uncertainties() -> None:
    """QHAResult maps property arrays to historical output names."""
    result = QHAResult(
        equilibrium_volume=np.array([[70.0, 69.0]]),
        isothermal_bulk_modulus=np.array([[160.0, 170.0]]),
        gibbs_free_energy=np.array([[-100.0, -99.0]]),
        uncertainties={"sigma_VT": np.array([[0.1, 0.2]])},
    )

    symbol = result.as_property_dict()

    assert set(symbol) == {"VT", "KT", "G"}
    np.testing.assert_allclose(result.uncertainty_dict()["sigma_VT"], [[0.1, 0.2]])
    assert result.has_thermodynamic_data()


def test_qha_result_stores_diagnostics_and_failed_points() -> None:
    """QHAResult stores local diagnostics and failed pressure-temperature points."""
    result = QHAResult(completed=False)
    result.add_fit_record(QHAFitRecord(quantity="F", method="poly"))
    result.add_failed_point(
        QHAFailedPoint(temperature=300.0, pressure=0.0, stage="fit", message="fail")
    )

    assert len(result.diagnostics_as_dict()) == 1
    assert len(result.failed_points_as_dict()) == 1
    assert not result.completed


def test_qha_result_partial_copy_filters_grid_arrays() -> None:
    """QHAResult can filter pressure-temperature arrays using a valid mask."""
    mask = np.array([[True, False], [False, True]])
    result = QHAResult(
        equilibrium_volume=np.array([[70.0, 71.0], [72.0, 73.0]]),
        isothermal_bulk_modulus=np.array([[160.0, 161.0], [162.0, 163.0]]),
        volume=np.array([68.0, 70.0, 72.0]),
        uncertainties={"sigma_VT": np.array([[0.1, 0.2], [0.3, 0.4]])},
    )

    partial = result.partial_copy(mask)

    np.testing.assert_allclose(partial.equilibrium_volume, [70.0, 73.0])
    np.testing.assert_allclose(partial.isothermal_bulk_modulus, [160.0, 163.0])
    np.testing.assert_allclose(partial.volume, [68.0, 70.0, 72.0])
    np.testing.assert_allclose(partial.uncertainties["sigma_VT"], [0.1, 0.4])


def test_diagnostics_to_dict_serializes_records() -> None:
    """diagnostics_to_dict converts records to serializable dictionaries."""
    records = [QHAFitRecord(quantity="KT", method="eos")]

    payload = diagnostics_to_dict(records)

    assert payload[0]["quantity"] == "KT"
    assert payload[0]["method"] == "eos"


def test_result_metadata_from_options_contains_workflow_information() -> None:
    """result_metadata_from_options returns scheme, method, units, and policies."""
    options = QHAOptions(
        scheme="td",
        minimization="eos",
        eos="Vinet",
        debug=True,
        max_consecutive_failures=3,
    )

    metadata = result_metadata_from_options(options)

    assert metadata["scheme"] == "td"
    assert metadata["minimization"] == "eos"
    assert metadata["eos"] == "Vinet"
    assert metadata["units"]["pressure"] == "GPa"
    assert metadata["diagnostics"]["debug"] is True
    assert metadata["diagnostics"]["max_consecutive_failures"] == 3


def test_merge_uncertainties_converts_values_to_arrays() -> None:
    """merge_uncertainties merges and converts uncertainty values."""
    merged = merge_uncertainties(
        {"sigma_VT": [0.1, 0.2]},
        {"sigma_KT": 1.5},
    )

    np.testing.assert_allclose(merged["sigma_VT"], [0.1, 0.2])
    assert merged["sigma_KT"].shape == ()
    assert merged["sigma_KT"] == pytest.approx(1.5)


def test_qha_input_reports_atoms_per_formula_unit() -> None:
    data = QHAInput(natoms=8, formula_units=4)

    assert data.natoms_per_formula_unit == 2.0
