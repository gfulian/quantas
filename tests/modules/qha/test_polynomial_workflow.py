from __future__ import annotations

import numpy as np
import pytest

from quantas.modules.qha.core.minimization import (
    evaluate_fitted_polynomial_at_pressure,
    fit_polynomial_free_energy_model,
    minimize_polynomial,
)
from quantas.modules.qha.models import QHAInput, QHAOptions
from quantas.modules.qha.workflow import run_volume_minimization_workflow


def test_scaled_polynomial_model_reproduces_narrow_volume_curve() -> None:
    volume = np.linspace(99.5, 100.5, 9)
    x = volume - 100.0
    free_energy = -1000.0 + 0.25 * x**2 + 0.02 * x**3

    fitted = fit_polynomial_free_energy_model(volume, free_energy, degree=3)

    assert fitted.success
    assert fitted.model is not None
    np.testing.assert_allclose(
        fitted.model.free_energy(volume),
        free_energy,
        atol=1.0e-12,
    )
    assert fitted.fit.condition_number is not None
    assert fitted.fit.condition_number < 10.0
    assert fitted.fit.metadata["coordinate"] == "x=(V-center)/scale"


def test_polynomial_minimum_is_selected_from_gibbs_objective() -> None:
    # Two physical minima with a small tilt. The lower F minimum must be
    # selected independently of the absolute energy offset.
    base = np.polynomial.Polynomial.fromroots([9.0, 9.0, 11.0, 11.0])
    polynomial = base + 0.08 * np.polynomial.Polynomial([-10.0, 1.0]) - 500.0
    volume = np.linspace(8.0, 12.0, 17)
    free_energy = polynomial(volume)

    result = minimize_polynomial(
        volume,
        free_energy,
        parameters=polynomial.coef,
        derivative_method="analytic",
    )

    roots = polynomial.deriv().roots()
    real = np.real(roots[np.isclose(np.imag(roots), 0.0)])
    minima = real[polynomial.deriv(2)(real) > 0.0]
    expected = float(minima[np.argmin(polynomial(minima))])
    assert result.success
    assert result.volume == pytest.approx(expected)
    assert result.volume < 10.0


def test_derivatives_agree_for_exact_polynomial() -> None:
    volume = np.linspace(8.0, 12.0, 9)
    free_energy = 0.5 * (volume - 10.0) ** 2 + 0.01 * (volume - 10.0) ** 4
    fitted = fit_polynomial_free_energy_model(volume, free_energy, degree=4)

    analytic = evaluate_fitted_polynomial_at_pressure(
        fitted,
        0.0,
        derivative_method="analytic",
    )
    local = evaluate_fitted_polynomial_at_pressure(
        fitted,
        0.0,
        derivative_method="local_grid",
        local_grid_points=5,
        local_grid_separation=0.05,
        local_degree=4,
    )

    assert analytic.success and local.success
    assert local.volume == pytest.approx(analytic.volume)
    assert local.bulk_modulus == pytest.approx(analytic.bulk_modulus, rel=1.0e-7)
    assert local.bulk_modulus_derivative == pytest.approx(
        analytic.bulk_modulus_derivative,
        rel=1.0e-5,
        abs=1.0e-7,
    )


def test_local_free_energy_controls_derivatives() -> None:
    volume = np.linspace(8.0, 12.0, 9)
    free_energy = (volume - 10.0) ** 2
    fitted = fit_polynomial_free_energy_model(volume, free_energy, degree=2)
    calls: list[np.ndarray] = []

    def evaluator(local_volume: np.ndarray) -> np.ndarray:
        calls.append(local_volume.copy())
        return 2.0 * (local_volume - 10.0) ** 2

    result = evaluate_fitted_polynomial_at_pressure(
        fitted,
        0.0,
        derivative_method="local_grid",
        local_free_energy=evaluator,
        local_grid_points=5,
        local_grid_separation=0.05,
        local_degree=2,
    )

    assert result.success
    assert len(calls) == 1
    assert calls[0].size == 5
    assert result.volume == pytest.approx(10.0)
    assert result.bulk_modulus == pytest.approx(40.0, rel=1.0e-8)
    assert result.bulk_modulus_derivative == pytest.approx(-1.0, abs=1.0e-8)
    assert result.metadata["analytic_bulk_modulus"] == pytest.approx(20.0)


def test_polynomial_workflow_fits_once_per_temperature() -> None:
    volume = np.linspace(9.0, 11.0, 7)
    energy = 0.5 * (volume - 10.0) ** 2 - 5.0
    input_data = QHAInput(volume=volume, energy=energy)
    options = QHAOptions(
        temperature_min=300.0,
        temperature_max=400.0,
        temperature_step=100.0,
        pressure_min=0.0,
        pressure_max=2.0,
        pressure_step=1.0,
        minimization="poly",
        free_energy_degree=2,
        energy_degree=2,
        polynomial_derivative_method="analytic",
        debug=True,
    )

    result = run_volume_minimization_workflow(input_data, options)

    assert result.completed
    assert result.metadata["polynomial_workflow"]["fit_count"] == 2
    assert result.metadata["polynomial_workflow"]["state_count"] == 6
    assert len(result.fit_records) == 2
    assert all(record.quantity == "F" for record in result.fit_records)
    assert all(record.pressure is None for record in result.fit_records)


def test_workflow_local_grid_uses_external_free_energy_evaluator() -> None:
    volume = np.linspace(9.0, 11.0, 7)
    energy = (volume - 10.0) ** 2
    input_data = QHAInput(volume=volume, energy=energy)
    options = QHAOptions(
        temperature_min=300.0,
        temperature_max=300.0,
        pressure_min=0.0,
        pressure_max=0.0,
        minimization="poly",
        free_energy_degree=2,
        energy_degree=2,
        polynomial_derivative_method="local_grid",
    )
    calls: list[tuple[float, int, int]] = []

    def evaluator(
        local_volume: np.ndarray, temperature: float, index: int
    ) -> np.ndarray:
        calls.append((temperature, index, local_volume.size))
        return 2.0 * (local_volume - 10.0) ** 2

    result = run_volume_minimization_workflow(
        input_data,
        options,
        local_free_energy_evaluator=evaluator,
    )

    assert result.completed
    assert calls == [(300.0, 0, 5)]
    assert result.isothermal_bulk_modulus[0, 0] > 0.0


def test_polynomial_options_validate_local_grid_settings() -> None:
    with pytest.raises(ValueError, match="odd integer"):
        QHAOptions(polynomial_grid_points=4).validate()
    with pytest.raises(ValueError, match="positive"):
        QHAOptions(polynomial_grid_separation=0.0).validate()
    with pytest.raises(ValueError, match="polynomial_derivative_method"):
        QHAOptions(polynomial_derivative_method="invalid").validate()  # type: ignore[arg-type]


def test_frequency_evaluator_reproduces_sampled_free_energy() -> None:
    from quantas.modules.qha.thermodynamics import (
        FrequencyThermodynamicEvaluator,
        calculate_sampled_thermodynamics,
        free_energy_grid,
    )

    volume = np.linspace(9.5, 10.5, 5)
    energy = 1.0e-3 * (volume - 10.0) ** 2
    base = 400.0 * (10.0 / volume) ** 2
    frequencies = np.stack([base, 1.2 * base, 1.4 * base], axis=0)[None, :, :]
    input_data = QHAInput(
        natoms=1,
        qpoints=1,
        volume=volume,
        energy=energy,
        frequencies=frequencies,
        weights=np.array([1.0]),
    )
    options = QHAOptions(
        temperature_min=300.0,
        temperature_max=300.0,
        temperature_step=1.0,
        scheme="freq",
        energy_degree=4,
        frequency_degree=4,
    )

    sampled = free_energy_grid(calculate_sampled_thermodynamics(input_data, options))
    evaluator = FrequencyThermodynamicEvaluator(input_data, options)
    evaluated = evaluator.free_energy_at(volume, 300.0)

    np.testing.assert_allclose(evaluated, sampled[0], rtol=1.0e-10, atol=1.0e-12)
