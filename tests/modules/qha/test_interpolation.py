from __future__ import annotations

import numpy as np

from quantas.core.math.fitting import FitStatus
from quantas.modules.qha.core.interpolation import (
    InterpolationStatus,
    evaluate_polynomial_series,
    fit_polynomial_series,
    interpolate_frequency_grid,
    interpolate_polynomial_series,
    interpolate_property_grid,
    interpolate_static_energy,
    validate_volume_grid,
)


def test_validate_volume_grid_rejects_duplicate_values():
    with np.testing.assert_raises(ValueError):
        validate_volume_grid([10.0, 11.0, 11.0])


def test_fit_polynomial_series_fits_scalar_dataset():
    volume = np.linspace(8.0, 12.0, 5)
    energy = 2.0 + 0.25 * (volume - 10.0) ** 2

    fit = fit_polynomial_series(volume, energy, degree=2)

    assert fit.success
    assert fit.status is InterpolationStatus.SUCCESS
    assert fit.coefficients is not None
    assert fit.coefficients.shape == (3,)
    assert fit.fits[()].status is FitStatus.SUCCESS
    np.testing.assert_allclose(
        np.polynomial.polynomial.polyval(volume, fit.coefficients),
        energy,
        atol=1.0e-12,
    )


def test_fit_polynomial_series_fits_frequency_grid_along_last_axis():
    volume = np.array([9.0, 10.0, 11.0, 12.0])
    frequencies = np.empty((2, 3, volume.size), dtype=float)
    for iq in range(2):
        for band in range(3):
            frequencies[iq, band] = 100.0 + 10.0 * iq + band + 0.5 * volume

    fit = fit_polynomial_series(volume, frequencies, degree=1, axis=-1)

    assert fit.success
    assert fit.coefficients is not None
    assert fit.coefficients.shape == (2, 3, 2)
    assert fit.failed_indices == []
    np.testing.assert_allclose(fit.coefficients[..., 1], 0.5, atol=1.0e-12)


def test_evaluate_polynomial_series_appends_point_axis():
    coefficients = np.array(
        [
            [1.0, 2.0],
            [3.0, 4.0],
        ]
    )

    result = evaluate_polynomial_series(np.array([10.0, 20.0]), coefficients)

    assert result.success
    assert result.values is not None
    assert result.values.shape == (2, 2)
    np.testing.assert_allclose(result.values, np.array([[21.0, 41.0], [43.0, 83.0]]))


def test_interpolate_polynomial_series_returns_fit_and_values():
    volume = np.linspace(8.0, 12.0, 5)
    values = np.vstack((volume**2, 2.0 * volume + 1.0))

    fit, result = interpolate_polynomial_series(
        volume, values, np.array([9.5, 10.5]), degree=2, axis=-1
    )

    assert fit.success
    assert result.success
    assert result.values is not None
    np.testing.assert_allclose(
        result.values[0], np.array([9.5**2, 10.5**2]), atol=1.0e-12
    )
    np.testing.assert_allclose(result.values[1], np.array([20.0, 22.0]), atol=1.0e-12)


def test_interpolate_static_energy_uses_single_fit_result():
    volume = np.linspace(8.0, 12.0, 5)
    static_energy = (volume - 10.0) ** 2

    fit, result = interpolate_static_energy(volume, static_energy, 10.5, degree=2)

    assert fit.success
    assert result.success
    assert result.values is not None
    np.testing.assert_allclose(float(result.values), 0.25, atol=1.0e-12)


def test_interpolate_frequency_grid_warns_on_negative_values():
    volume = np.linspace(8.0, 12.0, 5)
    frequencies = np.array([[[volume - 11.0]]])

    fit, result = interpolate_frequency_grid(volume, frequencies, 10.0, degree=1)

    assert fit.success
    assert result.success
    assert result.metadata["contains_negative_values"] is True
    assert result.warnings


def test_interpolate_property_grid_accepts_volume_axis_zero():
    volume = np.array([9.0, 10.0, 11.0, 12.0])
    properties = np.column_stack((volume, volume**2))

    fit, result = interpolate_property_grid(volume, properties, 10.5, degree=2, axis=0)

    assert fit.success
    assert result.success
    assert result.values is not None
    np.testing.assert_allclose(result.values, np.array([10.5, 10.5**2]), atol=1.0e-12)


def test_fit_polynomial_series_reports_invalid_axis_length():
    volume = np.array([1.0, 2.0, 3.0])
    values = np.ones((2, 2))

    fit = fit_polynomial_series(volume, values, degree=1, axis=-1)

    assert not fit.success
    assert fit.status is InterpolationStatus.INVALID_INPUT
