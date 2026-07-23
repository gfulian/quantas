from __future__ import annotations

import numpy as np

from quantas.core.math.fitting import (
    FitQuality,
    FitStatus,
    least_squares_fit,
    residual_metrics,
    validate_xy,
)


def _linear_model(x, a, b):
    return a * x + b


def test_validate_xy_accepts_matching_finite_arrays():
    x, y = validate_xy([1, 2, 3], [2, 4, 6])
    assert x.dtype == np.float64
    assert y.dtype == np.float64
    np.testing.assert_allclose(x, np.array([1.0, 2.0, 3.0]))


def test_validate_xy_accepts_multicoordinate_arrays():
    x, y = validate_xy([[1, 2, 3], [4, 5, 6]], [7, 8, 9])
    assert x.shape == (2, 3)
    assert y.shape == (3,)


def test_validate_xy_rejects_invalid_arrays():
    with np.testing.assert_raises(ValueError):
        validate_xy([[[1, 2]]], [1, 2])
    with np.testing.assert_raises(ValueError):
        validate_xy([1, 2], [1])
    with np.testing.assert_raises(ValueError):
        validate_xy([1, np.nan], [1, 2])


def test_residual_metrics_reports_expected_values():
    metrics = residual_metrics([1.0, 2.0, 3.0], [1.0, 2.5, 2.5])
    assert np.isclose(metrics["rmse"], np.sqrt((0.0**2 + (-0.5) ** 2 + 0.5**2) / 3.0))
    assert np.isclose(metrics["mae"], 1.0 / 3.0)
    assert np.isclose(metrics["max_abs_error"], 0.5)
    assert metrics["r_squared"] < 1.0


def test_least_squares_fit_returns_diagnostics_for_successful_fit():
    x = np.linspace(0.0, 10.0, 11)
    y = _linear_model(x, 2.0, -1.0)

    result = least_squares_fit(_linear_model, x, y, p0=[1.0, 0.0])

    assert result.success
    assert result.status is FitStatus.SUCCESS
    assert result.quality in {FitQuality.GOOD, FitQuality.POOR}
    np.testing.assert_allclose(result.parameters, np.array([2.0, -1.0]), atol=1.0e-12)
    assert result.residuals is not None
    np.testing.assert_allclose(result.residuals, np.zeros_like(x), atol=1.0e-12)
    assert result.rmse is not None and result.rmse < 1.0e-12
    assert result.n_points == x.size
    assert result.n_parameters == 2
    assert result.dof == x.size - 2
    assert result.as_dict()["status"] == "success"


def test_least_squares_fit_returns_failed_result_when_optimizer_fails():
    x = np.linspace(0.0, 1.0, 5)
    y = np.ones_like(x)

    def bad_model(values, a):
        raise RuntimeError("model evaluation failed")

    result = least_squares_fit(bad_model, x, y, p0=[1.0])

    assert not result.success
    assert result.status is FitStatus.FAILED
    assert result.quality is FitQuality.FAILED
    assert "model evaluation failed" in result.message


def test_least_squares_fit_marks_too_few_points_as_invalid_input():
    result = least_squares_fit(_linear_model, [1.0], [2.0], p0=[1.0, 0.0])

    assert not result.success
    assert result.status is FitStatus.INVALID_INPUT
    assert result.quality is FitQuality.FAILED
