"""Tests for model-oriented fitting infrastructure."""

from __future__ import annotations

import numpy as np

from quantas.core.math.fitting import (
    BaseFitModel,
    FitStatus,
    FittedModel,
    LeastSquaresFitter,
)


class LinearModel(BaseFitModel):
    """Simple linear model used to exercise the fitter contract."""

    @property
    def name(self) -> str:
        return "linear"

    @property
    def parameter_names(self) -> tuple[str, ...]:
        return ("slope", "intercept")

    def evaluate(self, x, parameters):
        slope, intercept = np.asarray(parameters, dtype=np.float64)
        return slope * np.asarray(x, dtype=np.float64) + intercept

    def initial_guess(self, x, y):
        return np.asarray([1.0, 0.0], dtype=np.float64)


def test_least_squares_fitter_operates_on_model_contract():
    x = np.linspace(-2.0, 2.0, 9)
    y = 3.5 * x - 0.25

    result = LeastSquaresFitter().fit(LinearModel(), x, y)

    assert result.success
    assert result.status is FitStatus.SUCCESS
    np.testing.assert_allclose(result.parameters, [3.5, -0.25], atol=1.0e-12)
    assert result.metadata["model"] == "linear"
    assert result.metadata["parameter_order"] == ["slope", "intercept"]


def test_fitted_model_reuses_optimized_parameters():
    x = np.linspace(0.0, 4.0, 5)
    y = 2.0 * x + 1.0
    model = LinearModel()
    result = LeastSquaresFitter().fit(model, x, y)

    fitted = FittedModel(model, result)

    np.testing.assert_allclose(fitted.evaluate([5.0, 6.0]), [11.0, 13.0])
