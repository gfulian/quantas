"""Tests for reusable scaled polynomial models."""

from __future__ import annotations

import numpy as np

from quantas.core.math.polynomials import (
    FittedPolynomial,
    fit_polynomial,
    fit_polynomial_result,
)


def test_scaled_polynomial_preserves_values_and_physical_derivatives():
    x = np.linspace(10.0, 12.0, 9)
    y = 4.0 + 2.0 * x + 3.0 * x**2

    result, model = fit_polynomial(x, y, 2, scale_coordinate=True)

    assert result.success
    assert model is not None
    np.testing.assert_allclose(model.evaluate(x), y, atol=1.0e-11)
    np.testing.assert_allclose(model.derivative(x, 1), 2.0 + 6.0 * x, atol=1.0e-10)
    np.testing.assert_allclose(model.derivative(x, 2), 6.0, atol=1.0e-10)


def test_scaled_polynomial_finds_minimum_with_derivative_addendum():
    model = FittedPolynomial(
        [0.0, 0.0, 1.0],
        center=5.0,
        scale=2.0,
        sampled_coordinates=[3.0, 4.0, 5.0, 6.0, 7.0],
    )

    minima = model.local_minima(derivative_addendum=0.5)

    # f = ((x - 5) / 2)^2, hence df/dx = (x - 5) / 2.
    np.testing.assert_allclose(minima, [4.0])


def test_fit_polynomial_result_returns_only_fit_diagnostics():
    x = np.linspace(-1.0, 1.0, 7)
    y = 1.0 + 2.0 * x + 3.0 * x**2

    result = fit_polynomial_result(x, y, 2)

    assert result.success
    np.testing.assert_allclose(result.parameters, [1.0, 2.0, 3.0], atol=1.0e-12)
    assert result.metadata["model"] == "polynomial"
