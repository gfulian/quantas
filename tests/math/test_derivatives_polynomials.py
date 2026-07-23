import numpy as np

from quantas.core.math.derivative import (
    derivative,
    polynomial_derivative,
    polynomial_derivative_from_coefficients,
    relative_derivative,
)
from quantas.core.math.polynomials import (
    R_squared,
    find_polynomial_minimum,
    interpolate,
    polyfit,
)


def test_derivative_of_quadratic_function():
    x = np.linspace(0.0, 10.0, 101)
    y = x**2

    dy = derivative(x, y)

    np.testing.assert_allclose(dy[1:-1], 2.0 * x[1:-1], rtol=1e-2)


def test_polyfit_single_dataset():
    x = np.linspace(-2.0, 2.0, 20)
    y = 1.0 + 2.0 * x + 3.0 * x**2

    pars = polyfit(x, y, degree=2)

    np.testing.assert_allclose(pars, [1.0, 2.0, 3.0], atol=1e-10)


def test_interpolate_single_polynomial():
    pars = np.array([1.0, 2.0, 3.0])
    result = interpolate(2.0, pars)

    assert result == 17.0


def test_r_squared_perfect_fit():
    x = np.linspace(-2.0, 2.0, 20)
    y = 1.0 + 2.0 * x + 3.0 * x**2
    pars = np.array([1.0, 2.0, 3.0])

    r2 = R_squared(x, y, pars)

    np.testing.assert_allclose(r2, 1.0)


def test_find_polynomial_minimum():
    pars = np.array([1.0, 0.0, 1.0])
    minima = find_polynomial_minimum(pars)

    np.testing.assert_allclose(minima, [0.0])


def test_relative_derivative_of_exponential_function():
    x = np.linspace(0.0, 2.0, 101)
    y = np.exp(3.0 * x)

    values = relative_derivative(x, y)

    np.testing.assert_allclose(values[2:-2], 3.0, rtol=2.0e-3)


def test_polynomial_derivative_fit_returns_analytical_derivative():
    x = np.linspace(-2.0, 2.0, 21)
    y = 1.0 - 4.0 * x + 3.0 * x**2
    x_eval = np.array([-1.5, 0.0, 1.5])

    values = polynomial_derivative(x, y, x_eval, degree=2)

    np.testing.assert_allclose(values, -4.0 + 6.0 * x_eval, atol=1.0e-12)


def test_polynomial_derivative_from_coefficients_uses_existing_fit():
    coefficients = np.array([1.0, -4.0, 3.0])
    x_eval = np.array([-1.5, 0.0, 1.5])

    values = polynomial_derivative_from_coefficients(coefficients, x_eval)

    np.testing.assert_allclose(values, -4.0 + 6.0 * x_eval, atol=1.0e-12)
