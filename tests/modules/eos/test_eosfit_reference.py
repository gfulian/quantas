"""External BM3 regression checks against EosFit7 reference output."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.math.fitting import (
    CovarianceScaling,
    EffectiveVarianceOptions,
)
from quantas.modules.eos import (
    EOSFitOptions,
    EOSFitRequest,
    EOSFitter,
    read_eos_input,
)

_DATA = Path(__file__).with_name("data")


def _parameter_errors(result) -> dict[str, float]:
    """Return named parameter standard errors from a successful fit."""
    assert result.fit.errors is not None
    return {
        name: float(error)
        for name, error in zip(
            result.fit.parameter_names,
            result.fit.errors,
            strict=True,
        )
    }


@pytest.mark.parametrize(
    ("filename", "expected", "expected_errors", "max_residual"),
    (
        (
            "PV_quartz.dat",
            {"V0": 112.96752, "K0": 37.28543, "KP": 5.93351, "KPP": -0.25642},
            {"V0": 0.02249, "K0": 0.21292, "KP": 0.06570},
            -0.02,
        ),
        (
            "PV_topaz.dat",
            {"V0": 346.97214, "K0": 135.74806, "KP": 4.38050, "KPP": -0.03252},
            {"V0": 0.68056, "K0": 7.20284, "KP": 0.41237},
            -1.26,
        ),
    ),
)
def test_bm3_ols_reproduces_eosfit7_reference(
    filename: str,
    expected: dict[str, float],
    expected_errors: dict[str, float],
    max_residual: float,
) -> None:
    """Validate unweighted BM3 parameters and standard errors."""
    dataset = read_eos_input(_DATA / filename)

    result = EOSFitter().fit(dataset, EOSFitRequest(model="BM3"))

    assert result.fit.success
    for name, value in expected.items():
        assert result.parameter_values[name] == pytest.approx(value, abs=3.0e-4)
    errors = _parameter_errors(result)
    for name, value in expected_errors.items():
        assert errors[name] == pytest.approx(value, abs=4.0e-4)
    assert result.fit.residuals is not None
    largest = float(result.fit.residuals[np.argmax(np.abs(result.fit.residuals))])
    assert largest == pytest.approx(max_residual, abs=6.0e-3)


@pytest.mark.parametrize(
    (
        "filename",
        "expected",
        "expected_errors",
        "reduced_chi_square",
        "max_residual",
    ),
    (
        (
            "PV_quartz.dat",
            {"V0": 112.98088, "K0": 37.12600, "KP": 5.98823, "KPP": -0.26478},
            {"V0": 0.00199, "K0": 0.09104, "KP": 0.04530},
            0.95,
            -0.03,
        ),
        (
            "PV_topaz.dat",
            {"V0": 345.50726, "K0": 161.99034, "KP": 2.97647, "KPP": -0.02416},
            {"V0": 0.08266, "K0": 2.55519, "KP": 0.16348},
            18.55,
            -2.04,
        ),
    ),
)
def test_bm3_effective_variance_reproduces_eosfit7_reference(
    filename: str,
    expected: dict[str, float],
    expected_errors: dict[str, float],
    reduced_chi_square: float,
    max_residual: float,
) -> None:
    """Validate EosFit-like effective variance and inflate-only covariance."""
    dataset = read_eos_input(_DATA / filename)
    request = EOSFitRequest(
        model="BM3",
        options=EOSFitOptions(solver_options=EffectiveVarianceOptions()),
    )

    result = EOSFitter().fit(dataset, request)

    assert result.fit.success
    for name, value in expected.items():
        tolerance = 2.0e-3 if filename == "PV_topaz.dat" else 3.0e-4
        assert result.parameter_values[name] == pytest.approx(value, abs=tolerance)
    errors = _parameter_errors(result)
    for name, value in expected_errors.items():
        tolerance = 3.0e-4 if name == "V0" else 5.0e-4
        assert errors[name] == pytest.approx(value, abs=tolerance)
    diagnostics = result.fit.diagnostics
    assert diagnostics is not None
    assert diagnostics.reduced_chi_square == pytest.approx(
        reduced_chi_square, abs=3.0e-2
    )
    assert diagnostics.metadata["covariance_scaling"] == (
        CovarianceScaling.INFLATE_ONLY.value
    )
    expected_scale = max(1.0, diagnostics.reduced_chi_square)
    assert diagnostics.metadata["covariance_scale_factor"] == pytest.approx(
        expected_scale
    )
    assert result.fit.residuals is not None
    largest = float(result.fit.residuals[np.argmax(np.abs(result.fit.residuals))])
    assert largest == pytest.approx(max_residual, abs=8.0e-3)


def test_topaz_effective_variance_inflates_only_parameter_uncertainty() -> None:
    """Check that covariance scaling does not change the fitted minimum."""
    dataset = read_eos_input(_DATA / "PV_topaz.dat")
    absolute = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="BM3",
            options=EOSFitOptions(
                solver_options=EffectiveVarianceOptions(
                    covariance_scaling=CovarianceScaling.ABSOLUTE,
                ),
            ),
        ),
    )
    inflated = EOSFitter().fit(
        dataset,
        EOSFitRequest(
            model="BM3",
            options=EOSFitOptions(solver_options=EffectiveVarianceOptions()),
        ),
    )

    np.testing.assert_allclose(inflated.fit.parameters, absolute.fit.parameters)
    assert absolute.fit.errors is not None
    assert inflated.fit.errors is not None
    diagnostics = inflated.fit.diagnostics
    assert diagnostics is not None
    factor = float(diagnostics.metadata["parameter_error_scale_factor"])
    np.testing.assert_allclose(inflated.fit.errors, absolute.fit.errors * factor)
    assert factor == pytest.approx(np.sqrt(18.5254927825), rel=2.0e-8)
