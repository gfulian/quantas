"""Tests for frontend-neutral Quantas EOS fitting requests and results."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.math.fitting import (
    FitMethod,
    FitQuality,
    FitResult,
    FitStatus,
    WLSOptions,
)
from quantas.modules.eos import (
    EOSDataset,
    EOSFitDomain,
    EOSFitOptions,
    EOSFitRequest,
    EOSFitResult,
    ParameterConstraint,
)


def _dataset() -> EOSDataset:
    return EOSDataset(
        jobname="synthetic",
        columns={
            "pressure": np.asarray([0.0, 1.0, 2.0]),
            "sigma_pressure": np.asarray([0.01, 0.02, 0.03]),
            "volume": np.asarray([100.0, 99.0, 98.2]),
            "sigma_volume": np.asarray([0.1, 0.1, 0.1]),
        },
        units={
            "pressure": "GPa",
            "sigma_pressure": "GPa",
            "volume": "angstrom^3",
            "sigma_volume": "angstrom^3",
        },
    )


def test_pressure_series_builds_pressure_residual_observations():
    series = _dataset().series("volume")

    observations = series.pressure_observations()

    assert observations.x_name == "volume"
    assert observations.y_name == "pressure"
    np.testing.assert_allclose(observations.x, [100.0, 99.0, 98.2])
    np.testing.assert_allclose(observations.y, [0.0, 1.0, 2.0])
    np.testing.assert_allclose(observations.sigma_x, [0.1, 0.1, 0.1])
    np.testing.assert_allclose(observations.sigma_y, [0.01, 0.02, 0.03])


def test_eos_fit_request_normalizes_model_method_and_constraints():
    request = EOSFitRequest(
        model="BM3",
        target="VOLUME",
        constraints=(
            ParameterConstraint.free("K0", 150.0, lower_bound=0.0),
            ParameterConstraint.fixed("KP", 4.0),
        ),
        options=EOSFitOptions(solver_options=WLSOptions()),
        request_id="volume-001",
    )

    assert request.model.tag == "BM3"
    assert request.target == "volume"
    assert request.domain is EOSFitDomain.PRESSURE_VOLUME
    assert request.options.method is FitMethod.WLS
    assert request.constraints[1].state.value == "fixed"
    assert request.as_dict()["model"]["tag"] == "BM3"


def test_eos_fit_request_rejects_angles_as_pressure_eos_targets():
    with pytest.raises(ValueError, match="volume or linear"):
        EOSFitRequest(model="BM3", target="alpha")


def test_eos_energy_request_requires_integrated_model_and_energy_target():
    with pytest.raises(ValueError, match="target='energy'"):
        EOSFitRequest(
            model="BM3",
            target="volume",
            domain=EOSFitDomain.ENERGY_VOLUME,
        )
    with pytest.raises(ValueError, match="no integrated"):
        EOSFitRequest(
            model="T3",
            target="energy",
            domain=EOSFitDomain.ENERGY_VOLUME,
        )


def test_eos_request_does_not_infer_weighting_from_dataset_uncertainties():
    request = EOSFitRequest(model="V3", target="volume")
    assert _dataset().has("sigma_pressure")
    assert request.options.method is FitMethod.OLS


def test_eos_fit_result_serializes_generic_result_and_predictions():
    request = EOSFitRequest(model="BM3", target="volume")
    fit = FitResult(
        success=True,
        status=FitStatus.SUCCESS,
        quality=FitQuality.GOOD,
        parameters=np.asarray([100.0, 150.0, 4.0]),
    )
    result = EOSFitResult(
        request=request,
        fit=fit,
        predictions={"pressure": np.asarray([0.0, 1.0])},
        derived={"compression_ratio": 0.98},
    )

    payload = result.as_dict()

    assert payload["request"]["model"]["tag"] == "BM3"
    assert payload["fit"]["status"] == "success"
    assert payload["predictions"]["pressure"] == [0.0, 1.0]
