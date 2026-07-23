"""Tests for the typed EOS solver-options API."""

from __future__ import annotations

from pathlib import Path

import pytest

from quantas.core.math.fitting import (
    EffectiveVarianceOptions,
    FitMethod,
    FitOptions,
    OLSOptions,
    ODROptions,
    OrthogonalDistanceOptions,
    WLSOptions,
)
from quantas.modules.eos import (
    EOSFitOptions,
    EOSFitRequest,
    EOSFitter,
    read_eos_input,
)

_DATA = Path(__file__).with_name("data")


@pytest.mark.parametrize(
    ("method", "expected_type"),
    (
        (FitMethod.OLS, OLSOptions),
        (FitMethod.WLS, WLSOptions),
        (FitMethod.EFFECTIVE_VARIANCE, EffectiveVarianceOptions),
        (FitMethod.ODR, OrthogonalDistanceOptions),
    ),
)
def test_method_convenience_creates_typed_solver_options(method, expected_type):
    options = EOSFitOptions(method=method)

    assert options.method is method
    assert isinstance(options.solver_options, expected_type)


def test_solver_options_and_method_must_agree():
    with pytest.raises(ValueError, match="conflicts with solver_options"):
        EOSFitOptions(
            method=FitMethod.WLS,
            solver_options=ODROptions(),
        )


def test_all_typed_solver_options_share_common_control_names():
    options = (
        OLSOptions(max_iterations=100, ftol=1.0e-9, xtol=1.0e-9),
        WLSOptions(max_iterations=100, ftol=1.0e-9, xtol=1.0e-9),
        EffectiveVarianceOptions(
            max_iterations=20,
            inner_max_iterations=100,
            ftol=1.0e-9,
            xtol=1.0e-9,
        ),
        ODROptions(max_iterations=100, ftol=1.0e-9, xtol=1.0e-9),
    )

    for item in options:
        assert item.max_iterations is not None
        assert item.ftol == 1.0e-9
        assert item.xtol == 1.0e-9
        assert hasattr(item, "gtol")
        assert hasattr(item, "covariance_scaling")
        assert not hasattr(item, "optimizer_options")


def test_untyped_optimizer_options_mapping_is_removed():
    with pytest.raises(TypeError, match="unexpected keyword"):
        FitOptions(optimizer_options={"maxfev": 100})
    with pytest.raises(TypeError, match="unexpected keyword"):
        EOSFitOptions(optimizer_options={"maxfev": 100})


def test_solver_options_are_serialized_as_a_nested_typed_contract():
    options = EOSFitOptions(
        solver_options=ODROptions(
            max_iterations=80,
            difference_scheme="central",
            ndigit=12,
        ),
        metadata={"origin": "api-test"},
    )

    payload = options.as_dict()

    assert payload["method"] == FitMethod.ODR.value
    assert payload["solver_options"]["type"] == "OrthogonalDistanceOptions"
    assert payload["solver_options"]["max_iterations"] == 80
    assert payload["solver_options"]["difference_scheme"] == "central"
    assert payload["solver_options"]["ndigit"] == 12
    assert payload["metadata"] == {"origin": "api-test"}


def test_eos_options_do_not_mutate_the_callers_solver_options():
    solver = WLSOptions(metadata={"solver": "original"})
    options = EOSFitOptions(
        solver_options=solver,
        metadata={"workflow": "eos"},
    )
    numerical = options.to_fit_options({"dataset": "quartz"})

    assert solver.metadata == {"solver": "original"}
    assert solver.covariance_policy.value == "absolute"
    assert options.solver_options is not solver
    assert options.solver_options is not None
    assert options.solver_options.covariance_policy.value == "inflate_only"
    assert numerical.metadata == {
        "solver": "original",
        "workflow": "eos",
        "dataset": "quartz",
    }


def test_original_odr_api_example_now_executes_unchanged():
    dataset = read_eos_input(_DATA / "PV_quartz.dat")
    request = EOSFitRequest(
        model="BM3",
        options=EOSFitOptions(
            method=FitMethod.ODR,
            solver_options=OrthogonalDistanceOptions(
                method=FitMethod.ODR,
            ),
        ),
    )

    result = EOSFitter().fit(dataset, request)

    assert result.fit.success
    assert result.fit.method is FitMethod.ODR
    assert result.parameter_values["K0"] == pytest.approx(37.12518, abs=3.0e-3)
    assert result.parameter_values["V0"] == pytest.approx(112.980885, abs=3.0e-4)


def test_new_concise_wls_api_uses_the_same_workflow_service():
    dataset = read_eos_input(_DATA / "PV_quartz.dat")
    request = EOSFitRequest(
        model="BM3",
        options=EOSFitOptions(
            solver_options=WLSOptions(),
        ),
    )

    result = EOSFitter().fit(dataset, request)

    assert result.fit.success
    assert result.fit.method is FitMethod.WLS
    assert result.fit.diagnostics is not None
    assert result.fit.diagnostics.weighted
