"""Tests for the backend-neutral orthogonal distance regression service."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from quantas.core.math.fitting import (
    BaseFitModel,
    CovarianceScaling,
    FitMethod,
    FitObservations,
    FitSolver,
    FitStatus,
    ODRBackendResult,
    ODRBackendUnavailableError,
    ODRPackBackend,
    OrthogonalDistanceFitter,
    ODROptions,
    OrthogonalDistanceOptions,
    ParameterDefinition,
    ParameterMap,
    ParameterState,
)


class AffineModel(BaseFitModel):
    """Affine model used to validate generic ODR behavior."""

    @property
    def name(self) -> str:
        return "affine"

    @property
    def parameter_names(self) -> tuple[str, ...]:
        return ("slope", "intercept")

    def evaluate(self, x, parameters):
        slope, intercept = np.asarray(parameters, dtype=np.float64)
        return slope * np.asarray(x, dtype=np.float64) + intercept

    def initial_guess(self, x, y):
        del x, y
        return np.asarray([1.0, 0.0], dtype=np.float64)

    def derivative_x(self, x, parameters):
        slope = float(np.asarray(parameters, dtype=np.float64)[0])
        return np.full(np.asarray(x).shape, slope, dtype=np.float64)


class AffineWithImpliedModel(AffineModel):
    """Affine model whose third physical parameter is implied and unused."""

    @property
    def parameter_names(self) -> tuple[str, ...]:
        return ("slope", "intercept", "twice_slope")

    def evaluate(self, x, parameters):
        slope, intercept, _ = np.asarray(parameters, dtype=np.float64)
        return slope * np.asarray(x, dtype=np.float64) + intercept


class DeterministicBackend:
    """Minimal fake backend proving that Quantas does not expose odrpack data."""

    @property
    def name(self) -> str:
        return "deterministic-odr"

    @property
    def version(self) -> str:
        return "1"

    def fit(
        self,
        model,
        x,
        y,
        initial_parameters,
        *,
        sigma_x,
        sigma_y,
        bounds,
        options,
    ):
        del initial_parameters, sigma_x, sigma_y, bounds, options
        parameters = np.asarray([3.0], dtype=np.float64)
        adjusted = np.asarray(x, dtype=np.float64) + 0.01
        fitted = np.asarray(model(adjusted, parameters), dtype=np.float64)
        return ODRBackendResult(
            parameters=parameters,
            x_corrections=adjusted - x,
            y_corrections=fitted - y,
            adjusted_x=adjusted,
            fitted_y=fitted,
            covariance=np.asarray([[0.25]], dtype=np.float64),
            residual_variance=2.0,
            sum_square=8.0,
            sum_square_x=2.0,
            sum_square_y=6.0,
            success=True,
            status_code=1,
            stop_reason="synthetic convergence",
            n_iterations=4,
            n_function_evaluations=9,
            n_jacobian_evaluations=3,
            rank_deficiency=0,
            inverse_condition_number=0.5,
            backend_name=self.name,
            backend_version=self.version,
        )


def _fixed_intercept_map() -> ParameterMap:
    return ParameterMap(
        (
            ParameterDefinition.free("slope", 1.0, lower_bound=0.0, upper_bound=10.0),
            ParameterDefinition.fixed("intercept", 2.0),
            ParameterDefinition.implied("twice_slope"),
        ),
        resolver=lambda values: {"twice_slope": 2.0 * values["slope"]},
    )


def test_base_quantas_import_keeps_runtime_odrpack_lazy():
    project_root = Path(__file__).resolve().parents[2]
    environment = os.environ.copy()
    current_pythonpath = environment.get("PYTHONPATH")
    source_path = str(project_root / "src")
    environment["PYTHONPATH"] = (
        source_path
        if not current_pythonpath
        else os.pathsep.join((source_path, current_pythonpath))
    )
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import sys; import quantas; import quantas.modules.eos; "
                "assert 'odrpack' not in sys.modules"
            ),
        ],
        check=False,
        capture_output=True,
        text=True,
        cwd=project_root,
        env=environment,
    )
    assert completed.returncode == 0, completed.stderr


def test_odr_options_fix_method_and_validate_backend_controls():
    assert ODROptions().method is FitMethod.ODR
    assert ODROptions is OrthogonalDistanceOptions
    with pytest.raises(ValueError, match="does not support gtol"):
        OrthogonalDistanceOptions(gtol=1.0e-8)
    with pytest.raises(ValueError, match="strictly positive"):
        OrthogonalDistanceOptions(
            parameter_scales=np.asarray([1.0, 0.0]),
        )


def test_orthogonal_fitter_implements_common_solver_protocol():
    assert isinstance(OrthogonalDistanceFitter(DeterministicBackend()), FitSolver)


def test_backend_neutral_result_maps_fixed_and_implied_parameters():
    observations = FitObservations(
        x=[0.0, 1.0, 2.0, 3.0],
        y=[2.0, 5.0, 8.0, 11.0],
        sigma_x=[0.1, 0.1, 0.1, 0.1],
        sigma_y=[0.2, 0.2, 0.2, 0.2],
    )
    result = OrthogonalDistanceFitter(DeterministicBackend()).fit_problem(
        AffineWithImpliedModel(),
        observations,
        _fixed_intercept_map(),
        ODROptions(covariance_scaling=CovarianceScaling.INFLATE_ONLY),
    )

    assert result.success
    assert result.method is FitMethod.ODR
    assert result.parameter_names == ("slope", "intercept", "twice_slope")
    assert result.parameter_states == (
        ParameterState.FREE,
        ParameterState.FIXED,
        ParameterState.IMPLIED,
    )
    np.testing.assert_allclose(result.parameters, [3.0, 2.0, 6.0])
    assert result.covariance is not None
    assert result.covariance[1, 1] == 0.0
    assert result.diagnostics is not None
    assert result.diagnostics.x_corrections is not None
    assert result.diagnostics.y_corrections is not None
    assert result.diagnostics.metadata["backend"] == "deterministic-odr"
    assert result.diagnostics.metadata["termination_category"] == "converged"
    assert result.diagnostics.metadata["free_parameter_names"] == ["slope"]
    assert (
        result.diagnostics.metadata["sigma_x_summary"]["positive_dynamic_range"] == 1.0
    )
    assert result.diagnostics.metadata["sum_square_x"] == 2.0
    assert result.diagnostics.metadata["sum_square_y"] == 6.0


def test_odr_requires_positive_uncertainties_on_both_coordinates():
    observations = FitObservations(
        x=[0.0, 1.0, 2.0],
        y=[2.0, 5.0, 8.0],
        sigma_y=[0.1, 0.1, 0.1],
    )
    result = OrthogonalDistanceFitter(DeterministicBackend()).fit_problem(
        AffineWithImpliedModel(),
        observations,
        _fixed_intercept_map(),
        ODROptions(),
    )

    assert not result.success
    assert result.status is FitStatus.INVALID_INPUT
    assert "sigma_x" in result.message


def test_missing_required_backend_returns_actionable_failed_result(monkeypatch):
    def unavailable(cls):
        raise ODRBackendUnavailableError(
            "required runtime dependency odrpack is unavailable"
        )

    monkeypatch.setattr(ODRPackBackend, "load", classmethod(unavailable))
    observations = FitObservations(
        x=[0.0, 1.0, 2.0],
        y=[2.0, 5.0, 8.0],
        sigma_x=[0.1, 0.1, 0.1],
        sigma_y=[0.1, 0.1, 0.1],
    )
    result = OrthogonalDistanceFitter().fit_problem(
        AffineWithImpliedModel(),
        observations,
        _fixed_intercept_map(),
        ODROptions(),
    )

    assert not result.success
    assert "odrpack" in result.message
    assert result.metadata["backend_available"] is False


def test_real_odrpack_backend_recovers_affine_parameters_and_corrections():
    x = np.linspace(0.0, 5.0, 15)
    x_observed = x + 0.02 * np.sin(np.arange(x.size))
    y_observed = 2.5 * x + 1.25 + 0.03 * np.cos(np.arange(x.size))
    mapping = ParameterMap(
        (
            ParameterDefinition.free("slope", 2.0, lower_bound=0.0, upper_bound=10.0),
            ParameterDefinition.free(
                "intercept", 1.0, lower_bound=-10.0, upper_bound=10.0
            ),
        )
    )
    result = OrthogonalDistanceFitter().fit_problem(
        AffineModel(),
        FitObservations(
            x=x_observed,
            y=y_observed,
            sigma_x=np.full_like(x, 0.02),
            sigma_y=np.full_like(x, 0.03),
        ),
        mapping,
        OrthogonalDistanceOptions(
            covariance_scaling=CovarianceScaling.ABSOLUTE,
        ),
    )

    assert result.success
    assert result.parameters is not None
    assert result.parameters[0] == pytest.approx(2.5, abs=0.03)
    assert result.parameters[1] == pytest.approx(1.25, abs=0.06)
    assert result.diagnostics is not None
    assert result.diagnostics.x_corrections is not None
    assert result.diagnostics.y_corrections is not None
    assert result.diagnostics.metadata["backend"] == "odrpack95"
    assert result.diagnostics.metadata["sum_square_x"] >= 0.0
    assert result.diagnostics.metadata["sum_square_y"] >= 0.0


def test_odr_options_reject_unknown_or_invalid_typed_controls():
    with pytest.raises(TypeError, match="unexpected keyword"):
        ODROptions(unsupported=1)
    with pytest.raises(ValueError):
        ODROptions(difference_scheme="invalid")


class NonConvergentBackend(DeterministicBackend):
    """Fake backend returning a valid final iterate without convergence."""

    def fit(
        self,
        model,
        x,
        y,
        initial_parameters,
        *,
        sigma_x,
        sigma_y,
        bounds,
        options,
    ):
        del model, sigma_x, sigma_y, bounds
        parameters = np.asarray(initial_parameters, dtype=np.float64) + 0.25
        adjusted = np.asarray(x, dtype=np.float64) + 0.02
        fitted = np.asarray(y, dtype=np.float64) + 0.03
        return ODRBackendResult(
            parameters=parameters,
            x_corrections=np.full_like(np.asarray(x, dtype=np.float64), 0.02),
            y_corrections=np.full_like(np.asarray(y, dtype=np.float64), 0.03),
            adjusted_x=adjusted,
            fitted_y=fitted,
            covariance=np.eye(parameters.size, dtype=np.float64),
            residual_variance=4.0,
            sum_square=12.0,
            sum_square_x=5.0,
            sum_square_y=7.0,
            success=False,
            status_code=4,
            stop_reason="maximum iterations reached",
            n_iterations=int(options.max_iterations or 50),
            n_function_evaluations=31,
            n_jacobian_evaluations=8,
            rank_deficiency=0,
            inverse_condition_number=0.25,
            backend_name=self.name,
            backend_version=self.version,
        )


def test_odr_nonconvergence_retains_last_backend_state_and_problem_ranges():
    observations = FitObservations(
        x=[0.0, 1.0, 2.0, 3.0],
        y=[2.0, 5.0, 8.0, 11.0],
        sigma_x=[0.1, 0.1, 0.1, 0.1],
        sigma_y=[0.2, 0.2, 0.2, 0.2],
    )
    options = ODROptions(max_iterations=3, metadata={"solver_debug": True})

    result = OrthogonalDistanceFitter(NonConvergentBackend()).fit_problem(
        AffineWithImpliedModel(),
        observations,
        _fixed_intercept_map(),
        options,
    )

    assert not result.success
    assert result.diagnostics is not None
    metadata = result.diagnostics.metadata
    assert metadata["termination_category"] == "iteration_limit"
    assert metadata["backend"] == "deterministic-odr"
    assert metadata["backend_status_code"] == 4
    assert metadata["free_parameter_names"] == ["slope"]
    assert metadata["last_valid_parameters"] == [1.25]
    assert metadata["sigma_x_summary"]["positive_dynamic_range"] == 1.0
    assert metadata["sigma_y_summary"]["positive_dynamic_range"] == 1.0
    assert result.diagnostics.n_iterations == 3
    assert result.diagnostics.n_evaluations == 31
    assert result.diagnostics.jacobian_rank == 1
