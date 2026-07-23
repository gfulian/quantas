"""Tests for backend-neutral fitting contracts and parameter mappings."""

from __future__ import annotations

import numpy as np
import pytest

from quantas.core.math.fitting import (
    BaseFitModel,
    CovarianceScaling,
    EffectiveVarianceFitter,
    EffectiveVarianceOptions,
    FitMethod,
    FitObservations,
    FitOptions,
    FitSolver,
    FitStatus,
    LeastSquaresFitter,
    LeastSquaresOptions,
    OLSOptions,
    WLSOptions,
    MappedFitModel,
    ParameterDefinition,
    ParameterMap,
    ParameterState,
    covariance_correlation,
    least_squares_fit,
    resolve_mapped_fit_result,
)


class AffineModel(BaseFitModel):
    """Affine model with complete slope and intercept parameters."""

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
        return np.asarray([1.0, 0.0], dtype=np.float64)


def _parameter_map() -> ParameterMap:
    return ParameterMap(
        (
            ParameterDefinition.free("slope", 1.0, lower_bound=0.0, upper_bound=10.0),
            ParameterDefinition.fixed("intercept", 2.0),
            ParameterDefinition.implied("twice_slope"),
            ParameterDefinition.derived("sum"),
        ),
        resolver=lambda values: {
            "twice_slope": 2.0 * values["slope"],
            "sum": values["slope"] + values["intercept"],
        },
    )


def test_parameter_definition_enforces_state_semantics():
    with pytest.raises(ValueError, match="requires an initial value"):
        ParameterDefinition(name="x", state=ParameterState.FREE)
    with pytest.raises(ValueError, match="requires a value"):
        ParameterDefinition(name="x", state=ParameterState.FIXED)
    with pytest.raises(ValueError, match="outside bounds"):
        ParameterDefinition.free("x", 3.0, lower_bound=0.0, upper_bound=2.0)


def test_parameter_map_expands_reduces_and_labels_complete_values():
    mapping = _parameter_map()

    resolved = mapping.expand([3.0])

    assert mapping.free_names == ("slope",)
    assert resolved.names == ("slope", "intercept", "twice_slope", "sum")
    assert resolved.states == (
        ParameterState.FREE,
        ParameterState.FIXED,
        ParameterState.IMPLIED,
        ParameterState.DERIVED,
    )
    np.testing.assert_allclose(resolved.values, [3.0, 2.0, 6.0, 5.0])
    np.testing.assert_allclose(resolved.model_values(), [3.0, 2.0, 6.0])
    np.testing.assert_allclose(mapping.reduce(resolved.values), [3.0])
    np.testing.assert_allclose(mapping.reduce(resolved.as_mapping()), [3.0])


def test_parameter_map_propagates_covariance_to_implied_and_derived_values():
    mapping = _parameter_map()

    jacobian = mapping.resolved_jacobian([3.0])
    covariance = mapping.propagate_covariance([[0.25]], [3.0])

    np.testing.assert_allclose(jacobian[:, 0], [1.0, 0.0, 2.0, 1.0], rtol=1e-6)
    expected_vector = np.asarray([1.0, 0.0, 2.0, 1.0])
    np.testing.assert_allclose(
        covariance, 0.25 * np.outer(expected_vector, expected_vector)
    )


def test_parameter_map_rejects_resolver_overwrite_and_unknown_values():
    overwrite = ParameterMap(
        (ParameterDefinition.free("x", 1.0),),
        resolver=lambda values: {"x": 2.0},
    )
    with pytest.raises(ValueError, match="cannot overwrite"):
        overwrite.expand([1.0])

    unknown = ParameterMap(
        (
            ParameterDefinition.free("x", 1.0),
            ParameterDefinition.implied("y"),
        ),
        resolver=lambda values: {"z": 2.0},
    )
    with pytest.raises(ValueError, match="unknown parameter"):
        unknown.expand([1.0])


def test_fit_observations_preserve_mask_and_select_uncertainties():
    observations = FitObservations(
        x=[1.0, 2.0, 3.0],
        y=[4.0, 5.0, 6.0],
        sigma_x=[0.1, 0.2, 0.3],
        sigma_y=[0.4, 0.5, 0.6],
        mask=[True, False, True],
        x_name="volume",
        y_name="pressure",
    )

    selected = observations.selected()

    assert observations.size == 3
    assert observations.selected_size == 2
    np.testing.assert_allclose(selected.x, [1.0, 3.0])
    np.testing.assert_allclose(selected.require_positive_sigma("y"), [0.4, 0.6])


def test_fit_observations_allow_zero_sigma_but_weighting_rejects_it():
    observations = FitObservations(x=[1.0, 2.0], y=[2.0, 3.0], sigma_y=[0.0, 0.1])
    with pytest.raises(ValueError, match="strictly positive"):
        observations.require_positive_sigma("y")


def test_weighted_least_squares_is_explicit_and_reports_chi_square():
    x = np.linspace(0.0, 4.0, 9)
    y = 2.0 * x + 1.0
    sigma = np.full_like(x, 0.2)

    result = least_squares_fit(
        lambda values, slope, intercept: slope * values + intercept,
        x,
        y,
        p0=[1.0, 0.0],
        sigma=sigma,
        method=FitMethod.WLS,
        absolute_sigma=True,
    )

    assert result.success
    assert result.method is FitMethod.WLS
    assert result.diagnostics is not None
    assert result.diagnostics.weighted
    assert result.diagnostics.chi_square is not None
    assert result.diagnostics.chi_square < 1.0e-20
    assert result.diagnostics.standardized_residuals is not None


def test_ordinary_least_squares_rejects_silent_sigma_use():
    result = least_squares_fit(
        lambda values, slope: slope * values,
        [1.0, 2.0],
        [2.0, 4.0],
        p0=[1.0],
        sigma=[0.1, 0.1],
        method=FitMethod.OLS,
    )
    assert not result.success
    assert result.status is FitStatus.INVALID_INPUT
    assert "select WLS explicitly" in result.message


def test_mapped_model_fits_only_free_values_and_resolves_complete_result():
    mapping = ParameterMap(
        (
            ParameterDefinition.free("slope", 1.0, lower_bound=0.0, upper_bound=10.0),
            ParameterDefinition.fixed("intercept", 2.0),
            ParameterDefinition.derived("twice_slope"),
        ),
        resolver=lambda values: {"twice_slope": 2.0 * values["slope"]},
    )
    model = MappedFitModel(AffineModel(), mapping)
    x = np.linspace(0.0, 5.0, 11)
    y = 3.0 * x + 2.0 + np.linspace(-0.01, 0.01, x.size)

    reduced = LeastSquaresFitter().fit(model, x, y)
    complete = resolve_mapped_fit_result(reduced, mapping)

    assert reduced.parameter_names == ("slope",)
    assert complete.parameter_names == ("slope", "intercept", "twice_slope")
    assert complete.parameter_states[1] is ParameterState.FIXED
    assert complete.parameters is not None
    assert complete.optimizer_parameters is not None
    assert complete.parameters[1] == 2.0
    assert complete.parameters[2] == pytest.approx(2.0 * complete.parameters[0])
    assert complete.covariance is not None
    assert complete.covariance[1, 1] == 0.0


def test_mapped_model_rejects_inconsistent_complete_parameter_order():
    mapping = ParameterMap(
        (
            ParameterDefinition.free("intercept", 0.0),
            ParameterDefinition.free("slope", 1.0),
        )
    )
    with pytest.raises(ValueError, match="must match"):
        MappedFitModel(AffineModel(), mapping)


def test_covariance_correlation_handles_fixed_zero_variance_rows():
    correlation = covariance_correlation([[4.0, 0.0], [0.0, 0.0]])
    assert correlation is not None
    np.testing.assert_allclose(correlation, [[1.0, 0.0], [0.0, 0.0]])


def test_fit_options_reject_invalid_controls():
    with pytest.raises(ValueError, match="max_iterations"):
        FitOptions(max_iterations=0)
    with pytest.raises(ValueError, match="ftol"):
        FitOptions(ftol=0.0)
    with pytest.raises(ValueError, match="supports only"):
        LeastSquaresOptions(method=FitMethod.ODR)


def test_covariance_scaling_policies_are_explicit_and_backward_compatible():
    assert (
        FitOptions(covariance_scaling=CovarianceScaling.INFLATE_ONLY).covariance_scaling
        is CovarianceScaling.INFLATE_ONLY
    )
    assert FitOptions(absolute_sigma=True).covariance_scaling is (
        CovarianceScaling.ABSOLUTE
    )
    assert FitOptions(absolute_sigma=False).covariance_scaling is (
        CovarianceScaling.REDUCED_CHI_SQUARE
    )
    with pytest.raises(ValueError, match="conflicts"):
        FitOptions(
            covariance_scaling=CovarianceScaling.INFLATE_ONLY,
            absolute_sigma=True,
        )


def test_inflate_only_covariance_never_shrinks_reported_uncertainty():
    x = np.linspace(0.0, 4.0, 9)
    sigma = np.full_like(x, 1.0)
    low_scatter = 2.0 * x + 1.0 + 0.01 * np.sin(np.arange(x.size))
    high_scatter = 2.0 * x + 1.0 + 2.0 * np.sin(np.arange(x.size))

    absolute = least_squares_fit(
        lambda values, slope, intercept: slope * values + intercept,
        x,
        low_scatter,
        p0=[1.0, 0.0],
        sigma=sigma,
        method=FitMethod.WLS,
        covariance_scaling=CovarianceScaling.ABSOLUTE,
    )
    low = least_squares_fit(
        lambda values, slope, intercept: slope * values + intercept,
        x,
        low_scatter,
        p0=[1.0, 0.0],
        sigma=sigma,
        method=FitMethod.WLS,
        covariance_scaling=CovarianceScaling.INFLATE_ONLY,
    )
    high = least_squares_fit(
        lambda values, slope, intercept: slope * values + intercept,
        x,
        high_scatter,
        p0=[1.0, 0.0],
        sigma=sigma,
        method=FitMethod.WLS,
        covariance_scaling=CovarianceScaling.INFLATE_ONLY,
    )

    assert absolute.errors is not None
    assert low.errors is not None
    assert high.errors is not None
    np.testing.assert_allclose(low.errors, absolute.errors)
    assert low.diagnostics is not None
    assert low.diagnostics.reduced_chi_square is not None
    assert low.diagnostics.reduced_chi_square < 1.0
    assert low.diagnostics.metadata["covariance_scale_factor"] == 1.0
    assert high.diagnostics is not None
    assert high.diagnostics.reduced_chi_square is not None
    assert high.diagnostics.reduced_chi_square > 1.0
    assert high.diagnostics.metadata["covariance_scale_factor"] == pytest.approx(
        high.diagnostics.reduced_chi_square
    )
    assert np.all(high.errors > absolute.errors)


def test_least_squares_fitter_implements_common_solver_protocol():
    assert isinstance(LeastSquaresFitter(), FitSolver)


def test_fit_problem_applies_mask_weights_and_complete_parameter_mapping():
    observations = FitObservations(
        x=[0.0, 1.0, 2.0, 3.0, 4.0],
        y=[2.0, 5.0, 8.0, 11.0, 14.0],
        sigma_y=[0.1, 0.1, 0.1, 0.1, 0.1],
        mask=[True, True, False, True, True],
    )
    parameters = ParameterMap(
        (
            ParameterDefinition.free("slope", 1.0, lower_bound=0.0, upper_bound=10.0),
            ParameterDefinition.fixed("intercept", 2.0),
        )
    )

    result = LeastSquaresFitter().fit_problem(
        AffineModel(),
        observations,
        parameters,
        WLSOptions(absolute_sigma=True),
    )

    assert result.success
    assert result.method is FitMethod.WLS
    assert result.n_points == 4
    assert result.n_parameters == 1
    assert result.parameter_names == ("slope", "intercept")
    assert result.parameter_states == (ParameterState.FREE, ParameterState.FIXED)
    assert result.parameters is not None
    np.testing.assert_allclose(result.parameters, [3.0, 2.0], atol=1.0e-10)
    assert result.diagnostics is not None
    assert result.diagnostics.weighted


def test_fit_problem_rejects_method_owned_by_another_solver():
    observations = FitObservations(x=[0.0, 1.0], y=[2.0, 5.0])
    parameters = ParameterMap(
        (
            ParameterDefinition.free("slope", 1.0),
            ParameterDefinition.fixed("intercept", 2.0),
        )
    )

    result = LeastSquaresFitter().fit_problem(
        AffineModel(),
        observations,
        parameters,
        FitOptions(method=FitMethod.ODR),
    )

    assert not result.success
    assert result.status is FitStatus.INVALID_INPUT
    assert "supports only" in result.message


def test_typed_solver_options_fix_method_and_validate_tolerances():
    assert OLSOptions().method is FitMethod.OLS
    assert WLSOptions().method is FitMethod.WLS
    assert EffectiveVarianceOptions().method is FitMethod.EFFECTIVE_VARIANCE
    with pytest.raises(ValueError, match="parameter_rtol"):
        EffectiveVarianceOptions(parameter_rtol=-1.0)
    with pytest.raises(ValueError, match="inner_max_iterations"):
        EffectiveVarianceOptions(inner_max_iterations=0)


def test_effective_variance_fitter_implements_common_solver_protocol():
    assert isinstance(EffectiveVarianceFitter(), FitSolver)


def test_effective_variance_combines_x_and_y_uncertainty_and_converges():
    x = np.linspace(0.0, 4.0, 9)
    y = 2.0 * x + 1.0 + 0.002 * np.sin(np.arange(x.size))
    sigma_x = np.full_like(x, 0.05)
    sigma_y = np.full_like(x, 0.20)
    observations = FitObservations(
        x=x,
        y=y,
        sigma_x=sigma_x,
        sigma_y=sigma_y,
    )
    parameters = ParameterMap(
        (
            ParameterDefinition.free("slope", 1.0),
            ParameterDefinition.free("intercept", 0.0),
        )
    )

    result = EffectiveVarianceFitter().fit_problem(
        AffineModel(),
        observations,
        parameters,
        EffectiveVarianceOptions(
            max_iterations=20,
        ),
    )

    assert result.success
    assert result.method is FitMethod.EFFECTIVE_VARIANCE
    assert result.parameters is not None
    np.testing.assert_allclose(result.parameters, [2.0, 1.0], atol=2.0e-3)
    assert result.diagnostics is not None
    metadata = result.diagnostics.metadata
    expected = np.sqrt(sigma_y**2 + (result.parameters[0] * sigma_x) ** 2)
    np.testing.assert_allclose(metadata["effective_sigma"], expected, rtol=1.0e-6)
    np.testing.assert_allclose(
        metadata["variance_x_component"],
        (result.parameters[0] * sigma_x) ** 2,
        rtol=1.0e-6,
    )
    assert result.diagnostics.n_iterations is not None
    assert result.diagnostics.n_iterations >= 2
    assert result.diagnostics.stop_reason == (
        "effective parameters and weights converged"
    )


def test_effective_variance_requires_at_least_one_uncertainty_component():
    observations = FitObservations(x=[0.0, 1.0, 2.0], y=[1.0, 3.0, 5.0])
    parameters = ParameterMap(
        (
            ParameterDefinition.free("slope", 1.0),
            ParameterDefinition.free("intercept", 0.0),
        )
    )

    result = EffectiveVarianceFitter().fit_problem(
        AffineModel(),
        observations,
        parameters,
        FitOptions(method=FitMethod.EFFECTIVE_VARIANCE),
    )

    assert not result.success
    assert result.status is FitStatus.INVALID_INPUT
    assert "requires sigma_x" in result.message


def test_effective_variance_respects_fixed_parameters():
    x = np.linspace(0.0, 4.0, 9)
    observations = FitObservations(
        x=x,
        y=3.0 * x + 2.0,
        sigma_x=np.full_like(x, 0.02),
        sigma_y=np.full_like(x, 0.10),
    )
    parameters = ParameterMap(
        (
            ParameterDefinition.free("slope", 1.0),
            ParameterDefinition.fixed("intercept", 2.0),
        )
    )

    result = EffectiveVarianceFitter().fit_problem(
        AffineModel(),
        observations,
        parameters,
        FitOptions(method=FitMethod.EFFECTIVE_VARIANCE),
    )

    assert result.success
    assert result.parameters is not None
    np.testing.assert_allclose(result.parameters, [3.0, 2.0], atol=1.0e-10)
    assert result.errors is not None
    assert result.errors[1] == 0.0


def test_least_squares_failure_retains_rich_debug_trace():
    x = np.linspace(0.0, 4.0, 9)
    y = 2.0 * x + 1.0

    result = least_squares_fit(
        lambda values, slope, intercept: slope * values + intercept,
        x,
        y,
        p0=[100.0, -100.0],
        max_iterations=1,
        method=FitMethod.OLS,
        metadata={
            "solver_debug": True,
            "parameter_order": ["slope", "intercept"],
        },
    )

    assert not result.success
    assert result.diagnostics is not None
    metadata = result.diagnostics.metadata
    assert metadata["termination_category"] == "iteration_limit"
    assert metadata["backend"] == "scipy.optimize.curve_fit"
    assert metadata["free_parameter_names"] == ["slope", "intercept"]
    assert metadata["recorded_model_evaluations"] >= 1
    assert metadata["first_evaluation"]["parameters"] == [100.0, -100.0]
    assert metadata["last_evaluation"] is not None
    assert metadata["evaluation_trace"]
    assert result.diagnostics.n_evaluations == metadata["recorded_model_evaluations"]


def test_effective_variance_failure_retains_last_valid_iteration():
    x = np.linspace(0.0, 4.0, 9)
    y = 2.0 * x + 1.0 + 0.1 * np.sin(np.arange(x.size))
    observations = FitObservations(
        x=x,
        y=y,
        sigma_x=np.full_like(x, 0.05),
        sigma_y=np.full_like(x, 0.20),
    )
    parameters = ParameterMap(
        (
            ParameterDefinition.free("slope", 10.0),
            ParameterDefinition.free("intercept", -10.0),
        )
    )

    result = EffectiveVarianceFitter().fit_problem(
        AffineModel(),
        observations,
        parameters,
        EffectiveVarianceOptions(
            max_iterations=1,
            metadata={"solver_debug": True},
        ),
    )

    assert not result.success
    assert result.diagnostics is not None
    metadata = result.diagnostics.metadata
    assert metadata["termination_category"] == "iteration_limit"
    assert metadata["effective_variance_converged"] is False
    assert metadata["last_valid_iteration"] == 1
    assert metadata["free_parameter_names"] == ["slope", "intercept"]
    assert len(metadata["last_valid_free_parameters"]) == 2
    assert len(metadata["last_valid_parameters"]) == 2
    assert metadata["total_inner_evaluations"] >= 1
    assert result.diagnostics.n_evaluations == metadata["total_inner_evaluations"]
