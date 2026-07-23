"""Workflow and event-contract tests for the elasticity calculator."""

from __future__ import annotations

import numpy as np

from quantas.core.events import EventLevel, ListObserver
from quantas.core.physics.elasticity import DirectionalExtrema
from quantas.modules.elasticity import calculator as calculator_module
from quantas.modules.elasticity.calculator import ElasticityCalculator
from quantas.modules.elasticity.models import ElasticityInput, ElasticityOptions


def _stable_input() -> ElasticityInput:
    return ElasticityInput(
        jobname="Synthetic orthorhombic",
        stiffness=np.diag([150.0, 160.0, 170.0, 50.0, 60.0, 70.0]),
    )


def _variation(value: float = 1.0) -> DirectionalExtrema:
    return DirectionalExtrema(
        minimum=value,
        maximum=value * 2.0,
        anisotropy=2.0,
        minimum_axis=[1.0, 0.0, 0.0],
        maximum_axis=[0.0, 1.0, 0.0],
    )


def test_calculator_emits_structured_result_events(monkeypatch) -> None:
    def fake_variations(tensor, result):
        result.add_variation("young_modulus", _variation())

    monkeypatch.setattr(
        calculator_module,
        "calculate_directional_variations",
        fake_variations,
    )
    observer = ListObserver()
    calculator = ElasticityCalculator(
        _stable_input(),
        options=ElasticityOptions(calculate_2d=False),
        observer=observer,
    )
    result = calculator.execute()
    kinds = [
        event.data.get("kind")
        for event in observer.events
        if event.level is EventLevel.RESULT and event.data.get("kind")
    ]
    assert kinds == [
        "settings",
        "input",
        "averages",
        "stability",
        "variations",
    ]
    elasticity = result.results["elasticity"]
    assert not hasattr(elasticity, "isotropic_velocities")
    assert set(elasticity.variations) == {"young_modulus"}
    assert elasticity.metadata["tensor_class"] == "OrthorhombicElasticTensor"


def test_calculator_converts_polar_callbacks_to_events(monkeypatch) -> None:
    monkeypatch.setattr(
        calculator_module,
        "calculate_directional_variations",
        lambda tensor, result: result.add_variation("young_modulus", _variation()),
    )

    def fake_polar(tensor, result, options, progress_callback=None, step_callback=None):
        step_callback("xy", "young_modulus")
        progress_callback("xy: Young modulus", 1, 4)
        progress_callback("xy: Young modulus", 4, 4)
        result.add_2d_data("xy", "theta", np.array([np.pi / 2.0]))
        result.add_2d_data("xy", "phi", np.array([0.0]))
        result.add_2d_data("xy", "young_modulus", np.array([100.0]))

    monkeypatch.setattr(calculator_module, "calculate_2d_properties", fake_polar)
    observer = ListObserver()
    result = ElasticityCalculator(
        _stable_input(),
        options=ElasticityOptions(
            calculate_2d=True,
            ntheta=4,
        ),
        observer=observer,
    ).execute()
    step_events = [
        event for event in observer.events if event.data.get("kind") == "2d_step"
    ]
    progress_events = [
        event for event in observer.events if event.level is EventLevel.PROGRESS
    ]
    assert step_events[0].data == {
        "kind": "2d_step",
        "plane": "xy",
        "property": "young_modulus",
    }
    assert [event.progress for event in progress_events] == [0.25, 1.0]
    assert result.results["elasticity"].has_2d_data() is True


def test_unstable_matrix_skips_directional_and_polar_analysis(monkeypatch) -> None:
    called = {"variations": False, "polar": False}

    def fake_variations(tensor, result):
        called["variations"] = True

    def fake_polar(*args, **kwargs):
        called["polar"] = True

    monkeypatch.setattr(
        calculator_module, "calculate_directional_variations", fake_variations
    )
    monkeypatch.setattr(calculator_module, "calculate_2d_properties", fake_polar)
    unstable = ElasticityInput(
        stiffness=np.diag([100.0, 110.0, 120.0, -1.0, 40.0, 50.0]),
    )
    result = ElasticityCalculator(
        unstable,
        options=ElasticityOptions(calculate_2d=True),
    ).execute()
    assert called == {"variations": False, "polar": False}
    elasticity = result.results["elasticity"]
    assert not hasattr(elasticity, "isotropic_velocities")
    assert elasticity.variations == {}
    assert elasticity.properties_2d == {}
    assert result.warnings
