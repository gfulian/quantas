"""Frontend-neutral sampled seismic workflow."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.events import EventLevel, ListObserver
from quantas.core.geometry import Hemisphere
from quantas.core.physics.elasticity import ElasticTensor
from quantas.core.physics.seismic import ElasticMedium, SamplingLevel
from quantas.modules.seismic import (
    SeismicCalculator,
    SeismicInput,
    SeismicOptions,
    SeismicResult,
    run_seismic,
)


DATA = (
    Path(__file__).parents[2]
    / "physics"
    / "seismic"
    / "data"
    / "hydroxylapatite.dat"
)


def _input() -> SeismicInput:
    lines = DATA.read_text(encoding="utf-8").splitlines()
    stiffness = np.zeros((6, 6), dtype=float)
    for row, line in enumerate(lines[1:7]):
        values = np.asarray([float(value) for value in line.split()])
        stiffness[row, row:] = values
    stiffness += np.triu(stiffness, 1).T
    return SeismicInput(
        jobname="Hydroxylapatite",
        stiffness=stiffness,
        density=float(lines[7]),
        source=DATA,
    )


def _options() -> SeismicOptions:
    return SeismicOptions(
        ntheta=4,
        nphi=7,
        hemisphere=Hemisphere.UPPER,
        level=SamplingLevel.GROUP,
        batch_size=5,
        track_polarization_axes=True,
    )


@pytest.mark.module
@pytest.mark.seismic
def test_calculator_builds_structured_result_and_progress_events() -> None:
    observer = ListObserver()
    result = SeismicCalculator(_input(), _options(), observer=observer).execute()

    assert result.metadata.module == "seismic"
    assert result.metadata.method == "christoffel"
    payload = result.results["seismic"]
    assert isinstance(payload, SeismicResult)
    assert payload.field.n_points == 28
    assert payload.field.group is not None
    assert payload.field.enhancement is None
    assert payload.field.tracking is not None
    assert payload.isotropic_velocities.shear > 0.0
    assert (
        payload.isotropic_velocities.compressional > payload.isotropic_velocities.shear
    )

    progress = [
        event for event in observer.events if event.level is EventLevel.PROGRESS
    ]
    assert [event.data["current"] for event in progress] == [5, 10, 15, 20, 25, 28]
    assert progress[-1].progress == pytest.approx(1.0)
    kinds = [
        event.data.get("kind")
        for event in observer.events
        if event.level is EventLevel.RESULT and event.data.get("kind")
    ]
    assert kinds == [
        "settings",
        "input",
        "isotropic_reference",
        "seismic_field",
    ]


@pytest.mark.module
@pytest.mark.seismic
def test_public_api_accepts_an_elastic_medium() -> None:
    source = _input()
    assert source.stiffness is not None
    medium = ElasticMedium(ElasticTensor(source.stiffness), source.density)

    result = run_seismic(medium, options=_options())

    assert result.metadata.module == "seismic"
    payload = result.results["seismic"]
    assert isinstance(payload, SeismicResult)
    assert payload.field.grid.hemisphere is Hemisphere.UPPER


@pytest.mark.module
@pytest.mark.seismic
@pytest.mark.parametrize("density", [0.0, -1.0, np.nan, np.inf])
def test_workflow_rejects_invalid_density(density: float) -> None:
    invalid = _input()
    invalid.density = density

    with pytest.raises(ValueError, match="finite and positive"):
        run_seismic(invalid, options=_options())


@pytest.mark.module
@pytest.mark.seismic
@pytest.mark.parametrize("hemisphere", list(Hemisphere))
def test_workflow_supports_each_hemisphere(hemisphere: Hemisphere) -> None:
    options = SeismicOptions(
        ntheta=4,
        nphi=7,
        hemisphere=hemisphere,
        level=SamplingLevel.PHASE,
        batch_size=8,
        track_polarization_axes=False,
    )

    result = run_seismic(_input(), options=options)
    payload = result.results["seismic"]

    assert isinstance(payload, SeismicResult)
    assert payload.grid.hemisphere is hemisphere
    assert payload.field.n_points == 28


@pytest.mark.module
@pytest.mark.seismic
def test_large_grid_progress_is_monotonic_and_not_persisted() -> None:
    observer = ListObserver()
    options = SeismicOptions(
        ntheta=20,
        nphi=30,
        hemisphere=Hemisphere.FULL,
        level=SamplingLevel.PHASE,
        batch_size=64,
        track_polarization_axes=False,
    )

    result = run_seismic(_input(), options=options, observer=observer)
    progress = [
        event for event in observer.events if event.level is EventLevel.PROGRESS
    ]

    values = [event.progress for event in progress]
    assert values == sorted(values)
    assert values[-1] == pytest.approx(1.0)
    assert progress[-1].data == {
        "kind": "sampling_progress",
        "current": 600,
        "total": 600,
    }
    assert all(event.level != EventLevel.PROGRESS.value for event in result.events)


@pytest.mark.module
@pytest.mark.seismic
def test_workflow_rejects_non_positive_definite_stiffness() -> None:
    unstable = SeismicInput(
        jobname="Unstable",
        stiffness=np.diag([100.0, 110.0, 120.0, -1.0, 40.0, 50.0]),
        density=3000.0,
    )
    with pytest.raises(ValueError, match="positive-definite"):
        run_seismic(unstable, options=_options())


@pytest.mark.module
@pytest.mark.seismic
def test_seismic_module_has_no_frontend_dependencies() -> None:
    package = Path(__file__).parents[3] / "src" / "quantas" / "modules" / "seismic"
    core = (
        Path(__file__).parents[3] / "src" / "quantas" / "core" / "physics" / "seismic"
    )
    source = "\n".join(
        path.read_text(encoding="utf-8")
        for root in (package, core)
        for path in root.glob("*.py")
    )
    for forbidden in (
        "import click",
        "import dash",
        "import matplotlib",
    ):
        assert forbidden not in source.lower()
    assert "print(" not in source


@pytest.mark.module
@pytest.mark.seismic
def test_options_normalize_string_enumerations() -> None:
    options = SeismicOptions(hemisphere="lower", level="phase")  # type: ignore[arg-type]
    assert options.hemisphere is Hemisphere.LOWER
    assert options.level is SamplingLevel.PHASE
