"""Package layout and dependency boundaries of the SEISMIC core."""

from __future__ import annotations

from pathlib import Path

import pytest

from quantas.core.physics.seismic import (
    MODE_ORDER,
    ChristoffelSolver,
    ElasticMedium,
    IsotropicSeismicVelocities,
    DirectionalEnhancementResult,
    DirectionalGroupResult,
    DirectionalPhaseResult,
    EnhancementFieldResult,
    EnhancementModeResult,
    GroupFieldResult,
    GroupModeResult,
    PhaseFieldResult,
    PhaseModeResult,
    PolarizationTrackingResult,
    WaveMode,
    build_enhancement_field,
    build_group_field,
    build_phase_field,
    isotropic_seismic_velocities,
    track_polarizations,
)


@pytest.mark.physics
@pytest.mark.seismic
def test_seismic_core_public_api_is_importable() -> None:
    assert ChristoffelSolver.__name__ == "ChristoffelSolver"
    assert ElasticMedium.__name__ == "ElasticMedium"
    assert IsotropicSeismicVelocities.__name__ == "IsotropicSeismicVelocities"
    assert isotropic_seismic_velocities.__name__ == "isotropic_seismic_velocities"
    assert DirectionalEnhancementResult.__name__ == ("DirectionalEnhancementResult")
    assert EnhancementModeResult.__name__ == "EnhancementModeResult"
    assert EnhancementFieldResult.__name__ == "EnhancementFieldResult"
    assert DirectionalGroupResult.__name__ == "DirectionalGroupResult"
    assert DirectionalPhaseResult.__name__ == "DirectionalPhaseResult"
    assert GroupModeResult.__name__ == "GroupModeResult"
    assert GroupFieldResult.__name__ == "GroupFieldResult"
    assert PhaseModeResult.__name__ == "PhaseModeResult"
    assert PhaseFieldResult.__name__ == "PhaseFieldResult"
    assert PolarizationTrackingResult.__name__ == "PolarizationTrackingResult"
    assert build_enhancement_field.__name__ == "build_enhancement_field"
    assert build_group_field.__name__ == "build_group_field"
    assert build_phase_field.__name__ == "build_phase_field"
    assert track_polarizations.__name__ == "track_polarizations"
    assert MODE_ORDER == (WaveMode.V_S2, WaveMode.V_S1, WaveMode.V_P)


@pytest.mark.physics
@pytest.mark.seismic
def test_seismic_core_has_no_frontend_or_workflow_dependencies() -> None:
    package = (
        Path(__file__).parents[3] / "src" / "quantas" / "core" / "physics" / "seismic"
    )
    source = "\n".join(
        path.read_text(encoding="utf-8") for path in sorted(package.glob("*.py"))
    ).lower()
    for forbidden in (
        "import click",
        "from click",
        "import dash",
        "from dash",
        "import matplotlib",
        "from matplotlib",
        "quantas.modules",
        "print(",
    ):
        assert forbidden not in source


@pytest.mark.physics
@pytest.mark.seismic
def test_seismic_core_uses_current_terminology() -> None:
    package = (
        Path(__file__).parents[3] / "src" / "quantas" / "core" / "physics" / "seismic"
    )
    source = "\n".join(
        path.read_text(encoding="utf-8") for path in sorted(package.glob("*.py"))
    ).lower()
    assert "legacy" not in source
