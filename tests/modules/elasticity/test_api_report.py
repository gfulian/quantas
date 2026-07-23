"""Public API and neutral-report characterization for elasticity."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.physics.elasticity import (
    DirectionalExtrema,
    ElasticAverages,
    IsotropicElasticProperties,
    StabilityResult,
)
from quantas.models import ReportTable, ResultData
from quantas.modules.elasticity.api import read_elasticity_input, run_elasticity
from quantas.modules.elasticity.models import (
    ElasticityInput,
    ElasticityOptions,
    ElasticityResult,
)
from quantas.modules.elasticity.report import (
    build_elasticity_report,
    stability_table,
)

DATA = Path(__file__).parent / "data" / "hydroxylapatite.dat"


def _properties(value: float) -> IsotropicElasticProperties:
    """Build deterministic isotropic properties for report tests."""
    return IsotropicElasticProperties(
        bulk_modulus=value,
        young_modulus=value + 10.0,
        shear_modulus=value / 2.0,
        poisson_ratio=0.25,
    )


def _averages() -> ElasticAverages:
    """Build deterministic Voigt-Reuss-Hill averages."""
    return ElasticAverages(
        voigt=_properties(100.0),
        reuss=_properties(90.0),
        hill=_properties(95.0),
    )


def test_read_elasticity_input_returns_public_input_model() -> None:
    """The library reader normalizes text input into ElasticityInput."""
    result = read_elasticity_input(DATA)

    assert isinstance(result, ElasticityInput)
    assert result.jobname == "Hydroxylapatite"
    assert result.source == DATA
    assert not hasattr(result, "density")
    assert result.stiffness.shape == (6, 6)


def test_read_elasticity_input_raises_value_error_for_invalid_content(tmp_path) -> None:
    """Malformed input is reported through the public API."""
    filename = tmp_path / "invalid.dat"
    filename.write_text("not enough data\n", encoding="utf-8")

    with pytest.raises(ValueError, match="too short"):
        read_elasticity_input(filename)


def test_run_elasticity_accepts_input_model_and_returns_result_data() -> None:
    """The public runner accepts a passive input object."""
    unstable = ElasticityInput(
        jobname="Fast API check",
        stiffness=np.diag([100.0, 110.0, 120.0, -1.0, 40.0, 50.0]),
    )
    result = run_elasticity(
        unstable,
        options=ElasticityOptions(calculate_2d=False),
    )

    assert isinstance(result, ResultData)
    assert result.metadata.module == "elasticity"
    assert result.metadata.method == "second_order"
    assert result.results["elasticity"].jobname == "Fast API check"
    assert result.warnings


def test_neutral_report_contains_optional_sections() -> None:
    """Elasticity results are converted into shared ReportTable objects."""
    variation = DirectionalExtrema(
        minimum=100.0,
        maximum=200.0,
        anisotropy=2.0,
        minimum_axis=[1.0, 0.0, 0.0],
        maximum_axis=[0.0, 1.0, 0.0],
    )
    result = ElasticityResult(
        jobname="Synthetic",
        crystal_system="orthorhombic",
        stiffness=np.eye(6) * 100.0,
        compliance=np.eye(6) * 0.01,
        averages=_averages(),
        stability=StabilityResult(True, np.ones(6) * 100.0),
        variations={"young_modulus": variation},
    )
    tables = build_elasticity_report(
        ElasticityInput(jobname="Synthetic"),
        ElasticityOptions(),
        result,
    )

    assert all(isinstance(table, ReportTable) for table in tables)
    assert [table.title for table in tables][-1:] == [
        "Directional elastic extrema",
    ]
    assert all("seismic" not in table.title.lower() for table in tables)
    assert stability_table(result).metadata["positive_definite"] is True


def test_stability_report_marks_instability() -> None:
    """The neutral table exposes positive-definiteness metadata."""
    table = stability_table(
        ElasticityResult(stability=StabilityResult(False, np.array([-1.0, 2.0])))
    )
    assert table.metadata["positive_definite"] is False
