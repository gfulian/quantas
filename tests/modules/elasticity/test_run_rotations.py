"""Run-time tensor transformation tests for Elasticity and SEISMIC."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quantas.core.geometry import TensorRotation
from quantas.core.physics.seismic import SamplingLevel
from quantas.modules.elasticity import ElasticityOptions, run_elasticity
from quantas.modules.elasticity.api import read_elasticity_input
from quantas.modules.elasticity.io import ElasticityInputCreator
from quantas.modules.seismic import SeismicOptions, run_seismic


CRYSTAL_CALCITE = (
    Path(__file__).parents[2]
    / "interfaces"
    / "data"
    / "calcite_crystal_elastcon_excerpt.out"
)


def _write_calcite_input(tmp_path: Path) -> Path:
    creator = ElasticityInputCreator("crystal")
    completed, error = creator.read(CRYSTAL_CALCITE)
    assert completed is True
    assert error is None
    filename = tmp_path / "calcite.dat"
    creator.write(filename, "Calcite source frame")
    return filename


@pytest.mark.module
@pytest.mark.elasticity
def test_inpgen_remains_a_source_frame_parser(tmp_path) -> None:
    """Input generation writes exactly the components reported by CRYSTAL."""
    filename = _write_calcite_input(tmp_path)
    text = filename.read_text(encoding="utf-8")
    assert "quantas_tensor_frame" not in text

    parsed = read_elasticity_input(filename)
    assert parsed.stiffness is not None
    assert parsed.stiffness[0, 3] == pytest.approx(0.0, abs=1.0e-6)
    assert parsed.stiffness[0, 4] == pytest.approx(20.670, abs=1.0e-6)


@pytest.mark.module
@pytest.mark.elasticity
def test_elasticity_default_analysis_preserves_source_frame(tmp_path) -> None:
    """Without a request, input and analysis stiffness are identical."""
    filename = _write_calcite_input(tmp_path)
    result = run_elasticity(filename, ElasticityOptions())
    payload = result.results["elasticity"]

    np.testing.assert_allclose(
        payload.stiffness,
        result.input_data.data["stiffness"],
        atol=0.0,
        rtol=0.0,
    )
    assert payload.metadata["tensor_frame"] == {
        "source_frame": "source",
        "analysis_frame": "source",
        "transformed": False,
    }


@pytest.mark.module
@pytest.mark.elasticity
def test_elasticity_rotation_is_applied_only_before_analysis(tmp_path) -> None:
    """A user rotation leaves input provenance intact and rotates results."""
    filename = _write_calcite_input(tmp_path)
    rotation = TensorRotation.from_xyz(0.0, 0.0, 30.0)
    source = run_elasticity(filename, ElasticityOptions())
    rotated = run_elasticity(
        filename,
        ElasticityOptions(rotation=rotation),
    )
    source_payload = source.results["elasticity"]
    rotated_payload = rotated.results["elasticity"]

    source_matrix = rotated.input_data.data["stiffness"]
    assert source_matrix[0, 4] == pytest.approx(20.670, abs=1.0e-6)
    assert rotated_payload.stiffness[0, 3] == pytest.approx(20.670, abs=1.0e-6)
    assert rotated_payload.stiffness[0, 4] == pytest.approx(0.0, abs=1.0e-10)
    assert rotated_payload.metadata["tensor_frame"]["transformed"] is True
    np.testing.assert_allclose(
        rotated_payload.metadata["tensor_frame"]["component_transform"],
        rotation.matrix,
    )
    np.testing.assert_allclose(
        rotated_payload.averages.as_array(),
        source_payload.averages.as_array(),
        atol=2.0e-10,
        rtol=2.0e-12,
    )


@pytest.mark.module
@pytest.mark.seismic
def test_seismic_uses_the_same_run_time_rotation(tmp_path) -> None:
    """SEISMIC applies the identical transformation before Christoffel sampling."""
    filename = _write_calcite_input(tmp_path)
    rotation = TensorRotation.from_xyz(0.0, 0.0, 30.0)
    result = run_seismic(
        filename,
        SeismicOptions(
            ntheta=3,
            nphi=5,
            level=SamplingLevel.PHASE,
            batch_size=8,
            rotation=rotation,
        ),
    )
    payload = result.results["seismic"]

    assert result.input_data.data["stiffness"][0, 4] == pytest.approx(20.670)
    assert payload.stiffness[0, 3] == pytest.approx(20.670, abs=1.0e-6)
    assert payload.stiffness[0, 4] == pytest.approx(0.0, abs=1.0e-10)
    assert payload.metadata["tensor_frame"]["transformed"] is True
    np.testing.assert_allclose(
        payload.metadata["tensor_frame"]["component_transform"],
        rotation.matrix,
    )


@pytest.mark.module
@pytest.mark.elasticity
def test_rotated_result_reports_and_hdf5_preserve_before_and_after(tmp_path) -> None:
    """Rotation provenance survives neutral reports and native persistence."""
    from quantas.modules.elasticity import (
        build_elasticity_report,
        read_elasticity_hdf5,
        write_elasticity_hdf5,
    )

    filename = _write_calcite_input(tmp_path)
    rotation = TensorRotation.from_xyz(0.0, 0.0, 30.0)
    result = run_elasticity(filename, ElasticityOptions(rotation=rotation))
    titles = [table.title for table in build_elasticity_report(result)]
    assert "Stiffness matrix before rotation (GPa)" in titles
    assert "Tensor component transformation" in titles
    assert "Stiffness matrix after rotation (GPa)" in titles

    output = write_elasticity_hdf5(result, tmp_path / "calcite_rotated")
    restored = read_elasticity_hdf5(output)
    payload = restored.results["elasticity"]
    assert payload.metadata["tensor_frame"]["transformed"] is True
    np.testing.assert_allclose(
        payload.metadata["tensor_frame"]["component_transform"],
        rotation.matrix,
    )
    np.testing.assert_allclose(
        restored.input_data.data["stiffness"],
        result.input_data.data["stiffness"],
    )
    np.testing.assert_allclose(
        payload.stiffness, result.results["elasticity"].stiffness
    )
