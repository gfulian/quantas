"""Input and native HDF5 I/O for sampled seismic results."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from quantas.core.geometry import Hemisphere
from quantas.core.physics.seismic import SamplingLevel
from quantas.models import ResultData, ResultMetadata
from quantas.modules.seismic import (
    SeismicOptions,
    SeismicResult,
    read_seismic_hdf5,
    read_seismic_input,
    run_seismic,
    write_seismic_hdf5,
)


DATA = (
    Path(__file__).parents[2]
    / "physics"
    / "seismic"
    / "data"
    / "hydroxylapatite.dat"
)


def _options(
    level: SamplingLevel = SamplingLevel.ENHANCEMENT,
    *,
    tracking: bool = True,
) -> SeismicOptions:
    return SeismicOptions(
        ntheta=3,
        nphi=5,
        hemisphere=Hemisphere.UPPER,
        level=level,
        batch_size=4,
        track_polarization_axes=tracking,
    )


def _write_matrix_input(
    path: Path,
    matrix: np.ndarray,
    form: str,
    *,
    jobname: bool = True,
    density: float = 3178.0,
) -> None:
    lines: list[str] = []
    if jobname:
        lines.append("Test crystal")
    if form == "full":
        rows = matrix
    elif form == "upper":
        rows = [matrix[index, index:] for index in range(6)]
    elif form == "lower":
        rows = [matrix[index, : index + 1] for index in range(6)]
    else:
        raise AssertionError(form)
    lines.extend(" ".join(f"{value:.10f}" for value in row) for row in rows)
    lines.append(f"{density:.10f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.mark.module
@pytest.mark.seismic
@pytest.mark.parametrize("form", ["full", "upper", "lower"])
def test_text_reader_supports_all_documented_matrix_forms(
    tmp_path: Path,
    form: str,
) -> None:
    reference = read_seismic_input(DATA)
    assert reference.stiffness is not None
    filename = tmp_path / f"{form}.dat"
    _write_matrix_input(filename, reference.stiffness, form)

    result = read_seismic_input(filename)

    assert result.jobname == "Test crystal"
    assert result.density == pytest.approx(3178.0)
    np.testing.assert_allclose(result.stiffness, reference.stiffness)
    assert result.source == filename
    assert result.raw == filename.read_text(encoding="utf-8")


@pytest.mark.module
@pytest.mark.seismic
def test_text_reader_accepts_a_full_matrix_without_jobname(tmp_path: Path) -> None:
    reference = read_seismic_input(DATA)
    assert reference.stiffness is not None
    filename = tmp_path / "no_job.dat"
    _write_matrix_input(filename, reference.stiffness, "full", jobname=False)

    result = read_seismic_input(filename)

    assert result.jobname == "Unknown"
    np.testing.assert_allclose(result.stiffness, reference.stiffness)


@pytest.mark.module
@pytest.mark.seismic
@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("Crystal\n1 2 3\n", "six stiffness rows"),
        (
            "Crystal\n"
            + "\n".join(
                [
                    "1 0 0 0 0 0",
                    "0 1 0 0 0 0",
                    "0 0 1 0 0 0",
                    "0 0 0 1 0 0",
                    "0 0 0 0 1 0",
                    "0 0 0 0 0 1",
                ]
            )
            + "\n0\n",
            "finite and positive",
        ),
        (
            "Crystal\n" + "\n".join(["1 2", "1", "1", "1", "1", "1"]) + "\n3000\n",
            "full, upper triangular, or lower triangular",
        ),
    ],
)
def test_text_reader_rejects_invalid_inputs(
    tmp_path: Path,
    content: str,
    message: str,
) -> None:
    filename = tmp_path / "invalid.dat"
    filename.write_text(content, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        read_seismic_input(filename)


@pytest.mark.module
@pytest.mark.seismic
def test_public_run_api_accepts_a_text_input_path() -> None:
    result = run_seismic(DATA, options=_options(SamplingLevel.PHASE, tracking=False))

    payload = result.results["seismic"]
    assert isinstance(payload, SeismicResult)
    assert payload.jobname == "Hydroxylapatite"
    assert result.input_data is not None
    assert result.input_data.source == DATA
    assert result.input_data.raw == DATA.read_text(encoding="utf-8")


@pytest.mark.module
@pytest.mark.seismic
def test_enhancement_hdf5_round_trip_is_complete_and_typed(tmp_path: Path) -> None:
    original = run_seismic(DATA, options=_options())
    output = write_seismic_hdf5(
        original,
        tmp_path / "seismic_result",
        report_text="Seismic report",
    )

    assert output.suffix == ".hdf5"
    with h5py.File(output, "r") as h5:
        assert set(h5) == {
            "metadata",
            "input",
            "options",
            "results",
            "diagnostics",
            "events",
        }
        assert h5["metadata"].attrs["module"] == "seismic"
        assert h5["results"].attrs["density_unit"] == "kg m^-3"
        assert h5["results/fields/phase"].attrs["mode_order"] == "v_s2,v_s1,v_p"
        assert h5["results/fields/phase/phase_speeds"].attrs["unit"] == "km s^-1"
        assert h5["results/fields/group/power_flow_angles"].attrs["unit"] == "rad"
        enhancement = h5["results/fields/enhancement"]
        assert "enhancement" not in enhancement
        assert "log10_enhancement" in enhancement
        assert enhancement["log10_enhancement"].attrs["representation"] == "log10(A)"
        assert enhancement["log10_enhancement"].compression == "gzip"
        assert enhancement["log10_enhancement"].chunks is not None
        assert "seismic" in h5["diagnostics"]
        assert h5["diagnostics/report_text"][()].decode() == "Seismic report"
        assert int(h5["events"].attrs["count"]) == len(original.events)

    restored = read_seismic_hdf5(output)
    assert restored.metadata.module == "seismic"
    assert restored.metadata.method == "christoffel"
    assert restored.options == original.options
    assert restored.warnings == original.warnings
    assert len(restored.events) == len(original.events)
    assert restored.input_data is not None
    assert restored.input_data.raw == DATA.read_text(encoding="utf-8")

    source = original.results["seismic"]
    target = restored.results["seismic"]
    assert isinstance(source, SeismicResult)
    assert isinstance(target, SeismicResult)
    _assert_result_equal(source, target)
    assert not target.field.phase.phase_speeds.flags.writeable
    assert target.field.group is not None
    assert not target.field.group.group_velocities.flags.writeable
    assert target.field.enhancement is not None
    assert not target.field.enhancement.log10_enhancement.flags.writeable
    assert target.field.tracking is not None
    assert not target.field.tracking.polarizations.flags.writeable


@pytest.mark.module
@pytest.mark.seismic
@pytest.mark.parametrize("level", [SamplingLevel.PHASE, SamplingLevel.GROUP])
def test_hdf5_preserves_optional_sampling_levels(
    tmp_path: Path,
    level: SamplingLevel,
) -> None:
    original = run_seismic(DATA, options=_options(level, tracking=False))
    output = write_seismic_hdf5(original, tmp_path / level.value)

    with h5py.File(output, "r") as h5:
        fields = h5["results/fields"]
        assert "phase" in fields
        assert ("group" in fields) is (level is SamplingLevel.GROUP)
        assert "enhancement" not in fields
        assert "tracking" not in fields

    restored = read_seismic_hdf5(output)
    payload = restored.results["seismic"]
    assert isinstance(payload, SeismicResult)
    assert payload.field.level is level
    assert (payload.field.group is not None) is (level is SamplingLevel.GROUP)
    assert payload.field.enhancement is None
    assert payload.field.tracking is None


@pytest.mark.module
@pytest.mark.seismic
def test_hdf5_reader_rejects_incompatible_units(tmp_path: Path) -> None:
    result = run_seismic(DATA, options=_options(SamplingLevel.PHASE, tracking=False))
    output = write_seismic_hdf5(result, tmp_path / "bad_units")
    with h5py.File(output, "r+") as h5:
        h5["results/fields/phase/phase_speeds"].attrs["unit"] = "m s^-1"

    with pytest.raises(ValueError, match="expected 'km s\\^-1'"):
        read_seismic_hdf5(output)


@pytest.mark.module
@pytest.mark.seismic
def test_hdf5_reader_rejects_another_module(tmp_path: Path) -> None:
    result = run_seismic(DATA, options=_options(SamplingLevel.PHASE, tracking=False))
    output = write_seismic_hdf5(result, tmp_path / "wrong_module")
    with h5py.File(output, "r+") as h5:
        h5["metadata"].attrs["module"] = "elasticity"

    with pytest.raises(ValueError, match="not 'seismic'"):
        read_seismic_hdf5(output)


@pytest.mark.module
@pytest.mark.seismic
def test_hdf5_writer_rejects_a_non_seismic_result(tmp_path: Path) -> None:
    result = ResultData(metadata=ResultMetadata(module="elasticity"))
    with pytest.raises(ValueError, match="valid seismic result"):
        write_seismic_hdf5(result, tmp_path / "invalid")


def _assert_result_equal(source: SeismicResult, target: SeismicResult) -> None:
    assert target.jobname == source.jobname
    assert target.density == pytest.approx(source.density)
    np.testing.assert_allclose(target.stiffness, source.stiffness)
    assert target.stability.is_stable is source.stability.is_stable
    np.testing.assert_allclose(
        target.stability.eigenvalues,
        source.stability.eigenvalues,
    )
    np.testing.assert_allclose(target.averages.as_array(), source.averages.as_array())
    np.testing.assert_allclose(
        target.isotropic_velocities.as_array(),
        source.isotropic_velocities.as_array(),
    )
    assert target.metadata == source.metadata
    assert target.grid.hemisphere is source.grid.hemisphere
    np.testing.assert_allclose(target.grid.theta, source.grid.theta)
    np.testing.assert_allclose(target.grid.phi, source.grid.phi)
    np.testing.assert_allclose(target.grid.directions, source.grid.directions)
    assert target.field.level is source.field.level
    assert target.field.batch_size == source.field.batch_size

    for name in (
        "eigenvalues",
        "phase_speeds",
        "polarizations",
        "mode_eigenvalue_gaps",
        "mode_relative_eigenvalue_gaps",
        "pair_eigenvalue_gaps",
        "pair_relative_eigenvalue_gaps",
        "eigenvalue_thresholds",
        "degeneracy_thresholds",
    ):
        np.testing.assert_allclose(
            getattr(target.field.phase, name),
            getattr(source.field.phase, name),
            equal_nan=True,
        )
    for name in (
        "valid_mask",
        "clamped_mask",
        "degeneracy_mask",
        "pair_degeneracy_mask",
    ):
        np.testing.assert_array_equal(
            getattr(target.field.phase, name),
            getattr(source.field.phase, name),
        )

    assert source.field.group is not None
    assert target.field.group is not None
    for name in (
        "eigenvalue_gradients",
        "group_velocities",
        "group_speeds",
        "ray_directions",
        "power_flow_angles",
    ):
        np.testing.assert_allclose(
            getattr(target.field.group, name),
            getattr(source.field.group, name),
            equal_nan=True,
        )
    for name in ("valid_mask", "resolved_mask"):
        np.testing.assert_array_equal(
            getattr(target.field.group, name),
            getattr(source.field.group, name),
        )

    assert source.field.enhancement is not None
    assert target.field.enhancement is not None
    for name in (
        "eigenvalue_hessians",
        "ray_direction_gradients",
        "area_factors",
        "caustic_thresholds",
        "log10_enhancement",
    ):
        np.testing.assert_allclose(
            getattr(target.field.enhancement, name),
            getattr(source.field.enhancement, name),
            equal_nan=True,
        )
    np.testing.assert_allclose(
        target.field.enhancement.enhancement,
        source.field.enhancement.enhancement,
        rtol=2.0e-15,
        atol=0.0,
        equal_nan=True,
    )
    for name in (
        "valid_mask",
        "resolved_mask",
        "finite_mask",
        "caustic_candidate_mask",
    ):
        np.testing.assert_array_equal(
            getattr(target.field.enhancement, name),
            getattr(source.field.enhancement, name),
        )

    assert source.field.tracking is not None
    assert target.field.tracking is not None
    for name in ("polarizations", "continuity_scores"):
        np.testing.assert_allclose(
            getattr(target.field.tracking, name),
            getattr(source.field.tracking, name),
            equal_nan=True,
        )
    for name in (
        "branch_mode_indices",
        "sign_flip_mask",
        "resolved_mask",
        "segment_start_mask",
        "shear_swap_mask",
        "shear_permutation_ambiguous_mask",
        "shear_subspace_rotation_mask",
    ):
        np.testing.assert_array_equal(
            getattr(target.field.tracking, name),
            getattr(source.field.tracking, name),
        )
