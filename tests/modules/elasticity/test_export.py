"""HDF5 and text-export characterization for elasticity results."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from quantas.core.physics.elasticity import (
    DirectionalExtrema,
    ElasticAverages,
    IsotropicElasticProperties,
    StabilityResult,
)
from quantas.models import InputData, ResultData, ResultMetadata
from quantas.modules.elasticity.io.export import (
    ElasticityHDF5Export,
    ElasticityTableExport,
)
from quantas.modules.elasticity.io.reader import (
    ElasticityHDF5Reader,
    read_elasticity_hdf5,
)
from quantas.modules.elasticity.models import ElasticityResult


def _properties(value: float) -> IsotropicElasticProperties:
    """Build deterministic isotropic properties for serialization tests."""
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


def _result_data(with_polar: bool = True) -> ResultData:
    """Create a generic Quantas result with structured elasticity data."""
    angle = np.linspace(0.0, 2.0 * np.pi, 5)
    variation = DirectionalExtrema(
        minimum=100.0,
        maximum=200.0,
        anisotropy=2.0,
        minimum_axis=[1.0, 0.0, 0.0],
        maximum_axis=[0.0, 1.0, 0.0],
        minimum_measurement_axis=[0.0, 0.0, 1.0],
        maximum_measurement_axis=[0.0, 0.0, -1.0],
    )
    polar: dict[str, dict[str, np.ndarray]] = {}
    if with_polar:
        polar = {
            "xy": {
                "theta": np.full(5, np.pi / 2.0),
                "phi": angle,
                "young_modulus": np.linspace(100.0, 120.0, 5),
                "shear_modulus": np.column_stack(
                    (np.linspace(40.0, 45.0, 5), np.linspace(50.0, 55.0, 5))
                ),
            }
        }
    elasticity = ElasticityResult(
        jobname="Synthetic",
        crystal_system="orthorhombic",
        stiffness=np.eye(6) * 100.0,
        compliance=np.eye(6) * 0.01,
        averages=_averages(),
        stability=StabilityResult(True, np.arange(1.0, 7.0)),
        variations={"young_modulus": variation},
        properties_2d=polar,
        metadata={"tensor_class": "OrthorhombicElasticTensor"},
    )
    return ResultData(
        metadata=ResultMetadata(module="elasticity", method="second_order"),
        input_data=InputData(
            source="synthetic.dat",
            data={"jobname": "Synthetic"},
        ),
        options={"calculate_2d": with_polar, "ntheta": 5},
        results={"elasticity": elasticity},
        warnings=["example warning"],
    )


def test_hdf5_roundtrip_preserves_elasticity_result(tmp_path) -> None:
    """Structured elasticity arrays and metadata survive an HDF5 roundtrip."""
    filename = tmp_path / "elasticity.hdf5"
    exporter = ElasticityHDF5Export()
    exporter.export(_result_data(), filename, report_text="Elasticity report")

    assert exporter.completed is True
    with h5py.File(filename, "r") as h5:
        assert h5["metadata"].attrs["module"] == "elasticity"
        assert h5["metadata"].attrs["method"] == "second_order"
        assert h5["options"].attrs["ntheta"] == 5
        assert h5["diagnostics/report_text"].asstr()[()] == "Elasticity report"
        assert "density" not in h5["results"].attrs
        assert "isotropic_velocities" not in h5["results"]

    result_data = read_elasticity_hdf5(filename)
    assert result_data.metadata.module == "elasticity"
    assert result_data.metadata.method == "second_order"
    assert result_data.input_data is not None
    assert result_data.input_data.source == "synthetic.dat"
    assert "density" not in result_data.input_data.data
    assert result_data.options["ntheta"] == 5
    assert result_data.warnings == ["example warning"]

    loaded = result_data.results["elasticity"]
    assert loaded.jobname == "Synthetic"
    assert loaded.crystal_system == "orthorhombic"
    np.testing.assert_allclose(loaded.stiffness, np.eye(6) * 100.0)
    assert loaded.averages is not None
    np.testing.assert_allclose(loaded.averages.as_array(), _averages().as_array())
    assert loaded.stability is not None
    assert loaded.stability.is_stable is True
    np.testing.assert_allclose(loaded.stability.eigenvalues, np.arange(1.0, 7.0))
    assert not hasattr(loaded, "isotropic_velocities")
    assert loaded.variations["young_modulus"].minimum_measurement_axis == [
        0.0,
        0.0,
        1.0,
    ]
    np.testing.assert_allclose(
        loaded.properties_2d["xy"]["shear_modulus"][:, 1],
        np.linspace(50.0, 55.0, 5),
    )
    assert loaded.metadata["tensor_class"] == "OrthorhombicElasticTensor"


def test_hdf5_reader_rejects_non_elasticity_file(tmp_path) -> None:
    """Metadata validation prevents loading unrelated HDF5 files."""
    filename = tmp_path / "other.hdf5"
    with h5py.File(filename, "w") as h5:
        metadata = h5.create_group("metadata")
        metadata.attrs["program"] = "quantas"
        metadata.attrs["module"] = "ha"
        metadata.attrs["schema_version"] = "2.0"
        for group_name in ("input", "options", "results", "diagnostics", "events"):
            h5.create_group(group_name)

    reader = ElasticityHDF5Reader()
    with pytest.raises(ValueError, match="not 'elasticity'"):
        reader.load(filename)
    assert reader.completed is False
    assert "not 'elasticity'" in reader.error


def test_table_export_uses_degrees_and_all_available_columns(tmp_path) -> None:
    """User-facing polar tables expose angular coordinates in degrees."""
    filename = tmp_path / "polar.dat"
    exporter = ElasticityTableExport()
    exporter.export(_result_data(), filename)

    text = filename.read_text(encoding="utf-8")
    assert exporter.completed is True
    assert "theta_deg" in text
    assert "phi_deg" in text
    assert "shear_modulus_1" in text
    assert "shear_modulus_2" in text
    assert "3.6000000000e+02" in text


def test_table_export_from_hdf5_honors_requested_output_path(tmp_path) -> None:
    """The HDF5 convenience path uses its output_file argument."""
    hdf5_file = tmp_path / "result.hdf5"
    requested = tmp_path / "custom_name.txt"
    ElasticityHDF5Export().export(_result_data(), hdf5_file)

    exporter = ElasticityTableExport()
    exporter.export_from_hdf5(hdf5_file, requested)

    output = requested.with_suffix(".dat")
    assert exporter.completed is True
    assert output.exists()
    assert "# Plane: xy" in output.read_text(encoding="utf-8")


def test_table_export_without_polar_data_is_a_completed_operation(tmp_path) -> None:
    """An explicit no-data table remains a successful export."""
    filename = tmp_path / "no_polar.dat"
    exporter = ElasticityTableExport()
    exporter.export(_result_data(with_polar=False), filename)

    assert exporter.completed is True
    assert "# No 2D directional data available." in filename.read_text(encoding="utf-8")
