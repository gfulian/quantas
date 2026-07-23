"""QHA integration tests for structural thermal expansion."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from quantas.models.structures import (
    CellNormalization,
    CrystalStructure,
    StructureVolumeSeries,
)
from quantas.modules.qha.calculator import QHACalculator
from quantas.modules.qha.io.hdf5_payload import read_qha_payload, write_qha_payload
from quantas.modules.qha.models import QHAInput, QHAOptions, QHAResult
from quantas.modules.qha.report import structural_property_tables
from quantas.modules.qha.structural import calculate_structural_thermal_expansion

pytestmark = [pytest.mark.qha, pytest.mark.scientific, pytest.mark.fast]


def _cubic_series() -> StructureVolumeSeries:
    edge = np.asarray([3.9, 4.0, 4.1, 4.2], dtype=np.float64)
    lattice = np.asarray([np.eye(3) * value for value in edge])
    volume = edge**3
    reference = CrystalStructure(
        lattice=lattice[1],
        fractional_positions=np.asarray([[0.0, 0.0, 0.0]]),
        atomic_numbers=np.asarray([12]),
    )
    return StructureVolumeSeries(
        reference=reference,
        lattices=lattice,
        fractional_positions=np.zeros((4, 1, 3), dtype=np.float64),
        volumes=volume,
        normalization=CellNormalization(
            basis="primitive",
            source_basis="primitive",
            expansion_matrix=np.eye(3, dtype=np.int64),
            repetitions=1,
            source_atoms=1,
            normalized_atoms=1,
        ),
        reference_index=1,
    )


def _result(series: StructureVolumeSeries) -> QHAResult:
    equilibrium = np.asarray(
        [[series.volumes[1]], [series.volumes[2]]],
        dtype=np.float64,
    )
    alpha_v = np.asarray([[1.2e-5], [1.8e-5]], dtype=np.float64)
    return QHAResult(
        jobname="structural qha",
        temperature=np.asarray([300.0, 600.0]),
        pressure=np.asarray([0.0]),
        volume=series.volumes,
        equilibrium_volume=equilibrium,
        thermal_expansion=alpha_v,
        isochoric_heat_capacity=np.asarray([[2.0e-5], [3.0e-5]]),
        uncertainties={
            "sigma_VT": np.asarray([[2.0e-3], [3.0e-3]], dtype=np.float64),
        },
        metadata={
            "units": {"temperature": "K", "pressure": "GPa", "volume": "A"},
            "thermal_expansion": {"selected_method": "mixed_derivative"},
        },
    )


def test_structural_analysis_populates_tensor_results() -> None:
    series = _cubic_series()
    result = _result(series)
    calculate_structural_thermal_expansion(result, series, QHAOptions())

    assert result.lattice_parameters is not None
    assert result.lattice_parameters.shape == (2, 1, 6)
    assert result.axial_thermal_expansion is not None
    np.testing.assert_allclose(
        result.axial_thermal_expansion[:, 0],
        np.broadcast_to(result.thermal_expansion[:, 0, None] / 3.0, (2, 3)),
        atol=1.0e-14,
    )
    assert result.thermal_expansion_tensor is not None
    np.testing.assert_allclose(
        np.trace(result.thermal_expansion_tensor, axis1=-2, axis2=-1),
        result.thermal_expansion,
        atol=1.0e-14,
    )
    metadata = result.metadata["structural_thermal_expansion"]
    assert metadata["equilibrium_volume_source"] == "qha_free_energy_minimization"
    assert metadata["full_anisotropic_qha"] is False
    assert metadata["calculation_branch"] == "cubic_exact"
    sigma_parameters = result.uncertainties["sigma_cell_parameters"]
    assert sigma_parameters.shape == (2, 1, 6)
    expected = (
        result.lattice_parameters[:, 0, 0]
        * result.uncertainties["sigma_VT"][:, 0]
        / (3.0 * result.equilibrium_volume[:, 0])
    )
    np.testing.assert_allclose(
        sigma_parameters[:, 0, :3],
        np.broadcast_to(expected[:, None], (2, 3)),
    )
    np.testing.assert_allclose(sigma_parameters[:, 0, 3:], 0.0)


def test_structural_qha_payload_round_trip(tmp_path) -> None:
    series = _cubic_series()
    result = _result(series)
    calculate_structural_thermal_expansion(result, series, QHAOptions())
    destination = tmp_path / "structural_qha.hdf5"

    with h5py.File(destination, "w") as handle:
        write_qha_payload(handle, result)
    with h5py.File(destination, "r") as handle:
        restored = read_qha_payload(handle["results"])

    np.testing.assert_allclose(restored.temperature, result.temperature)
    np.testing.assert_allclose(restored.pressure, result.pressure)
    np.testing.assert_allclose(
        restored.equilibrium_volume,
        result.equilibrium_volume,
    )
    np.testing.assert_allclose(
        restored.isochoric_heat_capacity,
        result.isochoric_heat_capacity,
    )
    np.testing.assert_allclose(
        restored.thermal_expansion,
        result.thermal_expansion,
    )
    np.testing.assert_allclose(
        restored.equilibrium_lattice,
        result.equilibrium_lattice,
    )
    np.testing.assert_allclose(restored.lattice_parameters, result.lattice_parameters)
    np.testing.assert_allclose(
        restored.lattice_parameter_derivatives,
        result.lattice_parameter_derivatives,
    )
    np.testing.assert_allclose(
        restored.axial_thermal_expansion,
        result.axial_thermal_expansion,
    )
    np.testing.assert_allclose(
        restored.thermal_expansion_tensor,
        result.thermal_expansion_tensor,
    )
    np.testing.assert_array_equal(
        restored.structural_extrapolation_mask,
        result.structural_extrapolation_mask,
    )
    np.testing.assert_allclose(
        restored.uncertainties["sigma_cell_parameters"],
        result.uncertainties["sigma_cell_parameters"],
    )


def test_structural_report_tables_are_created_per_pressure() -> None:
    series = _cubic_series()
    result = _result(series)
    calculate_structural_thermal_expansion(result, series, QHAOptions())

    tables = structural_property_tables(result)
    assert len(tables) == 3
    assert tables[0].columns == [
        "T",
        "V",
        "a",
        "b",
        "c",
        "alpha",
        "beta",
        "gamma",
    ]
    assert "alpha_a" in tables[1].columns[1]
    assert "sigma_a" in tables[2].columns


def test_qha_calculator_runs_structural_analysis_automatically() -> None:
    volume = np.linspace(9.5, 10.8, 7)
    edge = np.cbrt(volume)
    lattices = np.asarray([np.eye(3) * value for value in edge])
    reference = CrystalStructure(
        lattice=lattices[2],
        fractional_positions=np.asarray([[0.0, 0.0, 0.0]]),
        atomic_numbers=np.asarray([12]),
    )
    structure = StructureVolumeSeries(
        reference=reference,
        lattices=lattices,
        fractional_positions=np.zeros((volume.size, 1, 3), dtype=np.float64),
        volumes=volume,
        normalization=CellNormalization(
            basis="primitive",
            source_basis="primitive",
            expansion_matrix=np.eye(3, dtype=np.int64),
            repetitions=1,
            source_atoms=1,
            normalized_atoms=1,
        ),
        reference_index=2,
    )
    energy = 2.0e-3 * (volume - 10.0) ** 2
    base = 450.0 * (10.0 / volume) ** 2
    frequencies = np.stack([base, 1.15 * base, 1.4 * base], axis=0)[None, :, :]
    input_data = QHAInput(
        jobname="automatic structural QHA",
        natoms=1,
        qpoints=1,
        volume=volume,
        energy=energy,
        frequencies=frequencies,
        weights=np.asarray([1.0]),
        structure=structure,
    )
    options = QHAOptions(
        temperature_min=10.0,
        temperature_max=1000.0,
        temperature_step=990.0,
        pressure_min=0.0,
        pressure_max=0.0,
        pressure_step=1.0,
        scheme="td",
        minimization="poly",
        free_energy_degree=3,
    )

    result = QHACalculator(input_data, options).execute().results["qha"]

    assert result.equilibrium_lattice is not None
    assert result.lattice_parameters is not None
    assert result.axial_thermal_expansion is not None
    assert result.thermal_expansion_tensor is not None
    metadata = result.metadata["structural_thermal_expansion"]
    assert metadata["automatic"] is True
    assert metadata["full_anisotropic_qha"] is False
