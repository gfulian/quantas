from __future__ import annotations

import numpy as np

from quantas.models import InputData, ResultData, ResultMetadata
from quantas.core.math.fitting import FitQuality, FitResult, FitStatus
from quantas.modules.qha.io.export import QHAHDF5Export, QHATableExport
from quantas.modules.qha.io.hdf5 import QHAHDF5Reader, read_qha_hdf5
from quantas.modules.qha.models import QHAFitRecord, QHAResult


def _result_data() -> ResultData:
    temperature = np.array([0.0, 100.0])
    pressure = np.array([0.0, 1.0])
    values = np.array([[18.0, 17.9], [18.1, 18.0]])
    qha = QHAResult(
        jobname="MgO",
        temperature=temperature,
        pressure=pressure,
        volume=np.array([17.5, 18.0, 18.5]),
        static_energy=np.array([-10.0, -10.1, -10.0]),
        equilibrium_volume=values,
        zero_point_energy=np.full_like(values, 0.01),
        thermal_energy=np.full_like(values, 0.02),
        internal_energy=np.full_like(values, -10.0),
        entropy=np.full_like(values, 1.0e-4),
        isochoric_heat_capacity=np.full_like(values, 2.0e-4),
        isobaric_heat_capacity=np.full_like(values, 2.1e-4),
        heat_capacity_difference=np.full_like(values, 1.0e-5),
        vibrational_free_energy=np.full_like(values, -0.03),
        free_energy=np.full_like(values, -10.03),
        enthalpy=np.full_like(values, -10.02),
        gibbs_free_energy=np.full_like(values, -10.04),
        isothermal_bulk_modulus=np.array([[160.0, 170.0], [158.0, 168.0]]),
        adiabatic_bulk_modulus=np.array([[161.0, 171.0], [159.0, 169.0]]),
        bulk_modulus_derivative=np.full_like(values, 4.0),
        thermal_expansion=np.full_like(values, 2.0e-5),
        thermal_expansion_mixed=np.full_like(values, 2.0e-5),
        thermal_expansion_mode=np.full_like(values, 2.1e-5),
        thermal_expansion_numerical=np.full_like(values, 1.9e-5),
        thermal_expansion_source=np.ones_like(values, dtype=np.int8),
        gruneisen=np.full_like(values, 1.2),
        uncertainties={"sigma_VT": np.full_like(values, 0.001)},
        valid_mask=np.ones_like(values, dtype=bool),
        metadata={
            "units": {
                "energy": "Ha",
                "pressure": "GPa",
                "temperature": "K",
                "volume": "A^3",
                "entropy": "Ha cell^-1 K^-1",
                "heat_capacity": "Ha cell^-1 K^-1",
            },
            "thermodynamic_unit_convention": ("native_energy_per_cell_per_kelvin"),
            "normalization": {
                "native_basis": "cell",
                "formula_units_per_cell": 2,
                "natoms_per_cell": 4,
                "natoms_per_formula_unit": 2.0,
                "molar_basis": "formula_unit",
            },
            "state_variables": ["VT", "KT"],
        },
    )
    fit = FitResult(
        success=True,
        status=FitStatus.SUCCESS,
        quality=FitQuality.GOOD,
        parameters=np.array([1.0, 2.0]),
        errors=np.array([0.1, 0.2]),
        r_squared=0.999,
        rmse=1.0e-6,
        metadata={"parameter_names": ["E0", "K0"]},
    )
    qha.add_fit_record(
        QHAFitRecord(
            quantity="volume",
            method="eos",
            temperature=0.0,
            pressure=0.0,
            fit=fit,
            metadata={"parameter_names": ["E0", "K0"]},
        )
    )
    return ResultData(
        metadata=ResultMetadata(module="qha", method="freq-eos"),
        input_data=InputData(
            source="input.yaml",
            data={"jobname": "MgO", "natoms": 4, "formula_units": 2},
        ),
        options={"debug": True, "scheme": "freq"},
        results={"qha": qha},
        warnings=["test warning"],
    )


def test_hdf5_roundtrip_preserves_qha_fields(tmp_path) -> None:
    filename = tmp_path / "result.hdf5"
    QHAHDF5Export().export(_result_data(), filename, report_text="report")

    result = read_qha_hdf5(filename)
    qha = result.results["qha"]

    assert result.metadata.module == "qha"
    assert result.options["scheme"] == "freq"
    assert result.input_data is not None
    assert result.input_data.data["formula_units"] == 2
    assert result.warnings == ["test warning"]
    assert qha.jobname == "MgO"
    np.testing.assert_allclose(qha.equilibrium_volume, [[18.0, 17.9], [18.1, 18.0]])
    np.testing.assert_allclose(
        qha.isothermal_bulk_modulus, [[160.0, 170.0], [158.0, 168.0]]
    )
    np.testing.assert_allclose(qha.thermal_expansion_mixed, 2.0e-5)
    np.testing.assert_allclose(qha.thermal_expansion_mode, 2.1e-5)
    np.testing.assert_allclose(qha.thermal_expansion_numerical, 1.9e-5)
    np.testing.assert_array_equal(qha.thermal_expansion_source, 1)
    np.testing.assert_allclose(qha.uncertainties["sigma_VT"], 0.001)
    assert qha.metadata["units"]["pressure"] == "GPa"
    assert qha.metadata["normalization"]["formula_units_per_cell"] == 2
    assert (
        qha.metadata["thermodynamic_unit_convention"]
        == "native_energy_per_cell_per_kelvin"
    )
    assert qha.fit_records[0].fit is not None
    np.testing.assert_allclose(qha.fit_records[0].fit.parameters, [1.0, 2.0])


def test_hdf5_reader_sets_completed_flag(tmp_path) -> None:
    filename = tmp_path / "result.hdf5"
    QHAHDF5Export().export(_result_data(), filename)

    reader = QHAHDF5Reader()
    result = reader.load(filename)

    assert reader.completed is True
    assert result.results["qha"].completed is True


def test_table_export_from_loaded_hdf5_result(tmp_path) -> None:
    hdf5_file = tmp_path / "result.hdf5"
    table_file = tmp_path / "VT.dat"
    QHAHDF5Export().export(_result_data(), hdf5_file)

    result = read_qha_hdf5(hdf5_file)
    QHATableExport().export(result, table_file, property_name="VT")

    text = table_file.read_text()
    assert "(K)" in text
    assert "Pressure (GPa)" in text
    assert "(A^3)" in text
    assert "18.000000" in text


def test_reader_migrates_schema_1_0_implicit_milli_units(tmp_path) -> None:
    result = _result_data()
    result.metadata.schema_version = "1.0"
    qha = result.results["qha"]
    qha.metadata.pop("thermodynamic_unit_convention", None)
    qha.metadata.pop("normalization", None)
    qha.metadata["units"].pop("entropy", None)
    qha.metadata["units"].pop("heat_capacity", None)
    old_entropy = np.asarray(qha.entropy).copy()
    old_cv = np.asarray(qha.isochoric_heat_capacity).copy()

    filename = tmp_path / "schema_1_0.hdf5"
    QHAHDF5Export().export(result, filename)
    loaded = read_qha_hdf5(filename)
    migrated = loaded.results["qha"]

    assert loaded.metadata.schema_version == "1.1"
    np.testing.assert_allclose(migrated.entropy, old_entropy * 1.0e-3)
    np.testing.assert_allclose(migrated.isochoric_heat_capacity, old_cv * 1.0e-3)
    assert migrated.metadata["units"]["entropy"] == "Ha cell^-1 K^-1"
    assert migrated.metadata["unit_migration"]["scale_factor"] == 1.0e-3
    assert any("implicit milli-unit scale" in item for item in loaded.warnings)
