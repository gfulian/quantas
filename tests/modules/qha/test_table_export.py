"""Tests for QHA table export utilities."""

from __future__ import annotations

import csv
import re

import numpy as np

from quantas.models import ResultData, ResultMetadata
from quantas.modules.qha.io.export import QHATableExport
from quantas.modules.qha.models import QHAResult
from quantas.modules.qha.thermodynamics import _calculate_heat_capacity_correction
from quantas.modules.qha.models import QHAOptions


def _result_data() -> ResultData:
    qha = QHAResult(
        jobname="MgO",
        temperature=np.array([0.0, 10.0, 20.0]),
        pressure=np.array([0.0, 1.0]),
        equilibrium_volume=np.array([[10.0, 9.8], [10.1, 9.9], [10.2, 10.0]]),
        isothermal_bulk_modulus=np.array(
            [[100.0, 110.0], [99.0, 109.0], [98.0, 108.0]]
        ),
        adiabatic_bulk_modulus=np.array([[100.0, 110.0], [99.5, 109.5], [98.5, 108.5]]),
        bulk_modulus_derivative=np.full((3, 2), 4.0),
        thermal_expansion=np.full((3, 2), 1.0e-5),
        free_energy=np.array([[-1.0, -0.9], [-0.8, -0.7], [-0.6, -0.5]]),
        uncertainties={"sigma_VT": np.full((3, 2), 0.01)},
        metadata={
            "units": {
                "temperature": "K",
                "pressure": "GPa",
                "volume": "A^3",
                "energy": "Ha",
            }
        },
    )
    return ResultData(
        metadata=ResultMetadata(module="qha", method="qha"), results={"qha": qha}
    )


def _structural_result_data() -> ResultData:
    """Return one QHA result with complete structural thermal expansion."""
    result = _result_data()
    qha = result.results["qha"]
    shape = (3, 2)
    parameters = np.zeros(shape + (6,), dtype=np.float64)
    parameters[..., 0] = 4.8
    parameters[..., 1] = 4.8
    parameters[..., 2] = 16.2
    parameters[..., 3] = 90.0
    parameters[..., 4] = 90.0
    parameters[..., 5] = 120.0
    axial = np.zeros(shape + (3,), dtype=np.float64)
    axial[..., 0] = 5.0e-6
    axial[..., 1] = 5.0e-6
    axial[..., 2] = 20.0e-6
    tensor = np.zeros(shape + (3, 3), dtype=np.float64)
    tensor[..., 0, 0] = axial[..., 0]
    tensor[..., 1, 1] = axial[..., 1]
    tensor[..., 2, 2] = axial[..., 2]
    qha.lattice_parameters = parameters
    qha.axial_thermal_expansion = axial
    qha.thermal_expansion_tensor = tensor
    qha.thermal_expansion = np.trace(tensor, axis1=-2, axis2=-1)
    qha.structural_extrapolation_mask = np.asarray(
        [[False, True], [False, False], [False, False]],
        dtype=bool,
    )
    qha.uncertainties["sigma_cell_parameters"] = np.full(
        shape + (6,),
        1.0e-3,
        dtype=np.float64,
    )
    qha.metadata["structural_thermal_expansion"] = {
        "uncertainty_method": "first_order_delta_method",
    }
    return result


def test_single_property_export_is_grouped_by_pressure_then_temperature(
    tmp_path,
) -> None:
    filename = tmp_path / "volume.dat"

    QHATableExport().export(_result_data(), filename, property_name="VT")

    text = filename.read_text(encoding="utf8")
    assert "Pressure (GPa) = 0.00" in text
    assert "Pressure (GPa) = 1.00" in text
    assert "sigma_VT" in text
    rows = [
        line.split()
        for line in text.splitlines()
        if re.match(r"^\s+\d+\.\d{2}\s+\d+\.\d{6}", line)
    ]
    assert rows[:3] == [
        ["0.00", "10.000000", "0.010000"],
        ["10.00", "10.100000", "0.010000"],
        ["20.00", "10.200000", "0.010000"],
    ]
    assert rows[3] == ["0.00", "9.800000", "0.010000"]


def test_full_export_includes_available_pressure_temperature_properties(
    tmp_path,
) -> None:
    filename = tmp_path / "all.csv"

    QHATableExport().export(_result_data(), filename, property_name=None, delimiter=",")

    text = filename.read_text(encoding="utf8")
    assert not text.startswith("#")
    header = text.splitlines()[0]
    assert header.startswith("P (GPa),T (K),VT (A^3),sigma_VT (A^3)")
    assert "VT (A^3)" in header
    assert "KT (GPa)" in header
    assert "KS (GPa)" in header
    assert "F (Ha cell^-1)" in header

    with filename.open(newline="", encoding="utf8") as stream:
        rows = list(csv.reader(stream))
    assert all(len(row) == len(rows[0]) for row in rows)
    assert rows[1][:4] == ["0.00", "0.00", "10.000000", "0.010000"]
    assert rows[4][:4] == ["1.00", "0.00", "9.800000", "0.010000"]


def test_full_csv_export_includes_structural_thermal_expansion(tmp_path) -> None:
    filename = tmp_path / "structural_all.csv"

    QHATableExport().export(
        _structural_result_data(),
        filename,
        property_name=None,
        delimiter=",",
    )

    with filename.open(newline="", encoding="utf8") as stream:
        rows = list(csv.reader(stream))
    header = rows[0]
    assert header.index("a (A)") < header.index("KT (GPa)")
    assert "sigma_a (A)" in header
    assert "gamma (degree)" in header
    assert "alpha_a (K^-1)" in header
    assert "trace(alpha)-alphaV (K^-1)" in header
    assert "alpha_23 (K^-1)" in header
    assert "structural_extrapolated" in header
    assert all(len(row) == len(header) for row in rows)

    first = dict(zip(header, rows[1]))
    assert first["a (A)"] == "4.80000000"
    assert first["sigma_a (A)"] == "0.00100000"
    assert first["alpha_a (K^-1)"] == "5.000000E-06"
    assert first["trace(alpha)-alphaV (K^-1)"] == "0.000000E+00"
    assert first["structural_extrapolated"] == "0"


def test_human_export_splits_structural_data_into_readable_blocks(tmp_path) -> None:
    filename = tmp_path / "structural_all.dat"

    QHATableExport().export(_structural_result_data(), filename)

    text = filename.read_text(encoding="utf8")
    thermoelastic = text.index("Thermoelastic properties")
    lengths = text.index("Equilibrium lattice lengths")
    angles = text.index("Equilibrium lattice angles")
    axial = text.index("Axial thermal expansion")
    tensor = text.index("Cartesian thermal-expansion tensor")
    assert thermoelastic < lengths < angles < axial < tensor
    assert "sigma_a" in text
    assert "trace-alphaV" in text
    assert "alpha_12 x 10^5" in text
    assert "# Structural uncertainty method: first_order_delta_method" in text


def test_structure_property_alias_exports_only_structural_tables(tmp_path) -> None:
    text_file = tmp_path / "structure.dat"
    csv_file = tmp_path / "structure.csv"
    result = _structural_result_data()

    QHATableExport().export(result, text_file, property_name="structure")
    QHATableExport().export(
        result,
        csv_file,
        property_name="lattice",
        delimiter=",",
    )

    text = text_file.read_text(encoding="utf8")
    assert "Equilibrium lattice lengths" in text
    assert "Thermodynamic potentials" not in text
    with csv_file.open(newline="", encoding="utf8") as stream:
        header = next(csv.reader(stream))
    assert header[:3] == ["P (GPa)", "T (K)", "a (A)"]
    assert "KT (GPa)" not in header


def test_export_reports_thermal_expansion_method_per_point(tmp_path) -> None:
    result = _result_data()
    qha = result.results["qha"]
    qha.thermal_expansion_source = np.array([[1, 1], [1, 4], [1, 4]], dtype=np.int8)
    qha.metadata["thermal_expansion"] = {
        "requested_method": "mixed_derivative",
        "selected_method": "mixed_derivative+numerical_fallback",
        "fallback_method": "numerical",
        "source_codes": {
            "invalid": 0,
            "mixed_derivative": 1,
            "mode_gruneisen": 2,
            "numerical": 3,
            "numerical_fallback": 4,
        },
        "source_counts": {
            "invalid": 0,
            "mixed_derivative": 4,
            "mode_gruneisen": 0,
            "numerical": 0,
            "numerical_fallback": 2,
        },
    }
    csv_file = tmp_path / "alpha.csv"
    text_file = tmp_path / "alpha.dat"

    QHATableExport().export(
        result,
        csv_file,
        property_name="alphaV",
        delimiter=",",
    )
    QHATableExport().export(result, text_file, property_name="alphaV")

    with csv_file.open(newline="", encoding="utf8") as stream:
        rows = list(csv.reader(stream))
    assert rows[0] == ["P (GPa)", "T (K)", "alphaV (K^-1)", "alphaV method"]
    assert rows[1][-1] == "mixed_derivative"
    assert rows[5][-1] == "numerical_fallback"

    text = text_file.read_text(encoding="utf8")
    assert "# Thermal-expansion requested method: mixed_derivative" in text
    assert (
        "# Thermal-expansion effective method: mixed_derivative+numerical_fallback"
    ) in text
    assert "alphaV method" in text
    assert "numerical_fallback" in text


def test_bulk_moduli_match_at_zero_kelvin() -> None:
    qha = QHAResult(
        temperature=np.array([0.0, 10.0]),
        pressure=np.array([0.0, 1.0]),
        equilibrium_volume=np.ones((2, 2)),
        isothermal_bulk_modulus=np.array([[100.0, 110.0], [99.0, 109.0]]),
        thermal_expansion=np.full((2, 2), 1.0e-5),
        isochoric_heat_capacity=np.array([[0.0, 0.0], [1.0, 1.0]]),
    )

    _calculate_heat_capacity_correction(qha, QHAOptions())

    np.testing.assert_allclose(
        qha.adiabatic_bulk_modulus[0], qha.isothermal_bulk_modulus[0]
    )


def test_heat_capacity_correction_is_stored_in_native_energy_per_kelvin() -> None:
    from scipy import constants as cs

    qha = QHAResult(
        temperature=np.array([300.0]),
        pressure=np.array([0.0]),
        equilibrium_volume=np.array([[10.0]]),
        isothermal_bulk_modulus=np.array([[100.0]]),
        thermal_expansion=np.array([[1.0e-5]]),
        isochoric_heat_capacity=np.array([[1.0e-5]]),
    )

    _calculate_heat_capacity_correction(qha, QHAOptions())

    expected_j_mol_k = (1.0e-5**2) * 100.0e9 * 10.0e-30 * cs.Avogadro * 300.0
    expected_ha_k = expected_j_mol_k / (
        cs.Avogadro * cs.physical_constants["Hartree energy"][0]
    )
    np.testing.assert_allclose(qha.heat_capacity_difference, [[expected_ha_k]])


def test_full_export_writes_all_available_eos_uncertainties(tmp_path) -> None:
    result = _result_data()
    qha = result.results["qha"]
    qha.uncertainties.update(
        {
            "sigma_KT": np.full((3, 2), 0.5),
            "sigma_Kp": np.full((3, 2), 0.02),
        }
    )
    qha.metadata["eos_workflow"] = {"uncertainty_method": "covariance"}
    filename = tmp_path / "all_uncertainties.dat"

    QHATableExport().export(result, filename)

    text = filename.read_text(encoding="utf8")
    assert "sigma_VT" in text
    assert "sigma_KT" in text
    assert "sigma_Kp" in text
    assert "(A^3)" in text
    assert "(GPa)" in text
    assert "(-)" in text
    assert "# sigma_* columns contain one-standard-deviation uncertainties." in text
    assert "# Uncertainty method: covariance" in text


def test_human_export_uses_sign_aware_fixed_width_energy_columns(tmp_path) -> None:
    result = _result_data()
    qha = result.results["qha"]
    qha.free_energy = np.array([[1.0, -1.0], [2.0, -2.0], [3.0, -3.0]])
    filename = tmp_path / "energy.dat"

    QHATableExport().export(result, filename, property_name="F")

    text = filename.read_text(encoding="utf8")
    rows = [
        line
        for line in text.splitlines()
        if re.match(r"^\s+\d+\.\d{2}\s+[ +-]\d\.\d{12}E[+-]\d{2}", line)
    ]
    positive_field = rows[0][11:31]
    negative_field = rows[3][11:31]
    assert len(positive_field) == 20
    assert len(negative_field) == 20
    assert positive_field.index(".") == negative_field.index(".")
    assert positive_field.strip() == "1.000000000000E+00"
    assert negative_field.strip() == "-1.000000000000E+00"
