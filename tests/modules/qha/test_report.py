from __future__ import annotations

import numpy as np
import pytest

from quantas.modules.qha.inspect import (
    PressureEstimate,
    PressureVolumePreview,
    pressure_volume_preview,
)
from quantas.core.math.fitting import FitQuality, FitResult, FitStatus
from quantas.modules.qha.models import (
    QHAFailedPoint,
    QHAFitRecord,
    QHAInput,
    QHAOptions,
    QHAResult,
)
from quantas.modules.qha.report import (
    all_tables,
    diagnostics_table,
    failed_points_table,
    input_table,
    options_table,
    pressure_volume_preview_table,
    preview_diagnostics_table,
    preview_parameters_table,
    property_table,
    resolve_property_name,
    result_summary_table,
    selected_property_tables,
    thermal_expansion_provenance_table,
)


def _input() -> QHAInput:
    return QHAInput(
        jobname="MgO",
        natoms=2,
        qpoints=1,
        volume=np.array([18.0, 19.0, 20.0]),
        energy=np.array([-10.0, -10.2, -10.1]),
        frequencies=np.ones((1, 6, 3)),
        weights=np.array([1.0]),
    )


def _result() -> QHAResult:
    temperature = np.array([300.0, 400.0])
    pressure = np.array([0.0, 1.0])
    values = np.array([[18.2, 17.9], [18.5, 18.1]])
    return QHAResult(
        jobname="MgO",
        temperature=temperature,
        pressure=pressure,
        volume=np.array([18.0, 19.0, 20.0]),
        static_energy=np.array([-10.0, -10.2, -10.1]),
        equilibrium_volume=values,
        isothermal_bulk_modulus=np.array([[160.0, 165.0], [150.0, 155.0]]),
        uncertainties={"sigma_VT": np.full_like(values, 0.02)},
    )


def test_input_table_reports_dimensions() -> None:
    table = input_table(_input())

    assert table.title == "Input data"
    assert ["Job name", "MgO"] in table.rows
    assert ["Mode continuity", "assumed"] in table.rows
    assert ["Formula units per cell (Z)", 1] in table.rows
    assert ["Atoms per formula unit", 2.0] in table.rows


def test_options_table_reports_methods_and_units() -> None:
    options = QHAOptions(scheme="td", minimization="eos", eos="V")

    table = options_table(options)

    assert ["Scheme", "td"] in table.rows
    assert ["Minimization", "eos"] in table.rows
    assert ["EOS", "V"] in table.rows


def test_result_summary_marks_available_properties_and_uncertainties() -> None:
    table = result_summary_table(_result())

    row = next(row for row in table.rows if row[0] == "VT")
    assert row[2] == "yes"
    assert row[4] == "yes"


def test_property_table_accepts_compatibility_key_and_uncertainty() -> None:
    table = property_table(_result(), "VT")

    assert table.columns == ["T", "P", "VT", "sigma_VT"]
    assert table.rows[0] == [300.0, 0.0, 18.2, 0.02]
    assert table.metadata["has_uncertainty"] is True


def test_property_table_can_limit_rows() -> None:
    table = property_table(_result(), "equilibrium_volume", max_rows=2)

    assert len(table.rows) == 2
    assert table.metadata["truncated"] is True


def test_property_table_rejects_unknown_property() -> None:
    with pytest.raises(KeyError):
        resolve_property_name("bad_property")


def test_diagnostics_table_reports_fit_records() -> None:
    result = _result()
    result.fit_records.append(
        QHAFitRecord(
            quantity="F",
            method="poly",
            temperature=300.0,
            pressure=0.0,
            fit=FitResult(
                success=True,
                status=FitStatus.SUCCESS,
                quality=FitQuality.GOOD,
                parameters=np.array([1.0, 2.0]),
                residuals=np.array([0.0, 0.1]),
                rmse=0.05,
            ),
            message="ok",
        )
    )

    table = diagnostics_table(result)

    assert table.rows[0][0] == "F"
    assert table.rows[0][1] == "poly"
    assert table.rows[0][4] is True


def test_failed_points_table_reports_failed_points() -> None:
    result = _result()
    result.failed_points.append(
        QHAFailedPoint(
            temperature=300.0,
            pressure=5.0,
            stage="minimization",
            message="outside range",
        )
    )

    table = failed_points_table(result)

    assert table.rows == [[300.0, 5.0, "minimization", "outside range"]]


def test_pressure_volume_preview_table_reports_both_estimates() -> None:
    estimate = PressureEstimate(
        method="poly",
        pressure=np.array([1.0, 0.0]),
        fit=FitResult(
            success=True,
            status=FitStatus.SUCCESS,
            quality=FitQuality.GOOD,
            parameters=np.array([1.0]),
        ),
        unit="GPa",
    )
    preview = PressureVolumePreview(
        volume=np.array([18.0, 19.0]),
        energy=np.array([-10.0, -10.2]),
        pressure_unit="GPa",
        polynomial=estimate,
        eos=estimate,
    )

    table = pressure_volume_preview_table(preview)
    diag = preview_diagnostics_table(preview)

    assert table.rows[0] == ["18.000000", "-1.000000000000E+01", "1.00000", "1.00000"]
    assert diag.rows[0][0] == "poly"


def test_all_tables_contains_standard_report_sections() -> None:
    tables = all_tables(_input(), QHAOptions(), _result())

    assert [table.title for table in tables] == [
        "Input data",
        "Selected options",
        "QHA properties",
    ]


def test_thermal_expansion_provenance_reports_effective_sources() -> None:
    result = _result()
    result.thermal_expansion_source = np.array([[1, 1], [1, 4]], dtype=np.int8)
    result.metadata["thermal_expansion"] = {
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
            "mixed_derivative": 3,
            "mode_gruneisen": 0,
            "numerical": 0,
            "numerical_fallback": 1,
        },
    }

    table = thermal_expansion_provenance_table(result)

    assert table is not None
    assert table.title == "Thermal-expansion provenance"
    assert ["Requested method", "mixed_derivative", None, None] in table.rows
    assert [
        "Effective method",
        "mixed_derivative+numerical_fallback",
        None,
        None,
    ] in table.rows
    fallback = next(row for row in table.rows if row[1] == "numerical_fallback")
    assert fallback[2:] == [1, 25.0]


def test_preview_parameters_table_reports_values_and_errors() -> None:
    fit = FitResult(
        success=True,
        status=FitStatus.SUCCESS,
        quality=FitQuality.GOOD,
        parameters=np.array([1.0, 2.0]),
        errors=np.array([0.1, 0.2]),
        metadata={"parameter_order": ["a0", "a1"]},
    )
    estimate = PressureEstimate(
        method="polynomial",
        pressure=np.array([1.0, 0.0]),
        fit=fit,
        unit="GPa",
        metadata={"degree": 1},
    )
    preview = PressureVolumePreview(
        volume=np.array([18.0, 19.0]),
        energy=np.array([-10.0, -10.2]),
        pressure_unit="GPa",
        polynomial=estimate,
    )

    table = preview_parameters_table(preview)

    assert table.title == "Pressure-volume fit parameters"
    assert table.columns == [
        "Method",
        "Model",
        "Parameter",
        "Source",
        "Value",
        "Sigma",
        "Unit",
    ]
    assert table.rows[0] == [
        "polynomial",
        1,
        "a0",
        "fitted",
        1.0,
        0.1,
        "model unit",
    ]


def test_preview_parameters_table_reports_implied_bm2_parameters() -> None:
    volume = np.linspace(68.0, 76.0, 13)
    energy = -100.0 + 0.003 * (volume - 72.0) ** 2
    preview = pressure_volume_preview(
        QHAInput(jobname="bm2-preview", volume=volume, energy=energy),
        QHAOptions(
            eos="BM2",
            energy_unit="eV",
            volume_unit="A",
            pressure_unit="GPa",
        ),
        include_polynomial=False,
    )

    table = preview_parameters_table(preview)
    eos_rows = {row[2]: row for row in table.rows if row[0] == "eos"}

    assert list(eos_rows) == ["E0", "K0", "KP", "KPP", "V0"]
    assert eos_rows["KP"][3] == "implied"
    assert eos_rows["KP"][4] == pytest.approx(4.0)
    assert eos_rows["KP"][5] == pytest.approx(0.0)
    assert eos_rows["KPP"][3] == "implied"
    assert eos_rows["KPP"][5] is not None
    assert eos_rows["KPP"][5] > 0.0
    assert eos_rows["KPP"][6] == "GPa^-1"
    assert eos_rows["K0"][3] == "fitted"
    assert eos_rows["K0"][6] == "GPa"


def test_preview_parameters_table_reports_fitted_bm4_kpp() -> None:
    from quantas.core.physics.eos import EnergyEOS
    from quantas.core.physics.units import pressure_to_energy

    volume = np.linspace(68.0, 76.0, 17)
    k0 = float(pressure_to_energy(160.0, "eV", "A", "GPa"))
    kpp = -0.02 / float(pressure_to_energy(1.0, "eV", "A", "GPa"))
    parameters = np.array([-100.0, k0, 4.2, kpp, 72.0])
    energy = EnergyEOS().evaluate("BM4", volume, parameters)
    preview = pressure_volume_preview(
        QHAInput(jobname="bm4-preview", volume=volume, energy=energy),
        QHAOptions(
            eos="BM4",
            energy_unit="eV",
            volume_unit="A",
            pressure_unit="GPa",
        ),
        include_polynomial=False,
    )

    table = preview_parameters_table(preview)
    eos_rows = {row[2]: row for row in table.rows if row[0] == "eos"}

    assert eos_rows["KPP"][3] == "fitted"
    assert eos_rows["KPP"][4] == pytest.approx(-0.02, rel=1.0e-5)
    assert eos_rows["KPP"][5] is not None
    assert eos_rows["KPP"][6] == "GPa^-1"


def test_selected_property_tables_are_grouped_by_pressure() -> None:
    tables = selected_property_tables(_result())

    assert tables
    assert tables[0].title == "QHA results at P = 0.00"
    assert tables[0].columns[:4] == ["T", "V", "sigma_V", "KT"]
    assert tables[0].metadata["column_units"][:4] == [
        "K",
        "unknown",
        "unknown",
        "unknown",
    ]
    assert tables[0].metadata["uncertainty_columns"] == ["sigma_V"]
    assert tables[0].rows[0][0] == "300.00"


def test_selected_property_tables_include_all_available_eos_uncertainties() -> None:
    result = _result()
    result.bulk_modulus_derivative = np.full((2, 2), 4.0)
    result.uncertainties.update(
        {
            "sigma_KT": np.full((2, 2), 0.5),
            "sigma_Kp": np.full((2, 2), 0.02),
        }
    )

    table = selected_property_tables(result)[0]

    assert table.columns[:7] == [
        "T",
        "V",
        "sigma_V",
        "KT",
        "sigma_KT",
        "Kp",
        "sigma_Kp",
    ]
    assert table.metadata["uncertainty_columns"] == [
        "sigma_V",
        "sigma_KT",
        "sigma_Kp",
    ]


def test_selected_property_tables_do_not_add_empty_uncertainty_columns() -> None:
    result = _result()
    result.uncertainties.clear()

    table = selected_property_tables(result)[0]

    assert table.columns[:3] == ["T", "V", "KT"]
    assert table.metadata["uncertainty_columns"] == []


def test_options_table_includes_polynomial_derivative_settings() -> None:
    from quantas.modules.qha.report import options_table

    options = QHAOptions(
        polynomial_derivative_method="analytic",
        polynomial_grid_points=7,
        polynomial_grid_separation=0.1,
    )
    table = options_table(options)

    assert ["Polynomial derivative method", "analytic"] in table.rows
    assert ["Polynomial local-grid points", 7] in table.rows
    assert ["Polynomial local-grid separation (%)", 0.1] in table.rows
