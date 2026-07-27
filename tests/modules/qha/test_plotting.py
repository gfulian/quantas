"""Tests for neutral QHA plot builders and the Matplotlib renderer."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg", force=True)

import numpy as np
import pytest

from quantas.models import PlotCollection, ResultData, ResultMetadata
from quantas.modules.qha.models import QHAResult
from quantas.modules.qha.plot import (
    QHAPlotOptions,
    build_qha_plot_collection,
    list_available_plot_properties,
)
from quantas.renderers.plots import MatplotlibOptions, render_plot_collection


def _result_data() -> ResultData:
    temperature = np.array([0.0, 100.0, 200.0])
    pressure = np.array([0.0, 5.0])
    volume = np.array([[18.0, 17.7], [18.1, 17.8], [18.3, 18.0]])
    qha = QHAResult(
        jobname="MgO",
        temperature=temperature,
        pressure=pressure,
        equilibrium_volume=volume,
        isothermal_bulk_modulus=np.array(
            [[160.0, 180.0], [158.0, 178.0], [154.0, 174.0]]
        ),
        adiabatic_bulk_modulus=np.array(
            [[160.0, 180.0], [159.0, 179.0], [155.0, 175.0]]
        ),
        isochoric_heat_capacity=np.array([[0.0, 0.0], [10.0, 10.2], [12.0, 12.3]]),
        isobaric_heat_capacity=np.array([[0.0, 0.0], [10.1, 10.3], [12.2, 12.5]]),
        thermal_expansion=np.array([[0.0, 0.0], [1.1e-5, 1.0e-5], [1.5e-5, 1.3e-5]]),
        metadata={
            "units": {
                "temperature": "K",
                "pressure": "GPa",
                "volume": "A^3",
                "heat_capacity": "J mol^-1 K^-1",
            }
        },
    )
    return ResultData(
        metadata=ResultMetadata(module="qha", method="qha"), results={"qha": qha}
    )


def _render(
    data: ResultData,
    property_names: list[str],
    tmp_path,
    plot_options: QHAPlotOptions | None = None,
):
    collection = build_qha_plot_collection(
        data,
        property_names=property_names,
        options=plot_options,
    )
    return render_plot_collection(
        collection,
        MatplotlibOptions(output_dir=tmp_path, close=True),
    )


def test_build_qha_collection_is_frontend_neutral() -> None:
    collection = build_qha_plot_collection(_result_data(), property_names=["VT"])

    assert isinstance(collection, PlotCollection)
    assert collection.plots[0].key == "VT"
    assert collection.plots[0].filename_stem == "VT_1D"


def test_qha_renderer_writes_curve_and_contour_files(tmp_path) -> None:
    result = _render(
        _result_data(),
        ["VT"],
        tmp_path,
        QHAPlotOptions(include_contours=True),
    )

    assert not result.warnings
    paths = {
        artifact.path.name for artifact in result.artifacts if artifact.path is not None
    }
    assert "VT_1D.png" in paths
    assert "VT_2D.png" in paths


def test_qha_collection_warns_when_contour_grid_is_too_small(tmp_path) -> None:
    data = _result_data()
    qha = data.results["qha"]
    qha.pressure = np.array([0.0])
    qha.equilibrium_volume = qha.equilibrium_volume[:, :1]
    collection = build_qha_plot_collection(
        data,
        property_names=["VT"],
        options=QHAPlotOptions(include_contours=True),
    )
    result = render_plot_collection(
        collection,
        MatplotlibOptions(output_dir=tmp_path, close=True),
    )

    assert any(
        "requires at least two temperatures" in warning for warning in result.warnings
    )
    assert (tmp_path / "VT_1D.png").exists()
    assert not (tmp_path / "VT_2D.png").exists()


def test_heat_capacity_plot_combines_cp_and_cv(tmp_path) -> None:
    result = _render(_result_data(), ["heat_capacities"], tmp_path)

    assert result.artifacts[0].property_name == "heat_capacities"
    assert (tmp_path / "heat_capacities_1D.png").exists()


def test_list_available_plot_properties_includes_heat_capacities() -> None:
    rows = list_available_plot_properties(_result_data())

    keys = {row[0] for row in rows}
    assert "VT" in keys
    assert "KT" in keys
    assert "heat_capacities" in keys


def test_contour_plot_accepts_colormap_and_smooth_options(tmp_path) -> None:
    result = _render(
        _result_data(),
        ["alphaV"],
        tmp_path,
        QHAPlotOptions(
            include_contours=True,
            cmap="plasma",
            contour_mode="smooth",
            levels=8,
            isolines=True,
            isoline_labels=True,
        ),
    )

    assert not result.warnings
    assert (tmp_path / "alphaV_2D.png").exists()


def test_plot_unit_conversion_is_applied_to_axes_and_values(tmp_path) -> None:
    result = _render(
        _result_data(),
        ["Cv"],
        tmp_path,
        QHAPlotOptions(
            temperature_unit="Celsius",
            pressure_unit="bar",
            energy_unit="kJ/mol",
        ),
    )
    ax = result.artifacts[0].figure.axes[0]

    assert "Celsius" in ax.get_xlabel()
    assert "kJ/mol" in ax.get_ylabel()
    np.testing.assert_allclose(ax.lines[0].get_xdata(), [-273.15, -173.15, -73.15])
    np.testing.assert_allclose(ax.lines[0].get_ydata(), [0.0, 0.0100, 0.0120])
    assert "P = 50000 bar" in ax.lines[1].get_label()
    assert (tmp_path / "Cv_1D.png").exists()


def test_heat_capacity_conversion_accounts_for_fahrenheit_intervals(tmp_path) -> None:
    result = _render(
        _result_data(),
        ["Cv"],
        tmp_path,
        QHAPlotOptions(
            temperature_unit="Fahrenheit",
            energy_unit="J/mol",
        ),
    )
    ax = result.artifacts[0].figure.axes[0]

    np.testing.assert_allclose(ax.lines[0].get_xdata(), [-459.67, -279.67, -99.67])
    np.testing.assert_allclose(
        ax.lines[0].get_ydata(),
        np.array([0.0, 10.0, 12.0]) * (5.0 / 9.0),
    )
    assert "Fahrenheit" in ax.get_ylabel()


def test_heat_capacity_plot_draws_dulong_petit_when_natoms_are_available(
    tmp_path,
) -> None:
    data = _result_data()
    data.results["qha"].metadata["natoms"] = 2
    result = _render(data, ["heat_capacities"], tmp_path)
    labels = [line.get_label() for line in result.artifacts[0].figure.axes[0].lines]

    assert "Dulong-Petit limit" in labels


def test_dulong_petit_and_native_cv_use_formula_unit_normalization(tmp_path) -> None:
    from scipy import constants as cs

    natoms = 4
    formula_units = 2
    native_limit = (
        3.0 * natoms * cs.R / (cs.Avogadro * cs.physical_constants["Hartree energy"][0])
    )
    qha = QHAResult(
        temperature=np.array([300.0, 600.0]),
        pressure=np.array([0.0]),
        isochoric_heat_capacity=np.full((2, 1), native_limit),
        metadata={
            "units": {
                "temperature": "K",
                "pressure": "GPa",
                "energy": "Ha",
                "heat_capacity": "Ha cell^-1 K^-1",
            },
            "normalization": {
                "formula_units_per_cell": formula_units,
                "natoms_per_cell": natoms,
                "natoms_per_formula_unit": natoms / formula_units,
            },
            "input": {"natoms": natoms, "formula_units": formula_units},
        },
    )
    data = ResultData(
        metadata=ResultMetadata(module="qha", method="qha"),
        results={"qha": qha},
    )
    result = _render(
        data,
        ["Cv"],
        tmp_path,
        QHAPlotOptions(energy_unit="J/mol"),
    )
    ax = result.artifacts[0].figure.axes[0]
    curve = ax.lines[0].get_ydata()
    limit_line = next(
        line for line in ax.lines if line.get_label() == "Dulong-Petit limit"
    )
    expected = 3.0 * (natoms / formula_units) * cs.R

    np.testing.assert_allclose(curve, expected)
    np.testing.assert_allclose(limit_line.get_ydata(), expected)


def test_pressure_axis_sections_select_exact_native_temperatures() -> None:
    collection = build_qha_plot_collection(
        _result_data(),
        property_names=["VT"],
        options=QHAPlotOptions(
            curve_axis="pressure",
            selected_temperatures=(100.0,),
        ),
    )
    spec = collection.plots[0]

    assert spec.x_axis.key == "pressure"
    assert spec.filename_stem == "VT_P"
    assert len(spec.series) == 1
    np.testing.assert_allclose(spec.series[0].x, [0.0, 5.0])
    np.testing.assert_allclose(spec.series[0].y, [18.1, 17.8])
    assert spec.series[0].metadata["temperature_native"] == 100.0
    assert spec.metadata["curve_axis"] == "pressure"


def test_temperature_axis_sections_can_select_native_pressures() -> None:
    collection = build_qha_plot_collection(
        _result_data(),
        property_names=["KT"],
        options=QHAPlotOptions(
            curve_axis="temperature",
            selected_pressures=(5.0,),
        ),
    )
    spec = collection.plots[0]

    assert spec.x_axis.key == "temperature"
    assert len(spec.series) == 1
    np.testing.assert_allclose(spec.series[0].y, [180.0, 178.0, 174.0])
    assert spec.series[0].metadata["pressure_native"] == 5.0


def test_qha_sections_reject_coordinates_absent_from_native_grid() -> None:
    with pytest.raises(ValueError, match="not present in the native grid"):
        build_qha_plot_collection(
            _result_data(),
            property_names=["VT"],
            options=QHAPlotOptions(
                curve_axis="pressure",
                selected_temperatures=(150.0,),
            ),
        )


def test_pressure_axis_heat_capacity_comparison_uses_selected_temperatures() -> None:
    collection = build_qha_plot_collection(
        _result_data(),
        property_names=["heat_capacities"],
        options=QHAPlotOptions(
            curve_axis="pressure",
            selected_temperatures=(200.0,),
        ),
    )
    spec = collection.plots[0]

    assert spec.x_axis.key == "pressure"
    assert spec.filename_stem == "heat_capacities_P"
    assert len(spec.series) == 2
    np.testing.assert_allclose(spec.series[0].y, [12.2, 12.5])
    np.testing.assert_allclose(spec.series[1].y, [12.0, 12.3])


def test_qha_slice_builders_do_not_mutate_result_arrays() -> None:
    data = _result_data()
    payload = data.results["qha"]
    before = np.array(payload.equilibrium_volume, copy=True)

    build_qha_plot_collection(
        data,
        property_names=["VT"],
        options=QHAPlotOptions(
            curve_axis="pressure",
            selected_temperatures=(0.0, 200.0),
            include_contours=True,
        ),
    )

    np.testing.assert_array_equal(payload.equilibrium_volume, before)


def test_pressure_sections_require_multiple_stored_pressures() -> None:
    data = _result_data()
    qha = data.results["qha"]
    qha.pressure = np.asarray(qha.pressure)[:1]
    for attribute in (
        "equilibrium_volume",
        "isothermal_bulk_modulus",
        "isochoric_heat_capacity",
        "isobaric_heat_capacity",
        "thermal_expansion",
    ):
        values = getattr(qha, attribute, None)
        if values is not None and np.asarray(values).ndim == 2:
            setattr(qha, attribute, np.asarray(values)[:, :1])

    with pytest.raises(ValueError, match="at least two stored pressures"):
        build_qha_plot_collection(
            data,
            property_names=["VT"],
            options=QHAPlotOptions(curve_axis="pressure"),
        )


def test_qha_sections_reject_duplicate_native_pressures() -> None:
    data = _result_data()
    data.results["qha"].pressure = np.array([0.0, 0.0], dtype=np.float64)

    with pytest.raises(ValueError, match="duplicate coordinates"):
        build_qha_plot_collection(data, property_names=["VT"])
