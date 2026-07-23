"""Tests for neutral QHA plot builders and the Matplotlib renderer."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg", force=True)

import numpy as np

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
