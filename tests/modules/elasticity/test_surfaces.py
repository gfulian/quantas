from __future__ import annotations

import h5py
import numpy as np
import pytest

from quantas.core.physics.elasticity import ElasticTensor, sample_elasticity_surfaces
from quantas.models import SurfacePlotSpec
from quantas.modules.elasticity import (
    ElasticityOptions,
    ElasticitySurfaceOptions,
    build_elasticity_3d_plots,
    calculate_elasticity_surfaces,
    read_elasticity_hdf5,
    read_elasticity_input,
    run_elasticity,
    write_elasticity_hdf5,
)
from quantas.renderers.plots import MatplotlibOptions, render_plot_collection


def _isotropic_tensor() -> ElasticTensor:
    lam = 50.0
    mu = 30.0
    stiffness = np.array(
        [
            [lam + 2 * mu, lam, lam, 0, 0, 0],
            [lam, lam + 2 * mu, lam, 0, 0, 0],
            [lam, lam, lam + 2 * mu, 0, 0, 0],
            [0, 0, 0, mu, 0, 0],
            [0, 0, 0, 0, mu, 0],
            [0, 0, 0, 0, 0, mu],
        ],
        dtype=float,
    )
    return ElasticTensor(stiffness)


def test_isotropic_surfaces_and_transverse_branches() -> None:
    surfaces = sample_elasticity_surfaces(
        _isotropic_tensor(),
        ntheta=5,
        nphi=7,
    )
    np.testing.assert_allclose(surfaces.surfaces["young"].radius, 78.75)
    np.testing.assert_allclose(
        surfaces.surfaces["compressibility_positive"].radius,
        1000.0 / 210.0,
    )
    np.testing.assert_allclose(surfaces.surfaces["shear_minimum"].values, 30.0)
    np.testing.assert_allclose(surfaces.surfaces["shear_maximum"].values, 30.0)
    np.testing.assert_allclose(surfaces.surfaces["poisson_minimum"].values, 0.3125)
    np.testing.assert_allclose(surfaces.surfaces["poisson_maximum"].values, 0.3125)
    assert "compressibility_negative" not in surfaces.surfaces
    assert "poisson_negative" not in surfaces.surfaces


def test_3d_xy_section_matches_existing_2d_sampling() -> None:
    input_data = read_elasticity_input(
        "tests/modules/elasticity/data/hydroxylapatite.dat"
    )
    result_data = run_elasticity(
        input_data,
        ElasticityOptions(calculate_2d=True, ntheta=10),
    )
    surfaces = calculate_elasticity_surfaces(
        result_data,
        ElasticitySurfaceOptions(ntheta=3, nphi=9),
    )
    expected = result_data.results["elasticity"].properties_2d["xy"]
    row = 1
    np.testing.assert_allclose(
        surfaces.surfaces["young"].values[row],
        expected["young_modulus"][:-1],
        rtol=1.0e-10,
        atol=1.0e-10,
    )
    positive = surfaces.surfaces["compressibility_positive"].radius[row]
    negative_surface = surfaces.surfaces.get("compressibility_negative")
    negative = (
        np.zeros_like(positive)
        if negative_surface is None
        else negative_surface.radius[row]
    )
    np.testing.assert_allclose(
        np.column_stack((positive, negative)),
        expected["linear_compressibility"][:-1],
        rtol=1.0e-10,
        atol=1.0e-10,
    )
    np.testing.assert_allclose(
        np.column_stack(
            (
                surfaces.surfaces["shear_minimum"].values[row],
                surfaces.surfaces["shear_maximum"].values[row],
            )
        ),
        expected["shear_modulus"][:-1],
        rtol=1.0e-4,
        atol=3.0e-3,
    )
    poisson_negative = surfaces.surfaces.get("poisson_negative")
    negative_values = (
        np.zeros(9)
        if poisson_negative is None
        else np.nan_to_num(poisson_negative.values[row])
    )
    np.testing.assert_allclose(
        np.column_stack(
            (
                negative_values,
                np.nan_to_num(surfaces.surfaces["poisson_minimum"].values[row]),
                surfaces.surfaces["poisson_maximum"].values[row],
            )
        ),
        expected["poisson_ratio"][:-1],
        rtol=2.0e-4,
        atol=5.0e-5,
    )


def test_surface_meshes_are_transient_and_not_written_to_hdf5(tmp_path) -> None:
    result = run_elasticity(
        read_elasticity_input("tests/modules/elasticity/data/hydroxylapatite.dat")
    )
    calculate_elasticity_surfaces(
        result,
        ElasticitySurfaceOptions(ntheta=5, nphi=7, properties=("young",)),
    )
    filename = write_elasticity_hdf5(result, tmp_path / "elasticity.hdf5")
    with h5py.File(filename, "r") as h5:
        assert "surfaces" not in h5["results"]
        assert "surface" not in h5["results"]
        assert "properties_3d" not in h5["results"]


def test_matplotlib_renders_the_neutral_surface_spec(tmp_path) -> None:
    """The supported renderer consumes the frontend-neutral surface contract."""
    result = run_elasticity(
        read_elasticity_input("tests/modules/elasticity/data/hydroxylapatite.dat")
    )
    collection = build_elasticity_3d_plots(
        result,
        ElasticitySurfaceOptions(ntheta=5, nphi=7, properties=("young",)),
    )
    assert isinstance(collection.plots[0], SurfacePlotSpec)
    rendered = render_plot_collection(
        collection,
        MatplotlibOptions(output_dir=tmp_path, close=True),
    )
    assert rendered.artifacts[0].kind == "surface"
    assert rendered.artifacts[0].path == tmp_path / "3d_young.png"
    assert rendered.artifacts[0].path.is_file()


def test_requested_3d_surfaces_round_trip_through_hdf5(tmp_path) -> None:
    """Persisted 3D fields retain grid, physical values, and diagnostics."""
    options = ElasticitySurfaceOptions(
        ntheta=5,
        nphi=7,
        properties=("young", "shear"),
    )
    result = run_elasticity(
        read_elasticity_input("tests/modules/elasticity/data/hydroxylapatite.dat"),
        ElasticityOptions(calculate_3d=True, surface_options=options),
    )
    payload = result.results["elasticity"]
    assert payload.has_3d_data()
    filename = write_elasticity_hdf5(result, tmp_path / "elasticity.hdf5")
    restored = read_elasticity_hdf5(filename).results["elasticity"]
    assert restored.has_3d_data()
    assert restored.properties_3d is not None
    assert set(restored.properties_3d.surfaces) == {
        "young",
        "shear_minimum",
        "shear_maximum",
    }
    for key in restored.properties_3d.surfaces:
        expected = payload.properties_3d.surfaces[key]
        actual = restored.properties_3d.surfaces[key]
        np.testing.assert_allclose(actual.radius, expected.radius)
        np.testing.assert_allclose(actual.values, expected.values)
        np.testing.assert_allclose(actual.x, expected.x, atol=1.0e-12)
        np.testing.assert_allclose(actual.y, expected.y, atol=1.0e-12)
        np.testing.assert_allclose(actual.z, expected.z, atol=1.0e-12)


def test_3d_plot_geometries_preserve_physical_values() -> None:
    """Geometry changes only radial coordinates, never physical color data."""
    result = run_elasticity(
        read_elasticity_input("tests/modules/elasticity/data/hydroxylapatite.dat")
    )
    sampling = ElasticitySurfaceOptions(
        ntheta=5,
        nphi=7,
        properties=("young",),
    )
    collections = {
        geometry: build_elasticity_3d_plots(
            result,
            sampling,
            geometry=geometry,
        )
        for geometry in ("physical", "unit_sphere", "normalized")
    }
    physical_values = collections["physical"].plots[0].layers[0].values
    for collection in collections.values():
        np.testing.assert_allclose(
            collection.plots[0].layers[0].values,
            physical_values,
        )
    radii = {}
    for geometry, collection in collections.items():
        layer = collection.plots[0].layers[0]
        radii[geometry] = np.sqrt(layer.x**2 + layer.y**2 + layer.z**2)
    np.testing.assert_allclose(radii["normalized"], radii["physical"])
    np.testing.assert_allclose(radii["unit_sphere"], 1.0)
    assert np.nanmax(radii["physical"]) > 100.0


def test_surface_renderer_has_clean_axes_and_optional_title(tmp_path) -> None:
    """Elasticity surfaces hide the Matplotlib box and use X/Y/Z labels."""
    result = run_elasticity(
        read_elasticity_input("tests/modules/elasticity/data/hydroxylapatite.dat")
    )
    collection = build_elasticity_3d_plots(
        result,
        ElasticitySurfaceOptions(ntheta=5, nphi=7, properties=("young",)),
    )
    rendered = render_plot_collection(
        collection,
        MatplotlibOptions(
            output_dir=tmp_path,
            close=False,
            show_title=True,
        ),
    )
    figure = rendered.artifacts[0].figure
    axis = figure.axes[0]
    assert axis.axison is False
    assert figure._suptitle is not None
    assert figure._suptitle.get_text() == "Young's modulus"
    labels = {text.get_text() for text in axis.texts}
    assert {"X", "Y", "Z"}.issubset(labels)
    assert figure.axes[1].get_ylabel() == r"Young's modulus ($\mathrm{GPa}$)"


def test_surface_plot_collection_can_enable_mesh_rendering() -> None:
    """Mesh rendering options are propagated to the neutral surface style."""
    result = run_elasticity(
        read_elasticity_input("tests/modules/elasticity/data/hydroxylapatite.dat")
    )
    collection = build_elasticity_3d_plots(
        result,
        ElasticitySurfaceOptions(ntheta=5, nphi=7, properties=("young",)),
        show_mesh=True,
        mesh_color="black",
        mesh_line_width=0.5,
    )
    style = collection.plots[0].layers[0].style
    assert style.show_mesh is True
    assert style.mesh_color == "black"
    assert style.mesh_line_width == pytest.approx(0.5)


def test_build_3d_plots_reuses_persisted_surface_data(monkeypatch) -> None:
    """Plot preparation does not resample a complete stored property set."""
    result = run_elasticity(
        read_elasticity_input("tests/modules/elasticity/data/hydroxylapatite.dat"),
        ElasticityOptions(
            calculate_3d=True,
            surface_options=ElasticitySurfaceOptions(
                ntheta=5,
                nphi=7,
                properties=("young",),
            ),
        ),
    )

    def fail_calculation(*args, **kwargs):
        raise AssertionError("stored 3D data should have been reused")

    monkeypatch.setattr(
        "quantas.modules.elasticity.surface.calculate_elasticity_surfaces",
        fail_calculation,
    )
    collection = build_elasticity_3d_plots(
        result,
        properties=("young",),
    )
    assert [plot.key for plot in collection.plots] == ["young"]


def test_missing_stored_3d_property_is_calculated_transiently() -> None:
    """A missing requested property is sampled without mutating stored data."""
    result = run_elasticity(
        read_elasticity_input("tests/modules/elasticity/data/hydroxylapatite.dat"),
        ElasticityOptions(
            calculate_3d=True,
            surface_options=ElasticitySurfaceOptions(
                ntheta=5,
                nphi=7,
                properties=("young",),
            ),
        ),
    )
    payload = result.results["elasticity"]
    assert payload.properties_3d is not None
    assert set(payload.properties_3d.surfaces) == {"young"}
    collection = build_elasticity_3d_plots(
        result,
        properties=("shear",),
    )
    assert {plot.key for plot in collection.plots} == {
        "shear_minimum",
        "shear_maximum",
    }
    assert set(payload.properties_3d.surfaces) == {"young"}
