"""Tests for QSA component fitting, HDF5 diagnostics, and depth profiles."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from quantas.core.physics.elasticity import (
    cold_finite_strain_component,
    reconstruct_stiffness_from_components,
)
from quantas.core.physics.eos import EnergyEOS
from quantas.core.physics.units import energy_to_pressure
from quantas.models import ResultData, ResultMetadata
from quantas.models.structures import CrystalStructure, SymmetryMetadata
from quantas.modules.qha.models import QHAResult
from quantas.modules.thermoelasticity import (
    ElasticVolumePoint,
    ElasticVolumeSeries,
    MODULE_CONTRACT,
    ThermoelasticAnalysisEngine,
    ThermoelasticContext,
    ThermoelasticDepthProfile,
    ThermoelasticInput,
    ThermoelasticOptions,
    analyze_thermoelastic_result,
    build_thermoelastic_analysis_report,
    build_thermoelastic_profile_preset,
    build_thermoelastic_report,
    thermoelastic_profile_presets,
    read_thermoelastic_hdf5,
    run_thermoelastic_context,
    write_thermoelastic_hdf5,
)


def _synthetic_context() -> ThermoelasticContext:
    volumes = np.linspace(90.0, 110.0, 9)
    v0 = 100.0
    k0_gpa = 120.0
    kp = 4.4
    pressure_factor = float(energy_to_pressure(1.0, "Ha", "A", "GPa"))
    energy_parameters = [-200.0, k0_gpa / pressure_factor, kp, v0]
    energies = EnergyEOS().evaluate("BM3", volumes, energy_parameters)
    component_parameters = {
        "C11": (220.0, 6.0),
        "C12": (85.0, 3.5),
        "C44": (72.0, 1.8),
    }
    points: list[ElasticVolumePoint] = []
    mass = 4.0e-25
    for index, (volume, energy) in enumerate(zip(volumes, energies, strict=True)):
        components = {
            label: float(
                cold_finite_strain_component(
                    volume,
                    reference_volume=v0,
                    bulk_modulus=k0_gpa,
                    bulk_modulus_derivative=kp,
                    reference_component=c0,
                    component_pressure_derivative=cp,
                    wallace_delta={"C11": -3.0, "C12": -1.0, "C44": -1.0}[label],
                    order=3,
                )
            )
            for label, (c0, cp) in component_parameters.items()
        }
        matrix = reconstruct_stiffness_from_components(components, "cubic")
        a = volume ** (1.0 / 3.0)
        points.append(
            ElasticVolumePoint(
                source=f"point-{index}.out",
                pressure=float(index - 4),
                stress_pressure=float(index - 4),
                volume=float(volume),
                density=mass / (volume * 1.0e-30),
                energy=float(energy),
                stiffness=matrix,
                lattice=np.diag([a, a, a]),
            )
        )
    reference = CrystalStructure(
        lattice=points[4].lattice,
        fractional_positions=np.asarray([[0.0, 0.0, 0.0]]),
        atomic_numbers=np.asarray([12]),
    )
    series = ElasticVolumeSeries(
        points=tuple(points),
        reference_structure=reference,
        symmetry=SymmetryMetadata(
            space_group_number=221,
            international_symbol="Pm-3m",
            point_group="m-3m",
        ),
        elastic_symmetry="cubic",
        reference_index=4,
    )
    temperature = np.asarray([300.0, 500.0, 700.0])
    pressure = np.asarray([0.0, 5.0])
    equilibrium_volume = np.asarray(
        [[100.5, 97.0], [102.0, 98.5], [103.5, 100.0]], dtype=np.float64
    )
    alpha_tensor = np.zeros(equilibrium_volume.shape + (3, 3), dtype=np.float64)
    for index in range(3):
        alpha_tensor[..., index, index] = 1.0e-5
    qha = QHAResult(
        jobname="synthetic",
        temperature=temperature,
        pressure=pressure,
        volume=volumes,
        static_energy=energies,
        equilibrium_volume=equilibrium_volume,
        isochoric_heat_capacity=np.full_like(equilibrium_volume, 1.0e-4),
        thermal_expansion_tensor=alpha_tensor,
        uncertainties={"sigma_VT": np.full_like(equilibrium_volume, 0.02)},
    )
    qha_result = ResultData(
        metadata=ResultMetadata(module="qha", method="quasi-harmonic"),
        options={
            "temperature_unit": "K",
            "pressure_unit": "GPa",
            "volume_unit": "A",
            "energy_unit": "Ha",
        },
        results={"qha": qha},
    )
    return ThermoelasticContext(
        input_data=ThermoelasticInput(
            jobname="synthetic QSA",
            elastic_series=series,
        ),
        qha_result_data=qha_result,
        qha=qha,
        extrapolation_mask=np.zeros_like(equilibrium_volume, dtype=np.bool_),
    )


def test_qsa_recovers_independent_components_and_reconstructs_grid() -> None:
    """Exact synthetic data recover C0 and Cprime and full cubic tensors."""
    context = _synthetic_context()
    calibration = run_thermoelastic_context(context)
    calibrated = calibration.results["thermoelasticity"]
    assert calibrated.completed
    assert calibrated.stiffness_isothermal is None
    result_data = analyze_thermoelastic_result(calibration)
    result = result_data.results["thermoelasticity"]
    assert result.independent_labels == ("C11", "C12", "C44")
    expected = {"C11": (220.0, 6.0), "C12": (85.0, 3.5), "C44": (72.0, 1.8)}
    for label, parameters in expected.items():
        np.testing.assert_allclose(
            result.component_fits[label].parameters,
            parameters,
            rtol=2.0e-8,
            atol=2.0e-8,
        )
        assert result.component_fits[label].fit is not None
        assert result.component_fits[label].fit.success
    assert result.stiffness_isothermal is not None
    assert result.stiffness_isothermal.shape == (3, 2, 6, 6)
    np.testing.assert_allclose(
        result.stiffness_isothermal[..., 0, 0],
        result.stiffness_isothermal[..., 1, 1],
    )
    np.testing.assert_allclose(
        result.stiffness_isothermal[..., 3, 3],
        result.stiffness_isothermal[..., 5, 5],
    )
    assert result.sigma_stiffness_isothermal is not None
    assert np.all(np.isfinite(result.sigma_stiffness_isothermal))
    assert result.stiffness_adiabatic is not None
    assert result.adiabatic_correction is not None
    assert result.adiabatic_valid_mask is not None
    assert np.all(result.adiabatic_valid_mask)
    assert np.all(result.stiffness_adiabatic >= result.stiffness_isothermal - 1.0e-12)
    assert result.independent_stiffness_covariance is not None
    assert result.independent_stiffness_covariance.shape == (3, 2, 3, 3)
    np.testing.assert_allclose(
        result.independent_stiffness_covariance,
        np.swapaxes(result.independent_stiffness_covariance, -1, -2),
    )
    assert result.stability is not None
    assert np.all(result.stability.stable_mask)
    for record in result.component_fits.values():
        assert record.quality is not None
        assert record.quality.level == "supported"


def test_depth_profile_and_hdf5_round_trip_preserve_all_diagnostics(
    tmp_path: Path,
) -> None:
    """Grid, profile tensors, covariance, correlation and traces survive HDF5."""
    context = _synthetic_context()
    profile = ThermoelasticDepthProfile.linear(
        name="lithosphere",
        depth_min=0.0,
        depth_max=100.0,
        npoints=6,
        pressure_gradient=0.03,
        temperature_at_depth_min=300.0,
        temperature_gradient=2.0,
    )
    calibration = run_thermoelastic_context(
        context,
        options=ThermoelasticOptions(solver_debug=True),
    )
    result_data = analyze_thermoelastic_result(calibration, profiles=[profile])
    path = tmp_path / "thermoelastic.hdf5"
    write_thermoelastic_hdf5(result_data, path)
    loaded = read_thermoelastic_hdf5(path)
    result = loaded.results["thermoelasticity"]
    profile_result = result.profiles["lithosphere"]
    assert profile_result.stiffness_isothermal.shape == (6, 6, 6)
    assert profile_result.independent_stiffness.shape == (6, 3)
    assert profile_result.independent_stiffness_covariance.shape == (6, 3, 3)
    assert result.independent_stiffness_covariance is not None
    fit = result.component_fits["C11"].fit
    assert fit is not None
    assert fit.covariance is not None
    assert fit.diagnostics is not None
    assert fit.diagnostics.correlation is not None
    assert (
        "reference_eos_parameter_sensitivity" in result.component_fits["C11"].metadata
    )
    assert result.reference_eos.fit.residuals is not None
    assert result.stability is not None
    assert profile_result.stability is not None
    assert result.stiffness_adiabatic is not None
    assert profile_result.stiffness_adiabatic is not None
    assert result.isochoric_heat_capacity_cell is not None
    assert result.thermal_expansion_tensor is not None
    assert result.component_fits["C11"].quality is not None
    assert result.component_fits["C11"].quality.level == "supported"


def test_report_levels_add_point_covariance_and_debug_tables() -> None:
    """Report verbosity changes rendering only, not stored diagnostics."""
    result_data = run_thermoelastic_context(
        _synthetic_context(),
        options=ThermoelasticOptions(solver_debug=True),
    )
    result = result_data.results["thermoelasticity"]
    standard = build_thermoelastic_report(result, level="standard")
    extended = build_thermoelastic_report(result, level="extended")
    debug = build_thermoelastic_report(result, level="debug")
    assert len(standard) == 3
    assert len(extended) > len(standard)
    assert len(debug) > len(extended)
    assert any("observed and fitted" in table.title for table in extended)
    assert any("covariance" in table.title for table in extended)
    assert any("solver diagnostics" in table.title for table in debug)


def test_debug_report_level_automatically_records_solver_traces() -> None:
    """Debug reporting requests detailed solver traces during calculation."""
    options = ThermoelasticOptions(report_level="debug")
    assert options.solver_debug
    result_data = run_thermoelastic_context(_synthetic_context(), options=options)
    fit = result_data.results["thermoelasticity"].component_fits["C11"].fit
    assert fit is not None
    assert fit.diagnostics is not None
    assert fit.diagnostics.metadata["detailed_trace_requested"]
    assert fit.diagnostics.metadata["evaluation_trace"]


def test_failed_component_fit_stops_with_diagnostics(
    tmp_path: Path,
) -> None:
    """Toxic input remains inspectable in reports and HDF5 without tensors."""
    context = _synthetic_context()
    original = context.input_data.elastic_series
    point = original.points[original.reference_index]
    short_series = ElasticVolumeSeries(
        points=(point,),
        reference_structure=original.reference_structure,
        symmetry=original.symmetry,
        elastic_symmetry=original.elastic_symmetry,
        reference_index=0,
    )
    toxic_context = ThermoelasticContext(
        input_data=ThermoelasticInput(
            jobname="toxic single-point QSA",
            elastic_series=short_series,
        ),
        qha_result_data=context.qha_result_data,
        qha=context.qha,
        extrapolation_mask=context.extrapolation_mask,
    )
    result_data = run_thermoelastic_context(
        toxic_context,
        options=ThermoelasticOptions(
            report_level="debug",
            fit_failure_policy="stop",
        ),
    )
    payload = result_data.results["thermoelasticity"]
    assert not payload.completed
    assert payload.stiffness_isothermal is None
    fit = payload.component_fits["C11"].fit
    assert fit is not None
    assert not fit.success
    assert "at least two points" in fit.message
    assert fit.diagnostics is not None
    assert fit.diagnostics.metadata["termination_category"] == "invalid_input"
    assert fit.diagnostics.metadata["exception_type"] == "_InvalidFitProblem"

    debug = build_thermoelastic_report(payload, level="debug")
    assert any("solver diagnostics" in table.title for table in debug)

    path = tmp_path / "toxic.hdf5"
    write_thermoelastic_hdf5(result_data, path)
    loaded = read_thermoelastic_hdf5(path).results["thermoelasticity"]
    loaded_fit = loaded.component_fits["C11"].fit
    assert loaded_fit is not None
    assert not loaded_fit.success
    assert loaded_fit.diagnostics is not None
    assert loaded_fit.diagnostics.metadata["termination_category"] == "invalid_input"


def test_report_levels_cover_reconstruction_and_extrapolation() -> None:
    """Post-fit reports avoid repeating fits and debug extrapolated states."""
    fit_result = run_thermoelastic_context(_synthetic_context())
    analyzed = analyze_thermoelastic_result(
        fit_result,
        pressure=[8.0],
        temperature=[900.0],
        extrapolation_policy="allow",
    )
    payload = analyzed.results["thermoelasticity"]
    standard = build_thermoelastic_analysis_report(payload, level="standard")
    extended = build_thermoelastic_analysis_report(payload, level="extended")
    debug = build_thermoelastic_analysis_report(
        payload, level="debug", extrapolation_policy="allow"
    )
    standard_titles = {table.title for table in standard}
    assert "Static reference EOS" not in standard_titles
    assert "Quasi-static elastic-component fits" not in standard_titles
    assert "Pressure-temperature analysis summary" in standard_titles
    assert any("elastic-component ranges" in table.title for table in extended)
    assert any("Extrapolation debug context" == table.title for table in debug)
    states = next(
        table for table in debug if table.title == "Grid extrapolation states"
    )
    distances = next(
        table
        for table in debug
        if table.title == "Grid extrapolation distances and uncertainty"
    )
    assert states.rows[0][5] is True
    assert distances.rows[0][1] > 0.0
    assert distances.rows[0][2] > 0.0


def test_built_in_profile_presets_are_scientific_only() -> None:
    """Public presets expose only provenance-rich Earth references."""
    presets = thermoelastic_profile_presets()
    names = tuple(preset.name for preset in presets)
    assert names == (
        "continental-cratonic",
        "continental-reference",
        "continental-active",
        "oceanic-10ma",
        "oceanic-50ma",
        "oceanic-100ma",
        "mantle-katsura-2022",
    )
    assert not any("thin-crust" == name for name in names)

    scientific = build_thermoelastic_profile_preset("mantle-katsura-2022")
    assert scientific.depth[0] == 50.0
    assert scientific.depth[-1] == 2800.0
    assert scientific.temperature[0] == 1646.0
    assert scientific.temperature[-1] == 2587.0
    assert scientific.metadata["kind"] == "composed_earth_profile"
    assert scientific.metadata["pressure"]["model"] == "prem"
    assert scientific.metadata["temperature"]["model"] == "katsura-2022"
    assert len(scientific.metadata["references"]) == 3


def test_public_api_builds_adiabatic_seismic_input(tmp_path: Path) -> None:
    """The top-level facade selects C^S and QHA-consistent density."""
    from quantas.api.interop import thermoelastic_to_seismic

    result_data = run_thermoelastic_context(_synthetic_context())
    archive = tmp_path / "thermoelastic_for_seismic.hdf5"
    write_thermoelastic_hdf5(result_data, archive)
    seismic_input = thermoelastic_to_seismic(
        archive,
        pressure=0.0,
        temperature=300.0,
        tensor_condition="adiabatic",
    )
    assert seismic_input.stiffness is not None
    assert seismic_input.stiffness.shape == (6, 6)
    assert seismic_input.density > 0.0
    evaluated = analyze_thermoelastic_result(
        result_data, pressure=[0.0], temperature=[300.0]
    ).results["thermoelasticity"]
    assert evaluated.stiffness_adiabatic is not None
    np.testing.assert_allclose(
        seismic_input.stiffness,
        evaluated.stiffness_adiabatic[0, 0],
    )


def test_calibration_records_static_eos_and_single_wallace_convention() -> None:
    """Calibration provenance identifies static EOS and no duplicate correction."""
    result_data = run_thermoelastic_context(_synthetic_context())
    result = result_data.results["thermoelasticity"]
    assert result.metadata["reference_eos_state"] == "static 0 K, P=0 reference"
    assert (
        "no thermal or zero-point contribution"
        in result.metadata["reference_eos_source"]
    )
    assert "no second pressure correction" in result.metadata["wallace_convention"]
    for record in result.component_fits.values():
        np.testing.assert_allclose(
            record.observed,
            _synthetic_context().input_data.elastic_series.stiffness[
                :, record.entries[0][0], record.entries[0][1]
            ],
            atol=1.0e-12,
        )


def test_calibration_and_analysis_share_reconstruction_engine(tmp_path: Path) -> None:
    """Calibration is frontend-neutral and point/grid evaluation is identical."""
    context = _synthetic_context()
    calibration = run_thermoelastic_context(context)
    payload = calibration.results["thermoelasticity"]
    assert payload.stiffness_isothermal is None
    assert payload.metadata["analysis_stage"] == "calibration"

    engine = ThermoelasticAnalysisEngine(calibration)
    point = engine.evaluate_point(pressure=5.0, temperature=500.0)
    grid = engine.evaluate_grid(pressure=[5.0], temperature=[500.0])
    point_payload = point.results["thermoelasticity"]
    grid_payload = grid.results["thermoelasticity"]
    assert point_payload.stiffness_isothermal == pytest.approx(
        grid_payload.stiffness_isothermal
    )
    assert point_payload.stiffness_adiabatic == pytest.approx(
        grid_payload.stiffness_adiabatic
    )
    assert point_payload.equilibrium_volume == pytest.approx(
        grid_payload.equilibrium_volume
    )

    archive = write_thermoelastic_hdf5(calibration, tmp_path / "calibration")
    with h5py.File(archive, "r") as handle:
        assert handle["results"].attrs["thermoelastic_schema_version"] == "1.0"


def test_module_contract_is_public() -> None:
    """Thermoelasticity exposes the uniform workflow module contract."""
    assert MODULE_CONTRACT.name == "thermoelasticity"
    assert MODULE_CONTRACT.result_key == "thermoelasticity"
    assert callable(MODULE_CONTRACT.build_plots)
    assert callable(MODULE_CONTRACT.build_report)
