"""Frozen contracts for the organized Quantas 2 public API."""

from __future__ import annotations

import ast
from pathlib import Path

import numpy as np

from quantas import api
from quantas.api import elasticity, eos, registry
from quantas.api.registry import Capability


EXPECTED_PUBLIC_SYMBOLS = {
    "common": {
        "CallbackObserver",
        "CrystalStructure",
        "Event",
        "EventLevel",
        "EventRecord",
        "InputData",
        "ListObserver",
        "NullObserver",
        "Observer",
        "PhononInputData",
        "PhononInterface",
        "PlotCollection",
        "ReportTable",
        "ResultData",
        "ResultMetadata",
        "StructureVolumeSeries",
        "SymmetryMetadata",
        "TensorRotation",
        "TensorRotationKind",
        "get_result_payload",
    },
    "elasticity": {
        "Input",
        "InputInterface",
        "Options",
        "PlotProperty",
        "Result",
        "SurfaceGeometry",
        "SurfaceOptions",
        "SurfaceProperty",
        "TensorRotation",
        "TensorRotationKind",
        "build_2d_plots",
        "build_3d_plots",
        "build_plots",
        "build_report",
        "create_input",
        "describe_plots",
        "get_result",
        "normalize_input",
        "read_input",
        "read_result",
        "run",
        "write_result",
        "write_table",
    },
    "eos": {
        "Archive",
        "ArchiveInspection",
        "ArchivePlotInventory",
        "ArchiveSizeInfo",
        "BatchFailurePolicy",
        "BatchJob",
        "BatchPlan",
        "BatchResult",
        "CalculationResult",
        "CovarianceScaling",
        "CapabilityStatus",
        "DOMAIN_CAPABILITIES",
        "Dataset",
        "DatasetPlotDescriptor",
        "DiagnosticResult",
        "DomainCapability",
        "EOSFamily",
        "EOSModel",
        "EffectiveVarianceOptions",
        "FitDomain",
        "FitRecord",
        "FitMethod",
        "FitOptions",
        "FitRequest",
        "FitResult",
        "MGDNormalization",
        "MGDVolumeBasis",
        "ODRDifferenceScheme",
        "ODROptions",
        "OLSOptions",
        "ParameterConstraint",
        "PVTModel",
        "PVTCouplingFamily",
        "PlotOptions",
        "RecordDisposition",
        "RecordInspection",
        "RecordPlotDescriptor",
        "ReportDetail",
        "ReportOptions",
        "ResolvedSpec",
        "ResultSlot",
        "SlotInspection",
        "SlotPlotDescriptor",
        "SlotState",
        "SlotStatus",
        "StateEvent",
        "StateEventType",
        "TemperatureEOSFamily",
        "TemperatureEOSModel",
        "ThermalPressureFamily",
        "ThermalPressureModel",
        "SPEC_TEMPLATE_FILENAME",
        "Session",
        "SpecDocument",
        "SpecError",
        "SpecInputOptions",
        "SolverOptions",
        "available_plot_types",
        "available_eos_models",
        "available_eos_tags",
        "available_pvt_couplings",
        "available_temperature_eos_models",
        "build_batch_preamble",
        "build_batch_report",
        "build_plots",
        "calculate",
        "calculation_summary_table",
        "calculation_table",
        "describe_plots",
        "diagnose",
        "diagnostic_summary_table",
        "diagnostic_table",
        "domain_capability",
        "fit",
        "normalize_input",
        "open_archive",
        "read_input",
        "read_spec",
        "record_domain",
        "resolve_spec",
        "run_batch",
        "validate_request",
        "WLSOptions",
        "default_solver_options",
        "write_calculation_csv",
        "write_diagnostic_csv",
        "write_spec_template",
    },
    "ha": {
        "CurveAxis",
        "Input",
        "KiefferVolumeSeries",
        "Options",
        "PlotOptions",
        "PhononInterface",
        "Result",
        "StructureVolumeSeries",
        "build_plots",
        "build_report",
        "describe_plots",
        "create_input",
        "get_result",
        "normalize_input",
        "read_input",
        "read_result",
        "run",
        "write_result",
        "write_table",
    },
    "interop": {
        "load_qha_result",
        "qha_to_thermoelastic_context",
        "thermoelastic_to_seismic",
    },
    "profiles": {
        "DepthProfile",
        "Model",
        "Preset",
        "build_preset",
        "from_mapping",
        "preset_names",
        "presets",
        "read_spec",
    },
    "plotting": {
        "AxisFieldLayer",
        "AxisLocation",
        "AxisOrientation",
        "BandOrientation",
        "ColoredPathSeries",
        "ColoredPathStyle",
        "ContourPlotSpec",
        "LinePlotSpec",
        "LineStyle",
        "PanelPlotSpec",
        "PlotAxis",
        "PlotBand",
        "PlotBandStyle",
        "PlotCollection",
        "PlotRepresentationDescriptor",
        "PlotPropertyDescriptor",
        "PlotKind",
        "PlotInventory",
        "PlotContextValue",
        "PlotContextDescriptor",
        "PlotMask",
        "PlotSeries",
        "PlotSeriesStyle",
        "PlotSpan",
        "PlotSpec",
        "PolarPlotPanel",
        "PolarPlotSpec",
        "ScalarBackground",
        "SecondaryAxis",
        "SphericalMapSpec",
        "SphericalMarker",
        "SphericalProjection",
        "SphericalSummarySpec",
        "SurfaceLayer",
        "SurfacePlotSpec",
        "SurfaceStyle",
        "VectorFieldLayer",
        "VectorFieldStyle",
    },
    "qha": {
        "CurveAxis",
        "FitFailurePolicy",
        "Input",
        "KiefferVolumeSeries",
        "Minimization",
        "ModeContinuity",
        "Options",
        "PhononInterface",
        "PlotOptions",
        "PolynomialDerivativeMethod",
        "Preview",
        "PropertyDifference",
        "Result",
        "Scheme",
        "StructureVolumeSeries",
        "TableFileFormat",
        "TableFormat",
        "ThermalExpansionMethod",
        "ValidationSummary",
        "available_energy_eos",
        "build_inspection_plots",
        "build_inspection_report",
        "build_plots",
        "build_report",
        "compare_results",
        "create_input",
        "describe_plots",
        "get_result",
        "inspect",
        "list_plot_properties",
        "normalize_input",
        "read_input",
        "read_result",
        "run",
        "validate_result",
        "write_result",
        "write_table",
    },
    "rendering": {
        "PlotRenderResult",
        "RenderedPlot",
        "render_plots",
        "render_table",
        "render_tables",
    },
    "registry": {
        "Capability",
        "ModuleDescriptor",
        "OperationDescriptor",
        "get",
        "iter_modules",
        "list_modules",
        "module_from_result",
        "open_result",
    },
    "seismic": {
        "ElasticMedium",
        "Hemisphere",
        "Input",
        "InputInterface",
        "Options",
        "PlotOptions",
        "Result",
        "SamplingLevel",
        "SurfaceGeometry",
        "SurfaceOptions",
        "SurfaceType",
        "TensorRotation",
        "TensorRotationKind",
        "WaveMode",
        "build_plots",
        "build_report",
        "describe_plots",
        "build_summary",
        "build_surfaces",
        "create_input",
        "get_result",
        "normalize_input",
        "read_input",
        "read_result",
        "run",
        "write_csv",
        "write_result",
    },
    "thermoelasticity": {
        "AdiabaticMode",
        "ComparePlotOptions",
        "ComponentGroup",
        "Context",
        "CrystalStructure",
        "DepthProfile",
        "DomainPlotOptions",
        "ElasticVolumePoint",
        "ElasticVolumeSeries",
        "ExtrapolationPolicy",
        "FitFailurePolicy",
        "FitMethod",
        "FitPlotOptions",
        "Input",
        "InputInterface",
        "Options",
        "PTPlotOptions",
        "PTQuantity",
        "PlotLayout",
        "PlotPreset",
        "PlotStyleOptions",
        "ProfileBackground",
        "ProfileColor",
        "ProfileMode",
        "ProfilePlotOptions",
        "ProfilePreset",
        "ProfileResult",
        "QHAInput",
        "QHAOptions",
        "QHAThermoelasticPayload",
        "QualityPolicy",
        "ReportLevel",
        "Result",
        "StabilityPolicy",
        "SymmetryMetadata",
        "TensorCondition",
        "UncertaintyMode",
        "analyze_grid",
        "analyze_profiles",
        "build_compare_plots",
        "build_domain_plot",
        "build_fit_plots",
        "build_plots",
        "build_profile_plots",
        "build_profile_preset",
        "build_pt_plots",
        "build_report",
        "create_input",
        "describe_plots",
        "get_result",
        "grid_info_table",
        "normalize_input",
        "point_table",
        "prepare_context",
        "profile_presets",
        "read_depth_profile",
        "read_input",
        "read_profile_spec",
        "read_result",
        "regular_grid",
        "resolve_components",
        "run",
        "run_context",
        "select_stiffness_tensor",
        "write_grid_table",
        "write_profile_table",
        "write_profile_template",
        "write_result",
        "write_state_input",
        "write_tensor_export",
    },
}


def test_public_namespace_symbol_snapshot_is_explicit() -> None:
    """Every supported public name is reviewed rather than exported implicitly."""
    assert set(api.__all__) == set(EXPECTED_PUBLIC_SYMBOLS)
    assert {name for name in dir(api) if not name.startswith("_")} == set(
        EXPECTED_PUBLIC_SYMBOLS
    )
    for namespace, expected in EXPECTED_PUBLIC_SYMBOLS.items():
        module = getattr(api, namespace)
        assert set(module.__all__) == expected
        assert {name for name in dir(module) if not name.startswith("_")} == expected


def test_public_names_do_not_expose_implementation_roles() -> None:
    """Calculators, readers, exporters, workflows, and renderers stay internal."""
    forbidden = ("Calculator", "Reader", "Exporter", "Workflow", "Renderer", "Engine")
    for namespace, names in EXPECTED_PUBLIC_SYMBOLS.items():
        offenders = sorted(name for name in names if name.endswith(forbidden))
        assert not offenders, f"{namespace} exposes implementation roles: {offenders}"


def test_registry_operations_are_declared_public_symbols() -> None:
    """Every capability resolves to an intentionally exported public operation."""
    for descriptor in registry.list_modules():
        namespace = descriptor.load()
        for capability, operation_name in descriptor.operations:
            assert capability in descriptor.capabilities
            assert operation_name in namespace.__all__
            assert descriptor.operation(capability) is getattr(
                namespace, operation_name
            )


def test_registry_dispatches_native_result_and_eos_archive(tmp_path: Path) -> None:
    """HDF5 metadata selects the same public reader used by frontends."""
    source = Path(__file__).resolve().parents[2] / "examples/elasticity/calcite.dat"
    result = elasticity.run(source, options=elasticity.Options(calculate_2d=False))
    result_path = elasticity.write_result(result, tmp_path / "calcite")

    descriptor = registry.module_from_result(result_path)
    assert descriptor.name == "elasticity"
    reopened = registry.open_result(result_path)
    np.testing.assert_allclose(
        elasticity.get_result(reopened).stiffness,
        elasticity.get_result(result).stiffness,
    )

    dataset = eos.read_input(
        Path(__file__).resolve().parents[2] / "examples/eos/PV_quartz.dat"
    )
    archive_path = tmp_path / "quartz_eos.hdf5"
    with eos.Archive.create(archive_path, dataset=dataset, creator="API test"):
        pass
    assert registry.module_from_result(archive_path).name == "eos"
    opened = registry.open_result(archive_path)
    try:
        assert isinstance(opened, eos.Archive)
        assert opened.dataset_ids == (1,)
    finally:
        opened.close()


def test_cli_adapters_use_only_public_api() -> None:
    """CLI numerical dispatch depends on ``quantas.api`` rather than old facades."""
    source_root = Path(__file__).resolve().parents[2] / "src/quantas/cli"
    forbidden = {
        "quantas.modules.elasticity",
        "quantas.modules.seismic",
        "quantas.modules.ha",
        "quantas.modules.qha",
        "quantas.modules.eos",
        "quantas.modules.thermoelasticity",
    }
    offenders: list[str] = []
    for path in source_root.glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and (
                node.module in forbidden
                or any(node.module == f"{name}.api" for name in forbidden)
            ):
                offenders.append(f"{path.name}:{node.lineno}:{node.module}")
    assert not offenders, offenders


def test_cli_input_and_export_commands_use_public_facades() -> None:
    """Command adapters do not instantiate private input or export helpers."""
    source_root = Path(__file__).resolve().parents[2] / "src/quantas/cli"
    command_files = (
        "elasticity.py",
        "ha.py",
        "qha.py",
        "seismic.py",
        "phonon_input.py",
    )
    forbidden_fragments = (
        ".io.inpgen",
        ".io.export",
        "ElasticityInputCreator",
        "ElasticityTableExport",
        "HATableExport",
        "QHATableExport",
    )
    offenders: list[str] = []
    for filename in command_files:
        path = source_root / filename
        text = path.read_text(encoding="utf-8")
        for fragment in forbidden_fragments:
            if fragment in text:
                offenders.append(f"{filename}:{fragment}")
    assert offenders == []


def test_eos_registry_uses_domain_specific_capabilities() -> None:
    """EOS is discoverable without pretending to be a single-shot workflow."""
    descriptor = registry.get("eos")
    assert descriptor.has(Capability.FIT)
    assert descriptor.has(Capability.BATCH)
    assert descriptor.has(Capability.ARCHIVE)
    assert not descriptor.has(Capability.RUN)
