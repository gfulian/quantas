# -*- coding: utf-8 -*-

"""Public API for quasi-static thermoelastic calculations and analysis."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from quantas.core.events import Observer
from quantas.modules.qha.models import QHAInput, QHAOptions
from quantas.modules.thermoelasticity.api import (
    build_thermoelastic_plots as _build_plots,
    create_thermoelastic_input as _create_input,
    normalize_thermoelastic_input as _normalize_input,
    read_thermoelastic_hdf5 as _read_result,
    read_thermoelastic_input as _read_input,
    run_thermoelastic_context as _run_context,
    thermoelastic_report as _build_report,
    write_thermoelastic_hdf5 as _write_result,
)
from quantas.modules.thermoelasticity.models import (
    ThermoelasticAdiabaticMode as AdiabaticMode,
    ThermoelasticContext as Context,
    ThermoelasticDepthProfile as DepthProfile,
    ThermoelasticExtrapolationPolicy as ExtrapolationPolicy,
    ThermoelasticFitFailurePolicy as FitFailurePolicy,
    ThermoelasticInput as Input,
    ThermoelasticOptions as Options,
    ThermoelasticProfileResult as ProfileResult,
    ThermoelasticQualityPolicy as QualityPolicy,
    ThermoelasticReportLevel as ReportLevel,
    ThermoelasticResult as Result,
    ThermoelasticStabilityPolicy as StabilityPolicy,
    ThermoelasticTensorCondition as TensorCondition,
    select_stiffness_tensor,
)
from quantas.modules.thermoelasticity.plot import (
    ThermoelasticComparePlotOptions as ComparePlotOptions,
    ThermoelasticComponentGroup as ComponentGroup,
    ThermoelasticDomainPlotOptions as DomainPlotOptions,
    ThermoelasticFitPlotOptions as FitPlotOptions,
    ThermoelasticPTPlotOptions as PTPlotOptions,
    ThermoelasticPTQuantity as PTQuantity,
    ThermoelasticPlotLayout as PlotLayout,
    ThermoelasticPlotPreset as PlotPreset,
    ThermoelasticPlotStyleOptions as PlotStyleOptions,
    ThermoelasticProfileBackground as ProfileBackground,
    ThermoelasticProfileColor as ProfileColor,
    ThermoelasticProfileMode as ProfileMode,
    ThermoelasticProfilePlotOptions as ProfilePlotOptions,
    ThermoelasticUncertaintyMode as UncertaintyMode,
    build_thermoelastic_compare_plots as build_compare_plots,
    build_thermoelastic_domain_plot as build_domain_plot,
    build_thermoelastic_fit_plots as build_fit_plots,
    describe_thermoelastic_plots as _describe_plots,
    build_thermoelastic_profile_plots as build_profile_plots,
    build_thermoelastic_pt_plots as build_pt_plots,
    resolve_components,
)
from quantas.modules.thermoelasticity.postfit import (
    analyze_thermoelastic_profiles as analyze_profiles,
    analyze_thermoelastic_result as analyze_grid,
    regular_grid,
)
from quantas.modules.thermoelasticity.profiles import (
    build_thermoelastic_profile_preset as build_profile_preset,
    thermoelastic_profile_presets as profile_presets,
)
from quantas.modules.thermoelasticity.io.profile import (
    read_thermoelastic_depth_profile as read_depth_profile,
    read_thermoelastic_profile_spec as read_profile_spec,
    write_thermoelastic_profile_template as write_profile_template,
)
from quantas.modules.thermoelasticity.io.state_export import (
    write_thermoelastic_state_input as write_state_input,
)
from quantas.modules.thermoelasticity.io.table_export import (
    thermoelastic_point_table as point_table,
    write_thermoelastic_grid_table as write_grid_table,
    write_thermoelastic_profile_table as write_profile_table,
)
from quantas.modules.thermoelasticity.io.tensor_export import (
    thermoelastic_grid_info_table as grid_info_table,
    write_thermoelastic_tensor_export as write_tensor_export,
)

from .common import (
    PhononInputData,
    ReportTable,
    ResultData,
    _public_dir,
    get_result_payload,
)
from .plotting import PlotCollection, PlotInventory


def create_input(
    sources: str | Path | Sequence[str | Path],
    destination: str | Path,
    *,
    interface: str = "crystal",
    is_list: bool = False,
    jobname: str = "Quantas quasi-static thermoelastic input",
    reference: int | None = None,
    symprec: float = 1.0e-5,
    angle_tolerance: float = -1.0,
    elastic_tolerance: float = 1.0e-3,
    pressure_tolerance: float = 5.0e-2,
    structure_correspondence_tolerance: float = 5.0e-1,
) -> Path:
    """Create a thermoelastic YAML input from elastic-tensor outputs.

    Parameters
    ----------
    sources : str, Path, or sequence of str or Path
        Code-specific elastic-tensor outputs or a file-list path.
    destination : str or Path
        Destination YAML path.
    interface : str, optional
        Code-specific interface identifier.
    is_list : bool, optional
        Interpret a single source as a text file containing output paths.
    jobname : str, optional
        Human-readable workflow title.
    reference : int or None, optional
        Reference elastic-volume point.
    symprec : float, optional
        Cartesian symmetry tolerance.
    angle_tolerance : float, optional
        Symmetry angle tolerance in degrees.
    elastic_tolerance : float, optional
        Elastic symmetry tolerance in GPa.
    pressure_tolerance : float, optional
        Pressure consistency tolerance in GPa.
    structure_correspondence_tolerance : float, optional
        Atomic correspondence tolerance in angstrom.

    Returns
    -------
    Path
        Written thermoelastic YAML path.

    Raises
    ------
    ValueError
        If outputs are inconsistent or cannot be normalized to one elastic
        volume series.
    """
    return _create_input(
        sources,
        destination,
        interface=interface,
        is_list=is_list,
        jobname=jobname,
        reference=reference,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
        elastic_tolerance=elastic_tolerance,
        pressure_tolerance=pressure_tolerance,
        structure_correspondence_tolerance=structure_correspondence_tolerance,
    )


def read_input(source: str | Path) -> Input:
    """Read one Quantas thermoelastic input file.

    Parameters
    ----------
    source : str or Path
        Quantas thermoelastic YAML path.

    Returns
    -------
    Input
        Validated elastic-volume input contract.

    Raises
    ------
    ValueError
        If the YAML schema or elastic series is invalid.
    """
    return _read_input(source)


def normalize_input(source: Input | str | Path) -> Input:
    """Return a normalized thermoelastic input contract.

    Parameters
    ----------
    source : Input, str, or Path
        Existing input contract or YAML path.

    Returns
    -------
    Input
        Validated input suitable for coupling.

    Raises
    ------
    TypeError
        If the source type is unsupported.
    ValueError
        If a supplied YAML file is invalid.
    """
    return _normalize_input(source)


def prepare_context(
    input_data: Input | str | Path,
    qha_source: ResultData | QHAInput | PhononInputData | str | Path,
    *,
    qha_options: QHAOptions | None = None,
) -> Context:
    """Build the validated QHA--elasticity coupling context.

    Parameters
    ----------
    input_data : Input, str, or Path
        Elastic-volume series or thermoelastic YAML path.
    qha_source : ResultData, QHAInput, PhononInputData, str, or Path
        Existing QHA result, HDF5 path, or QHA input requiring calculation.
    qha_options : QHAOptions or None, optional
        Options used only when a QHA calculation is required.

    Returns
    -------
    Context
        Coupling context with volume coverage and data-completeness diagnostics.

    Raises
    ------
    ValueError
        If the QHA and elastic datasets cannot be coupled consistently.
    """
    from .interop import qha_to_thermoelastic_context

    return qha_to_thermoelastic_context(
        input_data,
        qha_source,
        qha_options=qha_options,
    )


def run_context(
    context: Context,
    *,
    options: Options | None = None,
    observer: Observer | None = None,
) -> ResultData:
    """Run a thermoelastic workflow from a prevalidated context.

    Parameters
    ----------
    context : Context
        Validated QHA--elasticity coupling context.
    options : Options or None, optional
        Fit, quality, stability, frame, and reconstruction controls.
    observer : Observer or None, optional
        Frontend-neutral workflow observer.

    Returns
    -------
    ResultData
        Calibrated thermoelastic fit result envelope.

    Raises
    ------
    ValueError
        If fitting, quality, or stability requirements are not satisfied.
    """
    return _run_context(context, options=options, observer=observer)


def run(
    input_data: Input | str | Path,
    qha_source: ResultData | QHAInput | PhononInputData | str | Path,
    *,
    options: Options | None = None,
    qha_options: QHAOptions | None = None,
    profiles: Sequence[DepthProfile] = (),
    observer: Observer | None = None,
) -> ResultData:
    """Run the complete QHA--elasticity thermoelastic workflow.

    Parameters
    ----------
    input_data : Input, str, or Path
        Elastic-volume series or thermoelastic YAML path.
    qha_source : ResultData, QHAInput, PhononInputData, str, or Path
        QHA result, native HDF5 path, or input requiring a QHA calculation.
    options : Options or None, optional
        Thermoelastic fitting and reconstruction controls.
    qha_options : QHAOptions or None, optional
        Options used only when ``qha_source`` must be calculated.
    profiles : sequence of DepthProfile, optional
        Depth paths evaluated after model calibration.
    observer : Observer or None, optional
        Frontend-neutral workflow observer.

    Returns
    -------
    ResultData
        Calibrated result envelope, optionally extended with profile analyses.

    Raises
    ------
    ValueError
        If coupling, fitting, reconstruction, or profile evaluation fails.

    Notes
    -----
    Optional depth profiles use the same post-fit engine as the CLI and future
    GUI, so all frontends share the same numerical path.
    """
    context = prepare_context(input_data, qha_source, qha_options=qha_options)
    calibrated = run_context(context, options=options, observer=observer)
    if profiles:
        return analyze_grid(calibrated, profiles=profiles)
    return calibrated


def get_result(result: ResultData) -> Result:
    """Return the typed thermoelastic payload from a result envelope.

    Parameters
    ----------
    result : ResultData
        Complete Quantas result envelope.

    Returns
    -------
    Result
        Module-specific thermoelastic result.

    Raises
    ------
    ValueError
        If the envelope is not a valid thermoelastic result.
    """
    return get_result_payload(
        result,
        module="thermoelasticity",
        key="thermoelasticity",
        expected_type=Result,
    )


def read_result(source: str | Path) -> ResultData:
    """Read a native Quantas thermoelastic HDF5 result.

    Parameters
    ----------
    source : str or Path
        Native Quantas HDF5 path.

    Returns
    -------
    ResultData
        Restored result envelope.

    Raises
    ------
    ValueError
        If the file is not a supported thermoelastic result.
    """
    return _read_result(source)


def write_result(
    result: ResultData,
    destination: str | Path,
    *,
    report_text: str | None = None,
) -> Path:
    """Write a native Quantas thermoelastic HDF5 result.

    Parameters
    ----------
    result : ResultData
        Complete thermoelastic result envelope.
    destination : str or Path
        Destination path.
    report_text : str or None, optional
        Deterministic report text to embed in diagnostics.

    Returns
    -------
    Path
        Final native HDF5 path.

    Raises
    ------
    ValueError
        If the result envelope is invalid.
    """
    get_result(result)
    return _write_result(result, destination, report_text=report_text)


def build_report(
    result: ResultData | Result,
    *,
    level: ReportLevel = "standard",
) -> list[ReportTable]:
    """Build frontend-neutral thermoelastic report tables.

    Parameters
    ----------
    result : ResultData or Result
        Complete envelope or module-specific thermoelastic result.
    level : {"standard", "extended", "debug"}, optional
        Scientific report detail.

    Returns
    -------
    list of ReportTable
        Ordered fit, reconstruction, stability, and provenance tables.

    Raises
    ------
    ValueError
        If the result is invalid.
    """
    return list(_build_report(result, level=level))


def describe_plots(result: ResultData) -> PlotInventory:
    """Return result-aware thermoelastic plot families and scientific context.

    Parameters
    ----------
    result : ResultData
        Complete thermoelastic result envelope.

    Returns
    -------
    PlotInventory
        Available calibration, P-T, profile, comparison, and domain families,
        together with exact stored grids, components, tensor conditions, and
        profile names.
    """
    get_result(result)
    return _describe_plots(result)


def build_plots(result: ResultData | Result) -> PlotCollection:
    """Build default plots appropriate for the result stage.

    Parameters
    ----------
    result : ResultData or Result
        Complete envelope or module-specific thermoelastic result.

    Returns
    -------
    PlotCollection
        Fit, pressure-temperature, or profile plots selected from the data
        available in the result.

    Raises
    ------
    ValueError
        If the result is invalid or lacks plottable data.
    """
    return _build_plots(result)


def __dir__() -> list[str]:
    """Return only names in the supported public contract."""
    return _public_dir(__all__)


__all__ = [
    "AdiabaticMode",
    "ComparePlotOptions",
    "ComponentGroup",
    "Context",
    "DepthProfile",
    "DomainPlotOptions",
    "ExtrapolationPolicy",
    "FitFailurePolicy",
    "FitPlotOptions",
    "Input",
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
    "ProfileResult",
    "QualityPolicy",
    "ReportLevel",
    "Result",
    "StabilityPolicy",
    "TensorCondition",
    "UncertaintyMode",
    "analyze_grid",
    "analyze_profiles",
    "build_compare_plots",
    "build_domain_plot",
    "describe_plots",
    "build_fit_plots",
    "build_plots",
    "build_profile_plots",
    "build_profile_preset",
    "build_pt_plots",
    "build_report",
    "create_input",
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
]
