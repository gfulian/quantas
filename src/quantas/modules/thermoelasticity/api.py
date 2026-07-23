# -*- coding: utf-8 -*-

"""Internal thermoelastic workflow facade used by :mod:`quantas.api.thermoelasticity`."""

from __future__ import annotations

from pathlib import Path
from collections.abc import Sequence

from quantas.core.events import Observer
from quantas.models import ModuleContract, PlotCollection, ResultData

from quantas.modules.thermoelasticity.analysis_engine import (
    ThermoelasticAnalysisEngine,
)
from quantas.modules.thermoelasticity.calculator import ThermoelasticityCalculator
from quantas.modules.thermoelasticity.io.export import ThermoelasticityHDF5Export
from quantas.modules.thermoelasticity.io.inpgen import (
    create_thermoelastic_input as _create_thermoelastic_input,
)
from quantas.modules.thermoelasticity.io.hdf5 import read_thermoelastic_hdf5
from quantas.modules.thermoelasticity.io.profile import (
    read_thermoelastic_depth_profile,
    read_thermoelastic_profile_spec,
    write_thermoelastic_profile_template,
)
from quantas.modules.thermoelasticity.io.state_export import (
    write_thermoelastic_state_input,
)
from quantas.modules.thermoelasticity.io.table_export import (
    thermoelastic_point_table,
    write_thermoelastic_grid_table,
    write_thermoelastic_profile_table,
)
from quantas.modules.thermoelasticity.io.tensor_export import (
    thermoelastic_grid_info_table,
    write_thermoelastic_tensor_export,
)
from quantas.modules.thermoelasticity.io.reader import read_thermoelastic_input
from quantas.modules.thermoelasticity.models import (
    ThermoelasticContext,
    ThermoelasticInput,
    ThermoelasticOptions,
    ThermoelasticReportLevel,
    ThermoelasticResult,
)
from quantas.modules.thermoelasticity.plot import (
    ThermoelasticComparePlotOptions,
    ThermoelasticDomainPlotOptions,
    ThermoelasticFitPlotOptions,
    ThermoelasticPTPlotOptions,
    ThermoelasticPlotStyleOptions,
    ThermoelasticProfilePlotOptions,
    build_thermoelastic_compare_plots,
    build_thermoelastic_domain_plot,
    build_thermoelastic_fit_plots,
    build_thermoelastic_profile_plots,
    build_thermoelastic_pt_plots,
)
from quantas.modules.thermoelasticity.postfit import (
    analyze_thermoelastic_profiles,
    analyze_thermoelastic_result,
)
from quantas.modules.thermoelasticity.profiles import (
    build_thermoelastic_profile_preset,
    thermoelastic_profile_presets,
)
from quantas.modules.thermoelasticity.report import (
    build_thermoelastic_analysis_report,
    build_thermoelastic_report,
)


def create_thermoelastic_input(
    sources: str | Path | Sequence[str | Path],
    outfile: str | Path,
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
    """Generate a thermoelastic YAML input from CRYSTAL SOEC outputs.

    Parameters
    ----------
    sources : str, Path, or sequence
        One CRYSTAL output, a text file listing outputs, or an explicit sequence
        of output paths.
    outfile : str or Path
        Destination YAML file.
    interface : str, optional
        Electronic-structure interface. Currently ``"crystal"`` only.
    is_list : bool, optional
        Interpret a scalar source as a text file containing one path per line.
    jobname : str, optional
        Human-readable input description.
    reference : int or None, optional
        Reference point index after sorting by increasing volume. By default,
        the point closest to zero pressure is used.
    symprec : float, optional
        spglib Cartesian tolerance in angstrom.
    angle_tolerance : float, optional
        spglib angular tolerance in degrees.
    elastic_tolerance : float, optional
        Elastic symmetry tolerance in GPa.
    pressure_tolerance : float, optional
        Maximum permitted difference, in GPa, between the CRYSTAL ``PRESSURE``
        keyword and the pressure reported for the corrected elastic tensor.
    structure_correspondence_tolerance : float, optional
        Maximum ordered-atom displacement in angstrom along the sampled
        structural path.

    Returns
    -------
    Path
        Written YAML path.
    """
    return _create_thermoelastic_input(
        sources,
        outfile,
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


def normalize_thermoelastic_input(
    input_data: ThermoelasticInput | str | Path,
) -> ThermoelasticInput:
    """Return a parsed thermoelastic input object.

    Parameters
    ----------
    input_data : ThermoelasticInput, str, or Path
        Existing passive input contract or YAML path.

    Returns
    -------
    ThermoelasticInput
        Normalized input object.

    Raises
    ------
    TypeError
        If the input type is unsupported.
    """
    if isinstance(input_data, ThermoelasticInput):
        return input_data
    if isinstance(input_data, (str, Path)):
        return read_thermoelastic_input(input_data)
    raise TypeError("input_data must be ThermoelasticInput or a YAML path")


def run_thermoelastic_context(
    context: ThermoelasticContext,
    *,
    options: ThermoelasticOptions | None = None,
    observer: Observer | None = None,
) -> ResultData:
    """Run quasi-static thermoelastic fitting for a validated context.

    Parameters
    ----------
    context : ThermoelasticContext
        Validated QHA and elastic-volume pairing.
    options : ThermoelasticOptions or None, optional
        Fitting and reconstruction controls.
    observer : Observer or None, optional
        Quantas workflow observer.

    Returns
    -------
    ResultData
        Result envelope containing a ``thermoelasticity`` payload.
    """
    return ThermoelasticityCalculator(
        context,
        options=options,
        observer=observer,
    ).execute()


def write_thermoelastic_hdf5(
    result: ResultData,
    filename: str | Path,
    *,
    report_text: str | None = None,
) -> Path:
    """Write a complete thermoelastic result and return the HDF5 path."""
    path = Path(filename).with_suffix(".hdf5")
    ThermoelasticityHDF5Export().export(
        result,
        path,
        report_text=report_text,
    )
    return path


def thermoelastic_report(
    result: ResultData | ThermoelasticResult,
    *,
    level: ThermoelasticReportLevel = "standard",
):
    """Return neutral report tables at the requested detail level."""
    payload = (
        result
        if isinstance(result, ThermoelasticResult)
        else result.results.get("thermoelasticity")
    )
    if not isinstance(payload, ThermoelasticResult):
        raise ValueError("result does not contain a thermoelasticity payload")
    return list(build_thermoelastic_report(payload, level=level))


def build_thermoelastic_plots(
    result: ResultData | ThermoelasticResult,
) -> PlotCollection:
    """Build default plots appropriate for the archived workflow stage."""
    payload = (
        result
        if isinstance(result, ThermoelasticResult)
        else result.results.get("thermoelasticity")
    )
    if not isinstance(payload, ThermoelasticResult):
        raise ValueError("result does not contain a thermoelasticity payload")
    if payload.profiles:
        return build_thermoelastic_profile_plots(payload)
    if payload.stiffness_isothermal is not None:
        return build_thermoelastic_pt_plots(payload)
    return build_thermoelastic_fit_plots(payload)


MODULE_CONTRACT = ModuleContract(
    name="thermoelasticity",
    result_key="thermoelasticity",
    read_input=read_thermoelastic_input,
    run=run_thermoelastic_context,
    read_hdf5=read_thermoelastic_hdf5,
    write_hdf5=write_thermoelastic_hdf5,
    build_report=thermoelastic_report,
    build_plots=build_thermoelastic_plots,
)


__all__ = [
    "MODULE_CONTRACT",
    "ThermoelasticAnalysisEngine",
    "ThermoelasticComparePlotOptions",
    "ThermoelasticDomainPlotOptions",
    "ThermoelasticFitPlotOptions",
    "ThermoelasticPTPlotOptions",
    "ThermoelasticPlotStyleOptions",
    "ThermoelasticProfilePlotOptions",
    "analyze_thermoelastic_profiles",
    "analyze_thermoelastic_result",
    "build_thermoelastic_compare_plots",
    "build_thermoelastic_plots",
    "build_thermoelastic_domain_plot",
    "build_thermoelastic_fit_plots",
    "build_thermoelastic_profile_plots",
    "build_thermoelastic_pt_plots",
    "build_thermoelastic_analysis_report",
    "build_thermoelastic_profile_preset",
    "build_thermoelastic_report",
    "create_thermoelastic_input",
    "normalize_thermoelastic_input",
    "read_thermoelastic_depth_profile",
    "read_thermoelastic_profile_spec",
    "read_thermoelastic_hdf5",
    "read_thermoelastic_input",
    "run_thermoelastic_context",
    "thermoelastic_grid_info_table",
    "thermoelastic_point_table",
    "thermoelastic_profile_presets",
    "thermoelastic_report",
    "write_thermoelastic_grid_table",
    "write_thermoelastic_profile_template",
    "write_thermoelastic_hdf5",
    "write_thermoelastic_profile_table",
    "write_thermoelastic_state_input",
    "write_thermoelastic_tensor_export",
]
