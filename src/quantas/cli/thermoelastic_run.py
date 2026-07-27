# -*- coding: utf-8 -*-

"""Calibration command for quasi-static thermoelasticity."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, TypedDict, cast

import click

from quantas.api.thermoelasticity import (
    prepare_context as prepare_thermoelastic_context,
)
from quantas.cli.contracts import (
    DOMAIN_GROUP,
    NUMERICAL_GROUP,
    PLOTTING_GROUP,
    SCIENTIFIC_GROUP,
    VALIDATION_GROUP,
    default_hdf5_path,
    default_report_path,
    figure_preset_option,
    force_option,
    output_option,
    parse_verbosity,
    progress_option,
    quiet_option,
    report_option,
    verbosity_option,
)
from quantas.cli.grouped_options import GroupedCommand, grouped_option
from quantas.cli.messages import quantas_finish, quantas_title
from quantas.cli.thermoelastic_common import (
    require_output_replacement,
    thermoelastic_payload,
)
from quantas.cli.thermoelastic_observer import ThermoelasticTextObserver
from quantas.cli.thermoelastic_plot_common import render_collection
from quantas.api.common import EventLevel
from quantas.api.qha import (
    Minimization as QHAMinimization,
    Options as QHAOptions,
    Scheme as QHAScheme,
)
from quantas.api.thermoelasticity import (
    FitMethod,
    FitPlotOptions as ThermoelasticFitPlotOptions,
    Options as ThermoelasticOptions,
    PlotStyleOptions as ThermoelasticPlotStyleOptions,
    build_fit_plots as build_thermoelastic_fit_plots,
    build_report as build_thermoelastic_report,
    run_context as run_thermoelastic_context,
    write_result as write_thermoelastic_hdf5,
)
from quantas.references import module_citation_keys, render_citation_notice
from quantas.api.thermoelasticity import (
    AdiabaticMode as ThermoelasticAdiabaticMode,
    ExtrapolationPolicy as ThermoelasticExtrapolationPolicy,
    FitFailurePolicy as ThermoelasticFitFailurePolicy,
    QualityPolicy as ThermoelasticQualityPolicy,
    ReportLevel as ThermoelasticReportLevel,
    StabilityPolicy as ThermoelasticStabilityPolicy,
)


@click.command(name="run", cls=GroupedCommand)
@click.argument(
    "input_file", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.argument(
    "qha_source", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@grouped_option(
    "--validation",
    group=VALIDATION_GROUP,
    type=click.Choice(["standard", "strict", "exploratory"], case_sensitive=False),
    default="standard",
    show_default=True,
    help=(
        "Scientific validation profile. Standard preserves balanced warnings; "
        "strict rejects weak support and extrapolation; exploratory retains "
        "all finite results with diagnostics."
    ),
)
@grouped_option(
    "--plot/--no-plot",
    group=PLOTTING_GROUP,
    default=False,
    show_default=True,
    help="Generate observed-versus-fitted diagnostic figures after calibration.",
)
@grouped_option(
    "--plot-output",
    group=PLOTTING_GROUP,
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Fit-figure directory. Default: fit archive base name + '_plots'.",
)
@grouped_option(
    "--plot-format",
    group=PLOTTING_GROUP,
    type=click.Choice(["png", "pdf", "svg"], case_sensitive=False),
    default="png",
    show_default=True,
    help="Static format for optional fit-diagnostic figures.",
)
@figure_preset_option(
    option_name="--plot-preset",
    parameter_name="plot_preset",
    group=PLOTTING_GROUP,
)
@grouped_option(
    "--reference-eos",
    group=SCIENTIFIC_GROUP,
    type=click.Choice(["BM2", "BM3", "BM4"], case_sensitive=True),
    default="BM3",
    show_default=True,
    help="Static energy EOS providing fixed V0, K0 and Kprime.",
)
@grouped_option(
    "--finite-strain-order",
    group=SCIENTIFIC_GROUP,
    type=click.IntRange(2, 3),
    default=3,
    show_default=True,
    help="Eulerian finite-strain order for Cij(V).",
)
@grouped_option(
    "--zero-tolerance",
    group=SCIENTIFIC_GROUP,
    hidden=True,
    type=click.FloatRange(min=0.0),
    default=1.0e-10,
    show_default=True,
    help="Maximum absolute GPa value retained as an exact zero component.",
)
@grouped_option(
    "--relative-error-floor",
    group=SCIENTIFIC_GROUP,
    hidden=True,
    type=click.FloatRange(min=0.0, min_open=True),
    default=1.0e-8,
    show_default=True,
    help="GPa denominator floor for relative-residual diagnostics.",
)
@grouped_option(
    "--max-iterations",
    group=NUMERICAL_GROUP,
    hidden=True,
    type=click.IntRange(min=1),
    default=None,
    help="Maximum component-fit model evaluations.",
)
@grouped_option(
    "--fit-failure",
    "fit_failure_policy",
    group=NUMERICAL_GROUP,
    hidden=True,
    type=click.Choice(["stop", "continue", "raise"], case_sensitive=False),
    default="stop",
    show_default=True,
    help="Action when one independent elastic-component fit fails.",
)
@grouped_option(
    "--extrapolation",
    "extrapolation_policy",
    group=VALIDATION_GROUP,
    hidden=True,
    type=click.Choice(["fail", "warn", "allow"], case_sensitive=False),
    default="warn",
    show_default=True,
    help="Policy for QHA volumes outside the sampled elastic-volume interval.",
)
@grouped_option(
    "--fit-quality",
    "quality_policy",
    group=VALIDATION_GROUP,
    hidden=True,
    type=click.Choice(["fail", "warn", "allow"], case_sensitive=False),
    default="warn",
    show_default=True,
    help=(
        "Policy for component fits classified as scientifically unsupported; "
        "cautions are reported but retained."
    ),
)
@grouped_option(
    "--minimum-fit-points",
    group=VALIDATION_GROUP,
    hidden=True,
    type=click.IntRange(min=3),
    default=4,
    show_default=True,
    help="Preferred minimum number of elastic-volume observations per component.",
)
@grouped_option(
    "--minimum-strain-span",
    "minimum_eulerian_strain_span",
    group=VALIDATION_GROUP,
    hidden=True,
    type=click.FloatRange(min=0.0),
    default=5.0e-3,
    show_default=True,
    help="Preferred minimum sampled Eulerian finite-strain span.",
)
@grouped_option(
    "--maximum-design-condition",
    "maximum_design_condition_number",
    group=VALIDATION_GROUP,
    hidden=True,
    type=click.FloatRange(min=0.0, min_open=True),
    default=1.0e6,
    show_default=True,
    help="Preferred maximum condition number of the normalized fit design matrix.",
)
@grouped_option(
    "--maximum-symmetry-spread",
    "maximum_relative_symmetry_spread",
    group=VALIDATION_GROUP,
    hidden=True,
    type=click.FloatRange(min=0.0),
    default=1.0e-2,
    show_default=True,
    help="Preferred maximum relative spread among symmetry-equivalent entries.",
)
@grouped_option(
    "--maximum-loo-change",
    "maximum_leave_one_out_parameter_change",
    group=VALIDATION_GROUP,
    hidden=True,
    type=click.FloatRange(min=0.0),
    default=5.0e-1,
    show_default=True,
    help="Preferred maximum scale-aware leave-one-out parameter change.",
)
@grouped_option(
    "--maximum-order-change",
    "maximum_order_parameter_change",
    group=VALIDATION_GROUP,
    hidden=True,
    type=click.FloatRange(min=0.0),
    default=5.0e-1,
    show_default=True,
    help="Preferred maximum parameter change between finite-strain orders.",
)
@grouped_option(
    "--stability",
    "stability_policy",
    group=VALIDATION_GROUP,
    hidden=True,
    type=click.Choice(["fail", "warn", "allow"], case_sensitive=False),
    default="warn",
    show_default=True,
    help="Policy for unstable or indeterminate reconstructed Wallace tensors.",
)
@grouped_option(
    "--stability-tolerance",
    group=VALIDATION_GROUP,
    hidden=True,
    type=float,
    default=0.0,
    show_default=True,
    help="Minimum accepted stiffness eigenvalue in GPa.",
)
@grouped_option(
    "--adiabatic",
    "adiabatic_mode",
    group=SCIENTIFIC_GROUP,
    type=click.Choice(["auto", "off", "require"], case_sensitive=False),
    default="auto",
    show_default=True,
    help=(
        "Convert C^T to C^S automatically when QHA C_V and the Cartesian "
        "thermal-expansion tensor are available, disable conversion, or "
        "require complete inputs."
    ),
)
@grouped_option(
    "--adiabatic-uncertainty/--no-adiabatic-uncertainty",
    group=SCIENTIFIC_GROUP,
    default=True,
    show_default=True,
    help="Propagate available QHA and isothermal-tensor uncertainties to C^S.",
)
@grouped_option(
    "--eos-covariance/--no-eos-covariance",
    group=SCIENTIFIC_GROUP,
    default=True,
    show_default=True,
    help="Propagate shared reference-EOS covariance into predicted tensors.",
)
@grouped_option(
    "--volume-uncertainty/--no-volume-uncertainty",
    group=SCIENTIFIC_GROUP,
    default=True,
    show_default=True,
    help="Propagate QHA equilibrium-volume uncertainty when available.",
)
@grouped_option(
    "--qha-temperature",
    group=DOMAIN_GROUP,
    nargs=3,
    type=float,
    metavar="MIN MAX STEP",
    default=(298.15, 298.15, 1.0),
    show_default="298.15 298.15 1.0",
    help="QHA temperature grid used only when QHA_SOURCE is YAML.",
)
@grouped_option(
    "--qha-pressure",
    group=DOMAIN_GROUP,
    nargs=3,
    type=float,
    metavar="MIN MAX STEP",
    default=(0.0, 0.0, 1.0),
    show_default="0 0 1",
    help="QHA pressure grid used only when QHA_SOURCE is YAML.",
)
@grouped_option(
    "--qha-scheme",
    group=DOMAIN_GROUP,
    type=click.Choice(["freq", "td"], case_sensitive=True),
    default="freq",
    show_default=True,
    help="QHA interpolation scheme used only for YAML sources.",
)
@grouped_option(
    "--qha-minimization",
    group=DOMAIN_GROUP,
    type=click.Choice(["poly", "eos"], case_sensitive=True),
    default="poly",
    show_default=True,
    help="QHA volume-minimization method used only for YAML sources.",
)
@grouped_option(
    "--qha-eos",
    group=DOMAIN_GROUP,
    default="BM3",
    show_default=True,
    help="QHA energy EOS used when --qha-minimization=eos.",
)
@grouped_option(
    "--qha-degree",
    group=DOMAIN_GROUP,
    type=click.IntRange(min=1),
    default=3,
    show_default=True,
    help="QHA polynomial degrees for energy, free energy, frequencies and structure.",
)
@output_option(
    help="Fit HDF5 archive. Default: thermoelastic input base name + '_FIT.hdf5'.",
)
@force_option()
@report_option()
@verbosity_option()
@quiet_option()
@progress_option()
def run(
    input_file: Path,
    qha_source: Path,
    output: Path | None,
    force: bool,
    report: Path | None,
    verbosity: str,
    quiet: bool,
    progress: bool,
    validation: str,
    plot: bool,
    plot_output: Path | None,
    plot_format: str,
    plot_preset: str,
    reference_eos: str,
    finite_strain_order: int,
    zero_tolerance: float,
    relative_error_floor: float,
    max_iterations: int | None,
    fit_failure_policy: str,
    extrapolation_policy: str,
    quality_policy: str,
    minimum_fit_points: int,
    minimum_eulerian_strain_span: float,
    maximum_design_condition_number: float,
    maximum_relative_symmetry_spread: float,
    maximum_leave_one_out_parameter_change: float,
    maximum_order_parameter_change: float,
    stability_policy: str,
    stability_tolerance: float,
    adiabatic_mode: str,
    adiabatic_uncertainty: bool,
    eos_covariance: bool,
    volume_uncertainty: bool,
    qha_temperature: tuple[float, float, float],
    qha_pressure: tuple[float, float, float],
    qha_scheme: str,
    qha_minimization: str,
    qha_eos: str,
    qha_degree: int,
) -> None:
    """Fit the static EOS and independent Cij(V) components into HDF5."""
    destination = default_hdf5_path(input_file, output, suffix="_FIT")
    report = default_report_path(input_file, report)
    report_verbosity = parse_verbosity(verbosity)
    require_output_replacement(destination, force)
    observer = ThermoelasticTextObserver(
        report,
        silent=quiet,
        show_progress=progress,
        verbosity=report_verbosity,
    )
    try:
        qha_options = _build_qha_options(
            temperature=qha_temperature,
            pressure=qha_pressure,
            scheme=qha_scheme,
            minimization=qha_minimization,
            eos=qha_eos,
            degree=qha_degree,
        )
        options = _build_thermoelastic_options(
            validation=validation,
            reference_eos=reference_eos,
            finite_strain_order=finite_strain_order,
            zero_tolerance=zero_tolerance,
            relative_error_floor=relative_error_floor,
            max_iterations=max_iterations,
            fit_failure_policy=fit_failure_policy,
            extrapolation_policy=extrapolation_policy,
            quality_policy=quality_policy,
            minimum_fit_points=minimum_fit_points,
            minimum_eulerian_strain_span=minimum_eulerian_strain_span,
            maximum_design_condition_number=maximum_design_condition_number,
            maximum_relative_symmetry_spread=maximum_relative_symmetry_spread,
            maximum_leave_one_out_parameter_change=(
                maximum_leave_one_out_parameter_change
            ),
            maximum_order_parameter_change=maximum_order_parameter_change,
            stability_policy=stability_policy,
            stability_tolerance=stability_tolerance,
            adiabatic_mode=adiabatic_mode,
            adiabatic_uncertainty=adiabatic_uncertainty,
            report_level=report_verbosity.value,
            eos_covariance=eos_covariance,
            volume_uncertainty=volume_uncertainty,
        )
        _execute_calibration(
            input_file=input_file,
            qha_source=qha_source,
            destination=destination,
            qha_options=qha_options,
            options=options,
            validation=validation,
            observer=observer,
            plot=plot,
            plot_output=plot_output,
            plot_format=plot_format,
            plot_preset=plot_preset,
        )
    except click.ClickException:
        raise
    except Exception as exc:
        raise click.ClickException(str(exc)) from exc
    finally:
        observer.close()


def _build_qha_options(
    *,
    temperature: tuple[float, float, float],
    pressure: tuple[float, float, float],
    scheme: str,
    minimization: str,
    eos: str,
    degree: int,
) -> QHAOptions:
    """Build QHA controls used only when the source is a YAML input."""
    return QHAOptions(
        temperature_min=temperature[0],
        temperature_max=temperature[1],
        temperature_step=temperature[2],
        pressure_min=pressure[0],
        pressure_max=pressure[1],
        pressure_step=pressure[2],
        scheme=cast(QHAScheme, scheme),
        minimization=cast(QHAMinimization, minimization),
        eos=eos,
        energy_degree=degree,
        free_energy_degree=degree,
        frequency_degree=degree,
        structural_degree=degree,
    )


def _build_thermoelastic_options(
    *,
    validation: str,
    reference_eos: str,
    finite_strain_order: int,
    zero_tolerance: float,
    relative_error_floor: float,
    max_iterations: int | None,
    fit_failure_policy: str,
    extrapolation_policy: str,
    quality_policy: str,
    minimum_fit_points: int,
    minimum_eulerian_strain_span: float,
    maximum_design_condition_number: float,
    maximum_relative_symmetry_spread: float,
    maximum_leave_one_out_parameter_change: float,
    maximum_order_parameter_change: float,
    stability_policy: str,
    stability_tolerance: float,
    adiabatic_mode: str,
    adiabatic_uncertainty: bool,
    report_level: str,
    eos_covariance: bool,
    volume_uncertainty: bool,
) -> ThermoelasticOptions:
    """Resolve validation presets and construct calibration options."""
    values = _validation_values(validation)
    if validation.lower() != "standard":
        extrapolation_policy = values["extrapolation_policy"]
        quality_policy = values["quality_policy"]
        stability_policy = values["stability_policy"]
        minimum_fit_points = values["minimum_fit_points"]
        minimum_eulerian_strain_span = values["minimum_eulerian_strain_span"]
        maximum_design_condition_number = values["maximum_design_condition_number"]
        maximum_relative_symmetry_spread = values["maximum_relative_symmetry_spread"]
        maximum_leave_one_out_parameter_change = values[
            "maximum_leave_one_out_parameter_change"
        ]
        maximum_order_parameter_change = values["maximum_order_parameter_change"]
    return ThermoelasticOptions(
        reference_eos=reference_eos,
        finite_strain_order=cast(Literal[2, 3], finite_strain_order),
        fit_method=FitMethod.OLS,
        max_iterations=max_iterations,
        zero_tolerance=zero_tolerance,
        relative_error_floor=relative_error_floor,
        fit_failure_policy=cast(ThermoelasticFitFailurePolicy, fit_failure_policy),
        extrapolation_policy=cast(
            ThermoelasticExtrapolationPolicy, extrapolation_policy.lower()
        ),
        quality_policy=cast(ThermoelasticQualityPolicy, quality_policy.lower()),
        minimum_fit_points=minimum_fit_points,
        minimum_eulerian_strain_span=minimum_eulerian_strain_span,
        maximum_design_condition_number=maximum_design_condition_number,
        maximum_relative_symmetry_spread=maximum_relative_symmetry_spread,
        maximum_leave_one_out_parameter_change=(maximum_leave_one_out_parameter_change),
        maximum_order_parameter_change=maximum_order_parameter_change,
        stability_policy=cast(ThermoelasticStabilityPolicy, stability_policy.lower()),
        stability_tolerance=stability_tolerance,
        adiabatic_mode=cast(ThermoelasticAdiabaticMode, adiabatic_mode.lower()),
        propagate_adiabatic_uncertainty=adiabatic_uncertainty,
        report_level=cast(ThermoelasticReportLevel, report_level.lower()),
        propagate_reference_eos_covariance=eos_covariance,
        propagate_volume_uncertainty=volume_uncertainty,
    )


def _execute_calibration(
    *,
    input_file: Path,
    qha_source: Path,
    destination: Path,
    qha_options: QHAOptions,
    options: ThermoelasticOptions,
    validation: str,
    observer: ThermoelasticTextObserver,
    plot: bool,
    plot_output: Path | None,
    plot_format: str,
    plot_preset: str,
) -> None:
    """Execute calibration, persistence, optional plots, and final guidance."""
    observer.output.message(quantas_title(), bold=True)
    context = prepare_thermoelastic_context(
        input_file, qha_source, qha_options=qha_options
    )
    result_data = run_thermoelastic_context(
        context,
        options=options,
        observer=observer,
    )
    result_data.options["validation_preset"] = validation.lower()
    payload = thermoelastic_payload(result_data)
    observer.output.tables(
        build_thermoelastic_report(payload, level=options.report_level)
    )
    for warning in result_data.warnings:
        observer.output.message(warning, level=EventLevel.WARNING)
    written = write_thermoelastic_hdf5(
        result_data,
        destination,
        report_text=observer.text(),
    )
    observer.output.message("Thermoelastic model calibration completed.")
    observer.output.message(f"Reusable fit archive: {written}")
    if plot:
        _render_fit_plots(
            result_data,
            archive=written,
            output_dir=plot_output,
            image_format=plot_format,
            preset=plot_preset,
        )
    _write_next_steps(observer, written)
    observer.output.text_block(
        render_citation_notice(module_citation_keys("thermoelasticity"))
    )
    observer.output.message(quantas_finish())
    observer.save()


def _render_fit_plots(
    result_data,
    *,
    archive: Path,
    output_dir: Path | None,
    image_format: str,
    preset: str,
) -> None:
    """Render the standard fit diagnostics requested by ``run --plot``."""
    cli_preset = "analysis" if preset.lower() == "screen" else preset.lower()
    collection = build_thermoelastic_fit_plots(
        result_data,
        options=ThermoelasticFitPlotOptions(
            style=ThermoelasticPlotStyleOptions(
                preset=cast(
                    Literal["analysis", "publication", "monochrome"], cli_preset
                )
            )
        ),
    )
    render_collection(
        archive,
        collection,
        family="fit",
        output_dir=(
            archive.with_name(archive.stem + "_plots")
            if output_dir is None
            else output_dir
        ),
        image_format=image_format.lower(),
        preset=preset,
        dpi=None,
        transparent=False,
        show=False,
        figure_size=(7.0, 6.5),
        axis_label_font_size=None,
        legend_font_size=None,
        title_font_size=None,
        tick_label_font_size=None,
    )


def _write_next_steps(observer: ThermoelasticTextObserver, archive: Path) -> None:
    """Write concise stage-aware follow-up commands."""
    observer.output.message("Next steps:")
    observer.output.message(f"  quantas thermoelasticity inspect {archive}")
    observer.output.message(
        f"  quantas thermoelasticity analysis point {archive} 0 300"
    )
    observer.output.message(
        f"  quantas thermoelasticity analysis grid {archive} -P 0 10 1 -T 300 1500 100"
    )
    observer.output.message(
        f"  quantas thermoelasticity analysis profile {archive} "
        "--preset mantle-katsura-2022"
    )


class _ValidationValues(TypedDict):
    """Resolved scientific-support controls for one validation preset."""

    extrapolation_policy: str
    quality_policy: str
    stability_policy: str
    minimum_fit_points: int
    minimum_eulerian_strain_span: float
    maximum_design_condition_number: float
    maximum_relative_symmetry_spread: float
    maximum_leave_one_out_parameter_change: float
    maximum_order_parameter_change: float


def _validation_values(name: str) -> _ValidationValues:
    """Return a complete named scientific-support policy."""
    profiles: dict[str, _ValidationValues] = {
        "standard": {
            "extrapolation_policy": "warn",
            "quality_policy": "warn",
            "stability_policy": "warn",
            "minimum_fit_points": 4,
            "minimum_eulerian_strain_span": 5.0e-3,
            "maximum_design_condition_number": 1.0e6,
            "maximum_relative_symmetry_spread": 1.0e-2,
            "maximum_leave_one_out_parameter_change": 5.0e-1,
            "maximum_order_parameter_change": 5.0e-1,
        },
        "strict": {
            "extrapolation_policy": "fail",
            "quality_policy": "fail",
            "stability_policy": "fail",
            "minimum_fit_points": 5,
            "minimum_eulerian_strain_span": 7.5e-3,
            "maximum_design_condition_number": 1.0e5,
            "maximum_relative_symmetry_spread": 5.0e-3,
            "maximum_leave_one_out_parameter_change": 2.5e-1,
            "maximum_order_parameter_change": 2.5e-1,
        },
        "exploratory": {
            "extrapolation_policy": "allow",
            "quality_policy": "allow",
            "stability_policy": "allow",
            "minimum_fit_points": 3,
            "minimum_eulerian_strain_span": 0.0,
            "maximum_design_condition_number": 1.0e12,
            "maximum_relative_symmetry_spread": 1.0,
            "maximum_leave_one_out_parameter_change": 10.0,
            "maximum_order_parameter_change": 10.0,
        },
    }
    try:
        return profiles[name.lower()]
    except KeyError as exc:
        raise ValueError(f"unknown validation preset: {name}") from exc


__all__ = ["run"]
