# -*- coding: utf-8 -*-

"""EOS batch fitting command."""

from __future__ import annotations

from pathlib import Path
import traceback
import click

from quantas.cli.contracts import (
    OUTPUT_GROUP,
    SCIENTIFIC_GROUP,
    UNITS_GROUP,
    VALIDATION_GROUP,
    default_hdf5_path,
    default_report_path,
    force_option,
    output_option,
    parse_verbosity,
    quiet_option,
    report_option,
    verbosity_option,
)
from quantas.cli.grouped_options import GroupedCommand, grouped_option
from quantas.cli.eos_helpers import (
    _EOSCLIEventObserver,
    _build_request,
    _build_solver_options,
    _configure_plan_debug,
    _merge_spec_report_options,
    _parse_constraint_specs,
    _reject_spec_cli_conflicts,
    _resolve_targets,
    _validate_resolved_plan,
)
from quantas.cli.messages import quantas_finish, quantas_title
from quantas.cli.output import CLIOutput
from quantas.core.events import EventLevel
from quantas.api.eos import (
    BatchFailurePolicy as EOSBatchFailurePolicy,
    BatchJob as EOSBatchJob,
    BatchPlan as EOSBatchPlan,
    FitDomain as EOSFitDomain,
    ReportDetail as EOSReportDetail,
    ReportOptions as EOSReportOptions,
    build_batch_preamble,
    build_batch_report,
    read_input as read_eos_input,
    read_spec as read_eos_spec,
    resolve_spec as resolve_eos_spec,
    run_batch as run_eos_batch,
)
from quantas.references import module_citation_keys, render_citation_notice


@click.command(name="run", cls=GroupedCommand)
@click.argument(
    "filename", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@grouped_option(
    "--spec",
    "spec_path",
    group=SCIENTIFIC_GROUP,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Strict QUANTAS EOS SPEC 1 batch specification. File suffix is irrelevant.",
)
@grouped_option(
    "--pressure-unit",
    "--punit",
    group=UNITS_GROUP,
    default=None,
    metavar="UNIT",
    help="Input pressure unit override. File declaration or GPa is used when omitted.",
)
@grouped_option(
    "--length-unit",
    "--lunit",
    group=UNITS_GROUP,
    default=None,
    metavar="UNIT",
    help="Input length unit override; the corresponding cubic unit is used for volume.",
)
@grouped_option(
    "--temperature-unit",
    "--tunit",
    group=UNITS_GROUP,
    type=click.Choice(["K", "C", "F"], case_sensitive=False),
    default=None,
    help="Input temperature scale override. Normalized calculations always use kelvin.",
)
@grouped_option(
    "--domain",
    group="Fit selection",
    type=click.Choice(["pv", "vt", "pvt"], case_sensitive=False),
    default="pv",
    show_default=True,
    help="Scientific relationship fitted by every generated job.",
)
@grouped_option(
    "--fit",
    "targets",
    group="Fit selection",
    multiple=True,
    type=click.Choice(["volume", "a", "b", "c", "all"], case_sensitive=False),
    default=("volume",),
    show_default=True,
    help="Target to fit. Repeat the option or use 'all' for all available targets.",
)
@grouped_option(
    "--pv-eos",
    group="P-V model",
    default="birch-murnaghan",
    show_default=True,
    metavar="FAMILY",
    help="P-V EOS family or alias: M, BM, NS/PT, Vinet, or Tait.",
)
@grouped_option(
    "--pv-order",
    group="P-V model",
    type=click.IntRange(2, 4),
    default=3,
    show_default=True,
    help="P-V EOS order. Ignored only by the orderless Murnaghan family.",
)
@grouped_option(
    "--vt-eos",
    group="V-T model",
    default="berman",
    show_default=True,
    metavar="FAMILY",
    help="V-T family or alias: Berman, Fei, MHP, Salje, or KHP.",
)
@grouped_option(
    "--vt-variant",
    group="V-T model",
    default=None,
    metavar="VARIANT",
    help="Optional V-T variant, for example linear, quadratic, general, or simplified.",
)
@grouped_option(
    "--pvt-coupling",
    group="P-V-T coupling",
    type=click.Choice(
        ["linear", "anderson-gruneisen", "thermal-pressure"], case_sensitive=False
    ),
    default="linear",
    show_default=True,
    help=(
        "P-V-T coupling prescription. Thermal-pressure coupling defines its "
        "own thermal model and does not use --vt-eos."
    ),
)
@grouped_option(
    "--thermal-pressure-model",
    group="P-V-T coupling",
    type=click.Choice(
        ["holland-powell-einstein", "mgd", "mgd:q-compromise"],
        case_sensitive=False,
    ),
    default="holland-powell-einstein",
    show_default=True,
    help=(
        "Oscillator model for thermal-pressure coupling. MGD requires "
        "--atoms-per-cell or --formula with --formula-units-per-cell."
    ),
)
@grouped_option(
    "--atoms-per-cell",
    group="MGD normalization",
    type=click.FloatRange(min=0.0, min_open=True),
    default=None,
    help="Number of atoms in the crystallographic cell for MGD cell volumes.",
)
@grouped_option(
    "--formula",
    group="MGD normalization",
    default=None,
    metavar="FORMULA",
    help="Chemical formula of one formula unit, for example NaF.",
)
@grouped_option(
    "--formula-units-per-cell",
    group="MGD normalization",
    type=click.FloatRange(min=0.0, min_open=True),
    default=None,
    metavar="Z",
    help="Formula units per crystallographic cell; used with --formula.",
)
@grouped_option(
    "--fix",
    "fixed_parameters",
    group="Parameter constraints",
    multiple=True,
    metavar="[TARGET:]NAME=VALUE",
    help="Fix a parameter. Repeat as needed; an optional target prefix limits its scope.",
)
@grouped_option(
    "--initial",
    "initial_parameters",
    group="Parameter constraints",
    multiple=True,
    metavar="[TARGET:]NAME=VALUE",
    help="Override an initial parameter value without fixing it.",
)
@grouped_option(
    "--bound",
    "parameter_bounds",
    group="Parameter constraints",
    multiple=True,
    metavar="[TARGET:]NAME=LOW:HIGH",
    help="Override parameter bounds. Use 'none' for an open side.",
)
@grouped_option(
    "--solver",
    group="Numerical solver",
    type=click.Choice(
        ["ols", "wls", "effective-variance", "odr"], case_sensitive=False
    ),
    default="ols",
    show_default=True,
    help="Regression strategy. Weighted methods use the required uncertainties explicitly.",
)
@grouped_option(
    "--max-iterations",
    group="Numerical solver",
    type=click.IntRange(min=1),
    default=None,
    help="Maximum solver iterations or model evaluations, depending on the selected solver.",
)
@grouped_option(
    "--ftol",
    group="Numerical solver",
    type=click.FloatRange(min=0.0, min_open=True),
    default=None,
    help="Function convergence tolerance.",
)
@grouped_option(
    "--xtol",
    group="Numerical solver",
    type=click.FloatRange(min=0.0, min_open=True),
    default=None,
    help="Parameter convergence tolerance.",
)
@grouped_option(
    "--gtol",
    group="Numerical solver",
    type=click.FloatRange(min=0.0, min_open=True),
    default=None,
    help="Gradient convergence tolerance where supported.",
)
@grouped_option(
    "--covariance-scaling",
    group="Numerical solver",
    type=click.Choice(
        ["absolute", "reduced-chi-square", "inflate-only"], case_sensitive=False
    ),
    default=None,
    help="Covariance policy. Weighted EOS fits default to EosFit-like inflate-only scaling.",
)
@grouped_option(
    "--inner-max-iterations",
    group="Effective variance",
    type=click.IntRange(min=1),
    default=None,
    help="Maximum WLS evaluations inside each effective-variance cycle.",
)
@grouped_option(
    "--odr-difference",
    group="ODR",
    type=click.Choice(["central", "forward"], case_sensitive=False),
    default="central",
    show_default=True,
    help="Finite-difference scheme used by the ODRPACK95 backend.",
)
@grouped_option(
    "--odr-ndigit",
    group="ODR",
    type=click.IntRange(min=1),
    default=None,
    help="Reliable decimal digits in ODR model evaluations.",
)
@grouped_option(
    "--dry-run",
    group="Batch execution",
    is_flag=True,
    default=False,
    help="Resolve and validate the complete plan without fitting or creating HDF5.",
)
@grouped_option(
    "--failure-policy",
    group="Batch execution",
    type=click.Choice(["stop", "continue"], case_sensitive=False),
    default="stop",
    show_default=True,
    help="Stop after a failed job or continue with later jobs.",
)
@grouped_option(
    "--show-uncertainties",
    group=OUTPUT_GROUP,
    is_flag=True,
    default=False,
    help="Print a separate table of input standard uncertainties.",
)
@grouped_option(
    "--max-data-rows",
    group=OUTPUT_GROUP,
    type=click.IntRange(min=1),
    default=None,
    help="Limit displayed input rows; all rows are shown when omitted.",
)
@grouped_option(
    "--traceback/--no-traceback",
    "show_traceback",
    group=VALIDATION_GROUP,
    default=False,
    show_default=True,
    help="Include a Python traceback in the report when an unexpected error occurs.",
)
@output_option(
    help="Native EOS HDF5 archive. Default: input base name + '_EOS.hdf5'.",
)
@force_option()
@report_option()
@verbosity_option()
@quiet_option()
@click.pass_context
def run(
    ctx: click.Context,
    filename: Path,
    spec_path: Path | None,
    pressure_unit: str | None,
    length_unit: str | None,
    temperature_unit: str | None,
    domain: str,
    targets: tuple[str, ...],
    pv_eos: str,
    pv_order: int,
    vt_eos: str,
    vt_variant: str | None,
    pvt_coupling: str,
    thermal_pressure_model: str,
    atoms_per_cell: float | None,
    formula: str | None,
    formula_units_per_cell: float | None,
    fixed_parameters: tuple[str, ...],
    initial_parameters: tuple[str, ...],
    parameter_bounds: tuple[str, ...],
    solver: str,
    max_iterations: int | None,
    ftol: float | None,
    xtol: float | None,
    gtol: float | None,
    covariance_scaling: str | None,
    inner_max_iterations: int | None,
    odr_difference: str,
    odr_ndigit: int | None,
    dry_run: bool,
    failure_policy: str,
    show_uncertainties: bool,
    max_data_rows: int | None,
    show_traceback: bool,
    output: Path | None,
    force: bool,
    report: Path | None,
    verbosity: str,
    quiet: bool,
) -> None:
    """Run one non-interactive EOS batch from a data file and optional spec."""
    destination = default_hdf5_path(filename, output, suffix="_EOS")
    report = default_report_path(filename, report)
    report_verbosity = parse_verbosity(verbosity)
    debug = report_verbosity.includes_debug
    if destination.exists() and not force and not dry_run:
        raise click.UsageError(
            f"output archive exists: {destination}; use --force to replace it"
        )
    output_service = CLIOutput(report_file=report, silent=quiet, show_progress=False)
    try:
        output_service.message(quantas_title(), bold=True)
        if spec_path is not None:
            _reject_spec_cli_conflicts(ctx)
            document = read_eos_spec(spec_path)
            input_options = document.input_options
            dataset = read_eos_input(
                filename,
                pressure_unit=input_options.pressure_unit,
                length_unit=input_options.length_unit,
                temperature_unit=input_options.temperature_unit,
            )
            resolved_spec = resolve_eos_spec(document, dataset)
            plan = resolved_spec.plan
            report_options = _merge_spec_report_options(
                ctx,
                resolved_spec.report_options,
                verbosity=report_verbosity,
                show_uncertainties=show_uncertainties,
                max_data_rows=max_data_rows,
            )
            output_service.message(f"EOS specification: {spec_path}")
            creator = "quantas eos run --spec"
        else:
            dataset = read_eos_input(
                filename,
                pressure_unit=pressure_unit,
                length_unit=length_unit,
                temperature_unit=temperature_unit,
            )
            resolved_targets = _resolve_targets(
                dataset, EOSFitDomain(domain.lower()), targets
            )
            solver_options = _build_solver_options(
                solver,
                max_iterations=max_iterations,
                ftol=ftol,
                xtol=xtol,
                gtol=gtol,
                covariance_scaling=covariance_scaling,
                inner_max_iterations=inner_max_iterations,
                odr_difference=odr_difference,
                odr_ndigit=odr_ndigit,
            )
            constraint_specs = _parse_constraint_specs(
                fixed_parameters,
                initial_parameters,
                parameter_bounds,
            )
            jobs = tuple(
                EOSBatchJob(
                    request=_build_request(
                        target,
                        EOSFitDomain(domain.lower()),
                        pv_eos=pv_eos,
                        pv_order=pv_order,
                        vt_eos=vt_eos,
                        vt_variant=vt_variant,
                        pvt_coupling=pvt_coupling,
                        thermal_pressure_model=thermal_pressure_model,
                        atoms_per_cell=atoms_per_cell,
                        formula=formula,
                        formula_units_per_cell=formula_units_per_cell,
                        solver_options=solver_options,
                        constraint_specs=constraint_specs,
                    ),
                    job_id=f"{domain.lower()}-{target}",
                )
                for target in resolved_targets
            )
            plan = EOSBatchPlan(
                jobs=jobs,
                failure_policy=EOSBatchFailurePolicy(failure_policy.lower()),
                metadata={"source": "cli", "command": "quantas eos run"},
            )
            report_options = EOSReportOptions(
                detail=(
                    EOSReportDetail.EXTENDED
                    if report_verbosity.includes_extended
                    else EOSReportDetail.SHORT
                ),
                show_uncertainties=show_uncertainties,
                max_data_rows=max_data_rows,
                debug=debug,
            )
            creator = "quantas eos run"

        _configure_plan_debug(plan, debug)
        _validate_resolved_plan(dataset, plan)
        output_service.tables(
            build_batch_preamble(dataset, plan, destination, report_options)
        )
        if dry_run:
            output_service.message(
                "EOS dry run completed: the resolved plan is valid; "
                "no fit was executed and no HDF5 archive was created."
            )
            output_service.text_block(
                render_citation_notice(module_citation_keys("eos"))
            )
            output_service.message(quantas_finish())
            output_service.save()
            return

        observer = _EOSCLIEventObserver(output_service)
        result = run_eos_batch(
            dataset,
            plan,
            destination,
            observer=observer,
            overwrite=force,
            creator=creator,
        )
        output_service.tables(build_batch_report(result, report_options))
        output_service.text_block(render_citation_notice(module_citation_keys("eos")))
        output_service.message(quantas_finish())
        output_service.save()
        if not result.successful:
            message = "one or more EOS batch jobs failed"
            if not debug:
                message += "; rerun with -v debug for detailed solver diagnostics"
            raise click.ClickException(message)
    except click.ClickException:
        output_service.save()
        raise
    except Exception as exc:
        output_service.message(
            f"{type(exc).__name__}: {exc}",
            level=EventLevel.ERROR,
        )
        if show_traceback:
            output_service.message(traceback.format_exc(), level=EventLevel.ERROR)
        else:
            output_service.message(
                "Rerun with --traceback for the Python traceback or with "
                "-v debug for detailed numerical diagnostics.",
                level=EventLevel.INFO,
            )
        output_service.save()
        raise click.ClickException(str(exc)) from exc
    finally:
        output_service.close()


__all__ = ["run"]
