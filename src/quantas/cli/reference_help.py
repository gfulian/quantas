# -*- coding: utf-8 -*-

"""Canonical long-form help used by the public Click reference.

The scientific commands keep concise descriptions next to their decorators so
that their implementation remains readable.  This module expands those
summaries at the assembled command-tree boundary.  The resulting text is used
by both 'quantas ... --help' and 'sphinx-click'; the installed CLI and the
published command reference therefore cannot drift apart.

Only frontend metadata is modified.  No callback, default, type conversion,
validation rule, or scientific option value is changed here.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Final

import click


_COMMAND_HELP: Final[Mapping[tuple[str, ...], str]] = {
    ("quantas",): (
        "Quantas command-line interface.\n\n"
        "Choose a scientific command group and then a subcommand such as 'run', "
        "'plot', 'export', or 'inspect'.  The CLI is a frontend over the same "
        "public workflow APIs used from Python; presentation choices never change "
        "the underlying float64 scientific calculation.\n\n"
        "Use 'quantas COMMAND --help' and 'quantas COMMAND SUBCOMMAND --help' "
        "to inspect the version-matched options installed on the current system."
    ),
    ("elasticity",): (
        "Analyze a second-order elastic stiffness tensor.\n\n"
        "The group validates the tensor, computes compliance, mechanical stability, "
        "Voigt–Reuss–Hill averages, and directional elastic properties.  Optional "
        "commands generate Quantas input, persist sampled 2D/3D fields, export tables, "
        "and render figures from a native HDF5 result."
    ),
    ("elasticity", "run"): (
        "Run an elasticity calculation from a Quantas text input file.\n\n"
        "FILENAME contains the stiffness tensor, density metadata when available, and "
        "the tensor convention.  Scalar summaries and global directional extrema are "
        "always evaluated; '--2d' and '--3d' request sampled fields that can be "
        "stored and plotted.  Rotation options transform the tensor into an analysis "
        "frame before any property is calculated."
    ),
    ("elasticity", "export"): (
        "Export stored two-dimensional elasticity fields from a Quantas HDF5 result.\n\n"
        "FILENAME must be an elasticity archive produced with '--2d'.  The command "
        "writes deterministic numerical tables suitable for external plotting and "
        "does not recalculate the tensor analysis."
    ),
    ("elasticity", "inpgen"): (
        "Generate a Quantas elasticity input from an external electronic-structure output.\n\n"
        "FILENAME is read through the selected CRYSTAL or VASP interface.  The generated "
        "text file should be inspected for units, pressure convention, density, and frame "
        "before it is passed to 'quantas elasticity run'."
    ),
    ("elasticity", "plot"): (
        "Create two- and three-dimensional elasticity figures from a native HDF5 result.\n\n"
        "FILENAME supplies the persisted tensor analysis.  Stored fields are reused when "
        "possible; explicit angular options can rebuild a plotted 3D surface at another "
        "resolution without modifying the archive.  All styling options affect rendering "
        "only."
    ),
    ("seismic", "inpgen"): (
        "Generate a Quantas SEISMIC input from an external electronic-structure output.\n\n"
        "FILENAME is read through the selected CRYSTAL or VASP interface.  Generation "
        "requires a stiffness tensor and finite positive density metadata; inspect the "
        "resulting shared elastic input before passing it to 'quantas seismic run'."
    ),
    ("ha",): (
        "Evaluate harmonic vibrational thermodynamics on one or more fixed volumes.\n\n"
        "The group converts supported phonon outputs to the Quantas YAML contract, runs "
        "temperature-dependent harmonic sums, exports numerical properties, and renders "
        "figures from native HDF5 results.  HA does not minimize volume or introduce a "
        "pressure coordinate."
    ),
    ("ha", "inpgen"): (
        "Generate a Quantas HA/QHA YAML input from supported phonon output.\n\n"
        "FILENAME may be one output, a CRYSTAL QHA output, or a text list of volume-dependent "
        "outputs when '--list' is selected.  The generated YAML preserves supplied q-point "
        "weights, normalization, and mode-continuity metadata and should be reviewed before "
        "calculation."
    ),
    ("ha", "run"): (
        "Run a harmonic-approximation calculation from a Quantas YAML input file.\n\n"
        "FILENAME supplies static energies, volumes, phonon frequencies, q-point weights, "
        "and normalization metadata.  The command evaluates every requested temperature at "
        "every stored volume, writes a native HDF5 result and deterministic report, and can "
        "optionally render a compact plot set."
    ),
    ("ha", "export"): (
        "Export one harmonic property from a Quantas HA HDF5 result.\n\n"
        "FILENAME is read without rerunning HA.  Select the stored property and an optional "
        "display unit; '--ask-unit' is intended for interactive terminal use, whereas an "
        "explicit '--unit' is preferable in scripts."
    ),
    ("ha", "plot"): (
        "Generate static figures from a Quantas HA HDF5 result.\n\n"
        "FILENAME provides the stored temperature and volume grids.  Property and unit "
        "choices affect only the rendered representation; the numerical HDF5 arrays are "
        "never rounded or rewritten."
    ),
    ("qha",): (
        "Evaluate quasi-harmonic thermodynamics over pressure and temperature.\n\n"
        "The group provides input generation, a static energy–volume preflight, the full "
        "QHA workflow, pressure–temperature exports, and plotting.  Scientific choices such "
        "as frequency versus thermodynamic interpolation, polynomial versus EOS minimization, "
        "and the thermal-expansion route are explicit command options."
    ),
    ("qha", "inpgen"): (
        "Generate a Quantas HA/QHA YAML input from supported phonon output.\n\n"
        "FILENAME may be one output, a CRYSTAL QHA output, or a text list of volume-dependent "
        "outputs when '--list' is selected.  The generated YAML preserves supplied q-point "
        "weights, normalization, and mode-continuity metadata and should be reviewed before "
        "calculation."
    ),
    ("qha", "inspect"): (
        "Inspect the sampled static energy–volume relation before a full QHA run.\n\n"
        "FILENAME is fitted independently with a polynomial and an energy EOS unless either "
        "preview is disabled.  The command reports the implied pressure coverage and fit "
        "quality, helping identify unsupported pressure ranges before the more expensive "
        "phonon calculation is started."
    ),
    ("qha", "run"): (
        "Run a quasi-harmonic calculation from a Quantas YAML input file.\n\n"
        "FILENAME supplies the multi-volume static and vibrational dataset.  Quantas builds "
        "the selected volume representation, minimizes the Gibbs energy at every pressure-"
        "temperature state, reconstructs thermodynamic and structural properties, records "
        "fit diagnostics and fallback provenance, and writes a native HDF5 archive."
    ),
    ("qha", "plot"): (
        "Generate line plots or pressure–temperature maps from a QHA HDF5 result.\n\n"
        "FILENAME is not refitted.  Repeating '--property' selects several products, while "
        "'--2d' requests filled maps only when the stored pressure and temperature axes "
        "contain enough points.  Unit and style options affect presentation only."
    ),
    ("qha", "export"): (
        "Export QHA pressure–temperature tables from a native HDF5 result.\n\n"
        "FILENAME supplies the archived grids, uncertainties, and structural path.  Export a "
        "single property, the structural subset, or the complete pressure–temperature table "
        "without rerunning the calculation."
    ),
    ("seismic",): (
        "Calculate, persist, export, and plot anisotropic acoustic-wave properties.\n\n"
        "The workflow solves the Christoffel problem on a spherical grid and can progressively "
        "add group velocity, power-flow, enhancement, caustic, and polarization-tracking data. "
        "The HDF5 archive is the common source for reports, CSV exports, and figures."
    ),
    ("seismic", "run"): (
        "Run the sampled SEISMIC workflow and save a native HDF5 result.\n\n"
        "FILENAME contains a stiffness tensor and density.  The angular grid determines spatial "
        "resolution, '--level' determines the highest derivative-dependent property computed, "
        "and '--batch-size' controls memory and throughput only.  A mechanically unstable "
        "tensor is rejected before wave propagation is evaluated."
    ),
    ("seismic", "plot"): (
        "Render spherical maps, acoustic surfaces, or a six-panel summary from HDF5.\n\n"
        "FILENAME must contain the property level required by each requested plot.  Map, surface, "
        "polarization, projection, and typography options affect visualization only; they do not "
        "change the stored wave solution."
    ),
    ("seismic", "export"): (
        "Export sampled acoustic fields from a Quantas SEISMIC HDF5 result.\n\n"
        "The CSV contains one row per sampled direction and wave mode, including every phase, "
        "group, enhancement, degeneracy, and tracking field available in FILENAME.  No new "
        "Christoffel calculation is performed."
    ),
    ("eos",): (
        "Fit equations of state and manage persistent EOS archives.\n\n"
        "Unlike a single-pass scientific module, EOS commonly compares several models, solvers, "
        "constraints, and data selections.  Every immutable fit attempt can be retained in one "
        "archive, while accepted and candidate records identify the interpretations selected for "
        "diagnosis, calculation, and plotting."
    ),
    ("eos", "spec-template"): (
        "Generate a complete commented 'QUANTAS EOS SPEC 1' batch template.\n\n"
        "OUTPUT defaults to 'eos.spec'; its suffix has no semantic meaning.  The template "
        "contains one active minimal job and commented examples for supported domains, model "
        "families, solvers, constraints, selection rules, and acceptance behavior."
    ),
    ("eos", "run"): (
        "Run one non-interactive EOS batch from a data file and optional specification.\n\n"
        "FILENAME supplies observations, uncertainties, groups, exclusions, and unit metadata. "
        "Direct CLI options generate one or more homogeneous jobs; '--spec' supports a fully "
        "declarative heterogeneous plan.  Every attempted fit is persisted, and acceptance is "
        "applied only after a successful record has been written."
    ),
    ("eos", "diagnose"): (
        "Inspect residuals and model-specific diagnostics from an EOS archive.\n\n"
        "ARCHIVE is queried by accepted slot or immutable record ID.  The command reports selected "
        "and excluded observations, physical and standardized residuals, group provenance, and "
        "finite-strain diagnostics when the chosen formulation defines them; it never refits data."
    ),
    ("eos", "plot"): (
        "Generate EOS figures from an immutable HDF5 fit record.\n\n"
        "ARCHIVE is queried by accepted slot or explicit record ID.  Available plot types depend "
        "on the fitted domain and model.  Data-selection, grouping, error-bar, curve-sampling, and "
        "style controls change only the figure representation."
    ),
    ("eos", "calculate"): (
        "Evaluate a stored EOS model at selected states or regular grids.\n\n"
        "ARCHIVE is queried by accepted slot or immutable record ID.  Coordinates can be repeated, "
        "expanded from inclusive ranges, combined as Cartesian products, or paired element-wise. "
        "Optional covariance propagation estimates parameter-derived uncertainty without rerunning "
        "the fit."
    ),
    ("thermoelasticity",): (
        "Calibrate and evaluate quasi-static thermoelastic stiffness tensors.\n\n"
        "The group normalizes static elastic calculations, couples them to QHA, fits a reusable "
        "cold finite-strain model, evaluates points, grids, or geological profiles, exports tensors, "
        "and renders calibration and pressure–temperature diagnostics."
    ),
    ("thermoelasticity", "inpgen"): (
        "Generate a readable thermoelastic YAML input from CRYSTAL SOEC outputs.\n\n"
        "FILENAME is one CRYSTAL output or, with '--list', a text file containing a volume/pressure "
        "series.  Quantas validates structural identity, co-rotates compatible tensors into a common "
        "frame, and records provenance; the generated YAML should be inspected before calibration."
    ),
    ("thermoelasticity", "profile-template"): (
        "Write an editable geothermobarometric profile YAML template.\n\n"
        "The generated file demonstrates layered pressure and temperature models, joining policies, "
        "and an explicit depth grid.  It is a template rather than a universal Earth model and must "
        "be adapted to the geological problem."
    ),
    ("thermoelasticity", "run"): (
        "Calibrate the static EOS and independent cold finite-strain 'Cij(V)' models.\n\n"
        "THERMOELASTIC_INPUT supplies the normalized elastic series and QHA_SOURCE supplies either "
        "a validated QHA HDF5 result or a QHA YAML input to evaluate first.  The command writes a "
        "reusable fit archive; it does not automatically construct a pressure–temperature grid."
    ),
    ("thermoelasticity", "analysis"): (
        "Evaluate a calibrated thermoelastic model at points, grids, or geological profiles.\n\n"
        "Choose the smallest subcommand matching the scientific task.  Point evaluation is useful "
        "for validation and interoperability, grid evaluation builds rectangular P–T archives, and "
        "profile evaluation follows only the requested depth path."
    ),
    ("thermoelasticity", "analysis", "point"): (
        "Evaluate one pressure–temperature state from a calibrated fit archive.\n\n"
        "ARCHIVE, PRESSURE, and TEMPERATURE identify the state.  The command reports volume, density, "
        "stiffness, uncertainty, stability, and support flags, and can write a shared Elasticity/"
        "SEISMIC text input for the selected tensor condition."
    ),
    ("thermoelasticity", "analysis", "grid"): (
        "Evaluate a rectangular pressure–temperature grid and save an analysis archive.\n\n"
        "ARCHIVE supplies the calibration.  Inclusive pressure and temperature ranges define all "
        "states; optional wide tables expose independent components, uncertainties, stability, and "
        "extrapolation flags for external analysis."
    ),
    ("thermoelasticity", "analysis", "profile"): (
        "Evaluate one geological pressure–temperature-depth path.\n\n"
        "ARCHIVE supplies the calibration, while a tabulated profile, profile specification, preset, "
        "or controlled linear test profile supplies the states.  Only points on the selected path "
        "are evaluated and persisted; no enclosing rectangular grid is created."
    ),
    ("thermoelasticity", "table"): (
        "Print or export compact thermoelastic tables containing independent components.\n\n"
        "The table commands read existing fit, grid, or profile archives and never recalibrate the "
        "model.  Use them for concise inspection or portable row-oriented output."
    ),
    ("thermoelasticity", "table", "point"): (
        "Print one pressure–temperature state as a compact independent-component table.\n\n"
        "ARCHIVE, PRESSURE, and TEMPERATURE select the state.  Tensor condition, support policy, and "
        "uncertainty choices mirror the corresponding point analysis without creating a new HDF5 file."
    ),
    ("thermoelasticity", "table", "grid"): (
        "Export an archived pressure–temperature grid in wide row-oriented form.\n\n"
        "ARCHIVE must contain a grid analysis.  Select tensor conditions and uncertainty columns, "
        "then choose deterministic text or CSV output without recomputing any state."
    ),
    ("thermoelasticity", "table", "profile"): (
        "Export one archived geological profile in wide row-oriented form.\n\n"
        "ARCHIVE must contain at least one profile.  Select the profile explicitly when several are "
        "stored; the table preserves depth, pressure, temperature, support flags, and requested tensors."
    ),
    ("thermoelasticity", "export"): (
        "Inspect coverage and export an explicit thermoelastic tensor selection.\n\n"
        "ARCHIVE may be a fit, grid, or profile result.  Choose one point, a generated grid, one "
        "archived profile, or the complete available collection.  Point export can produce the shared "
        "Elasticity/SEISMIC input format."
    ),
    ("thermoelasticity", "inspect"): (
        "Explain the stage, scientific coverage, and recommended next steps of an archive.\n\n"
        "ARCHIVE is read without evaluation or modification.  The summary distinguishes calibration, "
        "grid, and profile content; reports QHA and elastic support; identifies tensor conditions; and "
        "prints commands appropriate to the stored stage."
    ),
    ("thermoelasticity", "plot"): (
        "Generate thermoelastic calibration and scientific figures.\n\n"
        "Invoke a plot family to visualize the fitted cold model, a P–T field, a geological profile, "
        "the validity domain, or the difference between isothermal and adiabatic tensors.  Use "
        "'--list-plots' to display the available families."
    ),
    ("thermoelasticity", "plot", "compare"): (
        "Compare isothermal and adiabatic stiffness along a one-dimensional P–T slice.\n\n"
        "ARCHIVE must contain a suitable analysis grid.  Fix either pressure or temperature, choose "
        "components, and inspect the magnitude and temperature dependence of the adiabatic correction."
    ),
    ("thermoelasticity", "plot", "fit"): (
        "Plot observed and fitted independent elastic components as functions of volume.\n\n"
        "ARCHIVE must contain the calibration stage.  Curves, residuals, observed points, confidence "
        "bands, and component layout expose fit support and should be inspected before P–T analysis."
    ),
    ("thermoelasticity", "plot", "pt"): (
        "Plot stiffness or propagated uncertainty on an archived pressure–temperature grid.\n\n"
        "ARCHIVE supplies the evaluated grid.  Select components, tensor condition, scalar field, "
        "contour treatment, and extrapolation overlays; these controls do not modify the archive."
    ),
    ("thermoelasticity", "plot", "profile"): (
        "Plot elastic components along one archived geological profile.\n\n"
        "ARCHIVE supplies depth, pressure, temperature, tensors, and support masks.  The command can "
        "show absolute values or changes relative to a reference depth, with optional uncertainty, "
        "pressure axes, temperature encoding, and extrapolation intervals."
    ),
    ("thermoelasticity", "plot", "domain"): (
        "Plot the QHA coordinate domain and elastic-volume support with geological paths.\n\n"
        "ARCHIVE supplies the equilibrium-volume field and any stored profiles.  Extrapolation overlays "
        "distinguish states outside the QHA coordinate grid from states whose equilibrium volume lies "
        "outside the calibrated elastic interval."
    ),
}


_SHARED_OPTION_HELP: Final[Mapping[str, str]] = {
    "version": (
        "Show the installed Quantas version banner and exit before command dispatch. "
        "Use this value when reporting results or reproducing a workflow."
    ),
    "report": (
        "Write a deterministic plain-text scientific report.  When omitted, run commands use the "
        "primary input name with a '.log' suffix; reports contain no ANSI styling or live progress."
    ),
    "verbosity": (
        "Select report detail: 'standard' gives the normal scientific summary, 'extended' adds "
        "complete tables, and 'debug' adds numerical diagnostics.  The choice never changes stored "
        "numerical results."
    ),
    "quiet": (
        "Suppress terminal presentation while preserving calculation, HDF5 output, and the report. "
        "This is useful for scripts, redirected jobs, and automated validation."
    ),
    "progress": (
        "Enable or disable transient Rich progress on an interactive terminal.  Progress events are "
        "operational only and are never written to reports or HDF5 histories."
    ),
    "force": (
        "Replace existing generated files without an interactive confirmation.  Use this deliberately "
        "in reproducible scripts because overwritten results cannot be recovered by Quantas."
    ),
    "figure_preset": (
        "Choose a named figure style.  'screen' favors interactive viewing, 'publication' uses "
        "print-oriented geometry and resolution, and 'monochrome' adds grayscale-safe rendering; "
        "the preset affects presentation only."
    ),
    "dpi": (
        "Override the raster resolution supplied by the selected figure preset.  This affects PNG or "
        "other raster output only and does not change sampled scientific data."
    ),
    "image_format": (
        "Choose the static figure file format supported by the Matplotlib renderer.  The format changes "
        "serialization only; vector formats are preferable when downstream editing is required."
    ),
    "show": (
        "Display figures interactively after files are generated.  Leave this disabled on headless, "
        "batch, or continuous-integration systems."
    ),
    "show_title": (
        "Show or hide figure titles.  This is a presentation choice; publication figures often omit "
        "titles in favor of external captions."
    ),
    "title": (
        "Override the figure-preset title setting.  This affects only annotations in the rendered "
        "figure and never the numerical result."
    ),
    "grid": (
        "Override the figure-preset Cartesian-grid setting.  Grid lines are graphical guides and do "
        "not represent additional calculated samples."
    ),
    "transparent": (
        "Save figures with a transparent or opaque background.  Transparency affects the image canvas "
        "only and is useful when composing figures in external software."
    ),
}


_CONTEXT_OPTION_HELP: Final[Mapping[tuple[tuple[str, ...], str], str]] = {
    (("elasticity", "run"), "ntheta"): (
        "Set angular samples along the polar coordinate.  The default is 361 points for each 2D closed "
        "curve and 61 points for a 3D surface; increase this only after checking directional convergence."
    ),
    (("elasticity", "run"), "nphi"): (
        "Set azimuthal samples for stored 3D surfaces.  The default '2 * ntheta - 1' gives comparable "
        "angular spacing without duplicating the periodic endpoint."
    ),
    (("elasticity", "plot"), "ntheta"): (
        "Rebuild plotted 3D surfaces with this many polar samples instead of the stored/default value. "
        "The override changes the plotted sampling but does not rewrite FILENAME."
    ),
    (("elasticity", "plot"), "nphi"): (
        "Rebuild plotted 3D surfaces with this many azimuthal samples.  Use it together with '--ntheta' "
        "for a controlled visual convergence check; the HDF5 archive is not modified."
    ),
    (("seismic", "run"), "batch_size"): (
        "Set the number of directions solved in each vectorized numerical batch.  The default 512 balances "
        "temporary memory and overhead; changing it affects performance and progress granularity, not the "
        "angular grid or numerical accuracy."
    ),
    (("seismic", "run"), "ntheta"): (
        "Set the polar sampling count for the selected hemisphere.  The default 91 gives about one-degree "
        "spacing on the upper hemisphere; increase it to converge narrow extrema or acoustic-axis candidates."
    ),
    (("seismic", "run"), "nphi"): (
        "Set the azimuthal sampling count.  The default 181 avoids a duplicate periodic endpoint while "
        "matching the angular scale of the default upper-hemisphere polar grid."
    ),
    (("qha", "run"), "poly_grid_points"): (
        "Set the odd number of local volumes used to estimate 'K_T' and 'K'_T' after polynomial "
        "minimization.  Five is the default; compare 3, 5, and 7 points before changing production results."
    ),
    (("qha", "run"), "poly_grid_separation"): (
        "Set adjacent local derivative-grid spacing as a percentage of equilibrium volume.  The default "
        "0.05 percent balances cancellation and truncation error; test sensitivity together with the point count."
    ),
    (("eos", "plot"), "curve_points"): (
        "Set the number of samples used to draw each smooth calculated curve.  The default 300 affects only "
        "visual smoothness and never changes fitted parameters, residuals, or covariance."
    ),
    (("eos", "calculate"), "relative_step"): (
        "Set the relative finite-difference step used when propagating parameter covariance to calculated "
        "quantities.  Change it only after checking that the reported uncertainty is stable to smaller and "
        "larger steps."
    ),
    (("thermoelasticity", "plot", "fit"), "curve_points"): (
        "Set samples used to draw each smooth cold finite-strain fit curve.  This controls graphical smoothness "
        "only and does not refit the component or alter its covariance."
    ),
    (("elasticity", "run"), "rotate_xyz"): (
        "Rotate tensor components into a new analysis frame using right-handed fixed-axis angles about source "
        "X, then Y, then Z.  The material is unchanged; only components and reported directions are transformed."
    ),
    (("elasticity", "run"), "rotation_matrix"): (
        "Read a proper 3 x 3 source-to-analysis component-transformation matrix from a text file.  The rotation "
        "changes the reporting frame, not invariant elastic properties or the physical material."
    ),
    (("seismic", "run"), "rotate_xyz"): (
        "Rotate stiffness components into a new analysis frame using right-handed fixed-axis angles about source "
        "X, then Y, then Z.  Wave directions are reported in the rotated frame while invariant speeds are preserved."
    ),
    (("seismic", "run"), "rotation_matrix"): (
        "Read a proper 3 x 3 source-to-analysis component-transformation matrix.  This changes the coordinate "
        "frame used for directions and polarizations, not the physical crystal."
    ),
    (("ha", "run"), "benchmark"): (
        "Render backend timing events in addition to the scientific report.  Benchmarking is diagnostic only and "
        "does not change harmonic sums, stored arrays, or numerical precision."
    ),
    (("qha", "run"), "gruneisen"): (
        "Calculate and persist the macroscopic thermodynamic Gruneisen parameter when its inputs are well defined. "
        "Disabling it omits this derived output but does not change equilibrium-volume minimization."
    ),
    (("qha", "run"), "mode_gruneisen"): (
        "Calculate mode-resolved Gruneisen parameters for the frequency interpolation scheme.  Disable this output "
        "to reduce work when modal interpretation is unnecessary; it is unavailable for the thermodynamic scheme."
    ),
    (("eos", "run"), "spec_path"): (
        "Read a strict 'QUANTAS EOS SPEC 1' batch plan.  The file can define heterogeneous jobs, selections, "
        "constraints, acceptance, and failure policy; its filename suffix has no semantic meaning."
    ),
    (("eos", "run"), "fixed_parameters"): (
        "Fix one physical parameter at the supplied value and remove it from the free optimization vector.  Repeat "
        "as needed; a target prefix limits scope.  A fixed value is an external assumption, not a fitted result."
    ),
    (("eos", "run"), "initial_parameters"): (
        "Override the starting value of a free parameter without fixing it.  This can improve convergence but should "
        "not change a well-posed optimum; compare final records when sensitivity is suspected."
    ),
    (("eos", "run"), "parameter_bounds"): (
        "Override the admissible interval of a free parameter; use 'none' for an open side.  A solution on a bound "
        "is a diagnostic condition and must not be interpreted automatically as a measured parameter."
    ),
    (("eos", "run"), "dry_run"): (
        "Resolve units, defaults, jobs, models, constraints, and acceptance rules, then validate the complete plan "
        "without fitting or creating an HDF5 archive.  Use this before long or heterogeneous batches."
    ),
    (("eos", "run"), "failure_policy"): (
        "Choose whether a failed job stops the batch immediately or allows later independent jobs to run.  Continued "
        "execution does not make the overall batch successful and does not reuse a failed job automatically."
    ),
    (("eos", "run"), "show_traceback"): (
        "Include or suppress a Python traceback in the report for unexpected exceptions.  This changes diagnostic "
        "detail only and has no effect on expected Click validation errors or fit behavior."
    ),
    (("eos", "diagnose"), "no_normalized_pressure"): (
        "Omit finite-strain and normalized-pressure columns from the diagnostic table.  Physical and standardized "
        "residuals remain available, and the archived fit record is not modified."
    ),
    (("eos", "calculate"), "no_uncertainty"): (
        "Evaluate central fitted values without propagating parameter covariance.  This can reduce work for dense "
        "grids but does not imply that the calculated state is uncertainty-free."
    ),
    (("thermoelasticity", "plot", "fit"), "residuals"): (
        "Include or omit a residual panel beneath each observed-versus-fitted component.  This is a diagnostic "
        "presentation choice and does not alter the calibration record."
    ),
    (("thermoelasticity", "plot", "fit"), "confidence_band"): (
        "Display or hide the propagated parameter-confidence band around each fitted curve.  The band visualizes "
        "stored covariance and is not a refit or a total model-discrepancy estimate."
    ),
    (("thermoelasticity", "plot", "fit"), "confidence"): (
        "Set the central probability represented by the confidence band.  This changes only the plotted interval "
        "derived from stored covariance; it does not change fitted parameters."
    ),
    (("thermoelasticity", "plot", "fit"), "symmetry_spread"): (
        "Show or hide the numerical spread among symmetry-equivalent observed entries.  The overlay is diagnostic "
        "and does not alter the symmetry-reduced values used by the fit."
    ),
}


def apply_reference_help(
    command: click.Command,
    root: tuple[str, ...],
    *,
    recursive: bool = True,
) -> None:
    """Apply canonical long-form help to one assembled Click command tree.

    Parameters
    ----------
    command : click.Command
        Root command or command group to document.
    root : tuple of str
        Public path associated with 'command'.  Top-level scientific groups
        use e.g. '("qha",)'; the application root uses '("quantas",)'.
    recursive : bool, optional
        Whether to document registered descendants.
    """

    _document_command(command, root)
    if not recursive or not isinstance(command, click.Group):
        return
    for name, child in command.commands.items():
        apply_reference_help(child, (*root, name), recursive=True)


def _document_command(command: click.Command, path: tuple[str, ...]) -> None:
    """Apply command and option text without changing Click behavior."""
    if not hasattr(command, "_quantas_base_help"):
        setattr(command, "_quantas_base_help", command.help or "")
    command.help = _COMMAND_HELP.get(
        path,
        str(getattr(command, "_quantas_base_help", command.help or "")).strip(),
    )
    command.short_help = command.help.split("\n", maxsplit=1)[0].strip()

    for parameter in command.params:
        if not isinstance(parameter, click.Option):
            continue
        if not hasattr(parameter, "_quantas_base_help"):
            setattr(parameter, "_quantas_base_help", parameter.help or "")
        base = str(getattr(parameter, "_quantas_base_help", parameter.help or "")).strip()
        explicit = _CONTEXT_OPTION_HELP.get((path, parameter.name))
        shared = _SHARED_OPTION_HELP.get(parameter.name)
        parameter.help = explicit or shared or _expand_option_help(
            base,
            parameter=parameter,
            command_path=path,
        )


def _expand_option_help(
    base: str,
    *,
    parameter: click.Option,
    command_path: tuple[str, ...],
) -> str:
    """Return a safe second explanatory sentence for concise option help."""
    text = base.rstrip()
    if not text:
        text = "Configure this command option."
    if len(text) >= 115 and _has_multiple_sentences(text):
        return text

    group = str(getattr(parameter, "help_group", "")).strip().lower()
    command_name = command_path[-1] if command_path else ""

    if parameter.name.startswith("list_"):
        suffix = (
            " This is a discovery option: it prints the available choices and exits "
            "without calculation or file creation."
        )
    elif command_name in {"diagnose", "calculate"} and parameter.name in {
        "slot",
        "record_id",
        "pressure",
        "pressure_range",
        "coordinate",
        "coordinate_range",
        "temperature",
        "temperature_range",
        "pairwise",
    }:
        suffix = (
            " This selects an archived record or evaluation coordinates and does not "
            "refit or modify the EOS archive."
        )
    elif _is_plotting_option(group, command_path):
        suffix = (
            " This affects the generated figure only and does not modify the stored "
            "scientific result."
        )
    elif "unit" in group or _looks_like_unit_option(parameter):
        suffix = (
            " Input values are converted to Quantas internal canonical units before "
            "float64 calculation and persistence."
        )
    elif any(token in group for token in ("output", "report", "persistence")):
        suffix = (
            " This controls persistence or presentation and does not change the "
            "scientific model."
        )
    elif any(token in group for token in ("domain", "selection", "coordinate", "profile")):
        suffix = (
            " Changing it changes the states, data, or directions evaluated and may "
            "therefore affect runtime and the stored result."
        )
    elif any(token in group for token in ("model", "coupling", "normalization", "adiabatic")):
        suffix = (
            " This changes the physical model or interpretation and should be compared "
            "through the relevant scientific diagnostics."
        )
    elif any(token in group for token in ("numerical", "solver", "constraint", "odr", "variance")):
        suffix = (
            " Treat it as a numerical-method control and check sensitivity or convergence "
            "before changing the documented default."
        )
    elif any(token in group for token in ("validation", "diagnostic", "stability", "support", "batch")):
        suffix = (
            " This controls validation, diagnostics, or failure handling and should not "
            "be used to conceal unsupported input data."
        )
    elif "advanced" in group:
        suffix = (
            " Change this advanced control only after a targeted convergence, stability, "
            "or memory-use test."
        )
    elif command_name in {"plot"}:
        suffix = (
            " This is a rendering control and does not alter the archived numerical arrays."
        )
    elif "table" in command_path or command_name in {"export", "table", "diagnose", "inspect"}:
        suffix = (
            " The command reads an existing result and does not rerun the scientific workflow."
        )
    elif command_name in {"inpgen", "profile-template", "spec-template"}:
        suffix = (
            " Review the generated file before using it as scientific input."
        )
    else:
        suffix = (
            " Consult the corresponding implementation-and-workflow page before changing "
            "this setting in a production calculation."
        )

    if not text.endswith((".", "!", "?")):
        text += "."
    return text + suffix


def _has_multiple_sentences(text: str) -> bool:
    """Return whether help already contains more than one explanatory sentence."""
    return sum(text.count(mark) for mark in ".!?") >= 2


def _is_plotting_option(group: str, command_path: tuple[str, ...]) -> bool:
    """Return whether an option belongs exclusively to figure presentation."""
    return "plot" in command_path or any(
        token in group
        for token in (
            "plot",
            "figure",
            "style",
            "line",
            "marker",
            "axes",
            "contour",
            "presentation",
            "annotation",
            "geometry",
            "typography",
        )
    )


def _looks_like_unit_option(parameter: click.Option) -> bool:
    """Return whether the public option name denotes a measurement unit."""
    names = {
        value.lstrip("-").replace("-", "_")
        for value in (*parameter.opts, *parameter.secondary_opts)
    }
    return any(
        name.endswith("unit")
        or name in {"eunit", "vunit", "funit", "tunit", "punit", "lunit"}
        for name in names
    )


__all__ = ["apply_reference_help"]
