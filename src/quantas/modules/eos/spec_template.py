# -*- coding: utf-8 -*-

"""Generation of the commented EOS specification template.

The template is intentionally maintained as application-neutral text.  The
command-line interface only chooses a destination and delegates file creation
to :func:`write_eos_spec_template`.
"""

from __future__ import annotations

from pathlib import Path

EOS_SPEC_TEMPLATE_FILENAME = "eos.spec"
"""Default filename used by ``quantas eos spec-template``."""


_EOS_SPEC_TEMPLATE = """# QUANTAS EOS SPEC 1
#
# Quantas EOS batch specification template
# =========================================
#
# This file describes how Quantas should analyse a separate EOS data file.
# The filename extension is irrelevant; recognition is based on the exact
# signature above.  Comments begin with '#'.  Unknown sections, unknown keys,
# duplicated keys, empty values, and incompatible options are errors.
#
# Run the active example at the end with:
#
#   quantas eos run DATAFILE --spec eos.spec --dry-run
#   quantas eos run DATAFILE --spec eos.spec --output results.hdf5
#
# The dry run reads the data, resolves defaults, expands targets, and validates
# every fit without executing a solver or creating an HDF5 archive.
#
# Precedence for each job:
#
#   Quantas internal defaults
#       -> [defaults]
#       -> [defaults.pv], [defaults.vt], or [defaults.pvt]
#       -> [job NAME]
#
# Only the final [job volume-example] section is active by default.  The other
# examples are fully commented and can be copied or uncommented as required.


[metadata]
# Free-form provenance.  Metadata keys are not interpreted scientifically.
title = EOS batch analysis
# description = Volume and axial equation-of-state fits
# author = User name
# project = Project identifier


[input]
# Optional input-unit overrides.  When omitted, units declared in the data file
# are used; otherwise Quantas defaults to GPa, angstrom/angstrom^3, and kelvin.
# Internal calculations are normalized to GPa, angstrom, angstrom^3, and K.
#
# Supported temperature scales are K, C, and F.  Pressure and length units are
# interpreted by the shared Quantas unit-conversion infrastructure.
# pressure_unit = GPa
# length_unit = angstrom
# temperature_unit = K


[batch]
# Behaviour after a fit failure: stop or continue.
failure_policy = stop


[defaults]
# Common settings inherited by every job unless overridden later.
#
# Solver choices:
#   ols                 ordinary least squares; input uncertainties ignored
#   wls                 weighted least squares; response uncertainties required
#   effective-variance  includes uncertainties of explanatory coordinates
#   odr                 orthogonal distance regression; odrpack runtime dependency required
solver = ols
#
# Covariance scaling:
#   absolute
#   reduced-chi-square
#   inflate-only        EosFit-like: never reduce reported covariance
# covariance_scaling = inflate-only
#
# Generic positive solver controls.  Omit them to use typed solver defaults.
# max_iterations = 200
# ftol = 1.0e-10
# xtol = 1.0e-10
# gtol = 1.0e-10
#
# Effective-variance-only option.
# inner_max_iterations = 100
#
# ODR-only options.
# odr_difference = central
# odr_ndigit = 12
#
# Successful jobs are accepted by default.  Replacing a previously accepted
# result for the same domain/target must be requested explicitly.
accept = yes
replace_accepted = no
#
# Non-destructive data selection. The default base respects USE=0 and trailing
# '*' markers in the EOS data file. Use 'all' only when a job intentionally
# reconsiders every archived observation. GROUP values and row numbers are
# positive integers; row numbering is one-based in data-table order.
# selection_base = default
# groups = all
# include_groups = 1, 2
# exclude_groups = 3
# include_rows = 1, 2, 5, 8
# exclude_rows = 17, 23
#
# Parameter declarations may also be inherited.  Use canonical parameter names
# used by the Python API.  Use V0 for volume and L0 for a cell-axis length.
# The historical X0 and generic value_ref names are intentionally rejected.
#
# fix.KP = 4.0
# initial.K0 = 160.0
# bound.K0 = 1.0 : 500.0
#
# A fixed parameter cannot also declare initial or bound.  A bound for a free
# parameter requires an explicit initial.PARAMETER value.


[defaults.pv]
# Default isothermal P-V model for P-V jobs.
#
# Compact canonical models:
#   M
#   BM2, BM3, BM4
#   PT2, PT3, PT4       natural strain / Poirier-Tarantola
#   V2, V3
#   T2, T3, T4
#
# Full aliases such as birch-murnaghan, natural-strain, vinet, and tait are
# accepted; an unqualified ordered family uses its documented default order.
model = BM3
# solver = effective-variance
# covariance_scaling = inflate-only
# max_iterations = 200
# ftol = 1.0e-10
# xtol = 1.0e-10
# gtol = 1.0e-10
# inner_max_iterations = 100
# odr_difference = central
# odr_ndigit = 12
# accept = yes
# replace_accepted = no
# fix.KP = 4.0
# initial.K0 = 160.0
# bound.K0 = 1.0 : 500.0


[defaults.vt]
# Default thermal V-T model for V-T jobs.
#
# Supported model tags:
#   berman:linear
#   berman:quadratic
#   fei:linear
#   fei:inverse-square
#   mhp:simplified
#   mhp:general
#   salje
#   khp
model = berman:quadratic
# solver = effective-variance
# covariance_scaling = inflate-only
# max_iterations = 200
# ftol = 1.0e-10
# xtol = 1.0e-10
# gtol = 1.0e-10
# inner_max_iterations = 100
# odr_difference = central
# odr_ndigit = 12
# accept = yes
# replace_accepted = no
#
# Examples of V-T parameters.  Availability depends on the selected model.
# fix.temperature_ref = 298.15
# For volumetric V-T fits use V0.  For axial V-T fits use L0 in the job.
# initial.V0 = 1.0
# initial.alpha0 = 3.0e-5
# initial.alpha1 = 1.0e-8
# initial.alpha2 = 1.0e-3
# initial.theta_e = 500.0
# initial.theta_sat = 300.0
# initial.p1 = 1.0e-5
# initial.alpha_ref = 3.0e-5
# initial.kp = 4.0
# bound.alpha0 = none : none


[defaults.pvt]
# Default coupled P-V-T composition.
#
# coupling choices:
#   linear
#   anderson-gruneisen
#   thermal-pressure
#
# linear and anderson-gruneisen combine pv_model with vt_model.
# thermal-pressure defines its own thermal term, so vt_model must be omitted in
# the individual job that selects that coupling. Its oscillator model is selected
# with thermal_pressure_model = holland-powell-einstein, mgd, or
# mgd:q-compromise. The plain mgd tag is the full theta_d0/gamma0/q model;
# mgd:q-compromise is the distinct EosFit-style approximation and has no q
# parameter. MGD cell volumes require atoms_per_cell or formula +
# formula_units_per_cell.
# pv_model = BM3
# vt_model = berman:quadratic
# coupling = linear
# solver = effective-variance
# covariance_scaling = inflate-only
# max_iterations = 200
# ftol = 1.0e-10
# xtol = 1.0e-10
# gtol = 1.0e-10
# inner_max_iterations = 100
# odr_difference = central
# odr_ndigit = 12
# accept = yes
# replace_accepted = no
#
# Examples of P-V, thermal, and coupling parameters.  Availability depends on
# the selected pressure model, thermal model, and coupling.
# initial.K0 = 160.0
# initial.KP = 4.0
# initial.KPP = -0.02
# initial.V0 = 100.0
# initial.temperature_ref = 298.15
# initial.alpha0 = 3.0e-5
# initial.alpha1 = 1.0e-8
# initial.alpha2 = 1.0e-3
# initial.alpha_ref = 3.0e-5
# initial.p1 = 1.0e-5
# initial.theta_sat = 300.0
# initial.dK0_dT = -0.02
# bound.dK0_dT = none : 0.0
# initial.delta = 4.0
# initial.theta_e = 500.0
#
# MGD-only examples:
# thermal_pressure_model = mgd
# volume_basis = cell
# formula = NaF
# formula_units_per_cell = 4
# fix.temperature_ref = 295.0
# initial.theta_d0 = 459.0
# bound.theta_d0 = 1.0 : 5000.0
# initial.gamma0 = 1.5
# bound.gamma0 = 0.01 : 10.0
# initial.q = 1.0
# bound.q = -10.0 : 10.0


[presentation]
# Presentation does not change scientific results or the HDF5 archive.
# detail: short or extended.
detail = short
# Print a second input table containing available standard uncertainties.
show_uncertainties = no
# Positive integer, all, or none.  The complete dataset is always archived.
max_data_rows = all


# Example of a grouped sensitivity job:
#
# [job campaign-2]
# domain = pv
# targets = volume
# model = BM3
# groups = 2
# exclude_rows = 17
# accept = no
#
# ---------------------------------------------------------------------------
# ACTIVE MINIMAL EXAMPLE
# ---------------------------------------------------------------------------
# Fits the volume with the P-V defaults above: BM3 and OLS.
[job volume-example]
domain = pv
targets = volume
# model = BM2
# solver = effective-variance
# note = Optional human-readable note persisted with the job
# fix.KP = 4.0
# initial.K0 = 160.0
# bound.K0 = 1.0 : 500.0
# accept = yes
# replace_accepted = no


# ---------------------------------------------------------------------------
# COMMENTED EXAMPLE: all available P-V targets with one configuration
# ---------------------------------------------------------------------------
# [job all-pv]
# domain = pv
# targets = all
# model = BM3
# solver = effective-variance
# covariance_scaling = inflate-only


# ---------------------------------------------------------------------------
# COMMENTED EXAMPLE: volume and axes with different models and solvers
# ---------------------------------------------------------------------------
# [job volume-bm2]
# domain = pv
# targets = volume
# model = BM2
# solver = effective-variance
# fix.KP = 4.0
#
# [job axis-a-natural-strain]
# domain = pv
# targets = a
# model = PT3
# solver = odr
# odr_difference = central
# odr_ndigit = 12
# initial.M0 = 450.0
# initial.MP = 10.0
# bound.M0 = 1.0 : 2000.0
# bound.MP = 0.0 : 30.0
#
# [job axes-bc]
# domain = pv
# targets = b, c
# model = BM3
# solver = ols


# ---------------------------------------------------------------------------
# COMMENTED EXAMPLE: thermal volume and axial fits
# ---------------------------------------------------------------------------
# [job thermal-volume]
# domain = vt
# targets = volume
# model = berman:quadratic
# solver = effective-variance
#
# [job thermal-axes]
# domain = vt
# targets = a, c
# model = fei:inverse-square
# solver = odr
# odr_difference = central


# ---------------------------------------------------------------------------
# COMMENTED EXAMPLE: coupled P-V-T fit
# ---------------------------------------------------------------------------
# [job global-pvt]
# domain = pvt
# targets = volume
# pv_model = BM3
# vt_model = berman:quadratic
# coupling = linear
# solver = effective-variance
# initial.dK0_dT = -0.02
# bound.dK0_dT = none : 0.0


# ---------------------------------------------------------------------------
# COMMENTED EXAMPLE: Anderson-Gruneisen coupling
# ---------------------------------------------------------------------------
# [job anderson-gruneisen-pvt]
# domain = pvt
# targets = volume
# pv_model = BM3
# vt_model = berman:quadratic
# coupling = anderson-gruneisen
# solver = effective-variance
# initial.delta = 4.0
# bound.delta = 0.0 : 20.0


# ---------------------------------------------------------------------------
# COMMENTED EXAMPLE: thermal-pressure coupling
# ---------------------------------------------------------------------------
# Do not declare vt_model for this coupling.
# [job thermal-pressure-pvt]
# domain = pvt
# targets = volume
# pv_model = BM3
# coupling = thermal-pressure
# solver = effective-variance
# initial.theta_e = 500.0
# bound.theta_e = 1.0 : 5000.0


# ---------------------------------------------------------------------------
# COMMENTED EXAMPLE: volume-only Mie-Gruneisen-Debye P-V-T fit
# ---------------------------------------------------------------------------
# No adiabatic bulk-modulus observations or multi-objective residuals are used.
# [job mgd-pvt]
# domain = pvt
# targets = volume
# pv_model = BM4
# coupling = thermal-pressure
# thermal_pressure_model = mgd
# volume_basis = cell
# formula = NaF
# formula_units_per_cell = 4
# solver = effective-variance
# fix.temperature_ref = 295.0
# initial.theta_d0 = 459.0
# bound.theta_d0 = 1.0 : 5000.0
# initial.gamma0 = 1.5
# bound.gamma0 = 0.01 : 10.0
# initial.q = 1.0
# bound.q = -10.0 : 10.0


# ---------------------------------------------------------------------------
# COMMENTED EXAMPLE: two accepted fits for one slot
# ---------------------------------------------------------------------------
# Both records are preserved.  The second becomes current only because
# replace_accepted = yes is explicit.
# [job volume-first]
# domain = pv
# targets = volume
# model = BM2
# accept = yes
#
# [job volume-replacement]
# domain = pv
# targets = volume
# model = BM3
# accept = yes
# replace_accepted = yes
"""


def eos_spec_template() -> str:
    """Return the complete commented EOS specification template.

    Returns
    -------
    str
        UTF-8 text beginning with the mandatory ``QUANTAS EOS SPEC 1``
        signature and ending with a newline.
    """
    return _EOS_SPEC_TEMPLATE


def write_eos_spec_template(
    path: str | Path = EOS_SPEC_TEMPLATE_FILENAME,
    *,
    overwrite: bool = False,
) -> Path:
    """Write the complete EOS specification template to disk.

    Parameters
    ----------
    path : str or Path, optional
        Destination path.  Its suffix has no semantic meaning.
    overwrite : bool, optional
        Replace an existing regular file when ``True``.

    Returns
    -------
    pathlib.Path
        Resolved destination object used for writing.

    Raises
    ------
    FileExistsError
        If the destination exists and ``overwrite`` is ``False``.
    IsADirectoryError
        If ``path`` identifies a directory.
    OSError
        If the destination cannot be written.
    """
    destination = Path(path)
    if destination.exists():
        if destination.is_dir():
            raise IsADirectoryError(str(destination))
        if not overwrite:
            raise FileExistsError(str(destination))
    destination.write_text(eos_spec_template(), encoding="utf-8", newline="\n")
    return destination


__all__ = [
    "EOS_SPEC_TEMPLATE_FILENAME",
    "eos_spec_template",
    "write_eos_spec_template",
]
