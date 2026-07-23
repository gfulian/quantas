Command-line conventions
========================

This page explains notation and behavior shared by the generated command
reference.  The exact option names, accepted values, and defaults shown in the
module pages are imported from the installed Click command tree during the
Sphinx build.

Reading a command signature
---------------------------

A typical signature is written as:

.. code-block:: text

   quantas MODULE COMMAND [OPTIONS] ARGUMENT

Upper-case names such as ``FILENAME`` and ``ARCHIVE`` are positional
arguments.  Square brackets indicate optional syntax; they are not typed
literally.  An option followed by several metavariables consumes all of them:

.. code-block:: console

   --temperature MIN MAX STEP

Ranges are inclusive.  Quantas checks that the step has the correct sign and
that the requested axis can be constructed.  Options described as repeatable
must be written more than once rather than supplied as a comma-separated
string:

.. code-block:: console

   --property KT --property alphaV

Boolean pairs such as ``--progress/--no-progress`` state both the positive and
negative form.  The generated reference marks the active default.

Option groups
-------------

Option-rich commands organize their help under stable semantic headings:

``Scientific model``
   Selects a physical formulation or approximation.  Changing one of these
   options can change the scientific interpretation and normally deserves a
   sensitivity comparison.

``Calculation domain``
   Selects temperatures, pressures, volumes, directions, components, or data
   subsets.  These options often control both runtime and result size.

``Numerical method``
   Selects interpolation, minimization, derivatives, regression, sampling, or
   convergence controls.  Defaults should be changed only with an explicit
   convergence or robustness check.

``Validation and diagnostics``
   Controls support checks, tolerances, failure policies, and diagnostic
   reporting.  Such options do not repair physically inadequate data.

``Units``
   Tells the reader how to interpret input values or how to display an output.
   Native calculations and HDF5 floating-point arrays use ``float64`` in
   canonical internal units.

``Plotting``
   Controls figures and visual sampling.  Plot styling never changes persisted
   scientific arrays.  A plot command may rebuild a transient surface at a
   different angular resolution, but it does not rewrite the source archive.

``Output and reporting``
   Controls paths, plain-text reports, terminal presentation, and overwrite
   policy.  These options do not alter the scientific model.

Paths and overwrite policy
--------------------------

Most ``run`` commands derive HDF5 and report names from the primary input when
``--output`` or ``--report`` is omitted.  The generated help states the exact
rule for each command.

Quantas does not silently replace an existing scientific result.  Depending on
the command, an existing path either triggers a confirmation prompt or an
error.  ``--force`` performs the replacement without prompting and is intended
for controlled scripts where that behavior is deliberate.

Parent directories should normally exist before command execution.  Output
suffixes are normalized where the format has a required extension, such as
``.hdf5`` or ``.csv``.

Reports, terminal output, and progress
--------------------------------------

Scientific ``run`` commands distinguish three output channels:

* the native HDF5 result, containing numerical arrays and provenance;
* a deterministic plain-text report;
* transient terminal presentation.

``--quiet`` suppresses terminal output while preserving the HDF5 result and
report.  ``--no-progress`` disables the interactive progress bar but leaves
ordinary messages visible.  Progress events are never persisted in reports or
HDF5 histories.

The report verbosity levels are:

``standard``
   Normal scientific summary and principal diagnostics.

``extended``
   Additional complete tables where supported.

``debug``
   Numerical diagnostics, local fits, fallbacks, and state-level information.

Verbosity changes presentation and persisted meaningful messages, not numerical
precision or calculated arrays.

Units
-----

Unit options associated with a ``run`` command describe the source data unless
the option explicitly says *display unit*.  Readers normalize values before
calculation.  Plot and export units are presentation conversions applied to an
existing result.

Temperature ranges may accept kelvin or Celsius where documented, but internal
thermodynamic calculations use kelvin.  Pressure is normally expressed in GPa,
length in angstrom or bohr, and phonon frequency in reciprocal centimetres or a
supported frequency unit.

Plot presets
------------

Static plotting frontends share three presets:

``screen``
   Geometry and resolution suitable for interactive inspection.

``publication``
   Print-oriented dimensions, line weights, and raster resolution.

``monochrome``
   Grayscale-safe presentation in addition to publication-style geometry.

Explicit options such as ``--dpi``, font sizes, line widths, or transparency
override individual preset values.  These controls have no effect on the
calculation or HDF5 archive.

Errors and exit status
----------------------

Click rejects invalid syntax, missing paths, out-of-range option values, and
incompatible option combinations before numerical execution when possible.
Scientific validation failures and backend errors are reported as non-zero exit
statuses.  A command that writes a failed attempt to a persistent EOS archive
may still exit unsuccessfully; persistence of diagnostic evidence is not the
same as scientific success.

Use the shell exit status in scripts and inspect the report or archive before
using generated results downstream.

Related documentation
---------------------

* :doc:`../workflows/index` explains model choices, approximations, defaults,
  convergence, and failure policies.
* :doc:`../tutorials/index` provides complete reproducible command sequences.
* :doc:`../formats/index` defines input files, HDF5 envelopes, and portable
  exports.
* :doc:`../api/index` documents the equivalent public Python workflows.
