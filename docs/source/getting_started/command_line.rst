Command-line interface
======================

The Quantas command-line interface is a Click/Rich frontend over the same
scientific workflows exposed by :mod:`quantas.api`.  Click owns command
parsing, type conversion, validation, prompts, and dispatch.  Rich owns the
interactive terminal presentation.  Neither frontend changes the numerical
method.

Discover the command tree
-------------------------

Start with:

.. code-block:: console

   quantas --help
   quantas --version

The scientific groups are:

.. code-block:: console

   quantas ha --help
   quantas qha --help
   quantas elasticity --help
   quantas seismic --help
   quantas eos --help
   quantas thermoelasticity --help

Help is available at every level.  For example:

.. code-block:: console

   quantas qha run --help
   quantas thermoelasticity analysis profile --help

The generated :doc:`../cli/index` reference is version-matched to the installed
Click application.  It combines exact command signatures with maintained,
long-form option descriptions.

The common command pattern
--------------------------

Most modules separate calculation from later presentation or export:

.. code-block:: text

   input
   → run
   → native HDF5
   → inspect / plot / export

For example:

.. code-block:: console

   quantas elasticity run calcite.dat \
       --2d \
       --output calcite_elasticity.hdf5 \
       --report calcite_elasticity.log

   quantas elasticity plot calcite_elasticity.hdf5 \
       --2d \
       --property young \
       --output calcite_young

The calculation is not repeated by the plotting command.  The plot builder
reads the stored scientific result and creates frontend-neutral plot
specifications before the Matplotlib renderer writes concrete figures.

Command groups are scientifically different
--------------------------------------------

The broad pattern is shared, but each module exposes the operations its
science requires.

HA and QHA
   Run thermodynamic calculations, inspect or generate phonon input where
   appropriate, then export or plot stored temperature--pressure results.

Elasticity and SEISMIC
   Characterize an elastic tensor or acoustic field, persist the result, and
   later export directional data or render selected maps and surfaces.

Thermoelasticity
   Uses a staged workflow: input generation, cold finite-strain calibration,
   archive inspection, and Point/Grid/Profile analysis.  The calibration HDF5
   is not itself a P--T grid.

EOS
   Uses a fit/archive lifecycle.  One dataset may produce several immutable
   fit records under different models and solvers.  Accepted and last records
   are not necessarily the same.  This is why EOS post-fit commands differ
   from the ordinary ``run → plot/export`` sequence.

Common output controls
----------------------

The exact set of options depends on the command, but the following concepts are
shared.

``--output``
   Selects the principal native result or exported artifact.  Scientific run
   commands normally write HDF5; export commands may write text or CSV.

``--report``
   Writes deterministic plain text.  Reports contain no ANSI sequences, Rich
   box drawing, or transient progress artifacts.

``--force``
   Allows an existing destination to be replaced.  Without it, Quantas refuses
   silent data loss.

``--verbosity standard|extended|debug``
   Controls the amount of diagnostic information shown or written.  It does
   not change numerical precision or scientific formulas.

``--quiet``
   Suppresses ordinary terminal presentation where supported.  Calculation and
   persistence still occur.

``--no-progress``
   Disables the live progress display.  It does not disable warnings, change
   the result, or reduce the amount of work performed.

Read :doc:`../cli/conventions` for range syntax, repeatable options, unit
selection, file naming, plotting presets, and exit-status conventions.

Terminal output and redirection
-------------------------------

On an interactive terminal, Rich renders tables, warnings, errors, and live
progress directly to the standard output streams:

- ordinary text follows the terminal theme;
- warnings are highlighted in yellow;
- errors are highlighted in red;
- progress bars are transient.

When output is redirected or terminal control is unavailable, Quantas emits
plain text without ANSI escape sequences.  Styling can also be disabled with:

.. code-block:: console

   QUANTAS_NO_ANSI=1 quantas elasticity run calcite.dat

or with the standard ``NO_COLOR`` environment variable.

Terminal text is not the native scientific result.  Use ``--report`` for a
stable human-readable record and keep the HDF5 output for later calculations,
exports, and plots.

Units and numerical precision
-----------------------------

CLI unit options control input interpretation or displayed/exported units as
described by each command.  They never act as a runtime floating-point
precision selector.

All real calculations and native HDF5 floating values use ``float64``;
complex-valued calculations use ``complex128``.  Display precision belongs to
report and plot renderers and does not alter stored data.

Errors and exit status
----------------------

Quantas validates command syntax and scientific input before or during the
workflow.  A non-zero exit status indicates that the requested command did not
complete successfully.  Typical causes include:

- malformed input;
- incompatible units or dimensions;
- missing optional capabilities;
- an invalid scientific state, such as a non-positive density;
- numerical non-convergence;
- refusal to overwrite an existing file.

Warnings do not automatically imply command failure.  They may mark
extrapolation, marginal fit support, degeneracy, fallback behavior, or another
condition that requires scientific interpretation.

EOS batch execution
-------------------

EOS deserves special attention because ``quantas eos run`` controls one or
more fit jobs rather than producing only one immutable result envelope:

.. code-block:: console

   quantas eos run PV_quartz.dat \
       --pv-eos BM \
       --pv-order 3 \
       --solver effective-variance

A heterogeneous analysis can be described by a strict specification file:

.. code-block:: console

   quantas eos run PV_topaz.dat --spec topaz.spec --dry-run

``--dry-run`` resolves and validates the complete plan without fitting or
creating HDF5.  Use it before a large candidate-model comparison.  See
:doc:`../workflows/eos` for lifecycle semantics and
:doc:`../formats/eos_spec` for the grammar.

Where to continue
-----------------

- :doc:`../cli/index` provides the complete generated command reference.
- :doc:`../workflows/index` explains method choices, defaults, convergence,
  and limitations.
- :doc:`../tutorials/index` demonstrates end-to-end CLI calculations.
- :doc:`results` explains the files created by the commands.
