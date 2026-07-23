Results, reports, and plots
===========================

Quantas deliberately separates the scientific result from its presentation.
A normal calculation may produce several complementary artifacts, each with a
different purpose.

.. list-table:: Output roles
   :header-rows: 1
   :widths: 22 30 48

   * - Artifact
     - Intended use
     - Important limitation
   * - Terminal output
     - Immediate interactive feedback
     - Transient and affected by verbosity and terminal capabilities
   * - Plain-text report
     - Human inspection, archiving, and reproducible records
     - Not a structured replacement for numerical arrays
   * - Native HDF5
     - Complete scientific persistence and downstream workflows
     - Requires Quantas-aware readers for full semantics
   * - CSV or text export
     - Interchange with plotting, spreadsheets, and external tools
     - Usually contains a selected projection of the result
   * - Plot files
     - Visual interpretation and publication
     - Display choices do not constitute new scientific calculations

Result envelopes
----------------

Single-shot workflows return a
:class:`quantas.api.common.ResultData` envelope.  It contains:

- program, version, module, method, and schema metadata;
- normalized input and resolved options;
- one or more module-specific passive payloads;
- warnings and meaningful workflow events;
- numerical-precision metadata;
- optional diagnostics and report text.

Use the module's public ``get_result`` function to obtain the typed payload:

.. code-block:: python

   from quantas.api import elasticity

   envelope = elasticity.read_result("calcite_elasticity.hdf5")
   payload = elasticity.get_result(envelope)

   print(payload.averages.hill.bulk_modulus)
   print(payload.stability.is_stable)

EOS uses a related but different native archive.  It stores one normalized
dataset, immutable fit records, accepted/candidate state, and post-fit history.
That lifecycle is documented in :doc:`../formats/eos_hdf5` and
:doc:`../workflows/eos`.

Native HDF5 results
-------------------

HDF5 is the native scientific format.  A native file is intended to preserve
enough information to:

- reproduce reports without parsing terminal text;
- build new plots without repeating the calculation;
- export selected tables;
- inspect warnings and provenance;
- pass compatible results to another Quantas workflow;
- distinguish file type and schema from metadata rather than from its filename.

Native floating-point arrays use ``float64``.  Complex arrays use
``complex128``.  A renderer may show fewer digits, but display formatting never
rounds or rewrites the stored values.

General result-envelope structure is described in :doc:`../formats/hdf5`.
Module-specific formats are documented under :doc:`../formats/index`.

Reports
-------

Report builders transform structured results into frontend-neutral
:class:`quantas.api.common.ReportTable` objects.  The CLI renders them as
plain text; another frontend can present the same rows as graphical tables.

A report is deterministic:

- no ANSI color sequences;
- no Rich box-drawing characters;
- no live progress artifacts;
- stable units and display profiles;
- warnings and meaningful scientific messages retained where appropriate.

From Python:

.. code-block:: python

   from pathlib import Path
   from quantas.api import elasticity, rendering

   result = elasticity.read_result("calcite_elasticity.hdf5")
   tables = elasticity.build_report(result)
   text = rendering.render_tables(tables)
   Path("calcite_elasticity.log").write_text(text, encoding="utf-8")

A report is designed for reading, not for lossless machine reconstruction.
Keep the HDF5 file whenever the result may be analyzed again.

Tabular exports
---------------

Export commands and API functions write selected numerical views to plain text
or CSV.  Depending on the module, an export may contain:

- temperature or pressure series;
- directional elastic fields;
- one row per seismic direction and mode;
- P--T grids;
- geological profiles;
- EOS diagnostics or calculated curves.

Exports preserve raw numerical values according to the declared format and
units, but they do not necessarily contain the complete input, options,
covariance, warning history, or provenance stored in HDF5.  Treat them as
interchange products rather than native archives.

The shared conventions are described in :doc:`../formats/tabular_outputs`.

Plots
-----

Plot builders create frontend-neutral plot specifications.  A concrete
renderer then produces Matplotlib figures or files.  This two-stage design
keeps plotting outside the numerical calculators.

From the CLI:

.. code-block:: console

   quantas elasticity plot calcite_elasticity.hdf5 \
       --2d \
       --property young \
       --preset publication \
       --output calcite_young

From Python:

.. code-block:: python

   from quantas.api import elasticity, rendering

   options = elasticity.Options(calculate_2d=True)
   result = elasticity.run("calcite.dat", options=options)
   plots = elasticity.build_2d_plots(result)
   rendered = rendering.render_plots(
       plots,
       output_dir="figures",
       filename_prefix="calcite_",
       preset="publication",
       image_format="png",
       dpi=180,
       close=True,
   )

Plot options such as DPI, figure size, colormap, mesh lines, labels, and
annotation density control presentation.  Scientific sampling controls such as
angular grid resolution are documented separately because they can affect the
sampled data rather than merely the appearance.

Plotting requires the optional plotting dependencies when they are not already
installed:

.. code-block:: console

   python -m pip install "quantas[plot]"

Progress and workflow events
----------------------------

Progress events describe current operation and are intended for observers such
as a terminal progress bar or GUI.  They are not persisted as scientific
history.

Warnings, errors, meaningful messages, and structured results can be retained
in the result or archive.  This distinction prevents transient updates such as
``317/7381 directions`` from polluting a permanent scientific record.

What should be kept?
--------------------

For a calculation that may be cited, revisited, or used downstream, retain at
least:

- the original input;
- the native HDF5 result or EOS archive;
- the plain-text report;
- the exact Quantas version;
- any specification file or custom profile used;
- selected exports and figures needed by the analysis.

Figures and CSV files alone are not sufficient to reconstruct the complete
scientific calculation.

Where to continue
-----------------

- :doc:`../formats/index` documents all input and output formats.
- :doc:`../tutorials/index` shows reports, exports, and plots in complete
  workflows.
- :doc:`../api/rendering` lists the public rendering contract.
- :doc:`../cli/conventions` explains common CLI output behavior.
