Python API
==========

The supported Python contract begins at :mod:`quantas.api`.  Applications,
notebooks, services, the command-line interface, and Quantas GUI should import
from this namespace rather than from implementation packages under
:mod:`quantas.core`, :mod:`quantas.modules`, or :mod:`quantas.renderers`.

Public namespaces are organized by scientific domain:

.. code-block:: python

   from quantas.api import (
       elasticity,
       eos,
       ha,
       interop,
       profiles,
       plotting,
       qha,
       registry,
       rendering,
       seismic,
       thermoelasticity,
   )

The API is deliberately not one flat collection of functions.  Each namespace
exposes the passive data contracts and operations appropriate to that
scientific workflow.

The single-shot workflow pattern
--------------------------------

Elasticity, SEISMIC, HA, QHA, and the calibration stage of Thermoelasticity use
the same broad pattern:

.. code-block:: python

   from quantas.api import elasticity

   options = elasticity.Options(calculate_2d=True)
   result = elasticity.run("calcite.dat", options=options)
   payload = elasticity.get_result(result)

The three objects have different responsibilities.

``Options``
   A passive, validated description of scientific and numerical choices.  It
   contains no terminal state and performs no calculation.

``ResultData``
   The complete frontend-neutral envelope returned by the workflow.  It
   records program and schema metadata, normalized input, resolved options,
   warnings, persistent events, and module-specific results.

Typed payload
   The domain result returned by ``get_result``.  It provides named attributes
   and numerical arrays without requiring the caller to navigate a generic
   string-keyed mapping.

Use the module-specific ``get_result`` function whenever one is available.  It
checks that the envelope belongs to the expected module and that the payload
has the correct type.

Creating, reading, and normalizing input
----------------------------------------

A workflow exposes ``create_input`` when Quantas can convert output from a
supported external code into its standard input format.  The same function is
available to the CLI, notebooks, scripts, and graphical frontends.

.. code-block:: python

   from quantas.api import elasticity, qha, seismic

   elasticity_path = elasticity.create_input(
       "calcite-elastcon.out",
       "calcite.dat",
       interface="crystal",
       jobname="Calcite",
   )
   seismic_path = seismic.create_input(
       "OUTCAR",
       "material-seismic.dat",
       interface="vasp",
       jobname="Material",
   )
   qha_path = qha.create_input(
       "phonon-list.txt",
       "calcite-qha.yaml",
       interface="crystal-qha",
       is_list=True,
   )

HA and QHA use the same phonon input format, but both namespaces expose
``create_input`` so users can remain within the workflow they are running.
Elasticity and SEISMIC likewise expose workflow-owned entry points to their
shared elastic input format; SEISMIC generation additionally requires finite
positive density metadata. Thermoelasticity provides the corresponding
conversion for supported elastic-volume output series. EOS reads its
documented input formats or accepts public Python data objects directly.

Public modules commonly expose two related parsing operations:

``read_input(path)``
   Parse one supported file format and return a validated passive input object.

``normalize_input(value)``
   Accept either an existing input object or a path and return the normalized
   input contract used by the workflow.

For example:

.. code-block:: python

   source = elasticity.read_input("calcite.dat")
   normalized = elasticity.normalize_input(source)
   result = elasticity.run(normalized)

Reading the input explicitly is useful when a notebook needs to inspect or
modify passive options before running.  Passing the path directly to ``run``
is the shorter equivalent for ordinary calculations.

Reports and plots remain frontend-neutral
-----------------------------------------

Scientific modules build neutral report tables and plot descriptions.  A
separate public rendering namespace turns those descriptions into concrete
artifacts:

.. code-block:: python

   from pathlib import Path
   from quantas.api import elasticity, plotting, rendering

   options = elasticity.Options(calculate_2d=True)
   result = elasticity.run("calcite.dat", options=options)

   tables = elasticity.build_report(result)
   report_text = rendering.render_tables(tables)
   Path("calcite.log").write_text(report_text, encoding="utf-8")

   plots = elasticity.build_plots(result)
   assert all(isinstance(item, plotting.PlotSpec) for item in plots.plots)
   rendered = rendering.render_plots(
       plots,
       output_dir="figures",
       preset="publication",
       close=True,
   )

The calculator does not import Rich, Matplotlib, or a GUI toolkit.  This
separation is what allows CLI, notebooks, services, and graphical frontends to
use the same numerical result.

Public table and CSV exports
----------------------------

HDF5 remains the complete scientific record.  Text and CSV tables are
convenient derived exports, built by public functions so the CLI, notebooks,
and GUI use the same property selection, units, and data layout.

.. code-block:: python

   elasticity.write_table(result, "calcite-directions.dat")
   ha.write_table(ha_result, "free-energy.dat", property_name="F")
   qha.write_table(
       qha_result,
       "equilibrium-volume.csv",
       property_name="VT",
       file_format="csv",
   )

SEISMIC, Thermoelasticity, and EOS expose their corresponding public CSV,
tensor, grid, profile, or post-fit writers in their own namespaces.  HDF5
remains the persisted scientific source of truth.

Native result persistence
-------------------------

Single-shot modules expose matching read and write operations:

.. code-block:: python

   elasticity.write_result(
       result,
       "calcite_elasticity.hdf5",
       report_text=report_text,
   )

   restored = elasticity.read_result("calcite_elasticity.hdf5")
   restored_payload = elasticity.get_result(restored)

HDF5 is the scientific source of truth.  Display precision in a report or
figure never changes the stored ``float64`` values.

Domain-specific operations
--------------------------

A common vocabulary is used where the operation has the same meaning, but the
API does not force scientifically different workflows into one artificial
base class.  Examples include:

- :func:`quantas.api.qha.inspect` for preflight analysis of a QHA input;
- :func:`quantas.api.thermoelasticity.analyze_grid` for one state or a P--T grid;
- :func:`quantas.api.thermoelasticity.analyze_profiles` for geological paths;
- :func:`quantas.api.thermoelasticity.point_table` for one-state reporting;
- EOS fit, batch, archive, diagnostics, and calculation operations.

Consult the domain page in the :doc:`../api/index` reference rather than
assuming that every namespace exposes the same functions.

EOS uses an archive lifecycle
-----------------------------

EOS differs from the single-shot modules because one dataset may be fitted
many times with different models, solvers, constraints, and data selections.
Its HDF5 file therefore contains immutable fit records and mutable acceptance
state rather than one result envelope.

A direct fit can be constructed as:

.. code-block:: python

   from quantas.api import eos

   dataset = eos.read_input("PV_quartz.dat")
   request = eos.FitRequest(model="BM3")
   fit = eos.fit(dataset, request)

For persistent batch work:

.. code-block:: python

   plan = eos.BatchPlan(jobs=(eos.BatchJob(request=request),))
   summary = eos.run_batch(
       dataset,
       plan,
       "quartz_eos.hdf5",
       overwrite=True,
   )

Post-fit calculations operate on the stored archive:

.. code-block:: python

   diagnostics = eos.diagnose("quartz_eos.hdf5")
   properties = eos.calculate(
       "quartz_eos.hdf5",
       pressure=[0.0, 5.0, 10.0],
   )
   plots = eos.build_plots("quartz_eos.hdf5")

The different lifecycle is scientific, not merely syntactic.  It is explained
in :doc:`../workflows/eos`.

Cross-module transformations
----------------------------

Scientifically meaningful coupling is explicit in :mod:`quantas.api.interop`.
For example, one thermoelastic state can be transformed into a SEISMIC input:

.. code-block:: python

   from quantas.api import interop, seismic

   seismic_input = interop.thermoelastic_to_seismic(
       "dolomite_thermoelastic.hdf5",
       pressure=2.0,
       temperature=600.0,
       tensor_condition="adiabatic",
   )
   seismic_result = seismic.run(seismic_input)

The transformation centralizes unit conversion, coverage checks, tensor
condition, and provenance instead of duplicating them in user code.

Events and progress
-------------------

Public ``run`` functions accept frontend-neutral observers.  They never create
a terminal progress bar or GUI widget themselves:

.. code-block:: python

   from quantas.api import elasticity
   from quantas.api.common import ListObserver

   observer = ListObserver()
   result = elasticity.run("calcite.dat", observer=observer)

   for event in observer.events:
       print(event.level.value, event.message, event.progress)

A caller can map the same events to logging, a notebook display, a web service,
a GUI progress indicator, or no visible output.  Operational progress events
are not persisted as scientific history.

Capability discovery
--------------------

A frontend can discover supported modules and operations without importing
private classes or hard-coding a single inheritance hierarchy:

.. code-block:: python

   from quantas.api import registry
   from quantas.api.registry import Capability

   for descriptor in registry.list_modules():
       capabilities = sorted(item.value for item in descriptor.capabilities)
       print(descriptor.name, capabilities)

   qha = registry.get("qha")
   run_qha = qha.operation(Capability.RUN)
   options_type = qha.options_type

   elasticity = registry.get("elasticity")
   if elasticity.has(Capability.PLOT_INVENTORY):
       describe = elasticity.operation(Capability.PLOT_INVENTORY)
       inventory = describe(registry.open_result("calcite.hdf5"))
       print([item.key for item in inventory.properties])

   for operation in registry.get("thermoelasticity").list_operations(
       Capability.EXPORT
   ):
       print(operation.key, operation.name)

EOS declares ``FIT``, ``BATCH``, and ``ARCHIVE`` capabilities instead of
pretending to be a single-shot ``RUN`` module.  Native HDF5 files can also be
dispatched from their metadata:

.. code-block:: python

   result_or_archive = registry.open_result("quantas_result.hdf5")

Errors and validation
---------------------

Public operations raise Python exceptions for invalid inputs or unsupported
scientific states.  Do not suppress these errors indiscriminately.  A useful
application should:

- show the exception message to the user;
- preserve warnings and diagnostics from successful results;
- distinguish invalid input from numerical non-convergence;
- avoid converting a failed or extrapolated state into an accepted result
  without an explicit scientific decision.

Where to continue
-----------------

- :doc:`../api/index` lists the complete public contract.
- :doc:`../workflows/index` explains implementation choices and limitations.
- :doc:`../tutorials/index` provides complete parallel CLI/API examples.
- :doc:`results` explains persistence, reports, exports, and rendering.
