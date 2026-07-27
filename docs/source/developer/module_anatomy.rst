Anatomy of a scientific module
================================

Most Quantas modules use the same conceptual parts even when their exact
lifecycle differs. HA is a compact single-shot example; EOS and
Thermoelasticity intentionally use more specialized controllers.

Typical layout
--------------

A single-shot module commonly contains:

.. code-block:: text

   quantas/modules/example/
   ├── models.py          passive Input, Options, and Result contracts
   ├── analysis.py        frontend-neutral scientific workflow
   ├── calculator.py      active orchestration and event translation
   ├── api.py             internal module facade
   ├── report.py          ReportTable builders
   ├── io/
   │   ├── reader.py      normalized input and HDF5 readers
   │   ├── hdf5_payload.py
   │   ├── export.py
   │   └── inpgen.py      optional interface-to-Quantas input generation
   └── plot/
       ├── spec.py        neutral plot builders
       └── labels.py      module-specific scientific labels

The supported application facade is a separate file:

.. code-block:: text

   quantas/api/example.py

The Click adapter is also separate:

.. code-block:: text

   quantas/cli/example.py

Passive contracts
-----------------

``Input``
   Normalized scientific data. Arrays should have documented shapes, units, and
   copy/ownership semantics. Source provenance belongs here or in metadata.

``Options``
   Scientific and numerical choices. Options should validate themselves or be
   validated before the calculation starts. Presentation-only choices should
   not be mixed into scientific options.

``Result``
   Module-specific payload containing full-precision arrays, masks,
   uncertainties, and diagnostic metadata. It should not contain Matplotlib
   figures, Rich tables, open files, callbacks, or calculator instances.

Use ``slots=True`` where appropriate and copy mutable input arrays when the
contract promises ownership.

Analysis functions
------------------

The analysis layer performs numerical work without emitting Quantas events or
performing frontend I/O. Long loops may accept simple callbacks such as:

.. code-block:: python

   from collections.abc import Callable

   ProgressCallback = Callable[[int, int], None]

   def evaluate_grid(
       values: object,
       *,
       progress_callback: ProgressCallback | None = None,
   ) -> object:
       for index in range(total):
           # Numerical work only.
           if progress_callback is not None:
               progress_callback(index + 1, total)
       return result

The callback communicates only numerical progress. It does not receive a Rich
progress bar or a Click context.

Calculator
----------

A calculator is an active workflow object. It owns:

* normalized input and options;
* workflow stages;
* warning collection;
* event emission;
* conversion of low-level callbacks into Quantas events;
* final assembly of :class:`quantas.models.ResultData`.

For a straightforward workflow, deriving from
:class:`quantas.models.BasicCalculator` provides ``prepare``, ``execute``,
``emit``, warning handling, error events, and finalization.

The calculator should delegate formulas to ``core`` or ``analysis`` functions.
It should not become the only place where a numerical method can be called.

Result envelope
---------------

Single-shot modules return a generic envelope:

.. code-block:: text

   ResultData
   ├── metadata
   ├── input_data
   ├── options
   ├── results[module_key] -> typed Result
   ├── warnings
   └── events

The typed payload contains the science. The envelope contains provenance and
workflow context shared by every module.

Module contract
---------------

A conventional module can define a :class:`quantas.models.ModuleContract` that
connects the standard operations:

* read input;
* run;
* read/write HDF5;
* build report tables;
* build plot specifications.

The contract validates module metadata and the expected result payload key.
The public registry exposes a higher-level capability description for
applications and frontends.

Report builder
--------------

A report builder converts a result into ordered
:class:`quantas.api.common.ReportTable` objects. The table rows retain raw values.
Formatting profiles, alignment, and units are metadata interpreted by a
renderer.

Do not pre-format numerical values as strings merely to align a terminal table.
That destroys reuse by CSV and GUI renderers.

Plot builder
------------

A module plot builder prepares :class:`quantas.api.plotting.PlotCollection` objects.
It selects scientific quantities, masks, labels, uncertainty bands, and
geometries. It does not create a Matplotlib figure.

A new scientific visualization normally belongs in the module plot builder. A
new concrete drawing primitive belongs in the renderer only when the neutral
plot model cannot represent the required data.

Persistence adapter
-------------------

The module HDF5 adapter writes and reads only the scientific payload inside the
shared result envelope. Writers emit the current schema. Readers validate
metadata, units, shapes, and module identity before constructing the typed
result.

A round trip must preserve values and interpretation:

.. code-block:: text

   typed result -> HDF5 -> typed result

Equality tests should cover arrays, masks, metadata required for
interpretation, and optional/``None`` fields.

Internal and public facades
---------------------------

``quantas.modules.<name>.api`` may coordinate internal readers, calculators,
reports, plots, and persistence. The supported application facade under
``quantas.api.<name>`` exposes only reviewed aliases and wrapper functions.

This extra boundary is deliberate. It permits internal reorganization without
changing notebook, CLI, or GUI code.

CLI adapter
-----------

The Click command should:

#. parse strings and paths;
#. construct public ``Input``/``Options`` objects or call public readers;
#. attach a CLI observer;
#. call :mod:`quantas.api`;
#. render reports and plots;
#. choose output filenames and enforce overwrite policy.

It should not contain formulas, build scientific arrays independently, or
repair results returned by the API.

Non-standard module lifecycles
------------------------------

Uniform names are useful, but scientific structure takes priority.

EOS
   Keeps multiple immutable fit records in a persistent archive. It therefore
   exposes fitting, batch, session, calculation, and diagnostic operations
   rather than pretending to be one single-shot result.

Thermoelasticity
   Separates calibration from Point/Grid/Profile analysis. The calibration HDF5
   is reused by later operations.

When a new module genuinely needs a different lifecycle, document it and expose
capabilities explicitly. Do not force it through an untyped universal
``run(**kwargs)`` abstraction.
