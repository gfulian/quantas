Architecture and dependency rules
=================================

The central architectural rule is simple:

   **Scientific results must not depend on the frontend that requested them.**

The Python API, Click CLI, and a future GUI must call the same normalized
workflows and obtain the same ``float64`` arrays.

Layered package model
---------------------

.. code-block:: text

   external files
        |
        v
   quantas.interfaces ------+
        |                    |
        v                    |
   normalized contracts      |
        |                    |
        v                    |
   quantas.modules ----------+----> ResultData / HDF5
        |                                  |
        +----------+-----------------------+
                   |
        +----------+-----------+
        |                      |
        v                      v
   quantas.api            report / plot specs
        |                      |
        +----------+-----------+
                   |
          +--------+--------+
          |                 |
          v                 v
      quantas.cli      future GUI/notebook

``quantas.core``
   Reusable scientific algorithms and policies. Core code does not know which
   module or frontend called it.

``quantas.models``
   Shared passive containers and small active base contracts. This layer has no
   module-specific scientific policy.

``quantas.interfaces``
   Code-specific parsing. Interfaces convert external syntax to explicit
   physical quantities and provenance; they do not run Quantas workflows.

``quantas.modules``
   Scientific orchestration. A module validates normalized input, calls core
   routines, assembles typed results, emits events, persists payloads, builds
   reports and plot specifications, and exports data.

``quantas.api``
   Stable application-facing facades. Applications use this layer instead of
   concrete calculators, readers, or exporters.

``quantas.renderers``
   Concrete presentation of neutral tables and plots.

``quantas.cli``
   Click adapters and Rich terminal behaviour. The CLI parses and validates
   user input, constructs public options, calls :mod:`quantas.api`, observes
   events, and renders the returned neutral objects.

Allowed dependency direction
----------------------------

A lower layer must not import a higher layer. The intended direction is:

.. code-block:: text

   core <- models <- interfaces/modules <- api <- cli or external application
                      |
                      +-> neutral report/plot models -> renderers

Important consequences are:

* ``core`` never imports ``modules``, ``interfaces``, ``renderers``, or ``cli``;
* one scientific module does not import another module directly;
* cross-module transformations belong in :mod:`quantas.api.interop` or a shared
  lower-level contract, not in an accidental module-to-module import;
* interface parsers do not import workflows;
* module plot builders import neutral plot contracts, not Matplotlib;
* the CLI calls public API functions rather than concrete calculators.

The architecture tests enforce these boundaries by inspecting imports.

Active objects and passive data
-------------------------------

Use a normal class for an object that owns behaviour, state, or a lifecycle:

* calculators;
* readers and exporters;
* observers;
* workflow controllers;
* numerical solvers;
* persistent EOS sessions.

Use a dataclass for passive information:

* normalized inputs;
* options;
* metadata;
* scientific results;
* report and plot specifications;
* immutable requests and diagnostics.

A useful test is: if the object performs I/O, emits events, advances through
stages, or mutates internal state, it is probably not a passive dataclass.

Scientific module isolation
---------------------------

Modules share core physics and common contracts but remain independent. QHA
must not import Thermoelasticity merely because Thermoelasticity consumes QHA
results. SEISMIC must not be hidden inside Elasticity merely because both use
elastic tensors.

This separation preserves:

* independent scientific testing;
* explicit transformations between domains;
* smaller import surfaces;
* reuse by frontends;
* the ability to replace one workflow without destabilizing another.

Precision and units
-------------------

The authoritative numerical policy is defined in
``quantas.core.numerics.precision``:

* real calculations: ``numpy.float64``;
* complex calculations: ``numpy.complex128``;
* native HDF5 floating values: ``float64``;
* no runtime precision option.

Input units are normalized at workflow boundaries. Core functions should state
and test the units they consume and return. Display precision and display-unit
selection belong to reports, exporters, and renderers; they must never round or
rescale the stored scientific arrays in place.

Error and warning ownership
---------------------------

Use exceptions for states in which the requested result cannot be interpreted
safely. Use warnings for completed calculations that require scientific
attention. Do not silently clip, repair, reorder, or extrapolate data unless the
workflow explicitly defines that policy and records it in results or masks.

Frontend ownership
------------------

The CLI owns:

* command names and option grouping;
* string-to-type conversion;
* prompts and overwrite confirmation;
* output-path conventions;
* Rich warnings, errors, progress, and tables.

The scientific workflow owns:

* validation of scientific meaning;
* numerical defaults and approximations;
* warnings and failure policies;
* result fields and masks;
* report and plot data.

A future GUI may choose a different visual layout, but it must not reimplement
or bypass the scientific decisions in the workflow.
