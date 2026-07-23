Adding a new scientific module
==============================

A new module is a scientific feature, an architectural feature, a persistence
contract, and a public interface at the same time. Build it from the numerical
core outward. Do not begin by writing a Click command.

Phase 1: define the scientific contract
---------------------------------------

Write a short design note that states:

* physical quantities and conventions;
* input data and units;
* output fields, shapes, and units;
* approximations and domain of validity;
* error, warning, and extrapolation policies;
* uncertainty model;
* expected computational cost;
* validation references and datasets;
* whether the lifecycle is single-shot or staged.

Preserve historically meaningful Quantas names unless there is a strong
scientific reason to change them.

Phase 2: implement reusable numerical code
------------------------------------------

Place general mathematics or physics under ``quantas.core``. A core function:

* has no frontend imports;
* has a complete technical docstring;
* accepts explicit arrays and scalar parameters;
* states units and shapes;
* returns passive values;
* may accept a simple progress callback for a long loop;
* has analytical/manufactured and edge-case tests.

Keep module-specific orchestration out of core. If a function is meaningful
only inside one workflow and is unlikely to be reused, it may live in the
module analysis layer while remaining frontend-neutral.

Phase 3: add passive contracts
------------------------------

Create ``models.py`` with dataclasses for input, options, and result. A minimal
pattern is:

.. code-block:: python

   from __future__ import annotations

   from dataclasses import dataclass, field
   import numpy as np
   from numpy.typing import NDArray

   FloatArray = NDArray[np.float64]

   @dataclass(slots=True)
   class ExampleInput:
       values: FloatArray
       unit: str = "GPa"

   @dataclass(slots=True)
   class ExampleOptions:
       tolerance: float = 1.0e-10

       def validate(self) -> None:
           if self.tolerance < 0.0:
               raise ValueError("tolerance must be non-negative")

   @dataclass(slots=True)
   class ExampleResult:
       transformed: FloatArray
       valid_mask: NDArray[np.bool_]
       metadata: dict[str, object] = field(default_factory=dict)

Normalize arrays to the Quantas dtype and copy them where the contract should
own its data. Avoid embedding open handles, observers, or figures.

Phase 4: create the analysis layer
----------------------------------

Implement validation and the complete scientific operation in ``analysis.py``.
The analysis layer should be directly unit-testable without constructing a
calculator.

Separate stages into named functions rather than one very long procedure. For
example:

.. code-block:: text

   validate_input
   prepare_coordinates
   solve_states
   derive_properties
   classify_quality
   assemble_result

This makes characterization, scientific validation, and later optimization
much safer.

Phase 5: add the calculator
---------------------------

Use an active calculator to:

* own input/options/observer;
* translate progress callbacks into events;
* emit meaningful stage and result events;
* collect warnings;
* assemble ``ResultMetadata`` and ``ResultData``;
* preserve persistent event records;
* fail explicitly when a result cannot be interpreted.

Do not duplicate formulas in the calculator.

Phase 6: input and external interfaces
--------------------------------------

Create a normalized Quantas input reader under the module ``io`` package. If
source data come from CRYSTAL, VASP, Phonopy, or another code, implement the
code-specific parser under ``quantas.interfaces`` and convert its output into
the module input contract.

Document and test:

* file detection;
* completion markers;
* units and conversion;
* tensor or coordinate conventions;
* malformed and incomplete output;
* source provenance;
* representative fixtures.

See :doc:`interfaces`.

Phase 7: persistence
--------------------

Add module payload read/write helpers while reusing the shared HDF5 envelope.
Writers should attach units and descriptions to numerical datasets. Readers
validate:

* program and module metadata;
* schema version;
* required groups and fields;
* units;
* shapes;
* optional/``None`` representation;
* scientific interpretation.

Add a write/read round-trip test and at least one malformed-file test. See
:doc:`persistence`.

Phase 8: reports, plots, and exports
------------------------------------

Build raw-value ``ReportTable`` objects and neutral ``PlotCollection`` objects.
Add a concrete renderer primitive only when existing neutral contracts cannot
represent the scientific visualization.

Exports should inspect result metadata, not infer the module from a filename.
Exported values retain full available precision; display formatting is a
separate choice.

Phase 9: internal and public APIs
---------------------------------

Create the internal module facade and a reviewed ``quantas.api.<name>`` facade.
For a conventional module expose, where meaningful:

.. code-block:: text

   Input, Options, Result
   read_input, normalize_input
   run, get_result
   read_result, write_result
   build_report, build_plots

Add a registry descriptor and capabilities. A staged or persistent workflow may
use a different vocabulary, but that difference must be deliberate and
documented.

Phase 10: CLI adapter
---------------------

Only after the public API is stable should the Click adapter be added. Group
options by scientific purpose. Each option needs a help string that explains
what it controls, not merely its type.

The CLI calls public API functions, attaches an observer, renders neutral
objects, and manages paths. It does not contain scientific calculations.

Phase 11: documentation
-----------------------

A complete module normally needs:

* Scientific Background;
* Implementation and Workflows;
* tutorial with real or representative data;
* input/output format specification;
* command reference;
* curated API reference;
* validation page;
* developer notes when the architecture is unusual;
* canonical citations.

Update the docs and code together. Do not postpone all user-facing explanation
until the end of the implementation.

Phase 12: validation matrix
---------------------------

The minimum test matrix includes:

* core analytical/manufactured tests;
* input validation and parser fixtures;
* calculator/result tests;
* events and warning policy;
* HDF5 round trip;
* report and export semantics;
* plot-specification tests;
* public API surface and usage;
* CLI/API equivalence;
* curated example smoke test;
* scientific regression against an independent reference;
* architecture and packaging tests.

Completion checklist
--------------------

A new module is not complete until:

* it can be used without the CLI;
* the CLI and API give equivalent numerical results;
* native files can be read by metadata;
* no scientific module imports another scientific module;
* plotting is optional;
* public names are documented exactly once;
* the full staged test runner passes;
* a clean wheel and source distribution contain all required data.
