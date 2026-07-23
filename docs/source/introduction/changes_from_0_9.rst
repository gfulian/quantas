Changes from Quantas 0.9
========================

Quantas 2.0 is not a minor incremental update of Quantas 0.9.  It is a full
scientific and architectural refactoring that preserves the core scientific
identity of the project while substantially changing how workflows are exposed
and maintained.

Main differences
----------------

Architecture
^^^^^^^^^^^^

Quantas 0.9 was primarily a command-line scientific program with a more tightly
coupled internal organization.  Quantas 2.0 is structured as a reusable Python
library with explicit separation between:

- numerical core;
- passive data contracts;
- persistence;
- rendering;
- command-line adapters.

Public Python usage
^^^^^^^^^^^^^^^^^^^

Quantas 0.9 did not provide a stable, documented public Python surface.  In
Quantas 2.0, the supported entry point is :mod:`quantas.api`.

Persistence and reproducibility
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Quantas 2.0 adopts HDF5 as the native result format.  Result files contain
normalized inputs, options, workflow metadata, results, warnings, and relevant
methodological provenance.

Workflows and interoperability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Quantas 2.0 makes workflow boundaries and interoperability explicit.  For
example:

- QHA results can be transformed into thermoelastic contexts;
- thermoelastic results can be exported to SEISMIC inputs;
- reports and plots are represented as frontend-neutral contracts.

Documentation and validation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The documentation is being reorganized into a full scientific and operational
manual, with clear separation between theory, workflows, tutorials, API, CLI,
formats, and developer guidance.  Validation is also treated as a first-class
part of the public documentation.

What remained scientifically continuous
---------------------------------------

The refactoring preserves the scientific identity of Quantas:

- historical workflow names are retained where scientifically meaningful;
- numerical behavior is validated rather than casually rewritten;
- real examples remain central to the user-facing material;
- theory continues to guide implementation choices.

Quantas 2.0 should therefore be understood as a new public generation of the
software rather than as a one-to-one migration target from the 0.9 series.
