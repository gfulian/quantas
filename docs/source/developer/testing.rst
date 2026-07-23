Testing strategy and workflow
=============================

Quantas tests protect three different things:

* numerical implementation;
* scientific interpretation;
* architectural and public contracts.

A single end-to-end test cannot replace focused tests at each level.

Test organization
-----------------

Tests are grouped by scientific domain and frontend category. Registered
markers describe physics, modules, architecture, CLI, plotting, HDF5, exports,
examples, and scientific regression.

Named suites can be listed with:

.. code-block:: console

   pytest --list-suites

Examples include:

.. code-block:: console

   pytest --suite elasticity -q
   pytest --suite seismic -q
   pytest --suite ha -q
   pytest --suite qha -q
   pytest --suite cli -q

Complete staged run
-------------------

Use:

.. code-block:: console

   python tools/run_tests.py all -- -q

The runner executes independent stages:

.. code-block:: text

   core
   elasticity
   seismic
   ha
   qha
   eos
   thermoelasticity
   cli
   module-specific plotting stages
   examples

It disables third-party pytest plugin autoloading by default, isolates
Matplotlib configuration, fixes numerical library thread counts, applies a
per-stage timeout, and summarizes status.

Use ``--fail-fast`` only when a complete diagnostic summary is unnecessary.

Test pyramid
------------

Core tests
   Analytical identities, manufactured data, invariance, units, shapes,
   tolerances, and failure cases.

Module tests
   Input normalization, workflow stages, options, results, warnings, and
   failure policy.

Persistence tests
   HDF5 schema, units, round trip, historical migration, and malformed files.

Frontend tests
   CLI parsing/dispatch, report rendering, plot specifications, Matplotlib
   smoke tests, and redirected output.

Public API tests
   ``__all__``, lazy imports, facade operations, registry capabilities, and
   representative workflows.

Example tests
   Curated real datasets and tutorial commands.

Scientific regression
   Frozen comparisons with legacy, analytical, published, experimental, or
   independent external results.

Characterization versus validation
----------------------------------

A characterization test freezes current behaviour. It is useful before a
refactor, but it does not prove that the behaviour is scientifically correct.

A scientific validation test states:

* reference source and version;
* dataset provenance and checksum;
* compared quantities and normalization;
* numerical and scientific tolerances;
* expected limitations.

Keep these purposes explicit in test names and markers.

Tolerance selection
-------------------

Use exact equality for:

* strings, enums, shapes, axis ordering, and masks;
* deterministic serialization where byte-level identity is intended;
* CLI/API outputs that traverse the exact same numerical path and are expected
  to be identical.

Use tolerance-controlled comparison for floating calculations. Choose
``rtol`` and ``atol`` based on scale, conditioning, platform variation, and the
scientific claim. Do not copy a loose tolerance from an unrelated test.

Random and resampling methods
-----------------------------

Tests involving Monte Carlo or bootstrap propagation should:

* use an explicit seed;
* test statistical summaries rather than every sample when appropriate;
* keep sample counts small for unit tests;
* retain a larger scientific regression separately when needed.

Parser tests
------------

External-code parser tests need both positive and negative fixtures. Do not
mock away the source syntax when the parser itself is under test.

Plotting tests
--------------

Test in this order:

#. result arrays;
#. neutral plot specification and masks;
#. renderer dispatch and output creation;
#. selected artist/layout properties only where scientifically relevant.

Avoid brittle pixel snapshots for ordinary formatting changes.

Architecture tests
------------------

The infrastructure suite inspects imports and definitions to preserve:

* frontend-neutral core and modules;
* absence of module-to-module dependencies;
* Matplotlib confinement;
* unique authoritative report/plot contracts;
* public API boundaries;
* complete documentation navigation;
* package and distribution structure.

When adding a new directory or lifecycle, update architecture tests deliberately
rather than disabling them.

Static analysis
---------------

Run:

.. code-block:: console

   ruff check src tests tools
   mypy
   python -m compileall -q src tests tools

Mypy checks the complete ``src/quantas`` package. Do not solve type errors by
removing useful annotations or broadly introducing ``Any``.

Distribution validation
-----------------------

Before a checkpoint or release:

.. code-block:: console

   python -m build
   python -m twine check dist/*
   python tools/check_distribution.py dist

This catches missing examples, documentation assets, ``py.typed``, entry-point
problems, and installed-package import differences that editable tests may
miss.

Adding tests for a change
-------------------------

A cross-cutting change may need several focused tests. For a new result field,
for example, add:

* calculation test;
* dtype/shape/unit test;
* HDF5 round trip;
* report/export/plot test when applicable;
* public API access test;
* CLI/API equivalence;
* documentation structure test;
* scientific regression when the field makes a scientific claim.
