:orphan:

Development guide
=================

This guide is for contributors who need to understand, modify, test, or extend
Quantas without breaking its scientific behaviour or its frontend-independent
architecture.

Start here
----------

* :doc:`getting_started` — environment, repository orientation, and a safe first
  contribution.
* :doc:`architecture` — dependency direction and ownership of scientific versus
  frontend responsibilities.
* :doc:`module_anatomy` — input/options/result contracts, calculators, reports,
  plots, persistence, API, and CLI.
* :doc:`numerical_methods` — discipline for formulas, tolerances, uncertainty,
  callbacks, and optimization.

Implement changes
-----------------

* :doc:`change_recipes` — add an option, result field, HDF5 field, plot, parser,
  public symbol, default, or citation.
* :doc:`extending` — complete procedure for a new scientific module.
* :doc:`interfaces` — parsers and adapters for external scientific codes.
* :doc:`public_api` — supported facade, capabilities, stability, and promotion
  rules.

Cross-cutting infrastructure
----------------------------

* :doc:`events` — observers, persistent events, warnings, and progress.
* :doc:`persistence` — shared HDF5 envelope, payloads, migrations, and round
  trips.
* :doc:`rendering_frontends` — neutral reports/plots, Rich, Matplotlib, and GUI
  integration.
* :doc:`citation_registry` — canonical bibliographic records and generated
  references.
* :doc:`testing` — test pyramid, staged runner, regression, and distribution
  validation.
* :doc:`documentation` — docstrings, manual structure, CLI/API reference, and
  generated assets.
* :doc:`packaging` — wheel/sdist contents, clean installation, and platforms.

Worked case study and final review
----------------------------------

* :doc:`thermoelastic_architecture` — application of the architecture to a
  staged multi-module scientific workflow.
* :doc:`review_checklist` — scientific, architectural, documentation, testing,
  and packaging checks before a checkpoint.
