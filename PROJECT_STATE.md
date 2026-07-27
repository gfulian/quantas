# Quantas project state

This document records the current operational state of the Quantas scientific
backend.  Historical public changes belong in `CHANGELOG.md`, planned release
work belongs in `ROADMAP.md`, and durable coordinated GUI constraints are
recorded separately in the Quantas/Quantas GUI architectural decisions.

## State metadata

| Item | Current value |
|---|---|
| Last updated | 2026-07-27 |
| Current development version | `2.0.0b7` |
| Baseline | `2.0.0b6`, `dev-refactor` line |
| Working branch target | `dev/public-plot-api` |
| Development status | Public plotting API stabilization — incomplete |
| Numerical baseline | Unchanged from `2.0.0b6` |
| Persistence schemas | Unchanged |

## Source hierarchy

Use sources in this order when public behavior or architecture is unclear:

1. the current Quantas `dev-refactor` or approved feature-branch snapshot;
2. this `PROJECT_STATE.md` file;
3. the coordinated architectural decisions for Quantas and Quantas GUI;
4. current tests, documentation, roadmap, and changelog;
5. Quantas `0.9.1`, only as a legacy behavioral and format reference.

Legacy code is not an architectural model for Quantas 2.

## Active objective

Consolidate the public `quantas.api` plotting contract so the CLI, Quantas GUI,
notebooks, and advanced library clients can use the same frontend-neutral
scientific builders without importing implementation packages or maintaining a
second scientific property catalog.

EOS remains session-oriented and is not forced into the one-shot workflow used
by the other scientific modules.

## Completed in the current development cycle

- added the public `quantas.api.plotting` namespace;
- exported the existing authoritative plot specifications and primitives without
  wrappers or data copies;
- retained `quantas.api.common.PlotCollection` as the same class identity;
- routed public workflow facades through the new plotting namespace;
- routed immediately affected CLI plotting type imports through the public API;
- added frozen surface, identity, typed-dispatch, and renderer-independence tests;
- documented the complete plotting namespace and updated API navigation;
- reopened the roadmap freeze only for the scoped `2.0.0b7` public plotting pass;
- bumped package, citation, minimum-baseline comment, and version tests to
  `2.0.0b7`.
- added frozen descriptors for plot properties, representations, contexts, and
  complete result-aware inventories;
- added public ``describe_plots(result)`` discovery for Elasticity, SEISMIC,
  and HA without introducing a generic mapping-based build request;
- added the incremental registry capability ``PLOT_INVENTORY`` only for modules
  whose inventories are currently implemented and tested;
- exposed Elasticity branches, principal planes, and physical/unit-sphere
  geometry compatibility;
- exposed SEISMIC properties according to the calculation level actually stored
  in the result, including group, power-flow, enhancement, surface, and
  polarization availability;
- exposed HA properties, resolved units, and exact stored temperature and
  sampled-volume context;
- centralized the SEISMIC scalar-property catalogue used by inventory discovery
  and existing plot construction.

No plot builder signature, numerical array, scientific selection, unit
conversion, HDF5 schema, or rendering result was intentionally changed in
these increments.

## Verified evidence

- infrastructure, public API, documentation, source hygiene, contracts,
  registry, and packaging selection: `201 passed`;
- complete Elasticity, SEISMIC, and HA module selections: `166 passed`;
- unchanged QHA, Thermoelasticity, and EOS plotting selections: `29 passed`,
  with one pre-existing EOS `RankWarning` for a deliberately sparse polynomial
  fit;
- broader QHA, Thermoelasticity, and EOS selection in the available environment:
  `461 passed` and `29 failed`; every failure requires the unavailable real
  `odrpack` or `spglib` runtime dependency;
- architecture audit: no failures, `66` complexity warnings and `23`
  information items;
- source and test trees compile successfully with ``compileall``;
- an intermediate `2.0.0b7` wheel builds successfully, contains the new
  inventory modules and `py.typed`, and passes an installed-package public API
  smoke test.  This wheel is a development artifact, not a completed release.

## Environment limitations

The current execution environment does not provide Ruff, mypy, Sphinx, a
working `odrpack` backend, or `spglib`.  Therefore the following checks remain
mandatory in the project development environment before merge:

```powershell
ruff check src tests tools docs/tools
mypy
python -m sphinx -W --keep-going -b html docs/source docs/build/html
python -m pytest -q
python tools/check_architecture.py --root .
```

The full test suite failure observed here is environmental and occurs in the
real ODRPACK backend test; it is not a plotting-contract failure.

## Immediate next operation

Review and validate this first common inventory contract in the real project
environment.  Then extend the same discovery principles to QHA, including
exact-grid pressure-axis sections, before implementing Thermoelasticity family
and stage discovery.  A universal generic build request remains explicitly out
of scope.

## Open work for `2.0.0b7`

- final review of the common inventory contract and terminology;
- QHA pressure-axis sections on exact native grid points;
- QHA property, representation, pressure, and temperature inventory;
- Thermoelasticity public family/stage discovery;
- EOS session-oriented archive plot inventory;
- frontend-neutral symbol and unit metadata consolidation;
- PlotSpec serialization contract and round-trip tests;
- incremental CLI migration and CLI/API scientific-equivalence tests;
- Quantas GUI migration from class-name and manual-family dispatch to public
  typed contracts.
