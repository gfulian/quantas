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
| Development status | Public plotting API stabilization |
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

## Completed in the current increment

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

No plot builder, numerical array, scientific selection, unit conversion, HDF5
schema, or rendering result was intentionally changed in this increment.

## Verified evidence

- final infrastructure, API, documentation, model, contract, packaging, and
  repository selection: `118 passed`;
- affected Elasticity and SEISMIC plotting/CLI selection: `50 passed`;
- affected HA, QHA, Thermoelasticity, and EOS plotting/CLI selection:
  `50 passed`, with one pre-existing EOS `RankWarning` for a deliberately sparse
  polynomial fit;
- full suite in the available environment: `468 passed` before the first failure,
  caused by the unavailable real `odrpack` backend;
- architecture audit: no failures, `64` complexity warnings and `23` information
  items.

## Environment limitations

The current execution environment does not provide Ruff, mypy, Sphinx, the
`wheel` build package, or a working `odrpack` backend.  Therefore the following
checks remain mandatory in the project development environment before merge:

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

Define and characterize the smallest shared plot-inventory descriptors, then
apply them first to Elasticity, SEISMIC, and HA.  The next increment must not yet
introduce a universal generic build request.  QHA pressure-axis sections and
Thermoelasticity family discovery follow only after the inventory contract has
been validated on those three representative modules.

## Open work for `2.0.0b7`

- public property, representation, and context inventory contracts;
- Elasticity inventory including branch, plane, and geometry compatibility;
- SEISMIC inventory conditioned on result calculation level;
- HA property and sampled-volume inventory;
- QHA pressure-axis sections on exact native grid points;
- Thermoelasticity public family/stage discovery;
- EOS session-oriented archive plot inventory;
- frontend-neutral symbol and unit metadata consolidation;
- PlotSpec serialization contract and round-trip tests;
- incremental CLI migration and CLI/API scientific-equivalence tests;
- Quantas GUI migration from class-name and manual-family dispatch to public
  typed contracts.
