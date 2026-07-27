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
  HA, QHA, and Thermoelasticity without introducing a generic mapping-based
  build request;
- added the incremental registry capability ``PLOT_INVENTORY`` only for modules
  whose inventories are currently implemented and tested;
- exposed Elasticity branches, principal planes, and physical/unit-sphere
  geometry compatibility;
- exposed SEISMIC properties according to the calculation level actually stored
  in the result, including group, power-flow, enhancement, surface, and
  polarization availability;
- exposed HA properties, resolved units, and exact stored temperature and
  sampled-volume context from one authoritative HA plot-property catalogue;
- added exact-grid HA temperature and volume sections plus optional native
  V--T contour specifications while preserving temperature curves as default;
- added QHA result-aware property, representation, pressure, and temperature
  inventory;
- added exact-grid QHA pressure-axis sections at selected stored temperatures,
  plus pressure selection for the existing temperature-axis sections;
- exposed the same HA/QHA scientific section controls through the CLI without
  moving renderer-specific output handling into the public API;
- centralized the SEISMIC scalar-property catalogue used by inventory discovery
  and existing plot construction;
- added public result-aware Thermoelasticity discovery for cumulative
  calibration, P--T grid, profile, comparison, and domain plot families;
- exposed only successful fitted components, available tensor conditions,
  stored profiles, valid uncertainty quantities, and buildable comparison
  coordinates through the Thermoelasticity inventory;
- added archive-aware CLI family listing through the same public
  ``describe_plots`` contract while retaining the previous static help/listing;
- corrected automatic Thermoelasticity plotting for point and one-dimensional
  analyses so it falls back to archived fit diagnostics instead of attempting
  an invalid P--T contour;
- added a focused/full Windows PowerShell validation script for the public plot
  API stabilization branch.
- added a separate public EOS archive plot inventory covering embedded
  datasets, persistent result slots, immutable record history, acceptance state,
  and representation availability without forcing EOS into the one-shot module
  contract;
- exposed detailed common ``PlotInventory`` metadata only after an explicit
  record, accepted slot, or unique accepted result is selected;
- routed EOS CLI plot discovery through ``quantas.api.eos.describe_plots`` and
  removed its duplicate manual representation-description catalogue;
- exposed the archive history and inspection contracts returned by the public
  EOS lifecycle so Python clients and Quantas GUI do not need implementation
  imports.

Existing default plot-builder behavior for complete profile and two-dimensional
grid results, numerical arrays, unit conversion, HDF5 schemas, and scientific
calculations were not intentionally changed.  New HA/QHA section directions and
coordinate selections are opt-in and operate only on exact, unique native grid
points.  Alternative volume or pressure sections are advertised only when at
least two coordinates exist.  The point/one-dimensional Thermoelasticity
fallback is a plotting-dispatch correction and does not change calculation
results.

## Verified evidence

- focused public API, registry, module-contract, documentation,
  source-hygiene, and cross-module plotting/CLI selection, including the EOS
  archive inventory: `164 passed`, with the pre-existing sparse-grid EOS
  ``RankWarning``;
- complete infrastructure plus Thermoelasticity selection: `253 passed`; the
  four remaining failures require ``spglib``, which is unavailable in the
  current environment and is exercised by the CRYSTAL input generator rather
  than the plotting contract;
- complete Elasticity and SEISMIC module selections: `100 passed`;
- complete HA and QHA module selections: `261 passed`;
- EOS inventory and plotting selection: `17 passed`, covering P--V, V--T,
  P--V--T, dataset-only archives, ambiguous accepted slots, canonical CLI keys,
  and builder compatibility, with the pre-existing sparse-grid ``RankWarning``;
- complete EOS module selection: `246 passed`; the remaining `25` failures all
  require the unavailable real ``odrpack`` backend and are not inventory or
  plotting-contract failures;
- fail-fast full-suite run: `479 passed` before the known real-ODRPACK backend
  test failed because the runtime backend is unavailable in this environment;
- architecture audit: no failures, `74` complexity warnings and `24`
  information items;
- source and test trees compile successfully with ``compileall``;
- the intermediate `2.0.0b7` wheel builds successfully through the declared
  setuptools backend, contains the HA/QHA/Thermoelasticity inventory modules
  and ``py.typed``, and passes an installed-package public API smoke test.  This
  wheel is a development artifact, not a completed release.

## Environment limitations

The current execution environment does not provide Ruff, mypy, Sphinx,
PowerShell, a working `odrpack` backend, or `spglib`.  Therefore the following
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

Review and validate the HA/QHA exact-grid sections, cumulative
Thermoelasticity inventory, and separate EOS session/archive inventory in the
real Windows project environment.  Then complete PlotSpec serialization and
CLI/API/GUI equivalence work.  A universal generic build request remains
explicitly out of scope.

## Open work for `2.0.0b7`

- final review of the common inventory contract and terminology;
- Windows validation of the EOS session-oriented archive plot inventory;
- frontend-neutral symbol and unit metadata consolidation;
- PlotSpec serialization contract and round-trip tests;
- incremental CLI migration and CLI/API scientific-equivalence tests;
- Quantas GUI migration from class-name and manual-family dispatch to public
  typed contracts.
