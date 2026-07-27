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
| Development status | Public lifecycle API stabilization — validation pending |
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

Consolidate the complete public `quantas.api` lifecycle so CLI, Quantas GUI,
notebooks, and advanced library clients can create supported inputs, execute
workflows, persist and reopen native results, build reports and plots, and
produce derived exports without importing implementation packages or
maintaining duplicate scientific catalogs.

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
- added a focused/full Windows PowerShell validation script for the public
  lifecycle API stabilization branch.
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
  imports;
- added public Elasticity input generation and principal-plane table export,
  then routed the corresponding CLI commands through those facade operations;
- exposed the shared phonon input generator from both HA and QHA namespaces so
  each workflow has a complete discoverable input lifecycle;
- added public HA and QHA table writers and routed their CLI export commands
  through the same functions used by notebooks and future GUI downloads;
- closed public ``Input`` and ``Options`` annotations by re-exporting the
  required rotation, structure, symmetry, seismic selector, thermoelastic
  fitting, QHA coupling, and EOS model contracts;
- extended the registry with named operation descriptors and multiple-operation
  discovery for input generation, templates, tabular exports, tensor exports,
  profile exports, interoperability exports, and EOS post-fit CSV products.

Existing default plot-builder behavior for complete profile and two-dimensional
grid results, numerical arrays, unit conversion, HDF5 schemas, and scientific
calculations were not intentionally changed.  New HA/QHA section directions and
coordinate selections are opt-in and operate only on exact, unique native grid
points.  Alternative volume or pressure sections are advertised only when at
least two coordinates exist.  The point/one-dimensional Thermoelasticity
fallback is a plotting-dispatch correction and does not change calculation
results.

## Verified evidence

- complete infrastructure, public-surface, type-closure, documentation, source
  hygiene, registry, and lifecycle selection: `215 passed`;
- complete Elasticity, HA, QHA, and SEISMIC module selections: `361 passed`;
- focused CLI/API lifecycle tests verify byte-identical Elasticity input text and
  HA/QHA table exports produced through the public functions and CLI commands;
- complete Thermoelasticity selection: `49 passed`; the remaining `4` failures
  require ``spglib``, unavailable in the current environment, and belong to the
  CRYSTAL structural input generator;
- complete EOS selection: `248 passed`; the remaining `25` failures require the
  unavailable real ``odrpack`` backend;
- fail-fast full-suite run: `488 passed` before the known real-ODRPACK backend
  test failed because the runtime backend is unavailable in this environment;
- public callable type hints resolve successfully and require no unexported
  Quantas type identities; public ``Input``, ``Options``, requests, and contexts
  are closed against implementation-only imports;
- architecture audit: no failures, `74` complexity warnings and `24`
  information items;
- source and test trees compile successfully with ``compileall``.

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

Run the complete Windows validation against the consolidated public lifecycle,
including real input generators, table exports, CLI/API equivalence, Ruff,
mypy, Sphinx, packaging, ``odrpack``, and ``spglib``.  Then migrate Quantas GUI
to require the validated backend and consume registry, inventory, report, and
plot contracts without payload inspection or duplicate catalogs.

## Open work for `2.0.0b7`

- Windows validation of public input generation and every derived export;
- full CLI/API scientific-equivalence tests for representative lifecycle paths;
- frontend-neutral symbol and unit metadata consolidation;
- PlotSpec and descriptor serialization contract with round-trip tests;
- review of remaining CLI-only incremental report helpers, which may stay
  frontend-specific only when a complete public report equivalent exists;
- Quantas GUI migration from class-name, payload inspection, and manual-family
  dispatch to public typed contracts;
- final cross-platform Ruff, mypy, Sphinx, packaging, ``odrpack``, and
  ``spglib`` validation before declaring ``2.0.0b7`` complete.
