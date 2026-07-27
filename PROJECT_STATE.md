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
| Development status | Public lifecycle API stabilization — validated |
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

Finalize and merge the validated `quantas.api` lifecycle baseline so the CLI,
Quantas GUI, notebooks, and library clients can create supported inputs,
execute workflows, persist and reopen native results, build reports and plots,
and produce derived exports without importing implementation packages or
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
- reopened the roadmap freeze only for the scoped `2.0.0b7` public lifecycle pass;
- bumped package, citation, minimum-baseline comment, and version tests to
  `2.0.0b7`.
- added frozen descriptors for plot properties, representations, contexts, and
  complete result-aware inventories;
- added public ``describe_plots(result)`` discovery for Elasticity, SEISMIC,
  HA, QHA, and Thermoelasticity without introducing a generic mapping-based
  build request;
- introduced ``PLOT_INVENTORY`` incrementally and completed it for every
  supported scientific module;
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
  lifecycle API stabilization branch;
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

Default plot construction for complete profile and two-dimensional grid
results remains unchanged.  New HA/QHA section directions are opt-in and use
only exact native grid points; alternative volume or pressure sections are
offered only when at least two coordinates exist.  The point/one-dimensional
Thermoelasticity fallback corrects plot selection without changing calculation
results.  Numerical arrays, unit conversion, HDF5 schemas, and approved
scientific calculations were not changed.

## Verified evidence

The complete Windows validation passed in the project virtual environment with:

| Component | Validated version |
|---|---|
| Python | `3.10.11` |
| NumPy | `2.2.6` |
| spglib | `2.7.0` |
| odrpack | `0.6.1` |

Validated checks include:

- public lifecycle surface and documentation: `83 passed`;
- input generation and CLI/API equivalence: `51 passed`;
- public exports and persistence: `81 passed`;
- frontend-neutral plot contracts: `88 passed`;
- every staged scientific, CLI, plotting, and example suite;
- Ruff with no findings;
- mypy with no issues in `365` source files;
- Sphinx HTML build with warnings treated as errors;
- wheel and source-distribution builds;
- `twine check` for both distributions;
- clean installation and public-API smoke testing from the built wheel;
- architecture audit with no errors (`74` maintainability warnings and `24`
  information items).

The EOS tests emit a small number of expected NumPy ``RankWarning`` messages
for deliberately sparse P--V--T fits.  They do not fail validation and do not
indicate a change in scientific behavior.

## Immediate next operation

Push the validated `dev/public-plot-api` branch, review the complete diff, and
merge it into `dev/refactor` after the normal repository review.  Then update
Quantas GUI to require this backend baseline and consume the public registry,
input, report, plot, persistence, and export contracts without internal imports
or duplicate scientific catalogs.

## Remaining release steps for `2.0.0b7`

No functional backend blocker remains from the public-lifecycle stabilization.
Before merge or publication:

- confirm that `git status` contains no build directories, caches, generated
  documentation, local reports, or virtual environments;
- review the final branch diff and commit history;
- push the feature branch and run the repository CI checks;
- merge into `dev/refactor` only after review;
- update the release date and final distribution checksums when the beta is
  actually published.

Quantas GUI integration is the next coordinated project step, but it is not a
backend validation blocker for this beta.
