# Quantas project state

This document records the current operational state of the Quantas scientific
backend. Historical public changes belong in `CHANGELOG.md`, planned release
work belongs in `ROADMAP.md`, and durable coordinated GUI constraints are
recorded separately in the Quantas/Quantas GUI architectural decisions.

## State metadata

| Item | Current value |
|---|---|
| Last updated | 2026-08-02 |
| Current development version | `2.0.0b9` |
| Baseline | `2.0.0b9`, `dev/refactor` |
| Working branch target | `dev/refactor` release-candidate freeze |
| Development status | Pre-RC scientific validation and release preparation |
| Numerical baseline | Unchanged from `2.0.0b6` |
| Persistence schemas | Unchanged |

## Source hierarchy

Use sources in this order when public behavior or architecture is unclear:

1. the current Quantas `dev/refactor` snapshot;
2. this `PROJECT_STATE.md` file;
3. the coordinated architectural decisions for Quantas and Quantas GUI;
4. current tests, documentation, roadmap, and changelog;
5. Quantas `0.9.1`, only as a legacy behavioral and format reference.

Legacy code is not an architectural model for Quantas 2.

## Active objective

Keep the scientific source tree and public lifecycle frozen while completing
the evidence required for `2.0.0rc1`:

- formal scientific validation for every public workflow;
- publication-quality validation documentation;
- final review of the public API and native HDF5 schemas;
- strict documentation and cross-platform release checks;
- complete TestPyPI publication and clean-install verification.

Quantas GUI consumes this backend baseline through `quantas.api`, but GUI
milestones are independent and are not backend RC blockers.

## Completed backend baseline

The `2.0.0b9` baseline provides:

- one supported public Python surface under `quantas.api`;
- frontend-neutral typed inputs, options, results, reports, plot
  specifications, events, and capability discovery;
- public input generation, execution, native persistence, reopening, reports,
  plots, and derived exports for supported workflows;
- result-aware plot inventories for Elasticity, SEISMIC, HA, QHA,
  Thermoelasticity, and the separate EOS archive lifecycle;
- frontend-neutral QHA input-inspection reports and sampled energy-volume plot
  specifications built from the same fit preview used by the CLI;
- public CRYSTAL/VASP Elasticity and SEISMIC input generation, including VASP
  density extraction and finite-positive-density validation for SEISMIC;
- native SEISMIC HDF5 rejection of unstable or inconsistent elastic media;
- CLI adapters routed through the same public operations used by notebooks,
  services, and Quantas GUI;
- Python support constrained to 3.10 through 3.13 until the complete scientific
  dependency stack is validated on Python 3.14.

EOS remains session-oriented and is not forced into the single-shot lifecycle
used by the other scientific modules.

No approved numerical algorithm, thermodynamic formulation, unit convention,
stored scientific array, or HDF5 schema was intentionally changed by the public
lifecycle, SEISMIC API stabilization, and QHA inspection-presentation work.

## Verified evidence

The complete Windows public-lifecycle validation passed with Python 3.10.11,
NumPy 2.2.6, spglib 2.7.0, and odrpack 0.6.1.

Validated checks include:

- public API surface, lifecycle, registry, input generation, exports, reports,
  persistence, and plot inventories;
- CLI/API equivalence for supported lifecycle operations;
- all staged scientific, CLI, plotting, example, and infrastructure suites;
- Ruff with no findings;
- mypy across the complete `src/quantas` package;
- Sphinx HTML construction with warnings treated as errors;
- wheel and source-distribution builds;
- `twine check`;
- clean installation and public-API smoke testing from the built wheel;
- repository CI on the supported operating-system and Python matrix.

The remaining EOS `RankWarning` messages occur only in deliberately sparse
P--V--T regression cases and are documented expected warnings rather than test
failures.

## Immediate next operation

Build the formal RC validation matrix workflow by workflow. Each public claim
must identify its input, independent reference or analytical identity,
observable, units, comparison method, tolerance, automated test, evidence, and
known methodological limitations.

Begin with exact or controlled references before using literature values:

1. analytical identities and synthetic manufactured data;
2. independent implementations using identical input data and conventions;
3. intentional legacy comparisons;
4. real-data end-to-end benchmarks;
5. literature or experimental comparisons where methodological equivalence can
   be established.

## Remaining blockers for `2.0.0rc1`

- complete and publish the scientific validation matrix;
- finish the validation pages for HA/QHA, Elasticity, SEISMIC, EOS,
  Thermoelasticity, precision, tolerances, and strategy;
- close any correctness defect revealed by validation without expanding scope;
- perform the final public API and HDF5 compatibility review;
- run the complete clean-checkout release validation;
- publish and reinstall the candidate distribution through TestPyPI;
- synchronize the final `2.0.0rc1` version, citation metadata, changelog date,
  tag, GitHub release, and distribution hashes.

New scientific features, empirical mechanical correlations, additional GUI
capabilities, and post-2.0 roadmap work are not RC blockers.
