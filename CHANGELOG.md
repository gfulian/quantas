# Changelog

All notable public changes to Quantas are recorded here.  The project follows
Semantic Versioning after the first stable Quantas 2 release.  During the current
beta, breaking changes are permitted when they simplify and stabilize the final
public contract; they must still be documented and validated.

## [2.0.0b7] - Unreleased

### Added

- Added the public ``quantas.api.plotting`` namespace for the existing
  frontend-neutral plot specifications, axes, series, masks, overlays, surface
  layers, directional fields, composite panels, and plot collections.
- Added frozen frontend-neutral descriptors for plot properties,
  representations, scientific contexts, and complete result-aware inventories.
- Added public ``describe_plots(result)`` discovery for Elasticity, SEISMIC,
  HA, QHA, and Thermoelasticity, together with the incremental registry capability
  ``PLOT_INVENTORY``.
- Added frozen API-surface, renderer-independence, and typed-dispatch tests for
  the public plotting contracts.
- Added inventory consistency and builder-compatibility tests covering
  Elasticity branches and geometries, result-conditioned SEISMIC fields, HA
  temperature-volume grids, and QHA pressure-temperature grids.
- Added exact native-grid HA sections along temperature or sampled volume,
  optional V--T contour specifications, and public ``HAPlotOptions``.
- Added exact native-grid QHA sections along temperature or pressure, including
  selected pressure/temperature conditions and result-aware P--T inventory
  metadata.
- Added cumulative Thermoelasticity discovery for elastic-volume fits, P--T
  stiffness maps, depth profiles, isothermal--adiabatic comparisons, and
  calibration-domain diagnostics.  The inventory reports only families,
  components, tensor conditions, profiles, and coordinates that are buildable
  from the supplied archive stage.
- Added result-aware CLI discovery through
  ``quantas thermoelasticity plot --list-plots --archive RESULT.hdf5`` while
  retaining the existing static family overview when no archive is supplied.
- Added a Windows PowerShell validation entry point for the in-development
  public plotting API, with focused checks by default and the complete static,
  scientific, documentation, and distribution checks under ``-Full``.
- Added a separate session-aware EOS archive plotting inventory exposing
  lightweight dataset, result-slot, immutable-record, acceptance-state, and
  representation descriptors.  Detailed common plot metadata is returned only
  for an explicit record, accepted slot, or unique accepted result; ambiguous
  archives remain browseable without an arbitrary selection.
- Added public EOS archive-history and inspection types required by notebooks,
  CLI adapters, Quantas GUI, and other advanced library clients to consume the
  persistent workflow without importing implementation modules.

### Changed

- Public workflow facades and the immediately affected CLI plotting adapters now
  reference plot contracts through ``quantas.api.plotting`` instead of relying
  on the implementation namespace.
- Centralized SEISMIC scalar-property discovery so the public inventory and the
  existing plot builders derive labels and availability from one authoritative
  result-aware catalogue.
- Centralized HA plot keys, names, mathematical symbols, plain symbols, and
  scientific categories so inventory discovery and existing builders use one
  authoritative backend catalogue.
- Extended the HA and QHA CLI plot commands with the same public scientific
  axis and exact-coordinate selections used by the library API.
- Routed result-aware Thermoelasticity family listing through the same public
  ``describe_plots`` contract used by Python clients and Quantas GUI.
- Routed EOS ``--list-plots`` through the public session-aware inventory and
  removed the CLI-owned plot-description catalogue.  Ambiguous archives now
  list result slots and plottable record identifiers instead of guessing which
  accepted fit should be used.
- Reopened the pre-RC public-API freeze for the narrowly scoped plotting-contract
  stabilization required by CLI, GUI, notebooks, and scientific-library users.

### Fixed

- Prevented the default Thermoelasticity plot builder from attempting an
  invalid P--T contour for point or one-dimensional analysis archives.  The
  pre-existing priority remains profile, two-dimensional P--T grid, then
  calibration fit; point and one-dimensional archives now fall back to their
  archived calibration plots.

### Scientific compatibility

The public-plotting increments completed so far alias the existing authoritative
plot specifications, add passive discovery metadata, and add opt-in exact-grid
sections over already stored HA/QHA arrays.  Alternative volume or pressure
sections are exposed only when at least two native coordinates are present.
Thermoelasticity discovery inspects archived fit, tensor, uncertainty, profile,
mask, and grid availability without recalculating or mutating scientific data.
EOS discovery opens archives read-only and exposes only lightweight dataset,
slot, history, and representation metadata; it neither changes acceptance state
nor replaces the persistent session lifecycle.  Existing default construction
for profile and two-dimensional P--T archives is unchanged; the
point/one-dimensional fallback corrects a previously invalid automatic contour
request.  No numerical algorithm, thermodynamic calculation, stored array,
tensor convention, or HDF5 schema is modified.  The ``2.0.0b7`` development
cycle remains open until serialization, frontend equivalence, and full release
validation are complete.

## [2.0.0b6] - 2026-07-24

### Added

- Added a complete EOS tutorial sequence for P--V, V--T, and P--V--T analysis,
  including strict batch specifications, OLS/WLS/effective-variance/ODR guidance,
  public-API equivalents, diagnostics, post-fit calculations, publication plots,
  and guided model-selection exercises based on curated quartz, rutile, and NaF
  datasets.

### Fixed

- Constrained V--T plotting grids to non-negative temperatures so curve
  padding cannot evaluate thermal-expansion models below zero kelvin.
- Corrected HA reporting for native multi-volume calculations. Zero-point energy
  is now reported once per sampled volume, while temperature-dependent
  thermodynamic properties retain their full temperature-volume matrices. The
  same temperature-independent representation is handled correctly by HA plots
  without changing stored arrays or HDF5 payloads.
- Corrected the remaining NumPy shape annotations reported by the supported
  Python 3.10 mypy baseline without changing runtime array values or shapes.
- Passed the standard ``screen`` preset through the thermoelastic profile-analysis
  plotting shortcut and added CLI coverage for the complete ``--plot`` path.
- Made the real quartz effective-variance regression assert scientific convergence
  and diagnostic consistency rather than a platform-specific iteration count.
- Made the documentation asset generator satisfy Ruff while preserving direct
  execution from an unpacked source tree.

### Validation

- Added a fail-fast release-validation script that creates a fresh virtual
  environment, installs Quantas from the current source tree, runs all static,
  scientific, documentation, packaging, distribution, and repository checks,
  and records the complete environment and output.

### Scientific compatibility

No numerical algorithm, equation, physical convention, HDF5 schema, or public API
contract was intentionally changed.

## [2.0.0b5] - 2026-07-22

### Changed

- Matplotlib is now the sole supported graphical renderer in the scientific
  package; neutral plot specifications remain available for future GUI adapters.
- Source-hygiene tests inspect only first-party project directories and therefore
  work correctly when a virtual environment is created inside an extracted tree.
- Static typing uses a reproducible Python 3.10 baseline with NumPy stubs below
  2.3, while runtime CI continues to exercise current dependency versions.
- Documentation builds regenerate tutorial figures and downloadable HDF5 results
  from curated examples before Sphinx reads the source tree.
- Repository-only validation is handled by a Git-aware developer tool that skips
  cleanly for ZIP and source-distribution checkouts.
- Added portable `.gitattributes` and `.editorconfig` policies and removed the
  unsupported `uv.lock` development contract.

### Removed

- Removed the provisional GUI renderer module, its compatibility test, marker,
  example path, and public exports from the scientific package.

### Fixed

- Relaxed only the absolute tolerance of near-zero SEISMIC eigenvalue-Hessian
  reference terms to remain portable across supported NumPy and LAPACK builds.
- Removed Sphinx failures caused by missing generated Elasticity tutorial assets.

### Scientific compatibility

No scientific formula, physical convention, HDF5 payload, public API domain, or
workflow result was intentionally changed.

## [2.0.0b4] - 2026-07-22

### Changed

- Removed historical per-file copyright banners while retaining explicit UTF-8
  declarations in the installed source package.
- Simplified the root package namespace to expose only the authoritative version;
  CLI display metadata now lives in a private module.
- Reworked public and legacy interface docstrings so every public module, class,
  function, and method has a concise technical description.
- Reduced the repository README to a stable project entry point, leaving detailed
  scientific guidance to the Sphinx manual.
- Restored the academic citation notice in a concise form consistent with the
  original Quantas presentation.

### Fixed

- Plain-text tables now calculate widths from visible lines and render multiline
  cells without producing extremely long whitespace-filled rows.
- HA and QHA API reports summarize structural volume-series metadata instead of
  dumping nested arrays and mappings into one table cell.

### Scientific compatibility

This final source-cleanup pass does not intentionally change formulas, units,
array shapes, HDF5 numerical payloads, tensor conventions, or validated
tolerances.

## [2.0.0b3] - 2026-07-22

### Added

- Shared ``screen``, ``publication``, and ``monochrome`` static-figure presets
  across every plotting CLI, with explicit rendering options taking priority.
- Source-tree hygiene tests protecting readable test names, developer-tool names,
  and the import mode required by repeated semantic test filenames.

### Changed

- Test directories, files, functions, scientific reference bundles, and developer
  tools now use domain-oriented names rather than refactoring milestones or
  historical implementation labels.
- Pytest uses importlib collection so different scientific domains can safely use
  concise filenames such as ``test_plotting.py`` and ``test_hdf5.py``.
- Public maintenance tools use task-oriented names and write local reports under
  ``project_internal/checks``.
- Plotting commands use ``--preset`` and calculation commands with optional plots
  use ``--plot-preset``; ``--dpi`` remains an explicit override.

### Scientific compatibility

This source-freeze cleanup does not intentionally change formulas, units, array
shapes, HDF5 numerical payloads, tensor conventions, or approved tolerances.

## [2.0.0b2] - 2026-07-22

### Added

- A single organized public Python API under `quantas.api`, with dedicated
  namespaces for Elasticity, SEISMIC, HA, QHA, EOS, and quasi-static
  thermoelasticity.
- A lazy, frontend-neutral module registry with explicit capabilities, input,
  option, and result contracts for future GUI discovery.
- Cross-module interoperability helpers for QHA-to-thermoelastic and
  thermoelastic-to-SEISMIC workflows.
- Public terrestrial-profile helpers and metadata-driven result opening.
- Frozen API-surface tests that detect accidental additions, removals, and
  exposure of calculator, reader, exporter, renderer, or workflow internals.

### Changed

- The Click CLI now dispatches calculation workflows through the same public API
  used by Python applications and the future GUI.
- Public examples and API documentation now use domain namespaces such as
  `quantas.api.qha`, `quantas.api.eos`, and `quantas.api.thermoelasticity`.
- EOS keeps its scientifically appropriate fit/archive/calculate/diagnose
  lifecycle instead of being forced into an artificial generic `run` contract.
- Shared phonon input data are exposed once through `quantas.api.common`, avoiding
  duplicate public identities in HA and QHA documentation.
- Python API stability follows the Quantas package version; HDF5 and structured
  text formats retain independent schema versions.

### Removed

- The historical flat `quantas.api` module and its large unscoped symbol list.
- The root-level `quantas.eos` exception.
- Per-module public API version constants and unsupported calculator/reader/
  exporter re-exports.

### Compatibility

This is an intentional pre-release API break. No compatibility aliases are kept.
Scientific algorithms, numerical payloads, HDF5 schemas, units, tensor conventions,
and validated tolerances are unchanged.

## [2.0.0b1] - 2026-07-21

### Added

- Reusable Python-library architecture for HA, QHA, Elasticity, SEISMIC, EOS, and
  quasi-static thermoelasticity.
- Native HDF5 result envelopes with numerical precision, workflow provenance,
  warnings, events, diagnostics, and structured module payloads.
- Frontend-neutral reports and plot specifications with Click/Rich CLI adapters.
- Structured citation registry and repository-level `CITATION.cff`.
- Curated real-data examples and an isolated scientific-regression test suite.
- Deterministic multi-process validation runner and distribution smoke tests.

### Changed

- Unified calculation commands around grouped scientific, numerical, unit, and
  output/reporting options.
- Standardized `-o/--output`, `-r/--report`, `-v/--verbosity`, `-q/--quiet`,
  `-f/--force`, and progress controls.
- Significant calculations generate a deterministic `.log` report by default.
- ODRPACK is a required runtime dependency rather than a development-only backend.
- Large QHA report/export and EOS model/specification files were decomposed behind
  stable module facades.
- Real CRYSTAL and experimental datasets are shared between examples and integration
  tests rather than duplicated in module-specific test directories.

### Removed

- Pre-release compatibility aliases for superseded CLI option names.
- Internal phase reports, transient validation artifacts, and refactoring diaries
  from the public source tree.

### Scientific compatibility

No intentional changes were made to established formulas, array shapes, floating
precision, tensor conventions, HDF5 numerical payloads, or validated tolerances during
the Quantas 2 beta cleanup.  One EOS input enhancement recognizes absolute molar-volume
units declared through the historical `VSCALE` keyword.

[2.0.0b7]: https://github.com/gfulian/quantas/releases/tag/v2.0.0b7
[2.0.0b6]: https://github.com/gfulian/quantas/releases/tag/v2.0.0b6
[2.0.0b5]: https://github.com/gfulian/quantas/releases/tag/v2.0.0b5
[2.0.0b4]: https://github.com/gfulian/quantas/releases/tag/v2.0.0b4
[2.0.0b3]: https://github.com/gfulian/quantas/releases/tag/v2.0.0b3
[2.0.0b2]: https://github.com/gfulian/quantas/releases/tag/v2.0.0b2
[2.0.0b1]: https://github.com/gfulian/quantas/releases/tag/v2.0.0b1
