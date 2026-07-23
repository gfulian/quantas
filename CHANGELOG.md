# Changelog

All notable public changes to Quantas are recorded here.  The project follows
Semantic Versioning after the first stable Quantas 2 release.  During the current
beta, breaking changes are permitted when they simplify and stabilize the final
public contract; they must still be documented and validated.

## [2.0.0b6] - Unreleased

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

[2.0.0b6]: https://github.com/gfulian/quantas/releases/tag/v2.0.0b6
[2.0.0b5]: https://github.com/gfulian/quantas/releases/tag/v2.0.0b5
[2.0.0b4]: https://github.com/gfulian/quantas/releases/tag/v2.0.0b4
[2.0.0b3]: https://github.com/gfulian/quantas/releases/tag/v2.0.0b3
[2.0.0b2]: https://github.com/gfulian/quantas/releases/tag/v2.0.0b2
[2.0.0b1]: https://github.com/gfulian/quantas/releases/tag/v2.0.0b1
