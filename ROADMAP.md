# Quantas 2 roadmap

## Source-freeze status

The `2.0.0b9` public lifecycle baseline, SEISMIC input-generation work, and
frontend-neutral QHA inspection presentation are merged into `dev/refactor`.
The backend has returned to release-candidate freeze.

The completed stabilization covers:

- public input, execution, persistence, report, plot, and export lifecycles;
- typed public contracts and capability discovery;
- result-aware plot inventories;
- workflow-owned Elasticity and SEISMIC input generation;
- CLI/API equivalence;
- VASP density extraction for SEISMIC;
- rejection of unstable native SEISMIC results.
- public QHA inspection tables and sampled energy-volume plot specifications
  derived from the existing inspection result without repeating the fit.

It did not intentionally change approved numerical baselines, physical
conventions, scientific array layouts, or HDF5 schemas.

Before the first release candidate, changes are limited to:

- corrections demonstrated by final validation;
- validation and manual-style documentation;
- release metadata and publishing configuration;
- narrowly scoped compatibility fixes that preserve approved scientific
  behavior.

## Before 2.0.0rc1

- Complete the formal scientific validation matrix for every public workflow.
- Publish the corresponding validation pages, datasets, methods, units,
  tolerances, and limitations.
- Finish the manual-style documentation, CLI reference, API guide, and
  tutorials.
- Exercise the complete CI matrix on every supported operating system and
  Python version.
- Review the frozen `quantas.api` inventory and native HDF5 schemas one final
  time.
- Run the complete validation from a fresh checkout and from built
  distributions.
- Verify documentation hosting.
- Publish the candidate build to TestPyPI and reinstall it in clean Windows,
  Linux, and macOS environments.
- Resolve only release-blocking defects found by those checks.
- Synchronize version, citation, changelog, tag, release notes, and checksums
  for `2.0.0rc1`.

Quantas GUI integration continues on its independent roadmap. It may reveal a
real backend contract defect, but completion of GUI milestones is not a
prerequisite for the Quantas backend release candidate.

## Release-candidate policy

After `2.0.0rc1`, the public API, scientific defaults, units, and HDF5 schemas
are frozen. Further release candidates contain only:

- correctness fixes;
- portability and packaging fixes;
- documentation corrections;
- validation additions that do not silently redefine approved behavior.

Any newly proposed scientific capability is deferred unless omission would make
an existing public result incorrect or unusable.

## After Quantas 2.0

- Add Kieffer acoustic thermodynamics to HA/QHA after a dedicated formula audit.
- Extend the standalone EOS workflow, including coupled P--V--T diagnostics.
- Evaluate additional non-empirical Elasticity observables in a dedicated
  change.
- Develop additional code interfaces and scientific modules behind the same
  public API and capability registry.
- Continue Quantas GUI as an independent frontend over `quantas.api` without
  duplicating numerical logic.
