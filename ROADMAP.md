# Quantas 2 roadmap

## Current pre-RC development

The ``2.0.0b10`` Kieffer, CRYSTAL phonon-continuity, and generalized QSA
pressure/provenance tranche has been merged into ``dev/refactor`` with the CI
matrix green.  Development has moved to ``2.0.0b11`` on ``dev/energyeos``.

The remaining scientific feature before release-candidate closure is the
standalone Energy EOS workflow.  Work proceeds in small validated steps:

- complete parity between P(V) and volume-integrated E(V) model families;
- retain and validate the newly added SJEOS Energy EOS in its physical
  equilibrium-parameter form, including derived P(V) inspection;
- complete canonical energy-unit handling while preserving all existing unit
  aliases and conversions;
- promote ``ev/energy`` from shared core capability to public fitting,
  diagnostics, persistence, reporting, and plotting;
- expose derived ``P(V) = -dE/dV`` as an inspectable Energy EOS result;
- add ``eos inpgen`` from backend output series, beginning with CRYSTAL and
  then generalizing to VASP and future interfaces;
- keep QHA and Thermoelasticity coupled to the shared numerical Energy EOS
  service, not to the standalone EOS workflow layer.

The first b11 step completes the analytical integrated modified-Tait T2/T3/T4
forms and validates them against the existing pressure equations.

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

- Evaluate additional non-empirical Elasticity observables in a dedicated
  change.
- Develop additional code interfaces and scientific modules behind the same
  public API and capability registry.
- Continue Quantas GUI as an independent frontend over `quantas.api` without
  duplicating numerical logic.
