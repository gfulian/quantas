# Quantas 2 roadmap

## Current pre-RC development

The standalone Energy EOS work completed the ``2.0.0b11`` scientific tranche.
Development is now on ``2.0.0b12`` in ``dev/prerelease-hardening``.  This branch
adds no new scientific workflow: it is closing architecture, documentation,
validation, packaging, and release-readiness gaps before the next pre-release
checkpoint.

The b12 hardening tranche has:

- separated documented EOS request errors from unexpected failures;
- unified energy-derived pressure assignment and provenance across Kieffer and
  thermoelastic workflows;
- clarified CRYSTAL finite-prestress conversion versus the Wallace QSA
  finite-strain term;
- completed the public docstring audit and active-object cleanup;
- normalized citation/bibliography handling and simplified documentation
  information ownership;
- made the validation record explicit about completed and work-in-progress
  scopes;
- prepared synchronized release metadata and final distribution checks.

VASP/Phonopy cleanup and common MgO cross-backend characterization are planned
for a small b13 follow-up rather than being mixed into b12.

The completed b11 Energy EOS tranche:

- completed parity between the registered P(V) families and their integrated
  E(V) forms where scientifically defined, including modified Tait;
- added and validated SJEOS in its physical equilibrium-parameter form;
- centralized EOS model discovery through the shared resolver, compact
  historical tags, ``quantas eos show-models``, and shell completion including
  native PowerShell support;
- added canonical Hartree normalization for energy and ``sigma_energy`` while
  retaining raw units and provenance;
- promoted ``ev/energy`` to public fitting, diagnostics, HDF5 persistence,
  reporting, plotting, and post-fit calculation;
- exposed ``P(V) = -dE/dV`` together with ``K(V)``, ``K'(V)``, and ``K''(V)``
  from accepted Energy EOS records;
- added CRYSTAL ``eos inpgen`` support for static, optimized, and native
  multi-volume EOS outputs through ``StructureEnergySeries``;
- added explicit primitive/crystallographic reference-cell normalization,
  symmetry provenance, and a shared structural-path response for theoretical
  E--V datasets;
- reused the QHA ``StructuralPathModel`` to derive equilibrium axes,
  ``d ln(l_i)/d ln(V)``, and axial moduli ``M_i = K/eta_i`` with propagated
  uncertainty rather than introducing a second EOS-specific lattice model;
- added an optional, explicitly secondary pressure-form axial EOS refit from
  EnergyEOS-derived pressures, including retained pressure covariance and clear
  diagonal-WLS provenance;
- kept QHA and Thermoelasticity coupled to the shared numerical Energy EOS
  service rather than to the standalone workflow layer.

The remaining pre-RC decisions are deliberately narrow: complete the final
release-readiness gate, keep unfinished validation scopes labelled explicitly,
and freeze the public API and native archive contracts before ``2.0.0rc1``.

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
