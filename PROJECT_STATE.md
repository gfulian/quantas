# Quantas project state

This file is a practical snapshot of where Quantas 2 currently stands.  It is
not a changelog and it is not a long-term roadmap: it records the scientific
and engineering baseline that a developer should assume when starting work on
the current tree.

Historical public changes belong in `CHANGELOG.md`.  Planned work belongs in
`ROADMAP.md`.  When either document temporarily lags behind an active feature
branch, this file and the current tested source tree describe the operational
state.

## State metadata

| Item | Current value |
|---|---|
| Last updated | 2026-09-09 |
| Current development version | `2.0.0b10` |
| Stable development baseline | `2.0.0b9`, `dev/refactor` |
| Active scientific branch | `dev/kieffer` |
| Current focus | Kieffer acoustic thermodynamics, CRYSTAL phonon continuity, and quasi-static thermoelastic validation |
| Development status | Pre-RC scientific closure and validation |
| Numerical precision | `float64` for real calculations and native HDF5 values; `complex128` for complex quantities |
| Persistence | Native HDF5 envelope retained; HA/QHA and Thermoelasticity payloads have been extended additively with Kieffer and pressure/provenance data |

## Source hierarchy

When current behavior or architecture is unclear, use sources in this order:

1. the current tested Quantas source tree, including the active scientific
   branch before merge;
2. this `PROJECT_STATE.md` file;
3. current tests and scientific validation material;
4. `CHANGELOG.md`, the manual, and `ROADMAP.md`;
5. Quantas `0.9.1`, only as a legacy behavioral or format reference.

Legacy Quantas is useful for characterization and historical naming, but it is
not the architectural model for Quantas 2.

## What Quantas 2 already is

Quantas is no longer a command-line program with a library wrapped around it.
The scientific implementation is shared by all frontends:

- `quantas.core` contains physics, mathematics, chemistry, geometry, events,
  and numerical building blocks and has no frontend dependencies;
- `quantas.models` contains calculators, readers, exporters, and passive data
  contracts;
- `quantas.interfaces` owns code-specific parsing and generation for CRYSTAL,
  VASP, Phonopy, and future backends;
- `quantas.modules` contains the scientific workflows: Elasticity, SEISMIC,
  HA, QHA, Thermoelasticity, and EOS;
- `quantas.io` owns shared HDF5 and persistence infrastructure;
- `quantas.renderers` owns frontend-neutral tables and plotting
  specifications;
- `quantas.cli` is an adapter layer for Click and Rich.

The Python API and CLI therefore execute the same calculators and numerical
code.  A future GUI is expected to consume the same public API rather than
reimplement scientific logic.

The native scientific result format remains HDF5.  Results retain the inputs,
options, units, warnings, meaningful workflow events, numerical precision, and
scientific provenance needed to understand how the result was produced.

## Mature baseline before the current branch

The `2.0.0b9` baseline already provides:

- one supported Python surface under `quantas.api`;
- typed public inputs, options, results, reports, plot specifications, events,
  and capability discovery;
- public input generation, execution, reopening, reports, plots, and derived
  exports for supported workflows;
- frontend-neutral table and plot descriptions used by both API and CLI;
- mature Elasticity and SEISMIC workflows, including directional properties,
  rotations, stability, Christoffel analysis, phase/group velocities,
  polarization tracking, degeneracy handling, plotting, and export;
- HA and QHA workflows with native HDF5 persistence and CLI/API equivalence;
- quasi-static thermoelastic reconstruction from volume-dependent elastic
  tensors and QHA thermodynamics;
- PV, VT, and PVT EOS infrastructure with diagnostics and fitting;
- staged scientific test execution, static analysis, documentation checks,
  distribution builds, and cross-platform CI.

Python support remains 3.10 through 3.13 until the complete scientific
stack is validated on Python 3.14.

## What `2.0.0b10` / `dev/kieffer` adds

The current branch is substantially larger than the original
`dev/crystal-parser` task.  Kieffer acoustic thermodynamics was intentionally
implemented before the release candidate because it fills an important
scientific gap for large primitive cells where a converged phonon supercell or
full Brillouin-zone dispersion is prohibitively expensive.

### CRYSTAL phonons and QHA mode continuity

Quantas now parses CRYSTAL phonon eigenvectors for both real and complex
q-points and stores them as unit-norm mass-weighted vectors.  Independent
multi-volume phonon calculations can be reordered into physically continuous
branches using eigenvector overlaps rather than frequency sorting alone.

The continuity machinery includes:

- global one-to-one assignment for non-degenerate modes;
- subspace/SVD comparisons for degenerate manifolds;
- explicit ambiguity, caution, and unresolved states;
- leave-one-out frequency-path validation for weak overlaps;
- structured provenance and diagnostics in generated QHA inputs;
- support for CRYSTAL real `R(...)` q-points that legitimately print only
  `MODES IN PHASE`, while complex `C(...)` q-points still require the
  anti-phase component.

This work changes input normalization and branch identity, not the accepted
HA/QHA thermodynamic formulas.

### Kieffer acoustic thermodynamics

Quantas now implements a frontend-neutral Kieffer acoustic model for the three
long-wavelength acoustic branches.  The workflow derives direction-dependent
phase velocities from an incremental elastic tensor and density, averages
`v^-3` over the sphere, converts the effective velocities to sine-dispersion
cutoffs, and evaluates the acoustic zero-point and thermal contributions.

Kieffer is deliberately restricted to primitive Gamma-only phonon inputs:

- exactly one q-point;
- `q = (0, 0, 0)`;
- identity phonon supercell;
- no explicit phonon dispersion already representing the acoustic branches.

The contribution is opt-in.  `ha add-kieffer` / `qha add-kieffer` enrich a
normal HA/QHA input without making the base input depend on elasticity, and
`ha run --kieffer` / `qha run --kieffer` activate the stored acoustic
contribution explicitly.

The Kieffer contribution is kept separate in results and HDF5 even when it is
added to the total thermodynamics.  This makes the approximation inspectable
rather than hiding it inside the ordinary phonon sum.

### CRYSTAL energy and finite-pressure elastic semantics

The current branch also closes several issues that became visible while
building Kieffer and reusing the same elastic series for quasi-static
thermoelasticity.

CRYSTAL energy parsing now distinguishes the converged SCF energy from the
physical total energy.  When CRYSTAL prints total energies including `DISP`,
`GCP`, or both, Quantas treats the most complete printed total as authoritative
and keeps the components as provenance.  QHA, elastic E(V) reconstruction,
and future Energy-EOS work should reuse this same resolver.

For ELASTCON/ELAPIEZO outputs, structure, volume, density, energy, and output
stress are tied to the unstrained elastic reference state before the first
`STRAIN MATRIX`.  Later strained or internally relaxed geometries are not
allowed to masquerade as the reference lattice.

Raw CRYSTAL energy-strain tensors can be converted to the finite-pressure
incremental coefficients using the CRYSTAL/Erba hydrostatic transformation.
Tensors already corrected by `PRESSURE`/`PRESSEOS` are preserved.  Raw tensors
may instead obtain pressure from output stress, an explicit manual value, or a
static E(V) EOS/polynomial fit.  Quantas records the pressure source and the
correction decision and rejects an unsafe repeated correction.

The CRYSTAL-specific conversion is not assumed to be universal.  Other
backends must define the meaning of their printed elastic tensor before a
finite-pressure correction is applied.

## Scientific validation achieved on the current branch

The branch is covered by analytical, characterization, regression, invalid-
input, persistence, and end-to-end tests.  Two real-material validation tracks
are especially important.

### Hydroxylapatite: Gamma-only QHA + Kieffer + QSA

The OHAp benchmark exercises the intended Kieffer use case: a large primitive
cell with Gamma-only phonons and volume-dependent elastic tensors.  The
validation checks multi-volume cutoff construction, alternative pressure
sources, acoustic quadrature convergence, QHA execution with and without
Kieffer, HDF5 round trips, and coupling to quasi-static thermoelasticity.

The results are regular and physically interpretable.  Kieffer produces the
expected additional acoustic free-energy contribution, thermal expansion, and
thermoelastic feedback without introducing discontinuities, instability, or
spurious double pressure corrections.

### Periclase MgO: cubic control benchmark

MgO provides a deliberately simpler and more discriminating control because it
is cubic and a full 4x4x4 phonon dispersion is affordable.

The seven-volume elastic series was generated without CRYSTAL
`PRESSURE`/`PRESSEOS`, so Quantas receives raw energy-strain tensors and must
reconstruct the static pressure and apply the finite-pressure correction
itself.  With BM3 pressure from the same static E(V) series, the bulk modulus
of the corrected cubic tensor agrees with the EOS bulk modulus to roughly
0.3% over approximately -13 to +20 GPa.  The uncorrected raw tensor differs by
several percent at the ends of the range.  This is strong independent evidence
that the CRYSTAL pressure sign and finite-pressure transformation are being
applied correctly and once.

The full-dispersion phonons also provide a useful Kieffer control.  At fixed
volume the Kieffer acoustic contribution reproduces the three explicitly
sampled acoustic branches well at moderate and high temperature; larger
differences between complete Gamma+Kieffer and full-dispersion QHA are mainly
caused by the neglected dispersion of the optical branches, not by a failure
of the acoustic model itself.

MgO also clarifies the QSA/QHA comparison.  The QSA bulk modulus follows the
cold static EOS evaluated at the QHA equilibrium volume very closely, while
the full QHA bulk modulus increasingly separates at high temperature.  The
remaining difference is therefore largely the intrinsic vibrational-elastic
physics omitted by the quasi-static approximation, rather than evidence of a
broken cubic finite-strain reconstruction.

These validation datasets are valuable enough to become documented test cases
or tutorials later, but that packaging is not a prerequisite for finishing the
current scientific closure.

## Known scientific questions still open

The remaining questions are narrow and should be characterized before changing
numerical formulas.

### Gamma translations when Kieffer is active

For an exact periodic crystal the three acoustic modes at Gamma are rigid
translations and have exactly zero frequency.  Real calculations may print
small positive or negative residual frequencies because translational
invariance is not numerically exact.

The ordinary harmonic core currently excludes non-positive frequencies but
accepts every positive frequency.  Therefore a small *positive* numerical
Gamma translation could be treated as a finite harmonic oscillator.  If
Kieffer is then activated, the three continuous acoustic branches are added as
well and that translational degree of freedom would be counted twice.

This is a Kieffer-application issue, not a reason to delete modes from the
phonon dataset.  Quantas should continue to preserve `3N` branches at every
q-point for parsing, tracking, and provenance.  The next characterization step
is to identify the three translational Gamma modes from their eigenvectors and
ensure that, only when Kieffer supplies the acoustic branches, numerical
residuals of those Gamma translations do not enter the discrete harmonic sum.
Frequency magnitude alone must not be used, because a genuine soft optical
mode can also lie close to zero or become imaginary.

### Acoustic branch identity in strongly anisotropic media

The current spherical average uses the locally ordered Christoffel phase
velocities.  This is well behaved for the validated materials, but a severe
anisotropic characterization case should still compare speed ordering with
polarization-aware branch tracking before the Kieffer implementation is frozen.

### Energy EOS

The main scientific workflow still missing before the planned Quantas 2
feature set is a complete standalone Energy-EOS path.  It should not create a
new CRYSTAL energy parser or a separate pressure engine.  The work completed
for QHA, Kieffer, and QSA already provides the common ingredients: resolved
static total energies, volume/state provenance, EOS fitting, and P(V)
reconstruction.

## Engineering and quality state

The project keeps the scientific core independent from Click, Rich, plotting
libraries, and GUI objects.  Calculators may emit Quantas events; pure
numerical objects do not.  Long numerical loops use callbacks when progress is
needed.  Display precision remains a renderer concern and never changes stored
scientific values.

The repository has extensive automated coverage spanning scientific kernels,
interfaces, workflows, CLI/API equivalence, HDF5 round trips, architecture,
public API contracts, source hygiene, documentation, and distribution builds.
The canonical complete validation remains:

```text
python tools/run_tests.py all -- -q
```

with Ruff, mypy, Sphinx warning-as-error builds, package builds, `twine check`,
and clean-install smoke tests also required before release.

The active `dev/kieffer` branch must pass the complete repository CI before
merge.  Local typing and targeted validation are already part of the branch
closure; CI remains the final cross-platform gate rather than a substitute for
scientific validation.

## Immediate next operations

The short sequence from the current state is:

1. finish the Gamma-translation characterization for Kieffer without changing
   the shape or branch count of the stored phonon dataset;
2. add the MgO and OHAp evidence to the formal validation matrix and update the
   corresponding validation documentation;
3. perform the final strongly-anisotropic acoustic-branch characterization;
4. merge the validated Kieffer/CRYSTAL/QSA tranche back into `dev/refactor`;
5. implement the standalone Energy-EOS workflow using the shared static-energy
   and pressure-resolution infrastructure rather than duplicating it;
6. synchronize `ROADMAP.md`, schema/version notes, validation pages, and release
   documentation with the actual source tree;
7. run the complete clean-checkout and built-distribution release validation;
8. freeze the public API and native schemas for `2.0.0rc1`.

## Remaining blockers for `2.0.0rc1`

- close the narrow Kieffer characterization items listed above;
- complete the standalone Energy-EOS workflow and its validation;
- publish the scientific validation matrix and the still-incomplete workflow
  validation pages;
- synchronize project state, roadmap, changelog, schema/version notes, and
  release documentation;
- perform the final public API and HDF5 compatibility review;
- run complete cross-platform CI from a clean checkout and from built
  distributions;
- publish and reinstall the candidate through TestPyPI;
- synchronize the final `2.0.0rc1` version, citation metadata, changelog date,
  tag, GitHub release, and distribution hashes.

New scientific features beyond this list should be deferred unless validation
shows that an existing public result is incorrect or unusable.
