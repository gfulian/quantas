# Quantas roadmap

Quantas 2 is currently in beta. The main scientific architecture and public
workflows are in place, and current development is increasingly focused on
scientific validation, interoperability, documentation, and preparation for the
first stable 2.0 release.

This roadmap describes the direction of the project rather than acting as a
record of completed development. Implemented changes are documented in
[CHANGELOG.md](CHANGELOG.md).

## Towards Quantas 2.0

The main objective before the first stable Quantas 2 release is to complete the
scientific and release validation of the existing public workflows.

This includes completing the validation matrix for Elasticity, SEISMIC,
harmonic and quasi-harmonic thermodynamics, Kieffer acoustic thermodynamics,
Thermoelasticity, and the EOS workflows; documenting the corresponding
reference datasets, methods, units, tolerances, and known limitations; and
finishing the remaining user and developer documentation.

The public `quantas.api` surface, native HDF5 schemas, scientific defaults, and
interoperability contracts will be reviewed before the release candidate. Once
the release-candidate series begins, changes to these contracts should be
limited to corrections required for scientific or technical reliability.

Release candidates will also be validated from clean source archives and built
distributions on the supported Python and operating-system combinations.

## Scientific development

Further development will continue behind the common library/API architecture
rather than introducing frontend-specific scientific implementations.

Important areas for continued work include more complete thermoelastic
treatments beyond the quasi-static approximation, further development and
validation of equation-of-state workflows, and broader interoperability with
external electronic-structure and lattice-dynamics codes.

Phonon interoperability is expected to expand beyond the current direct
CRYSTAL and VASP paths, including workflows in which phonon data are produced
by tools such as Phonopy. New interfaces should reuse the backend-neutral
structure, phonon, elasticity, and provenance contracts already present in
Quantas.

## Validation and reference data

Scientific validation remains a continuing part of development rather than a
one-time release task.

Where possible, Quantas results will be compared with analytical limits,
independent implementations, experimental data, published reference values,
and representative real-material calculations.

Validation datasets that are sufficiently small, reproducible, and useful for
future regression testing may become curated examples or documented reference
cases.

## Documentation and usability

The scientific manual will continue to grow alongside the software, with an
emphasis on reproducible workflows rather than command listings alone.

Tutorials and examples should show both the command-line and Python API routes
when this helps demonstrate that they use the same scientific backend.

Quantas GUI remains a separate frontend project. It is expected to consume the
public Quantas API and frontend-neutral result, table, and plotting contracts
without introducing a second numerical implementation.

## Beyond 2.0

After the first stable release, new scientific capabilities can resume without
the constraints of the pre-release freeze.

Likely areas include additional thermoelastic and equation-of-state methods,
new computational-code interfaces, expanded phonon workflows, and scientific
observables that can be built cleanly on the existing core architecture.

The exact order will be determined by scientific need, validation requirements,
and use in ongoing research rather than by a fixed release calendar.