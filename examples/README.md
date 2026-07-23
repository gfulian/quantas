# Quantas examples

This directory contains curated, real-world inputs used both as user examples and as
integration/scientific-regression fixtures.  The files exercise complete public
workflows and external-code readers; small synthetic arrays remain in unit tests when
an isolated formula, boundary condition, or failure mode must be tested directly.

## Layout

- `elasticity/`: second-order elastic tensors for calcite and hydroxylapatite.
- `seismic/`: public-API tutorial scripts using the curated hydroxylapatite tensor.
- `eos/`: experimental P-V, V-T, and P-V-T datasets, strict batch
  specifications, and public-API tutorial scripts.
- `ha/`: harmonic-approximation API example.
- `qha/`: CRYSTAL phonon-series and native CRYSTAL QHA examples.
- `thermoelasticity/`: quasi-static thermoelastic inputs, CRYSTAL SOEC outputs,
  geotherm/profile data, and API scripts.
- `io/`: format-inspection examples for traversing native Quantas HDF5 files
  and extracting selected datasets to NumPy archives without private imports.
- `interoperability/`: equivalent CLI and public-API chains for the supported
  QHA -> Thermoelasticity -> Elasticity/SEISMIC transformation.

All text inputs distributed by Quantas are UTF-8 with LF line endings.  Run
`python tools/update_examples_manifest.py --check` to verify the checked-in
manifest, file sizes, and SHA-256 hashes.

The examples are intentionally not imported by the Quantas runtime package.  They are
included in the source distribution and repository, while the wheel remains limited
to executable package code.

See `PROVENANCE.md` for scientific origin and usage notes.
