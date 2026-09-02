# QHA examples

- `crystal-qha/` demonstrates direct ingestion of a native CRYSTAL QHA
  calculation.
- `crystal-phonons/` demonstrates a volume series assembled from independent
  CRYSTAL phonon outputs.
- `run_api.py` is a compact public-API example.
- `tutorial_api.py` reproduces the detailed MgO QHA documentation tutorial,
  including input inspection, HDF5 output, reports, and selected plots.

Both source families include normalized Quantas YAML inputs.  The two routes
represent different continuity provenance:

- `crystal-qha/` preserves the mode-continuity verdict reported by the native
  CRYSTAL QHA workflow (`method: crystal-qha`).
- `crystal-phonons/` exercises Quantas' independent adjacent-volume
  eigenvector-overlap tracker, including degenerate-subspace and
  leave-one-out diagnostics.

The scientific details of input normalization and mode continuity are documented
in `docs/source/workflows/phonon_input_generation.rst`.
