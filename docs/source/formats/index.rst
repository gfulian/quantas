Input and output formats
========================

This section is the normative data-contract reference for Quantas.  It defines
what each workflow reads, which values are normalized before calculation, what
is stored in native HDF5 files, and which information is preserved in reports
and tabular exports.

The pages are organized in three groups.

Input contracts
---------------

- :doc:`elasticity_input` describes the shared Elasticity/SEISMIC stiffness
  text format.
- :doc:`phonon_yaml` defines the normalized HA/QHA phonon YAML contract.
- :doc:`eos_input` and :doc:`eos_spec` separate EOS measurements from the
  fitting plan.
- :doc:`thermoelastic_input` defines the normalized volume-dependent elastic
  series used by the quasi-static thermoelastic workflow.
- :doc:`earth_profile_spec` defines tabulated and composed depth--pressure--
  temperature paths.

Native result contracts
-----------------------

- :doc:`hdf5` defines the common Quantas envelope.
- :doc:`ha_qha_hdf5` documents HA and QHA payloads.
- :doc:`elasticity_seismic_hdf5` documents Elasticity and SEISMIC payloads.
- :doc:`eos_hdf5` documents the session-capable EOS archive.
- :doc:`thermoelastic_hdf5` documents calibration, grid, and profile payloads.
- :doc:`hdf5_inspection` shows how to inspect and extract native files from
  Python without relying on terminal output.

Human-readable outputs
----------------------

:doc:`tabular_outputs` explains deterministic reports, CSV and text exports,
plot files, units, missing values, and the distinction between stored precision
and displayed precision.

General rules
-------------

- File extensions are conveniences, not scientific type identifiers. Native
  HDF5 readers inspect metadata; text readers inspect the declared format.
- Numerical input is normalized to ``float64`` before scientific calculation.
- Unit labels and normalization bases are part of the data contract and must be
  read together with numerical arrays.
- Native HDF5 is the authoritative machine-readable result. Reports, CSV files,
  and figures are derived views.
- Unknown schema versions or incompatible units must be rejected rather than
  guessed.
- Historical pre-release files are read only through explicit migration paths
  implemented by the responsible module.
