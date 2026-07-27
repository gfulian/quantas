# Quantas 2 roadmap

## Source-freeze status

The ``2.0.0b7`` public-lifecycle stabilization has completed full Windows
validation.  The scientific source tree, CLI contracts, persistence layer,
curated examples, and staged test runner now return to the release-candidate
freeze while the feature branch is reviewed and merged.

The completed exception covered public input generation, public type closure,
operation discovery, derived exports, plot specifications and builders, and
their tests and documentation.  It did not change approved numerical baselines
or HDF5 schemas.

Other changes before the first release candidate are limited to:

- corrections required by final validation;
- manual-style documentation and tutorials;
- release metadata and external publishing configuration;
- narrowly scoped fixes that preserve the approved numerical baselines.

## Before 2.0.0rc1

- Complete the formal scientific validation matrix for every public workflow.
- Finish the manual-style documentation, CLI reference, API guide, and tutorials.
- Exercise CI on all supported operating systems and Python versions.
- Verify documentation hosting and the complete TestPyPI installation workflow.
- Merge the validated ``2.0.0b7`` public lifecycle-contract pass.
- Integrate the validated HA V--T and QHA P--T section controls into Quantas
  GUI without changing stored scientific arrays.
- Integrate cumulative Thermoelasticity discovery for calibration, P--T,
  point/grid, profile, comparison, and domain archives through the public API.
- Integrate the separate EOS session/archive inventory into Quantas GUI without
  forcing EOS into the one-shot module contract.
- Review the resulting public API inventory and HDF5 schema compatibility one final time.

Backend CLI/API validation is complete.  Quantas GUI adapter validation follows
on the independent GUI roadmap and does not reopen the backend scientific
implementation.

## After Quantas 2.0

- Add Kieffer acoustic thermodynamics to HA/QHA after a dedicated formula audit.
- Extend the standalone EOS workflow, including coupled P-V-T diagnostics.
- Develop additional code interfaces and scientific modules behind the same public
  API and capability registry.
- Continue Quantas GUI as an independent frontend over ``quantas.api`` without
  duplicating numerical logic; backend and GUI release milestones remain separate.
