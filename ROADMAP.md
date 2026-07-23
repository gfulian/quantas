# Quantas 2 roadmap

## Source-freeze status

The scientific source tree, public Python API, CLI contracts, persistence layer,
curated examples, and staged test runner are frozen for the Quantas 2.0 release
candidate line. Changes before the first release candidate are limited to:

- corrections required by final validation;
- manual-style documentation and tutorials;
- release metadata and external publishing configuration;
- narrowly scoped fixes that preserve the approved numerical baselines.

## Before 2.0.0rc1

- Complete the formal scientific validation matrix for every public workflow.
- Finish the manual-style documentation, CLI reference, API guide, and tutorials.
- Exercise CI on all supported operating systems and Python versions.
- Verify documentation hosting and the complete TestPyPI installation workflow.
- Review the frozen public API inventory and HDF5 schema compatibility one final time.

## After Quantas 2.0

- Add Kieffer acoustic thermodynamics to HA/QHA after a dedicated formula audit.
- Extend the standalone EOS workflow, including coupled P-V-T diagnostics.
- Develop additional code interfaces and scientific modules behind the same public
  API and capability registry.
- Build a GUI as a frontend over ``quantas.api`` without duplicating numerical logic.
