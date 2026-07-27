# Quantas 2 roadmap

## Source-freeze status

The scientific source tree, CLI contracts, persistence layer, curated examples,
and staged test runner remain frozen for the Quantas 2.0 release-candidate
line.  The public Python API is temporarily reopened for the ``2.0.0b7``
plotting-contract stabilization required to make the already accepted
frontend-neutral architecture usable by the CLI, Quantas GUI, notebooks, and
advanced library clients.  This exception is limited to public discovery,
plot specifications, scientifically necessary plot builders, and their tests
and documentation; it must not change approved numerical baselines or HDF5
schemas.

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
- Complete and validate the ``2.0.0b7`` public plotting-contract pass.
- Validate exact-grid HA V--T and QHA P--T section construction through API,
  CLI, and Quantas GUI adapters without changing stored scientific arrays;
  alternative axes require at least two unique native coordinates.
- Validate cumulative Thermoelasticity plot discovery for calibration, P--T,
  point/grid, profile, comparison, and domain archives through the same public
  API used by CLI and Quantas GUI.
- Validate the separate EOS session/archive plotting inventory through Python,
  CLI, and Quantas GUI without forcing EOS into the one-shot module contract.
- Review the resulting public API inventory and HDF5 schema compatibility one final time.

## After Quantas 2.0

- Add Kieffer acoustic thermodynamics to HA/QHA after a dedicated formula audit.
- Extend the standalone EOS workflow, including coupled P-V-T diagnostics.
- Develop additional code interfaces and scientific modules behind the same public
  API and capability registry.
- Continue Quantas GUI as an independent frontend over ``quantas.api`` without
  duplicating numerical logic; backend and GUI release milestones remain separate.
