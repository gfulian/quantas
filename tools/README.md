# Quantas developer tools

This directory contains the small command-line utilities used to maintain and
verify the Quantas source tree. They are development tools, not part of the
public Python API or the installed ``quantas`` command.

## Routine checks

- ``validate_release.sh`` recreates a dedicated virtual environment and runs
  the complete source, test, documentation, packaging, distribution, and Git
  validation sequence. It is the preferred final check for an unpacked source
  tree or release candidate.
- ``run_tests.py`` runs named pytest suites. The ``all`` target isolates the
  shared core, each scientific module, the CLI, plotting, and curated examples
  in separate processes.
- ``check_architecture.py`` checks package boundaries, module contracts, source
  size, and other maintainability signals.
- ``check_distribution.py`` builds and inspects wheel and source archives and
  can verify clean installation.

## Project inspection

- ``project_stats.py`` writes a concise JSON inventory of source files, tests,
  public modules, and cross-module dependencies.

## Reference maintenance

- ``update_examples_manifest.py`` refreshes the deterministic manifest for the
  curated example data.
- ``update_scientific_reference.py`` regenerates the approved cross-module
  numerical reference bundle.
- ``update_seismic_reference.py`` regenerates the independent SEISMIC formula
  reference data.
- ``compare_qha_workflows.py`` compares selected QHA workflows during numerical
  maintenance.

## Benchmarks

- ``benchmark_elasticity_sampling.py`` measures directional Elasticity sampling.
- ``benchmark_seismic_vectorization.py`` compares vectorized and pointwise
  SEISMIC kernels.

Generated reports are written to ``project_internal/checks/`` by default. That
directory is local developer state and is excluded from public distributions.
