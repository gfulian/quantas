# Quantas documentation

This directory contains the public Sphinx documentation source intended for the Quantas repository and Read the Docs.

## Install the documentation dependencies

From the repository root:

```bash
python -m pip install -e ".[docs]"
```

## Build on Linux or macOS

```bash
make -C docs html
```

## Build on Windows

From Command Prompt or PowerShell, while located in the repository root:

```bat
.\docs\make.bat html
```

The `html` argument may be omitted because it is the default target:

```bat
.\docs\make.bat
```

Other supported Sphinx targets include:

```bat
.\docs\make.bat clean
.\docs\make.bat linkcheck
```

The script uses the active Python interpreter. A different interpreter may be selected through the `PYTHON` environment variable, for example:

```bat
set PYTHON=py -3.13
.\docs\make.bat html
```

The build regenerates the Elasticity and SEISMIC tutorial figures from curated examples before Sphinx reads the sources. The generated HTML is written to `docs/_build/html`.

Maintainer-only refactoring notes, audits, and historical snapshots are stored in a separate development archive. Documentation tools may create an ignored local `project_internal/` workspace for generated audit output.
