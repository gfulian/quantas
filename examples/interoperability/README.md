# Interoperability example

This directory contains equivalent public CLI and Python API workflows for the
currently supported cross-module path:

```
QHA -> Thermoelasticity -> Elasticity / SEISMIC
```

Run from the project root:

```bash
bash examples/interoperability/workflow_cli.sh
python examples/interoperability/workflow_api.py
```

Both examples use the curated dolomite QHA and elastic-volume datasets.  The
CLI communicates through native HDF5 and shared text inputs; the Python example
keeps typed result contracts in memory while also writing the same reusable
artifacts.
