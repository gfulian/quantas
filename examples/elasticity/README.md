# Elasticity examples

The curated calcite and hydroxylapatite tensors can be processed through the
Python API or the CLI.

- `tutorial_api.py` reproduces the detailed calcite Elasticity tutorial,
  including the HDF5 result, deterministic report, and selected 2D/3D plots.
- `run_api.py` is a compact calculation and HDF5 example.
- `plot_3d_api.py` is a compact surface-construction example.

`hydroxylapatite.dat` also supplies the density required by the SEISMIC
tutorial, allowing both modules to use one physical input without duplicating
scientific data.
