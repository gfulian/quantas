# Quantas native HDF5 inspection examples

- `inspect_hdf5.py` prints groups, datasets, shapes, dtypes, units, descriptions,
  compression, and attributes without loading large arrays.
- `extract_hdf5.py` extracts selected dataset paths (or all datasets) to a
  compressed NumPy NPZ file and writes a JSON manifest preserving source paths,
  units, shapes, dtypes, root attributes, and the complete group-attribute tree.

Both scripts use only public file semantics and `h5py`; they do not import
Quantas implementation modules.
