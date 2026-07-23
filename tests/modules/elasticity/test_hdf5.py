"""HDF5 precision tests for elasticity results."""

from __future__ import annotations

import h5py
import numpy as np


def test_hdf5_stores_floating_results_in_double_precision(tmp_path) -> None:
    """Elasticity calculations and native persistence both use float64."""
    from quantas.modules.elasticity.api import (
        read_elasticity_input,
        run_elasticity,
        write_elasticity_hdf5,
    )

    result = run_elasticity(
        read_elasticity_input("tests/modules/elasticity/data/hydroxylapatite.dat")
    )
    payload = result.results["elasticity"]
    assert payload.stiffness.dtype == np.dtype("float64")

    filename = tmp_path / "double_precision.hdf5"
    write_elasticity_hdf5(result, filename)

    with h5py.File(filename, "r") as h5:
        assert h5["results/stiffness"].dtype == np.dtype("float64")
        assert h5["results/compliance"].dtype == np.dtype("float64")
        assert h5["metadata/numerics"].attrs["working"] == "double"
        assert h5["metadata/numerics"].attrs["storage"] == "double"
