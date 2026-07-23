# -*- coding: utf-8 -*-

"""Reader for native quasi-static thermoelastic HDF5 files."""

from __future__ import annotations

from pathlib import Path

import h5py

from quantas.io.hdf5 import (
    read_events,
    read_input_data,
    read_options,
    read_result_metadata,
    read_warnings,
)
from quantas.models import ResultData
from quantas.modules.thermoelasticity.io.hdf5_payload import (
    read_thermoelastic_payload,
)


def read_thermoelastic_hdf5(filename: str | Path) -> ResultData:
    """Read a native thermoelasticity result file."""
    path = Path(filename)
    with h5py.File(path, "r") as h5:
        metadata = read_result_metadata(h5)
        if metadata.module != "thermoelasticity":
            raise ValueError("HDF5 file is not a thermoelasticity result")
        if "results" not in h5:
            raise ValueError("HDF5 file does not contain a results group")
        payload = read_thermoelastic_payload(h5["results"])
        return ResultData(
            metadata=metadata,
            input_data=read_input_data(h5, path),
            options=read_options(h5),
            results={"thermoelasticity": payload},
            warnings=read_warnings(h5),
            events=read_events(h5),
        )


__all__ = ["read_thermoelastic_hdf5"]
