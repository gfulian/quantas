# -*- coding: utf-8 -*-

"""Native HDF5 exporter for quasi-static thermoelasticity."""

from __future__ import annotations

from pathlib import Path

import h5py

from quantas.io.hdf5 import (
    write_diagnostics,
    write_events,
    write_input_data,
    write_options,
    write_precision_metadata,
    write_result_metadata,
)
from quantas.io.path import ensure_suffix
from quantas.models import BasicHDF5Export, ResultData
from quantas.modules.thermoelasticity.io.hdf5_payload import (
    write_thermoelastic_payload,
)
from quantas.modules.thermoelasticity.models import ThermoelasticResult


class ThermoelasticityHDF5Export(BasicHDF5Export):
    """Write a complete native thermoelastic HDF5 result."""

    def export(
        self,
        result: ResultData,
        filename: str | Path,
        report_text: str | None = None,
    ) -> None:
        """Write generic envelope, scientific arrays, and all fit diagnostics."""
        payload = result.results.get("thermoelasticity")
        if not isinstance(payload, ThermoelasticResult):
            raise ValueError("result does not contain a thermoelasticity payload")
        path = ensure_suffix(filename, ".hdf5")
        with h5py.File(path, "w") as h5:
            write_result_metadata(h5, result)
            write_precision_metadata(h5)
            write_input_data(h5, result.input_data)
            write_options(h5, result.options)
            write_thermoelastic_payload(h5, payload)
            write_diagnostics(h5, result, report_text=report_text)
            write_events(h5, result.events)
        self.completed = True


__all__ = ["ThermoelasticityHDF5Export"]
