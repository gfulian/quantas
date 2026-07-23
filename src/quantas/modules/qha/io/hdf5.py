# -*- coding: utf-8 -*-

"""HDF5 readers for quasi-harmonic approximation results.

The reader defined here reconstructs the generic Quantas result container and
its module-specific :class:`quantas.modules.qha.models.QHAResult` payload from
the native QHA HDF5 layout.  The same object can then be used by API code, CLI
exporters, graphical viewers, and plotting functions.
"""

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
from quantas.models import (
    BasicReader,
    ResultData,
    RESULT_SCHEMA_VERSION,
    validate_result_schema,
)
from quantas.modules.qha.io.hdf5_payload import (
    migrate_schema_1_0_thermodynamic_units,
    read_historical_qha_warnings,
    read_qha_payload,
)


class QHAHDF5Reader(BasicReader):
    """Reader for Quantas QHA HDF5 result files."""

    def load(self, filename: str | Path) -> ResultData:
        """Load a QHA HDF5 result file.

        Parameters
        ----------
        filename : str or Path
            Path to the QHA HDF5 result file.

        Returns
        -------
        ResultData
            Generic Quantas result containing a :class:`QHAResult` under the
            ``"qha"`` key.

        Raises
        ------
        OSError
            If the file cannot be opened by HDF5.
        ValueError
            If the file does not contain a QHA result layout.
        """
        path = Path(filename)
        try:
            with h5py.File(path, "r") as h5:
                result = _read_result_data(h5, path)
        except Exception as exc:
            self.completed = False
            self.error = str(exc)
            raise

        self.completed = True
        self.error = None
        return result


class QHAHDF5ResultReader(QHAHDF5Reader):
    """Alias reader kept for explicit result-loading code paths."""


def read_qha_hdf5(filename: str | Path) -> ResultData:
    """Read a QHA HDF5 result file.

    Parameters
    ----------
    filename : str or Path
        Path to the HDF5 result file.

    Returns
    -------
    ResultData
        Generic Quantas result containing the QHA payload.
    """
    return QHAHDF5Reader().load(filename)


def _read_result_data(h5: h5py.File, source: Path) -> ResultData:
    """Read a generic result container from an open HDF5 file.

    Parameters
    ----------
    h5 : h5py.File
        Open HDF5 file.
    source : Path
        Source file path.

    Returns
    -------
    ResultData
        Reconstructed generic result.

    Raises
    ------
    ValueError
        If required QHA groups are missing.
    """
    if "results" not in h5:
        raise ValueError("HDF5 file does not contain a 'results' group")

    metadata = read_result_metadata(h5)
    if metadata.schema_version == RESULT_SCHEMA_VERSION:
        validate_result_schema(h5, expected_module="qha")
    input_data = read_input_data(h5, source)
    options = read_options(h5)
    qha_result = read_qha_payload(h5["results"])
    warnings = read_warnings(h5) or read_historical_qha_warnings(h5)
    migration_warning = migrate_schema_1_0_thermodynamic_units(
        qha_result,
        metadata.schema_version,
        input_data,
    )
    if migration_warning is not None:
        warnings.append(migration_warning)
        metadata.schema_version = "1.1"

    return ResultData(
        metadata=metadata,
        input_data=input_data,
        options=options,
        results={"qha": qha_result},
        warnings=warnings,
        events=read_events(h5),
    )
