"""Package-layout tests for the Quantas EOS workflow foundation."""

from __future__ import annotations

from quantas.modules import eos


def test_public_eos_workflow_foundation_is_available() -> None:
    assert eos.EOSDataset is not None
    assert eos.EOSSeries is not None
    assert eos.EOSInputFileReader is not None
    assert eos.EOSFitter is not None
    assert eos.EOSArchive is not None
    assert eos.EOSArchiveInspection is not None
    assert eos.EOSFitRecord is not None
    assert eos.EOSSession is not None
    assert eos.EOSPlotOptions is not None
    assert eos.EOSPlotter is not None
    assert eos.PressureEOSFitModel is not None
    assert eos.build_pressure_parameter_map is not None
    assert eos.read_eos_input is not None
    assert callable(eos.odr_backend_available)
