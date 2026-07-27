EOS contracts and model specifications
======================================

This page documents the passive data, model, parameter, solver, report, and
plot contracts used by the EOS operations.  Construct these objects through
:mod:`quantas.api.eos`; numerical evaluators and solver implementations remain
internal.

Dataset and domain capabilities
-------------------------------

.. autoclass:: quantas.api.eos.Dataset
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.FitDomain
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.CapabilityStatus
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.DomainCapability
   :members:
   :show-inheritance:

.. autodata:: quantas.api.eos.DOMAIN_CAPABILITIES

Model and parameter contracts
-----------------------------

.. autoclass:: quantas.api.eos.PVTModel
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.MGDNormalization
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.ParameterConstraint
   :members:
   :show-inheritance:

Fit request and result
----------------------

.. autoclass:: quantas.api.eos.FitRequest
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.FitOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.FitResult
   :members:
   :show-inheritance:

Regression selectors and solver options
---------------------------------------

``SolverOptions`` is the public union of the four supported passive solver
option contracts.  These objects select and configure a backend; they are not
numerical solver instances.

.. autoclass:: quantas.api.eos.FitMethod
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.CovarianceScaling
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.ODRDifferenceScheme
   :members:
   :show-inheritance:

.. autodata:: quantas.api.eos.SolverOptions

.. autoclass:: quantas.api.eos.OLSOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.WLSOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.EffectiveVarianceOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.ODROptions
   :members:
   :show-inheritance:

Archive slot, report, and plot contracts
----------------------------------------

.. autoclass:: quantas.api.eos.ResultSlot
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.ReportDetail
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.ReportOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.PlotOptions
   :members:
   :show-inheritance:

Session-aware plot discovery contracts
--------------------------------------

These lightweight immutable descriptors expose EOS archive, dataset, result-slot,
and record availability without retaining HDF5 handles or numerical fit arrays.
A detailed :class:`~quantas.api.plotting.PlotInventory` is attached only after
an explicit record, accepted slot, or unique accepted record has been selected.

.. autoclass:: quantas.api.eos.DatasetPlotDescriptor
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.RecordPlotDescriptor
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.SlotPlotDescriptor
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.ArchivePlotInventory
   :members:
   :show-inheritance:
