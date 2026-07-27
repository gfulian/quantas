Thermoelasticity contracts and selectors
========================================

This page documents the passive calibration, analysis, profile, and plotting
contracts.  Literal selectors are listed explicitly so type-checkers and
frontends can expose the same accepted values.

Scientific selectors
--------------------

.. autodata:: quantas.api.thermoelasticity.InputInterface

.. autoclass:: quantas.api.thermoelasticity.FitMethod
   :members:
   :show-inheritance:

.. autodata:: quantas.api.thermoelasticity.AdiabaticMode

.. autodata:: quantas.api.thermoelasticity.ExtrapolationPolicy

.. autodata:: quantas.api.thermoelasticity.FitFailurePolicy

.. autodata:: quantas.api.thermoelasticity.QualityPolicy

.. autodata:: quantas.api.thermoelasticity.StabilityPolicy

.. autodata:: quantas.api.thermoelasticity.TensorCondition

.. autodata:: quantas.api.thermoelasticity.ReportLevel

Reusable input structures
-------------------------

These aliases close the public annotations of :class:`Input` and
:class:`Context`.  They are passive scientific contracts, not calculator or
parser implementations.

.. autoclass:: quantas.api.thermoelasticity.CrystalStructure
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.SymmetryMetadata
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.ElasticVolumePoint
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.ElasticVolumeSeries
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.QHAInput
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.QHAOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.QHAThermoelasticPayload
   :members:
   :show-inheritance:

Calibration and result contracts
--------------------------------

.. autoclass:: quantas.api.thermoelasticity.Input
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.Context
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.Options
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.Result
   :members:
   :show-inheritance:

Profile contracts
-----------------

.. autoclass:: quantas.api.thermoelasticity.ProfilePreset
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.DepthProfile
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.ProfileResult
   :members:
   :show-inheritance:

Plot selectors
--------------

.. autodata:: quantas.api.thermoelasticity.ComponentGroup

.. autodata:: quantas.api.thermoelasticity.PTQuantity

.. autodata:: quantas.api.thermoelasticity.PlotLayout

.. autodata:: quantas.api.thermoelasticity.PlotPreset

.. autodata:: quantas.api.thermoelasticity.ProfileBackground

.. autodata:: quantas.api.thermoelasticity.ProfileColor

.. autodata:: quantas.api.thermoelasticity.ProfileMode

.. autodata:: quantas.api.thermoelasticity.UncertaintyMode

Plot option contracts
---------------------

.. autoclass:: quantas.api.thermoelasticity.PlotStyleOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.FitPlotOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.PTPlotOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.ProfilePlotOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.ComparePlotOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.thermoelasticity.DomainPlotOptions
   :members:
   :show-inheritance:
