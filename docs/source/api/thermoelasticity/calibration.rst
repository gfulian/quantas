Thermoelasticity calibration API
================================

Calibration combines a normalized elastic-volume series with a QHA result,
fits the reference EOS and independent stiffness components, and returns a
reusable :class:`quantas.api.common.ResultData` envelope.

Input construction and normalization
------------------------------------

.. autofunction:: quantas.api.thermoelasticity.create_input

.. autofunction:: quantas.api.thermoelasticity.read_input

.. autofunction:: quantas.api.thermoelasticity.normalize_input

QHA coupling context
--------------------

.. autofunction:: quantas.api.thermoelasticity.prepare_context

Calculation
-----------

.. autofunction:: quantas.api.thermoelasticity.run

.. autofunction:: quantas.api.thermoelasticity.run_context

Typed result and persistence
----------------------------

.. autofunction:: quantas.api.thermoelasticity.get_result

.. autofunction:: quantas.api.thermoelasticity.write_result

.. autofunction:: quantas.api.thermoelasticity.read_result

Calibration reporting and default plots
---------------------------------------

.. autofunction:: quantas.api.thermoelasticity.build_report

.. autofunction:: quantas.api.thermoelasticity.build_fit_plots

.. autofunction:: quantas.api.thermoelasticity.build_plots
