Quasi-Harmonic Approximation API
================================

:mod:`quantas.api.qha` exposes preflight inspection, quasi-harmonic
calculation, result validation, method comparison, plotting, reporting, and
native HDF5 persistence.

Recommended lifecycle
---------------------

.. code-block:: python

   from quantas.api import qha

   options = qha.Options(
       scheme="freq",
       minimization="poly",
       thermal_expansion_method="mixed_derivative",
       energy_degree=3,
       frequency_degree=3,
   )
   preview = qha.inspect(
       "mgo_b3lyp.yaml",
       options=options,
       polynomial_degree=3,
       eos="BM3",
   )
   result_data = qha.run("mgo_b3lyp.yaml", options=options)
   summary = qha.validate_result(result_data)
   result = qha.get_result(result_data)

Use :func:`inspect` before a production run to examine static volume support,
fit quality, and the pressure interval represented by the sampled energies.

Scientific option types
-----------------------

The public selectors accept the following literal values:

- ``Scheme``: ``freq`` or ``td``;
- ``Minimization``: ``poly`` or ``eos``;
- ``ThermalExpansionMethod``: ``mixed_derivative``, ``mode_gruneisen``, or ``numerical``;
- ``PolynomialDerivativeMethod``: ``local_grid`` or ``analytic``;
- ``ModeContinuity``: ``verified``, ``assumed``, ``unknown``, or ``unreliable``;
- ``FitFailurePolicy``: ``continue``, ``stop``, or ``raise``.


.. autodata:: quantas.api.qha.Scheme

.. autodata:: quantas.api.qha.Minimization

.. autodata:: quantas.api.qha.ThermalExpansionMethod

.. autodata:: quantas.api.qha.PolynomialDerivativeMethod

.. autodata:: quantas.api.qha.ModeContinuity

.. autodata:: quantas.api.qha.FitFailurePolicy

Passive contracts
-----------------

.. autoclass:: quantas.api.qha.Input
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.qha.Options
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.qha.PlotOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.qha.Result
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.qha.Preview
   :members:
   :show-inheritance:

Input, inspection, and calculation
----------------------------------

.. autofunction:: quantas.api.qha.read_input

.. autofunction:: quantas.api.qha.normalize_input

.. autofunction:: quantas.api.qha.inspect

.. autofunction:: quantas.api.qha.run

.. autofunction:: quantas.api.qha.get_result

Validation and comparison
-------------------------

``validate_result`` evaluates completeness and internal consistency of one
result.  ``compare_results`` compares two results on common P--T states and is
appropriate for method-sensitivity studies such as ``freq`` versus ``td`` or
polynomial versus EOS minimization.

.. autoclass:: quantas.api.qha.ValidationSummary
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.qha.PropertyDifference
   :members:
   :show-inheritance:

.. autofunction:: quantas.api.qha.validate_result

.. autofunction:: quantas.api.qha.compare_results

Reporting, plotting, and persistence
------------------------------------

.. autofunction:: quantas.api.qha.list_plot_properties

.. autofunction:: quantas.api.qha.build_report

.. autofunction:: quantas.api.qha.build_plots

.. autofunction:: quantas.api.qha.write_result

.. autofunction:: quantas.api.qha.read_result

See also
--------

- :doc:`../workflows/qha`
- :doc:`../tutorials/qha`
- :doc:`../formats/phonon_yaml`
- :doc:`../formats/ha_qha_hdf5`
