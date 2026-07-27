EOS reports, diagnostics, calculations, and plots
=================================================

Post-fit operations consume an immutable archive record.  They never refit the
data.  Unless a record identifier is supplied, operations resolve the accepted
record in the selected result slot.

Post-fit result contracts
-------------------------

.. autoclass:: quantas.api.eos.DiagnosticResult
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.CalculationResult
   :members:
   :show-inheritance:

Diagnostics
-----------

.. autofunction:: quantas.api.eos.diagnose

.. autofunction:: quantas.api.eos.diagnostic_summary_table

.. autofunction:: quantas.api.eos.diagnostic_table

.. autofunction:: quantas.api.eos.write_diagnostic_csv

Calculations
------------

.. autofunction:: quantas.api.eos.calculate

.. autofunction:: quantas.api.eos.calculation_summary_table

.. autofunction:: quantas.api.eos.calculation_table

.. autofunction:: quantas.api.eos.write_calculation_csv

Batch reporting
---------------

.. autofunction:: quantas.api.eos.build_batch_preamble

.. autofunction:: quantas.api.eos.build_batch_report

Plot preparation
----------------

.. autofunction:: quantas.api.eos.describe_plots

.. autofunction:: quantas.api.eos.available_plot_types

.. autofunction:: quantas.api.eos.build_plots
