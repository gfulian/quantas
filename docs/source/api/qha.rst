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
   inspection_tables = qha.build_inspection_report(preview)
   inspection_plots = qha.build_inspection_plots(preview)
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

.. autodata:: quantas.api.qha.CurveAxis

.. autodata:: quantas.api.qha.PhononInterface

.. autodata:: quantas.api.qha.TableFileFormat

.. autoclass:: quantas.api.qha.TableFormat
   :members:
   :show-inheritance:

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

.. autoclass:: quantas.api.qha.KiefferVolumeSeries
   :members:

.. autoclass:: quantas.api.qha.Preview
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.qha.StructureVolumeSeries
   :no-index:
   :members:
   :show-inheritance:

Input, inspection, and calculation
----------------------------------

``create_input`` shares the HA/QHA phonon generator.  For independent
multi-volume outputs, the normalized YAML may contain Quantas eigenvector-based
mode-continuity diagnostics; for native QHA sources it preserves source-managed
continuity provenance.  See :doc:`../workflows/phonon_input_generation`.

The ``kieffer_cutoffs`` argument accepts one direct cutoff state for every
sampled primitive-cell volume.  This stage supports ``Options(scheme="td")``;
``scheme="freq"`` is rejected until cutoff-frequency evaluation at arbitrary
equilibrium volumes is available.

.. autofunction:: quantas.api.qha.create_input

.. autofunction:: quantas.api.qha.available_energy_eos

.. autofunction:: quantas.api.qha.read_input

.. autofunction:: quantas.api.qha.normalize_input

.. autofunction:: quantas.api.qha.inspect

.. autofunction:: quantas.api.qha.build_inspection_report

.. autofunction:: quantas.api.qha.build_inspection_plots

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

Scalar QHA properties are stored on a native temperature-pressure grid.  The
default line representation shows each property as a function of temperature
at every selected pressure.  ``PlotOptions(curve_axis="pressure")`` produces
the complementary pressure sections at exact stored temperatures when at
least two pressure coordinates exist.  Filled P--T maps continue to use the
complete native grid.

Selections are expressed in the native units reported by
:func:`describe_plots`.  They must match stored coordinates exactly; the first
public implementation does not interpolate or snap nearby values.

.. code-block:: python

   inventory = qha.describe_plots(result_data)
   temperatures = inventory.context_by_key("temperature_grid").values

   pressure_sections = qha.build_plots(
       result_data,
       properties=("equilibrium_volume", "isothermal_bulk_modulus"),
       options=qha.PlotOptions(
           curve_axis="pressure",
           selected_temperatures=(temperatures[0], temperatures[-1]),
       ),
   )

.. autofunction:: quantas.api.qha.list_plot_properties

.. autofunction:: quantas.api.qha.describe_plots

.. autofunction:: quantas.api.qha.build_report

.. autofunction:: quantas.api.qha.build_plots

.. autofunction:: quantas.api.qha.write_result

.. autofunction:: quantas.api.qha.read_result

.. autofunction:: quantas.api.qha.write_table

See also
--------

- :doc:`../workflows/phonon_input_generation`
- :doc:`../workflows/qha`
- :doc:`../tutorials/qha`
- :doc:`../formats/phonon_yaml`
- :doc:`../formats/ha_qha_hdf5`
