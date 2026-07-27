Elasticity API
==============

:mod:`quantas.api.elasticity` exposes second-order elastic analysis, exact
transverse extrema, optional two-dimensional sections, sampled
three-dimensional surfaces, neutral reports and plots, and native HDF5
persistence.

Minimal lifecycle
-----------------

.. code-block:: python

   from quantas.api import elasticity

   options = elasticity.Options(calculate_2d=True)
   result_data = elasticity.run("calcite.dat", options=options)
   result = elasticity.get_result(result_data)
   print(result.averages.hill.bulk_modulus)

Use ``calculate_3d`` and ``SurfaceOptions`` only when three-dimensional fields
must be persisted in the result.  A transient 3-D surface can also be built
later from the stored stiffness tensor with :func:`build_3d_plots`.

Passive contracts and selectors
-------------------------------

.. autoclass:: quantas.api.elasticity.Input
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.elasticity.Options
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.elasticity.SurfaceOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.elasticity.Result
   :members:
   :show-inheritance:

.. autodata:: quantas.api.elasticity.InputInterface

.. autoclass:: quantas.api.elasticity.TensorRotation
   :no-index:
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.elasticity.TensorRotationKind
   :no-index:
   :members:
   :show-inheritance:

``PlotProperty`` and ``SurfaceProperty`` accept ``young``, ``compressibility``,
``shear``, and ``poisson``.  ``SurfaceGeometry`` accepts ``physical`` or
``unit_sphere``.

.. autodata:: quantas.api.elasticity.PlotProperty

.. autodata:: quantas.api.elasticity.SurfaceProperty

.. autodata:: quantas.api.elasticity.SurfaceGeometry

Input and calculation
---------------------

.. autofunction:: quantas.api.elasticity.create_input

.. autofunction:: quantas.api.elasticity.read_input

.. autofunction:: quantas.api.elasticity.normalize_input

.. autofunction:: quantas.api.elasticity.run

.. autofunction:: quantas.api.elasticity.get_result

Reporting and plotting
----------------------

.. autofunction:: quantas.api.elasticity.build_report

.. autofunction:: quantas.api.elasticity.describe_plots

.. autofunction:: quantas.api.elasticity.build_plots

.. autofunction:: quantas.api.elasticity.build_2d_plots

.. autofunction:: quantas.api.elasticity.build_3d_plots

Persistence
-----------

.. autofunction:: quantas.api.elasticity.write_result

.. autofunction:: quantas.api.elasticity.read_result

.. autofunction:: quantas.api.elasticity.write_table

See also
--------

- :doc:`../workflows/elasticity`
- :doc:`../tutorials/elasticity`
- :doc:`../formats/elasticity_input`
- :doc:`../formats/elasticity_seismic_hdf5`
