SEISMIC API
===========

:mod:`quantas.api.seismic` exposes directional Christoffel analysis from phase
velocity through group velocity, polarization tracking, acoustic enhancement,
caustic candidates, neutral summaries and surfaces, CSV export, and native
HDF5 persistence.

Minimal lifecycle
-----------------

.. code-block:: python

   from quantas.api import seismic

   options = seismic.Options(
       level="group",
       ntheta=61,
       nphi=121,
       track_polarization_axes=True,
   )
   result_data = seismic.run("hydroxylapatite.dat", options=options)
   result = seismic.get_result(result_data)
   print(result.field.phase.phase_speeds.shape)

The sampling resolution is controlled by ``ntheta`` and ``nphi``.  Batch size
controls memory and throughput only; it does not change the sampled directions
or numerical accuracy.

Passive contracts and selectors
-------------------------------

.. autoclass:: quantas.api.seismic.ElasticMedium
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.seismic.Input
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.seismic.Options
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.seismic.PlotOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.seismic.SurfaceOptions
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.seismic.Result
   :members:
   :show-inheritance:

The public enums below are the same types used by :class:`Options` and
:class:`SurfaceOptions`; CLI and GUI clients do not need imports from
implementation packages.

.. autoclass:: quantas.api.seismic.Hemisphere
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.seismic.SamplingLevel
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.seismic.WaveMode
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.seismic.TensorRotation
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.seismic.TensorRotationKind
   :members:
   :show-inheritance:

``SurfaceType`` accepts ``phase``, ``slowness``, or ``group``.
``SurfaceGeometry`` accepts ``physical`` or ``unit_sphere``.

.. autodata:: quantas.api.seismic.SurfaceType

.. autodata:: quantas.api.seismic.SurfaceGeometry

Input and calculation
---------------------

.. autofunction:: quantas.api.seismic.read_input

.. autofunction:: quantas.api.seismic.normalize_input

.. autofunction:: quantas.api.seismic.run

.. autofunction:: quantas.api.seismic.get_result

Reports, summaries, and surfaces
--------------------------------

.. autofunction:: quantas.api.seismic.build_report

.. autofunction:: quantas.api.seismic.describe_plots

.. autofunction:: quantas.api.seismic.build_summary

.. autofunction:: quantas.api.seismic.build_plots

.. autofunction:: quantas.api.seismic.build_surfaces

Export and persistence
----------------------

.. autofunction:: quantas.api.seismic.write_csv

.. autofunction:: quantas.api.seismic.write_result

.. autofunction:: quantas.api.seismic.read_result

See also
--------

- :doc:`../workflows/seismic`
- :doc:`../tutorials/seismic`
- :doc:`../formats/elasticity_input`
- :doc:`../formats/elasticity_seismic_hdf5`
