Common public contracts
=======================

:mod:`quantas.api.common` contains the frontend-neutral contracts shared by
all public workflows.  They define how normalized input, results, reports,
plots, and events move between numerical code and frontends.

Data lifecycle
--------------

``InputData`` and ``PhononInputData`` preserve normalized input and provenance.
``ResultData`` is the common HDF5-compatible result envelope.  A module result
is stored as a typed payload inside that envelope and should be retrieved with
the module-specific ``get_result`` function or, for generic frontend code,
``get_result_payload``.

.. code-block:: python

   from quantas.api import common, elasticity

   result_data = elasticity.run("calcite.dat")
   result = common.get_result_payload(
       result_data,
       module="elasticity",
       key="elasticity",
       expected_type=elasticity.Result,
   )

Input and result contracts
--------------------------

.. autoclass:: quantas.api.common.InputData
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.PhononInputData
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.ResultData
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.ResultMetadata
   :members:
   :show-inheritance:

Reusable scientific input structures
------------------------------------

.. autodata:: quantas.api.common.PhononInterface

.. autoclass:: quantas.api.common.StructureVolumeSeries
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.CrystalStructure
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.SymmetryMetadata
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.TensorRotation
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.TensorRotationKind
   :members:
   :show-inheritance:

.. autofunction:: quantas.api.common.get_result_payload

Neutral reporting and plotting
------------------------------

``ReportTable`` stores raw values, labels, units, alignment, and display
metadata without embedding Rich or terminal behavior.  ``PlotCollection``
stores frontend-neutral plot specifications.  Their concrete public types
are documented in :mod:`quantas.api.plotting`; rendering is delegated to
:mod:`quantas.api.rendering`.

.. autoclass:: quantas.api.common.ReportTable
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.PlotCollection
   :members:
   :show-inheritance:

Events and observers
--------------------

Observers are optional sinks for operational events.  ``ListObserver`` is
useful in notebooks and tests, ``CallbackObserver`` adapts an arbitrary
callable, and ``NullObserver`` explicitly discards events.  Progress events
are operational and are not persisted as scientific history.

.. code-block:: python

   from quantas.api import qha
   from quantas.api.common import ListObserver

   observer = ListObserver()
   result = qha.run("phonons.yaml", observer=observer)
   for event in observer.events:
       print(event.level.value, event.message, event.progress)

.. autoclass:: quantas.api.common.EventLevel
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.Event
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.EventRecord
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.Observer
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.CallbackObserver
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.ListObserver
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.common.NullObserver
   :members:
   :show-inheritance:

See also
--------

- :doc:`../getting_started/results`
- :doc:`../workflows/concepts/reporting_and_events`
- :doc:`plotting`
- :doc:`rendering`
- :doc:`../formats/hdf5`
