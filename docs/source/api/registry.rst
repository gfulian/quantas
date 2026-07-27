Capability registry API
=======================

:mod:`quantas.api.registry` allows a frontend to discover scientific modules,
passive types, and supported operations without importing implementation
classes or forcing every workflow into one inheritance hierarchy.

.. code-block:: python

   from quantas.api import registry
   from quantas.api.registry import Capability

   for descriptor in registry.list_modules():
       print(descriptor.name, sorted(item.value for item in descriptor.capabilities))

   qha = registry.get("qha")
   run_qha = qha.operation(Capability.RUN)

   elasticity = registry.get("elasticity")
   if elasticity.has(Capability.PLOT_INVENTORY):
       describe = elasticity.operation(Capability.PLOT_INVENTORY)

``PLOT_INVENTORY`` is introduced incrementally.  In the current development
snapshot it is declared by Elasticity, SEISMIC, and HA.  A module must not
advertise it until its result-aware inventory is public and covered by builder
compatibility tests.  EOS will use a separate session-oriented inventory rather
than pretending to be a one-shot result workflow.

Capability and descriptor contracts
-----------------------------------

.. autoclass:: quantas.api.registry.Capability
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.registry.ModuleDescriptor
   :members:
   :show-inheritance:

Discovery
---------

.. autofunction:: quantas.api.registry.get

.. autofunction:: quantas.api.registry.iter_modules

.. autofunction:: quantas.api.registry.list_modules

Result dispatch
---------------

``module_from_result`` inspects native metadata rather than filename
conventions.  ``open_result`` dispatches to the public module reader or EOS
archive as appropriate.

.. autofunction:: quantas.api.registry.module_from_result

.. autofunction:: quantas.api.registry.open_result

See also
--------

- :doc:`../getting_started/python_api`
- :doc:`../formats/hdf5`
- :doc:`../formats/hdf5_inspection`
