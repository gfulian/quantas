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

All supported scientific modules now advertise ``PLOT_INVENTORY``.  The
one-shot modules describe a typed result, while EOS describes an archive plus
an optional slot or immutable fit record.  The capability therefore remains
semantically uniform without forcing EOS into the one-shot lifecycle.

Capability and descriptor contracts
-----------------------------------

.. autoclass:: quantas.api.registry.Capability
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.registry.ModuleDescriptor
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.registry.OperationDescriptor
   :members:
   :show-inheritance:

A capability may have one canonical operation and additional named operations.
For example, Thermoelasticity and EOS expose several scientifically distinct
exports.  Use ``operations_for`` to enumerate them and ``named_operation`` to
resolve one stable operation key.

.. code-block:: python

   eos = registry.get("eos")
   for operation in eos.list_operations(Capability.EXPORT):
       print(operation.key, operation.name)

   write_diagnostics = eos.named_operation("export_diagnostics_csv")

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
