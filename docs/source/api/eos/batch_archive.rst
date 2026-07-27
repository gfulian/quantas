EOS batch, archive, and session API
===================================

Batch plans describe a deterministic sequence of fit requests.  Native archives
preserve every immutable attempt together with accepted/candidate state and
history.  ``Session`` provides the stateful public surface appropriate for a
notebook, service, or Quantas GUI.

Batch contracts
---------------

.. autoclass:: quantas.api.eos.BatchFailurePolicy
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.BatchJob
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.BatchPlan
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.eos.BatchResult
   :members:
   :show-inheritance:

Specification files
-------------------

.. autodata:: quantas.api.eos.SPEC_TEMPLATE_FILENAME

.. autofunction:: quantas.api.eos.read_spec

.. autofunction:: quantas.api.eos.resolve_spec

.. autofunction:: quantas.api.eos.write_spec_template

Batch execution
---------------

.. autofunction:: quantas.api.eos.run_batch

Native archive
--------------

.. autoclass:: quantas.api.eos.Archive
   :members:
   :show-inheritance:

.. autofunction:: quantas.api.eos.open_archive

Persistent session
------------------

.. autoclass:: quantas.api.eos.Session
   :members:
   :show-inheritance:
