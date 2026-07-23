Interoperability API
====================

:mod:`quantas.api.interop` contains explicit, scientifically meaningful
transformations between public workflows.  These functions centralize unit
conversion, tensor condition, coverage checks, normalization, and provenance;
they are preferable to manually moving arrays between modules.

The complete scientific contracts, CLI equivalents, failure policies, and
end-to-end examples are described in :doc:`../workflows/interoperability`.
This page defines the callable contract.

QHA loading and coupling
------------------------

.. autofunction:: quantas.api.interop.load_qha_result

.. autofunction:: quantas.api.interop.qha_to_thermoelastic_context

Thermoelastic state to SEISMIC
------------------------------

.. autofunction:: quantas.api.interop.thermoelastic_to_seismic

Complete public-API example
---------------------------

The downloadable example performs the complete supported chain and writes
reusable HDF5 checkpoints:

:download:`workflow_api.py <../_downloads/interoperability/workflow_api.py>`

Run it from the project root with::

   python examples/interoperability/workflow_api.py \
       --output-dir interoperability_api

See also
--------

- :doc:`../workflows/interoperability`
- :doc:`../workflows/thermoelasticity`
- :doc:`../tutorials/thermoelasticity`
- :doc:`qha`
- :doc:`thermoelasticity`
- :doc:`seismic`
