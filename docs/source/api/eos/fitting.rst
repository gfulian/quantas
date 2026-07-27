EOS input and fitting operations
================================

The direct fitting layer parses and normalizes datasets, validates one
scientific request, executes the selected regression strategy, and returns an
in-memory fit result.  It does not mutate an archive unless the caller later
adds the result through a batch or session.

Input
-----

.. autofunction:: quantas.api.eos.read_input

.. autofunction:: quantas.api.eos.normalize_input

Capabilities and request validation
-----------------------------------

.. autofunction:: quantas.api.eos.available_eos_models

.. autofunction:: quantas.api.eos.available_eos_tags

.. autofunction:: quantas.api.eos.available_pvt_couplings

.. autofunction:: quantas.api.eos.available_temperature_eos_models

.. autofunction:: quantas.api.eos.domain_capability

.. autofunction:: quantas.api.eos.default_solver_options

.. autofunction:: quantas.api.eos.validate_request

Direct fitting
--------------

.. autofunction:: quantas.api.eos.fit

Record metadata helper
----------------------

.. autofunction:: quantas.api.eos.record_domain
