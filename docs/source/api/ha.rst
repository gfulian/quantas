Harmonic Approximation API
==========================

:mod:`quantas.api.ha` exposes the complete supported lifecycle for harmonic
thermodynamics: create or read a normalized phonon input, run the calculation,
retrieve the typed payload, build reports and plots, and persist the result.

Minimal lifecycle
-----------------

.. code-block:: python

   from pathlib import Path
   from quantas.api import ha, rendering

   options = ha.Options(
       temperature_min=0.0,
       temperature_max=1000.0,
       temperature_step=100.0,
   )
   result_data = ha.run("mgo_b3lyp.yaml", options=options)
   result = ha.get_result(result_data)

   report = rendering.render_tables(ha.build_report(result_data))
   Path("mgo_ha.log").write_text(report, encoding="utf-8")
   ha.write_result(result_data, "mgo_ha.hdf5", report_text=report)

The returned payload retains the volume axis even when only one volume is
present.  Temperature-dependent properties normally have shape ``(nT, nV)``;
zero-point energy remains temperature independent.

Passive contracts
-----------------

.. autoclass:: quantas.api.ha.Input
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.ha.Options
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.ha.Result
   :members:
   :show-inheritance:

Input preparation
-----------------

``create_input`` converts supported code-specific outputs to the normalized
phonon YAML contract.  ``read_input`` parses that YAML.  ``normalize_input``
accepts either the public HA input, the shared phonon input contract, or a
path.

.. autofunction:: quantas.api.ha.create_input

.. autofunction:: quantas.api.ha.read_input

.. autofunction:: quantas.api.ha.normalize_input

Calculation and typed results
-----------------------------

.. autofunction:: quantas.api.ha.run

.. autofunction:: quantas.api.ha.get_result

Reports, plots, and persistence
-------------------------------

.. autofunction:: quantas.api.ha.build_report

.. autofunction:: quantas.api.ha.build_plots

.. autofunction:: quantas.api.ha.write_result

.. autofunction:: quantas.api.ha.read_result

See also
--------

- :doc:`../workflows/ha`
- :doc:`../tutorials/ha`
- :doc:`../formats/phonon_yaml`
- :doc:`../formats/ha_qha_hdf5`
