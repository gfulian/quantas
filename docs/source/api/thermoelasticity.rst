Thermoelasticity API
====================

:mod:`quantas.api.thermoelasticity` exposes the complete quasi-static
thermoelastic lifecycle: normalize an elastic-volume series, couple it to QHA,
calibrate a reusable cold finite-strain model, reconstruct P--T grids or depth
profiles, select isothermal or adiabatic tensors, and export frontend-neutral
reports, tables, inputs, and plots.

The API has two distinct stages:

.. code-block:: text

   calibration: elastic Cij(V) + QHA -> reusable fit ResultData
   analysis:    fit ResultData -> point / grid / profile ResultData

Do not assume that the output of :func:`run` already represents a requested
P--T analysis grid.

Reference sections
------------------

.. toctree::
   :maxdepth: 1

   thermoelasticity/contracts
   thermoelasticity/calibration
   thermoelasticity/analysis
   thermoelasticity/plotting

See also
--------

- :doc:`../workflows/thermoelasticity`
- :doc:`../tutorials/thermoelasticity`
- :doc:`../formats/thermoelastic_input`
- :doc:`../formats/thermoelastic_hdf5`
- :doc:`../formats/earth_profile_spec`
