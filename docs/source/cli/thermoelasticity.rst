``quantas thermoelasticity``
============================

Thermoelasticity is a staged workflow.  ``run`` calibrates a reusable cold
finite-strain model from static elastic calculations and QHA data.  Subsequent
``analysis`` commands evaluate points, rectangular P--T grids, or geological
profiles.  Table, export, inspection, and plotting commands consume those
archives without repeating calibration.

Recommended sequence
--------------------

.. code-block:: console

   quantas thermoelasticity inpgen soec-files.txt --list \
      --output material_thermoelastic.yaml
   quantas thermoelasticity run material_thermoelastic.yaml material_QHA.hdf5 \
      --output material_FIT.hdf5
   quantas thermoelasticity inspect material_FIT.hdf5
   quantas thermoelasticity analysis point material_FIT.hdf5 5 800
   quantas thermoelasticity analysis grid material_FIT.hdf5 \
      --pressure 0 10 1 --temperature 300 1200 50
   quantas thermoelasticity analysis profile material_FIT.hdf5 \
      --profile-spec profile.yaml

Command families
----------------

``inpgen`` and ``profile-template``
   Prepare normalized elastic input and editable geological profile
   specifications.

``run``
   Fit the reference static EOS and independent ``Cij(V)`` models.  The output
   is a calibration archive, not a P--T grid.

``analysis point/grid/profile``
   Evaluate calibrated tensors, density, uncertainty, stability, and support
   flags at the requested states.

``table`` and ``export``
   Produce compact portable selections from existing archives.  Point export
   can create the shared Elasticity/SEISMIC tensor input.

``inspect``
   Explain archive stage, coverage, tensor conditions, and useful next
   commands.

``plot fit/pt/profile/domain/compare``
   Separate calibration diagnostics from P--T fields, geological paths,
   validity domains, and isothermal/adiabatic comparisons.  Use
   ``quantas thermoelasticity plot --list-plots --archive RESULT.hdf5`` to query
   the public result-aware inventory and list only families buildable from that
   archive.

Scientific cautions
-------------------

* QHA-coordinate extrapolation and elastic-volume extrapolation are distinct
  and remain separate in command options and result masks.
* ``isothermal`` and ``adiabatic`` identify tensor conditions, not alternative
  cold finite-strain fits.
* validation presets describe support and policy; they do not improve weak
  source data.
* densifying an analysis grid does not add information to the underlying QHA
  grid or elastic-volume calibration.

See :doc:`../workflows/thermoelasticity` for implementation details,
:doc:`../tutorials/thermoelasticity` for the complete dolomite workflow,
:doc:`../formats/thermoelastic_input`, :doc:`../formats/earth_profile_spec`, and
:doc:`../formats/thermoelastic_hdf5` for data contracts.

Generated command reference
---------------------------

.. click:: quantas.cli.thermoelastic:thermoelasticity
   :prog: quantas thermoelasticity
   :nested: full
