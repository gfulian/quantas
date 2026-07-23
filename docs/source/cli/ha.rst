``quantas ha``
==============

The HA frontend evaluates vibrational thermodynamics independently at every
sampled volume.  It does not minimize the free energy with respect to volume
and does not create a pressure axis.

Recommended sequence
--------------------

.. code-block:: console

   quantas ha inpgen phonon-output.out --output material.yaml
   quantas ha run material.yaml --temperature 0 1000 10
   quantas ha plot material_HA.hdf5 --property Cv --property Fvib
   quantas ha export material_HA.hdf5 --property Cv --unit J/mol

``inpgen`` is optional when a valid Quantas phonon YAML already exists.  ``run``
creates the scientific archive and report.  ``plot`` and ``export`` read that
archive without repeating harmonic sums.

Important distinctions
----------------------

* ``--temperature`` controls the calculated grid and stored arrays.
* unit options on ``run`` describe the YAML values; unit options on ``plot`` or
  ``export`` convert an existing result for presentation.
* plotting during ``run`` is a convenience.  The standalone ``plot`` command is
  preferable when several figure variants are required.
* q-point weights are taken from the input and normalized; the CLI does not
  derive missing symmetry multiplicities.

See :doc:`../workflows/ha` for implementation choices,
:doc:`../tutorials/ha` for a complete MgO calculation, and
:doc:`../formats/phonon_yaml` for the input contract.

Generated command reference
---------------------------

.. click:: quantas.cli.ha:ha
   :prog: quantas ha
   :nested: full
