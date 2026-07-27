``quantas qha``
===============

The QHA frontend builds a volume-dependent free-energy representation and
minimizes it at every requested pressure-temperature state.  Its command line
therefore exposes scientific choices that are absent from HA: interpolation
scheme, minimization model, polynomial degrees, derivative strategy,
thermal-expansion route, and local-fit failure policy.

Recommended sequence
--------------------

.. code-block:: console

   quantas qha inpgen qha-outputs.txt --list --output material.yaml
   quantas qha inspect material.yaml --eos BM3 --degree 3
   quantas qha run material.yaml --scheme freq --minimization poly \
      --temperature 0 1000 25 --pressure 0 10 1
   quantas qha plot material_QHA.hdf5 --property VT --property alphaV --2d
   quantas qha plot material_QHA.hdf5 --property VT --axis pressure \
      --temperature 300 --temperature 1000
   quantas qha export material_QHA.hdf5 --property VT --format csv

Use ``inspect`` before a production run.  It compares the sampled static
energy-volume data with polynomial and EOS previews and reports the implied
pressure support.  A dense requested P--T grid does not extend the support of
the sampled volumes.

Choosing options
----------------

* ``--scheme=freq`` retains mode-resolved information but requires defensible
  mode continuity.  ``--scheme=td`` interpolates integrated harmonic
  properties and is less dependent on branch tracking.
* ``--minimization=poly`` is flexible near a well-sampled minimum;
  ``--minimization=eos`` imposes a selected physical EOS form.
* ``--thermal-expansion`` chooses among mixed-derivative, mode-Gruneisen, and
  numerical-volume routes.  Their agreement is a diagnostic, not an identity
  guaranteed for every dataset.
* ``--poly-grid-points`` and ``--poly-grid-separation`` control local
  derivatives after polynomial minimization.  They do not change the
  equilibrium volume itself.
* ``--failure-policy`` determines whether failed local states terminate,
  accumulate, or raise immediately; it does not convert an unsupported state
  into a valid one.
* ``plot --axis temperature`` selects exact native pressures with
  ``--pressure``; ``plot --axis pressure`` selects exact native temperatures
  with ``--temperature``.  Plot construction never interpolates or snaps the
  requested coordinate.

See :doc:`../workflows/qha` for the decision guide,
:doc:`../tutorials/qha` for reproducible calculations and method comparisons,
and :doc:`../formats/phonon_yaml` for input details.

Generated command reference
---------------------------

.. click:: quantas.cli.qha:qha
   :prog: quantas qha
   :nested: full
