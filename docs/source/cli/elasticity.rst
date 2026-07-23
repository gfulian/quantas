``quantas elasticity``
======================

The Elasticity frontend validates a second-order stiffness tensor, evaluates
mechanical stability and isotropic aggregate bounds, searches global
directional extrema, and optionally persists principal-plane or full-sphere
sampled fields.

Recommended sequence
--------------------

.. code-block:: console

   quantas elasticity inpgen crystal-elastic.out --interface crystal
   quantas elasticity run crystal-elastic_elasticity_input.dat --2d --3d
   quantas elasticity plot crystal-elastic_elasticity.hdf5 --2d --3d \
      --property young --preset publication
   quantas elasticity export crystal-elastic_elasticity.hdf5

``inpgen`` converts supported external outputs but does not remove the need to
verify pressure convention, units, density metadata, and tensor frame.  ``run``
always computes scalar summaries and global extrema; ``--2d`` and ``--3d`` add
sampled data products.

Sampling and plotting
---------------------

* ``--ntheta`` has different documented defaults for 2D and 3D calculation.
  ``--nphi`` applies only to 3D surfaces.
* global extrema use their own search-and-refinement algorithm and are not
  controlled by the plotting grid.
* ``physical`` geometry uses the property as radius; ``unit-sphere`` uses a
  sphere and encodes the property by color.
* angular overrides on the standalone ``plot`` command can rebuild a transient
  surface at another visual resolution without modifying the HDF5 archive.
* rotation options transform tensor components and reported directions into a
  new analysis frame; invariant scalar properties should remain unchanged.

See :doc:`../workflows/elasticity` for sampling and stability details,
:doc:`../tutorials/elasticity` for the calcite example, and
:doc:`../formats/elasticity_input` for the text-input contract.

Generated command reference
---------------------------

.. click:: quantas.cli.elasticity:elasticity
   :prog: quantas elasticity
   :nested: full
