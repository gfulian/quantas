``quantas seismic``
===================

The SEISMIC frontend solves the Christoffel eigenproblem over a spherical grid.
The selected calculation level determines whether the archive contains only
phase fields, first-derivative group properties, or second-derivative
enhancement and caustic diagnostics.

Recommended sequence
--------------------

.. code-block:: console

   quantas seismic inpgen OUTCAR --interface vasp --output material.dat
   quantas seismic run material.dat --level phase --ntheta 61 --nphi 121
   quantas seismic run material.dat --level enhancement --output material_final.hdf5
   quantas seismic plot material_final.hdf5 --summary
   quantas seismic plot material_final.hdf5 --2d --property phase_v_s1 \
      --polarizations
   quantas seismic export material_final.hdf5 --output material_fields.csv

``inpgen`` exposes the same CRYSTAL/VASP generator as
``quantas.api.seismic.create_input``.  Unlike Elasticity generation, SEISMIC
requires finite positive density metadata and refuses to create a runnable input
when the external output contains only a stiffness tensor.

A staged resolution study is normally faster than beginning with the densest
``enhancement`` calculation.  Converge phase velocities and splitting first,
then group properties, and only then curvature-dependent enhancement or
caustic candidates.

Key controls
------------

* ``--hemisphere`` changes the sampled angular interval.  The same ``ntheta``
  gives a coarser spacing on a full sphere than on an upper hemisphere.
* ``--level`` controls scientific content and archive size.
* ``--ntheta`` and ``--nphi`` control angular resolution and therefore sampled
  extrema and candidate locations.
* ``--batch-size`` controls vectorized throughput and temporary memory only;
  its default is 512 directions and it is not an accuracy parameter.
* eigenvalue and degeneracy tolerances classify near-zero numerical cases.
  They should not be enlarged to hide an unstable stiffness tensor.
* polarization tracking changes branch continuity and visualization, not phase
  velocities.

See :doc:`../workflows/seismic` for the detailed numerical workflow,
:doc:`../tutorials/seismic` for the hydroxylapatite example, and
:doc:`../formats/elasticity_input` for the shared tensor input.

Generated command reference
---------------------------

.. click:: quantas.cli.seismic:seismic
   :prog: quantas seismic
   :nested: full
