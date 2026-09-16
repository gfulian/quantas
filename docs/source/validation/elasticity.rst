Elasticity validation
=====================

.. admonition:: Work in progress

   The Elasticity implementation is already covered by analytical and
   regression tests, but the public validation record has not yet been assembled
   into the release-candidate evidence format used by this section.

Planned validation record
-------------------------

The final page will collect the reference datasets, tolerances, and test
traceability for:

- stiffness/compliance inversion and tensor rotations;
- symmetry specialization and mechanical-stability criteria;
- Voigt, Reuss, and Hill aggregate bounds;
- analytical directional properties and exact extrema;
- 2D/3D sampling and frame invariance;
- pressure/prestress characterization where it affects the interpreted tensor.

Until that record is completed, :doc:`../theory/elasticity`,
:doc:`../workflows/elasticity`, and the automated test suite describe the
implemented behavior but should not be read as a substitute for a consolidated
validation matrix.
