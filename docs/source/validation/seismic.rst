SEISMIC validation
==================

.. admonition:: Work in progress

   The SEISMIC implementation already has analytical, isotropic-limit,
   degeneracy, tracking, and field regression tests. The public validation
   record still needs to consolidate those checks into the release-candidate
   evidence format used by this section.

Planned validation record
-------------------------

The final page will collect the reference datasets, tolerances, and test
traceability for:

- Christoffel eigenvalues and isotropic analytical limits;
- phase and analytical group velocities;
- polarization tracking and degenerate eigenspaces;
- enhancement, area factors, and caustic-candidate diagnostics;
- rotational invariance and directional-grid convergence;
- HDF5/export round trips and field-level reproducibility.

Until that record is completed, :doc:`../theory/seismic`,
:doc:`../workflows/seismic`, and the automated test suite document the
implemented behavior but are not a consolidated validation report.
