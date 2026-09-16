Validation matrix
=================

The matrix below summarizes the public release-candidate evidence currently
assembled in this manual. ``validated`` refers only to the scope described in
the linked page; it is not a claim that every material, model choice, or input
dataset has been physically validated. ``work in progress`` means that the
implementation can already be well tested while the public evidence record is
still incomplete. See :doc:`strategy` for the status definitions.

.. list-table:: Current scientific validation matrix
   :header-rows: 1
   :widths: 17 22 20 22 19

   * - Capability
     - Reference / evidence
     - Observable
     - Traceability
     - Status
   * - Energy EOS E--V
     - Analytical E(V)/P(V) identities, synthetic recovery, and a seven-volume
       CRYSTAL/PBE MgO series
     - ``E0``, ``V0``, ``K0``, pressure derivatives, energy residuals, and
       structural response
     - :doc:`eos`; Energy EOS core/workflow tests; curated MgO example
     - **validated**
   * - Experimental P--V BM3 reference scope
     - Quartz and topaz reference fits frozen from EosFit7-compatible OLS and
       effective-variance analyses
     - ``V0``, ``K0``, ``KP``, implied ``KPP``, parameter errors, residuals,
       and reduced chi-square where defined
     - :doc:`eos`; ``tests/modules/eos/test_eosfit_reference.py``
     - **validated** for the documented BM3/solver scope
   * - V--T and P--V--T EOS
     - Existing analytical tests, synthetic datasets, and real tutorial data
     - Fitted thermal/coupled parameters, residuals, uncertainties,
       diagnostics, and derived properties
     - :doc:`eos`; EOS physics/module tests
     - **work in progress** -- consolidated external/reference matrix still to
       close
   * - Kieffer acoustic thermodynamics
     - Analytical limits and identities, historical characterization, isotropic
       acoustics, and multi-volume OHAp data
     - ``F``, ``S``, ``C_V``, ZPE, effective acoustic velocities, cutoffs, and
       HA/QHA composition
     - :doc:`ha_qha`; Kieffer core, HA, QHA, and OHAp reference tests
     - **validated**
   * - CRYSTAL phonon normalization and QHA mode continuity
     - Native MgO QHA output, seven-volume dolomite dispersion series, synthetic
       permutations, phase changes, and degenerate subspaces
     - Primitive reconstruction, eigenvector normalization, assignments,
       cautions, and unresolved links
     - :doc:`ha_qha`; CRYSTAL phonon-mode and QHA tracking tests
     - **validated**
   * - HA/QHA thermodynamic workflows
     - Analytical thermodynamic identities and frozen scientific-reference
       arrays across HA and four QHA scheme/minimization combinations
     - Free energies, heat capacities, equilibrium volumes, bulk moduli, masks,
       and stored result shapes
     - :doc:`ha_qha`; ``tests/baselines/scientific_reference.*`` and HA/QHA
       module tests
     - **validated** within the documented numerical workflow scope
   * - Thermoelastic QSA
     - Real CRYSTAL/QHA MgO and dolomite reference sets plus independent
       least-squares and adiabatic reconstructions
     - Fitted elastic coefficients, C(P,T), stability, frame normalization, and
       isothermal-to-adiabatic correction
     - :doc:`thermoelasticity`; frozen thermoelastic reference fixture
     - **validated** within the documented QSA scope
   * - Elasticity
     - Analytical tensor identities, isotropic limits, formula-oracle checks,
       hydroxylapatite directional baselines, rotations, stability, and HDF5
       round trips already exist
     - Stiffness/compliance, VRH bounds, directional properties, extrema,
       symmetry, frames, and persisted outputs
     - :doc:`elasticity`; elasticity physics/module tests and dedicated baseline
       files
     - **work in progress** -- public reference/tolerance matrix not yet
       consolidated
   * - SEISMIC
     - Isotropic analytical limits, an independent formula reference, frozen
       hydroxylapatite directional baselines, derivatives, tracking,
       degeneracies, enhancement, and rotation tests already exist
     - Phase/group velocities, polarization, power-flow angle, enhancement,
       caustic diagnostics, and sampled fields
     - :doc:`seismic`; seismic physics/module tests and dedicated baseline files
     - **work in progress** -- public reference/tolerance matrix not yet
       consolidated
   * - Numerical precision and tolerance policy
     - ``float64``/``complex128`` storage and module-specific numerical tests
     - dtype, solver/regression tolerances, conditioning, and uncertainty-aware
       comparisons
     - :doc:`precision`; infrastructure and numerical tests
     - **work in progress** -- tolerance catalogue still to consolidate

The matrix is deliberately conservative. A capability moves from *work in
progress* to *validated* only when its public page contains enough information
to reproduce and interpret the evidence, not merely because the implementation
has a large test count.
