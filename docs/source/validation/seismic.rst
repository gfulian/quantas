SEISMIC validation
==================

.. admonition:: Work in progress

   SEISMIC already has analytical-limit, independent-formula, derivative,
   degeneracy, polarization-tracking, enhancement, rotation, sampling, and
   persistence tests. The remaining task is to consolidate that evidence into
   the public release-candidate matrix with explicit reference values,
   tolerances, convergence criteria, and scope limits.

Evidence already protected
--------------------------

The present test suite covers the principal layers of the Christoffel analysis:

- isotropic analytical phase velocities and their density scaling;
- Christoffel eigenpairs against an independent formula reference;
- gradient and Hessian expressions checked against centred finite differences;
- analytical group velocity, the radial group/phase identity, power-flow angle,
  and isotropic radial propagation;
- the public mode order ``V_S2, V_S1, V_P`` and its separation from tracked
  polarization-branch continuity;
- shear and triple degeneracies, polarization-sign continuity, local shear-mode
  exchange, and degenerate-subspace alignment;
- enhancement/area-factor calculations, antipodal symmetry, density-scaling
  invariance, and caustic-candidate thresholding;
- Christoffel rotational covariance;
- spherical-grid geometry, periodic presentation seams, hemisphere selection,
  and batched/pointwise equivalence;
- frozen hydroxylapatite directional results in
  ``tests/baselines/seismic_reference.*``;
- HDF5 round trips, sampling-level persistence, CSV export, neutral plot
  specifications, CLI/API equivalence, and frontend-independent workflow
  execution.

The frozen SEISMIC baseline deliberately records known limitations of the
historical Quantas 0.9 reference. Current tests therefore use it as a numerical
characterization source where appropriate rather than treating every historical
convention as scientifically authoritative. New invariants and corrected
conventions are tested independently.

What still has to be assembled
------------------------------

Before this page can move to ``validated`` under :doc:`strategy`, the public
record still needs:

- a compact table of analytical and independent-reference observables with
  numerical tolerances;
- explicit directional-grid convergence targets for the fields that depend on
  angular resolution;
- material-level reference checkpoints for phase and group velocities;
- a documented acceptance criterion for polarization tracking around difficult
  near-degenerate regions;
- a reproducible interpretation of enhancement and caustic diagnostics that
  separates numerical field convergence from any material-specific physical
  conclusion;
- direct test/fixture traceability for each row of the final matrix.

These additions are primarily a validation-documentation task. They should not
be manufactured by copying current implementation output into the manual; the
reference must remain analytically or independently motivated.

Current traceability
--------------------

The main automated evidence is distributed across:

- ``tests/physics/seismic/test_solver_reference.py`` and
  ``tests/physics/seismic/test_invariants.py`` -- Christoffel eigenpairs,
  analytical derivatives, isotropic limits, and phase/group/enhancement
  invariants.
- ``tests/physics/seismic/test_group_velocity.py`` and
  ``tests/physics/seismic/test_derivatives.py`` -- group-velocity identities and
  finite-difference verification of derivatives.
- ``tests/physics/seismic/test_polarization.py`` and
  ``tests/physics/seismic/test_acoustic_axes.py`` -- sign continuity, local
  branch exchange, degeneracies, subspace alignment, and acoustic-axis
  behaviour.
- ``tests/physics/seismic/test_enhancement.py`` -- enhancement, curvature,
  density and antipodal invariance, and caustic diagnostics.
- ``tests/physics/seismic/test_reference_baseline.py`` -- frozen
  hydroxylapatite directional regression data and logarithmic enhancement
  convention.
- ``tests/modules/seismic/`` -- end-to-end workflow, reports, HDF5, export,
  rendering, CLI, and API parity.

Until the remaining public matrix is assembled, these tests establish strong
numerical and architectural characterization but do not by themselves establish
universal predictive accuracy for anisotropic seismic observables. See
:doc:`../theory/seismic` and :doc:`../workflows/seismic` for the scientific
model and implemented workflow.
