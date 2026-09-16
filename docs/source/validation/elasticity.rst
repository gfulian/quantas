Elasticity validation
=====================

.. admonition:: Work in progress

   The Elasticity implementation already has analytical, independent-formula,
   real-tensor regression, rotation, stability, sampling, and persistence
   coverage. The remaining work is to assemble those checks into a public
   release-candidate record with explicit reference values, tolerances, and
   traceability for each scientific observable.

Evidence already protected
--------------------------

The current automated suite covers substantially more than a placeholder
validation page would suggest. In particular it includes:

- exact Voigt--Cartesian stiffness and compliance round trips, including the
  engineering-shear factors;
- minor and major tensor symmetries and stiffness/compliance inversion;
- rotational covariance, rotation round trips, fixed-axis ``xyz`` rotations,
  and the trigonal ``C14``/``C15`` convention exchange;
- isotropic rotational invariance;
- symmetry detection and specialization, including hexagonal and orthorhombic
  cases;
- Voigt, Reuss, and Hill aggregate properties against characterized reference
  values;
- positive-definiteness and mechanical-stability reporting;
- pointwise Young modulus, linear compressibility, shear modulus, and Poisson
  ratio against an independent formula reference;
- exact transverse extrema and consistency between returned values and angles;
- batched and pointwise directional sampling equivalence, periodic plane
  closure, smooth isotropic degeneracy, and progress contracts;
- a frozen hydroxylapatite directional-surface baseline in
  ``tests/baselines/elasticity_reference.*``;
- source-frame preservation and explicit analysis-frame rotations across both
  Elasticity and SEISMIC consumers;
- HDF5 round trips, ``float64`` persistence, neutral plot specifications, and
  CSV/table export.

The shared elasticity physics layer also has independent coverage for the
finite-prestress helpers, cold finite-strain QSA relations, and the
isothermal-to-adiabatic correction. Their thermoelastic end-to-end validation is
reported separately in :doc:`thermoelasticity`.

What still has to be assembled
------------------------------

Before this page can be marked ``validated`` under :doc:`strategy`, the public
record still needs to bring the existing evidence together in a compact matrix
that states, for each major observable:

- the analytical or independent reference;
- the material/tensor used for anisotropic regression;
- the numerical tolerance and why it is appropriate;
- the exact implementing test or frozen fixture;
- the convergence criterion for directional extrema and sampled surfaces;
- the scientific limits of the comparison.

A small set of material-level comparisons should also be chosen so that the
release record does not rely only on self-consistency identities. Candidate
systems should exercise at least a high-symmetry tensor and a genuinely
anisotropic tensor while keeping the reference provenance reproducible.

Current traceability
--------------------

The main automated evidence is distributed across:

``tests/physics/elasticity/test_tensor_formulas.py``
   Tensor conversion, symmetry, directional formulas, isotropic rotation
   invariance, and trigonal frame conventions.

``tests/physics/elasticity/test_directional.py`` and
``tests/physics/elasticity/test_surface_formulas.py``
   Independent directional formulas, transverse extrema, surface sampling, and
   the frozen hydroxylapatite baseline.

``tests/physics/elasticity/test_frames.py``
   Rotation covariance, provenance, and reference-frame validation.

``tests/physics/elasticity/test_averages_stability.py``
   VRH aggregates and positive-definiteness diagnostics.

``tests/modules/elasticity/``
   End-to-end analysis, reports, rotations, persistence, export, plotting, CLI,
   and public-API behaviour.

Until the remaining public reference/tolerance matrix is assembled, the
existing tests establish implementation consistency but should not be read as a
claim of universal material-level accuracy. For the scientific definitions and
workflow semantics, see :doc:`../theory/elasticity` and
:doc:`../workflows/elasticity`.
