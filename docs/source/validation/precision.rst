Numerical precision and tolerances
==================================

Quantas calculations and native HDF5 values use ``float64`` for real quantities
and ``complex128`` for complex quantities. Display precision belongs to
renderers and does not modify stored data.

.. admonition:: Work in progress

   The storage policy is fixed and tested, but the release-candidate catalogue
   that explains every module-specific regression tolerance, conditioning rule,
   and uncertainty-aware comparison is still being consolidated. Until that
   catalogue is complete, the relevant validation page and test remain the
   authoritative source for a particular tolerance.

Three different kinds of numerical criterion appear in the project and should
not be confused.

Regression tolerances
---------------------

A regression tolerance protects a frozen numerical result across supported
Python, NumPy, SciPy, BLAS/LAPACK, and platform combinations. It should be tight
enough to detect a scientific change but not so tight that harmless differences
in nonlinear-solver termination become release failures.

Analytical identities and exactly constructed synthetic problems should use
much tighter tolerances than real nonlinear fits. When a real-data parameter is
known to be solver-sensitive, the stable observables should carry the CI gate
and the sensitive quantity should be validated independently by analytical or
synthetic tests. The Energy EOS validation in :doc:`eos` follows this pattern
for higher pressure derivatives.

Solver and convergence criteria
-------------------------------

A solver tolerance controls an iterative or discretized numerical method. It is
part of the scientific numerical method when changing it can alter a reported
result. Examples include optimization termination, quadrature convergence,
directional-grid resolution, degeneracy thresholds, and mode-tracking overlap
criteria.

Such values require method-specific characterization; they are not universal
``float64`` tolerances.

Physical acceptance criteria
----------------------------

A physical acceptance criterion concerns the scientific adequacy of the result,
not floating-point agreement. Mechanical stability, absence of imaginary modes,
fit residual quality, experimental agreement, or convergence with respect to an
upstream electronic-structure parameter belong in this category.

Passing a numerical regression does not imply passing a physical acceptance
criterion, and vice versa. :doc:`strategy` describes how this distinction is
used in the validation record.

Still to consolidate
--------------------

The final release-candidate precision page will add a module-by-module table
covering:

- frozen-regression ``rtol``/``atol`` values and their rationale;
- nonlinear fitting and inversion tolerances;
- directional sampling and quadrature convergence criteria;
- degeneracy, overlap, and ambiguity thresholds;
- uncertainty-aware comparisons and covariance conditioning;
- rules for changing a tolerance without silently changing the scientific
  contract.
