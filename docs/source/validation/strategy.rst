Validation strategy
===================

Scientific validation in Quantas is kept separate from tutorials and from
ordinary regression testing. A tutorial shows how to run a calculation; a
regression test protects a known behaviour; a validation record explains why a
numerical or scientific result can be trusted within a stated scope.

Evidence hierarchy
------------------

The strongest validation normally combines several independent forms of
evidence:

- analytical limits, identities, and invariants;
- synthetic parameter-recovery problems with known answers;
- comparison with the historical Quantas implementation when that behaviour is
  still scientifically meaningful;
- comparison with independent or established external implementations;
- comparison with published experimental or theoretical reference values;
- end-to-end checks on curated real datasets;
- persistence and frontend-equivalence checks when a scientific result crosses
  API, CLI, HDF5, or export boundaries.

No single level is sufficient for every capability. An analytical identity is
usually the strongest test of a formula, while a real material is more useful
for exposing interactions between parsing, units, symmetry, numerical fitting,
and workflow orchestration.

Validation status
-----------------

The public validation pages use two statuses.

``validated``
   The documented scope has a public evidence record that identifies the
   reference or invariant, observable, comparison method, acceptance criterion
   or numerical tolerance, test or fixture traceability, and the scientific
   limits of the conclusion. This status does **not** claim universal physical
   accuracy outside that scope.

``work in progress``
   The implementation may already have extensive analytical, characterization,
   and regression tests, but the release-candidate evidence record has not yet
   been consolidated to the standard above. A work-in-progress page states what
   is already protected and what still has to be documented or compared before
   the status can change.

This distinction is intentional. Test coverage is necessary for validation,
but test coverage alone is not a validation claim.

What a validation record contains
---------------------------------

For a major scientific capability, the release-candidate record should make the
following information recoverable without reading the implementation:

#. **scope** -- the formula, workflow, material class, or numerical domain being
   validated;
#. **reference** -- an analytical result, independent implementation, published
   datum, frozen historical result, or curated real dataset;
#. **observable** -- the quantity that is actually compared;
#. **comparison** -- how the independent and Quantas values are constructed;
#. **acceptance criterion** -- an analytical equality, invariant, physical
   criterion, or explicit numerical tolerance;
#. **traceability** -- the tests, fixtures, source digests, or validation
   driver that reproduce the evidence;
#. **limits** -- what the comparison does not establish.

The :doc:`matrix` is the compact index of these records. The individual pages
remain authoritative for the detailed evidence.

Regression tolerance is not physical acceptance
------------------------------------------------

A numerical regression tolerance answers a narrow question: *has this
implementation changed beyond the variation expected from supported numerical
environments?* A physical acceptance criterion answers a different question:
*is the model or result adequate for the scientific purpose?*

The two must not be conflated. Examples include:

- a nonlinear EOS fit can need a cross-platform regression tolerance even when
  an analytical pressure derivative is exact;
- a mode-tracking overlap threshold is an operational criterion, not a
  universal constant of lattice dynamics;
- convergence of a directional grid is a numerical requirement, while agreement
  with experiment is a physical validation question;
- a QSA thermoelastic regression can be internally exact while the
  quasi-static approximation itself remains an approximation whose material
  accuracy must be assessed separately.

The numerical precision policy and the still-developing catalogue of tolerances
are collected in :doc:`precision`.

Scientific changes and revalidation
-----------------------------------

Changes to formulas, conventions, units, strain measures, branch ordering,
normalization, tolerances, fitting objectives, or numerical solvers are treated
as scientific changes when they can alter reported values or acceptance
criteria. Such changes require characterization before refactoring, focused
validation after the change, and an update to the public validation record when
its evidence or limits have changed.

A display-only change does not require scientific revalidation when raw values
and stored results are unchanged. This is one reason Quantas keeps rendering
precision separate from numerical precision.
