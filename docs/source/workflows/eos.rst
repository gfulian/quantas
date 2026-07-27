EOS fitting: implementation and workflow
========================================

Purpose and scope
-----------------

The standalone EOS subsystem is different from most Quantas modules.  A
harmonic, elastic, seismic, or thermoelastic calculation normally transforms
one normalized scientific input into one result envelope.  Equation-of-state
analysis is instead an iterative model-selection problem: the same dataset may
be fitted repeatedly with different domains, formulations, regressions,
constraints, data selections, and initial conditions before one result is
accepted for later use.

Quantas therefore treats EOS analysis as a sequence of immutable fit attempts
stored in a persistent archive.  The workflow is designed to answer questions
such as:

- which observations and physical target are being fitted?;
- which P--V, V--T, or P--V--T formulation is being tested?;
- which statistical interpretation is appropriate for the available
  uncertainties?;
- are fixed parameters, bounds, or initial values scientifically justified?;
- does a numerically converged fit have acceptable residuals and parameter
  correlations?;
- which successful attempt is the currently accepted result?;
- can the accepted model be used safely for diagnostics, plots, or property
  calculation?;

The equations and their physical assumptions are described in
:doc:`../theory/eos`.  Worked examples are provided in
:doc:`../tutorials/eos`.  This page focuses on the implementation choices,
workflow semantics, diagnostics, and reasons why the EOS command line differs
from the other Quantas modules.

Why EOS uses a different command-line workflow
----------------------------------------------

Most Quantas workflows follow a nearly linear pattern:

.. code-block:: text

   scientific input → run → one HDF5 result → plot / export

EOS fitting is not naturally linear.  A realistic analysis more often follows:

.. code-block:: text

   dataset
      ├─ fit candidate A
      ├─ fit candidate B
      ├─ refit after excluding or regrouping observations
      ├─ test another regression or parameter constraint
      ├─ compare diagnostics
      └─ accept one record for each scientific slot
            ↓
         diagnose / calculate / plot

This difference explains several deliberate CLI choices.

``quantas eos run`` is a batch fitting command
   It may execute one job assembled from command-line options or an ordered
   collection of jobs defined in a strict specification file.  Every attempted
   fit is persisted before the command decides whether it should become the
   accepted result.

The native archive is a history, not only a final result
   It contains the normalized dataset, immutable fit records, solver
   diagnostics, selection metadata, accepted-result views, and meaningful
   workflow events.  Failed attempts are retained because they document what
   was tried and why it failed.

Post-fit commands select a record
   ``diagnose``, ``plot``, and ``calculate`` operate on an accepted
   domain/target slot or on an explicit immutable record identifier.  They do
   not silently rerun the fit.

No long interactive Click prompt tree is provided
   Reproducible advanced CLI use is expressed by the declarative specification
   format.  A frontend-neutral persistent :class:`quantas.api.eos.Session`
   remains available for Quantas GUI, where interactive acceptance, rejection,
   comparison, and reuse of previous guesses are more natural.

The purpose of this design is not to make EOS artificially complicated.  It is
to prevent a trial model, failed fit, or later refit from silently overwriting
the scientific history of the analysis.

Computational pipeline
----------------------

The complete non-interactive path is:

.. code-block:: text

   keyword-directed EOS data file
       │
       ├─ parse columns, units, groups, uncertainties, and row markers
       ├─ retain original text and raw numerical values
       └─ normalize all active arrays to float64 Quantas units

   command options or QUANTAS EOS SPEC 1
       │
       ├─ resolve domain and target slots
       ├─ resolve models, regression, constraints, and data selections
       ├─ validate the complete ordered plan
       └─ optionally stop after preflight (--dry-run)

   ordered batch execution
       │
       ├─ construct one immutable FitRequest per job
       ├─ fit the selected observations
       ├─ build common diagnostics and parameter metadata
       ├─ append every attempt to the HDF5 archive
       ├─ accept successful jobs only when explicitly allowed
       └─ stop or continue after a failure according to policy

   post-fit use
       │
       ├─ inspect residual and finite-strain diagnostics
       ├─ calculate properties from an accepted or explicit record
       ├─ generate frontend-neutral plot specifications
       └─ export diagnostic or calculated tables

Scientific domains and result slots
-----------------------------------

The public fitting domains are:

``pv``
   Isothermal pressure--volume analysis for volume or a linear cell parameter.

``vt``
   Volume--temperature or length--temperature analysis.

``pvt``
   Coupled pressure--volume--temperature analysis.  At the current checkpoint,
   the public P--V--T fitting workflow accepts the volume target only.

Energy--volume formulations are shared numerical core functionality used by
QHA, but there is no public standalone E--V fitting workflow in the EOS
subsystem.

A result slot is identified by domain and target, for example:

.. code-block:: text

   pv/volume
   pv/a
   vt/volume
   vt/c
   pvt/volume

Several immutable records may exist in the same slot, but only one may be the
current accepted result.  Different targets occupy independent slots, so a
single archive may legitimately contain accepted fits for volume and several
cell axes.

Input normalization and preserved provenance
--------------------------------------------

The reader is content-directed: the filename suffix does not select the
parser.  Column declarations determine which coordinates, targets,
uncertainties, and groups are present.

Unit precedence is:

#. explicit command or specification override;
#. declaration in the data file;
#. EOS default units.

Internally, Quantas uses GPa, angstrom or cubic angstrom where appropriate, and
kelvin.  Relative quantities such as ``V/V0`` and ``L/L0`` remain
dimensionless.  Standard uncertainties are converted with the same scale as
their corresponding observations.

The archive retains both:

- original text, raw arrays, unit labels, grouping, and exclusion markers;
- normalized float64 arrays used by the calculation.

This separation is important for auditability: normalization never destroys
the representation supplied by the user.

Targets, groups, and non-destructive data selection
---------------------------------------------------

``volume``, ``a``, ``b``, and ``c`` are separate scientific targets.  ``all``
expands only to targets that are valid and non-constant in the dataset.
Missing columns, constant independent coordinates, or unsupported target/domain
combinations are rejected during preflight.

EOS data selection is deliberately non-destructive.  The original dataset is
always archived, while each fit request carries its own boolean mask and
selection provenance.

A specification job may select data through:

- the default input selection or all rows;
- included or excluded groups;
- included or excluded one-based data rows;
- combinations of the above in a defined order.

This permits controlled experiments such as fitting independent experimental
groups, testing the influence of one high-pressure subset, or repeating a fit
without an obvious outlier.  Quantas does not decide that an observation is an
outlier merely because its residual is large.  Exclusion is always an explicit
scientific action.

Model construction
------------------

Model selection is separated into three layers.

Scientific domain
   P--V, V--T, or P--V--T determines the coordinates and response used by the
   fit.

Analytical formulation
   The request identifies one implemented family, order, and optional variant.
   P--V--T jobs additionally identify the coupling strategy and, for
   thermal-pressure models, the oscillator model.

Parameter state
   Each physical parameter is classified as free, fixed, implied by the model,
   derived, or unavailable.  Initial values and bounds apply only where they
   are meaningful.

Quantas validates model/order combinations before fitting.  An orderless model
cannot acquire an arbitrary order merely because a generic option was supplied,
and a model-implied parameter is not converted into a free parameter by an
initial-value override.

The equations themselves are not repeated here; see :doc:`../theory/eos` for
scientific definitions and :doc:`../tutorials/eos_pv`,
:doc:`../tutorials/eos_vt`, and :doc:`../tutorials/eos_pvt` for comparisons.

P--V implementation choices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

P--V fitting supports volumetric and linear targets.  For a linear target,
Quantas applies the established EOSFit convention internally: the length is
mapped to an auxiliary cubic quantity for evaluation by a volumetric EOS, then
reported using linear physical parameters.  The auxiliary cube is not treated
as a crystallographic volume.

The principal implementation decisions are:

- family and order are explicit;
- constraints use stable physical parameter names rather than positional
  coefficients;
- normalized-pressure quantities are produced as diagnostics, not as a
  different fitting objective;
- pressure is the dependent response for the public P--V fit;
- uncertainties in pressure and volume enter only through the selected
  regression method.

V--T implementation choices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

V--T fitting accepts volumetric or linear targets and keeps the reference
state, temperature scale, and model variant explicit.  Quantas evaluates the
chosen thermal-expansion formulation directly; it does not first smooth the
observations with an unrelated generic polynomial.

Temperature uncertainties affect the objective only in methods that account
for independent-coordinate error.  This is important because similar fitted
volumes may still imply noticeably different thermal-expansion derivatives.

Some V--T formulations have restricted physical domains.  Model evaluation
therefore distinguishes an invalid state from ordinary non-convergence.  A
finite last iterate outside the valid temperature domain is diagnostic evidence,
not a successful fit.

P--V--T implementation choices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A P--V--T request composes:

- one P--V reference model;
- either one V--T model and a coupling prescription, or a thermal-pressure
  model that supplies its own thermal contribution;
- a common set of reference-state and normalization parameters.

The current objective fits pressure residuals from one volume dataset.  It does
not combine pressure, volume, heat capacity, and adiabatic bulk modulus into a
multi-observable objective.

For Mie--Gruneisen--Debye fitting, normalization is part of the scientific
request rather than a hidden conversion.  Depending on the volume basis, the
request must state the number of atoms per cell or provide a chemical formula
and number of formula units.  An incorrect normalization can yield a
numerically excellent but physically meaningless thermal pressure; Quantas
therefore rejects incomplete MGD normalization before fitting.

The full MGD and the ``q-compromise`` variant are distinct models.  The latter
omits ``q`` according to its specific approximation; it is not represented by
silently fixing a hidden arbitrary value.

Regression methods and their interpretation
-------------------------------------------

All regressions return a common passive result and diagnostic contract, but
they do not represent the same statistical assumptions.

Ordinary least squares (OLS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OLS minimizes unweighted residuals in the dependent response.  Available
uncertainty columns are preserved and reported but do not enter the objective.

Use OLS when:

- uncertainties are unavailable or not comparable;
- a reproducible unweighted baseline is needed;
- the purpose is an initial model comparison.

Do not interpret OLS parameter standard errors as though externally calibrated
measurement variances had been supplied.

Weighted least squares (WLS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WLS uses the standard uncertainty of the dependent response.  Every selected
observation must therefore provide a valid positive uncertainty for that
response.

Use WLS when response uncertainties are meaningful and independent-coordinate
uncertainties are negligible relative to them.

Iterative effective variance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Effective variance augments the response variance by projecting uncertainty in
independent coordinates through local model derivatives.  Because those
derivatives depend on the current model parameters, the weights are updated
iteratively around repeated weighted fits.

This method is useful when both coordinate and response uncertainties matter,
but it remains a local variance projection rather than a full orthogonal
probability model.  Strong non-linearity, very large coordinate uncertainty, or
poor initial parameters can make the iteration difficult.

Orthogonal distance regression (ODR)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ODR treats corrections along all fitted coordinates and requires positive
uncertainties for them.  It is provided through the required ``odrpack``
runtime dependency installed with the base Quantas package.

ODR may be appropriate when coordinate uncertainties are scientifically
important and the orthogonal error model is defensible.  It is not
automatically superior to effective variance; the result depends on the
meaning of the supplied uncertainties and on the geometry of the model.

Choosing a regression
~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Initial regression guide
   :header-rows: 1
   :widths: 25 25 50

   * - Available information
     - Suggested first fit
     - Interpretation
   * - No reliable uncertainties
     - OLS
     - Unweighted scientific baseline.
   * - Reliable response uncertainty only
     - WLS
     - Heteroscedastic response errors enter explicitly.
   * - Reliable response and coordinate uncertainties
     - Effective variance
     - Practical derivative-projected treatment.
   * - Full orthogonal-error interpretation justified
     - ODR
     - Explicit coordinate corrections through the required ODRPACK95 backend.
   * - Unsure whether uncertainties are meaningful
     - Compare OLS and weighted fits
     - Large changes indicate that the uncertainty model controls the result.

A solver is never selected implicitly because uncertainty columns happen to be
present.

Initial values, fixed parameters, and bounds
--------------------------------------------

Non-linear EOS fitting often requires a physically reasonable starting state.
Quantas separates three concepts that are sometimes conflated.

Initial value
   Starting point for a free parameter.  It does not constrain the final
   result.

Fixed value
   Removes the parameter from the free numerical vector.  It changes the model
   hypothesis and the effective degrees of freedom.

Bound
   Defines an allowed interval for a free parameter.  A bound hit is a
   diagnostic condition and may indicate weak support, a poor starting state,
   or an inappropriate model.

Constraints use stable physical names and may be global or target-specific in a
specification.  Model-implied or derived parameters cannot be made free through
a contradictory declaration.

Good practice is to begin with as few constraints as the data support, inspect
correlations and residuals, and fix a parameter only when independent
scientific evidence justifies doing so.  A fixed parameter can stabilize a fit,
but it also transfers its assumed value directly into the remaining
parameters.

Covariance scaling
------------------

Quantas separates the numerical covariance from the policy used to scale it.
The available policies are:

``absolute``
   Treat supplied uncertainties as absolute and do not rescale by the observed
   residual variance.

``reduced-chi-square``
   Scale covariance by the reduced chi-square.

``inflate-only``
   Apply reduced-chi-square scaling only when it increases the covariance.
   Weighted EOS fits use this EosFit-like policy by default.

``inflate-only`` avoids shrinking parameter uncertainty merely because a small
dataset happens to produce a reduced chi-square below unity.  It does not prove
that the measurement uncertainties are correct.

The batch specification
-----------------------

Simple one-off fits can be described with command options.  Advanced EOS work
is better expressed with ``QUANTAS EOS SPEC 1`` because it records the complete
ordered experiment.

The specification supports:

- input-unit overrides;
- global defaults;
- domain defaults;
- multiple named jobs;
- model and coupling choices;
- solver-specific options;
- fixed values, initial values, and bounds;
- observation and group selections;
- acceptance and replacement policy;
- report detail;
- batch failure policy.

Resolution order is:

.. code-block:: text

   built-in defaults
       → [defaults]
       → [defaults.<domain>]
       → [job <name>]

The resolved plan, not merely the original text, is stored in the archive.  A
``--dry-run`` resolves and validates the complete plan without executing a fit
or creating HDF5.  It should be used whenever a specification contains many
jobs, selections, or constraints.

Batch ordering, acceptance, and replacement
-------------------------------------------

Jobs are executed in declaration order.  For each job Quantas:

#. fits the selected dataset;
#. appends the complete attempt to HDF5;
#. accepts it only if the fit succeeded and ``accept`` is true;
#. continues or stops according to the batch failure policy.

Acceptance is slot-specific.  A plan cannot accept two results into the same
slot unless the later job explicitly sets ``replace_accepted``.  Replacement is
never inferred from job order.

Useful patterns include:

Candidate comparison
   Run several models with ``accept = no`` and inspect their records before
   choosing one through the session API or a later controlled workflow.

Ordered replacement
   Accept a baseline, then explicitly replace it with a later job only when the
   specification itself encodes that decision.

Independent targets
   Accept volume, ``a``, ``b``, and ``c`` fits in separate slots within one
   archive.

A failed job can never be accepted.  Its last finite iterate, if available, is
stored only as diagnostic information.

Failure policy
--------------

``stop`` is the default.  The first failed job terminates later execution after
persisting the failed record.

``continue`` executes remaining independent jobs.  This is useful for a broad
comparison plan, but the overall batch is still unsuccessful if any executed
job failed.

``continue`` should not be used to hide a required dependency between jobs.
Batch jobs do not automatically initialize themselves from earlier jobs; such
stateful reuse belongs to the persistent Session workflow.

Persistent archives and immutable history
-----------------------------------------

An EOS HDF5 archive contains:

- normalized and raw dataset representations;
- one append-only record for every fit attempt;
- request, model, solver, selection, and constraint metadata;
- fitted parameters, covariance, residuals, and diagnostics;
- current state for every domain/target slot;
- materialized accepted-result views;
- acceptance, rejection, candidate, and note events;
- batch manifest and creator provenance.

The accepted view is mutable, but the underlying fit records are immutable.
Revoking or replacing acceptance changes the slot state and appends an event;
it does not rewrite the original record.

The archive therefore distinguishes:

``last record``
   Most recently attempted fit in a slot.

``accepted record``
   Fit currently selected for ordinary post-fit operations.

``candidate record``
   Bookmarked alternative that does not change acceptance.

This distinction is central to Quantas GUI operation and to reproducible
scientific review.

The persistent Session API
--------------------------

:class:`quantas.api.eos.Session` exposes the mutable decision workflow without
embedding Click or terminal rendering.  It can:

- create or resume an archive;
- fit another request;
- accept or reject a record;
- revoke acceptance;
- mark candidates;
- attach notes;
- reuse compatible values from an earlier record as an initial guess.

Initial-guess reuse is matched by stable physical parameter name.  Explicit
constraints in the new request take precedence, and implied, fixed, derived,
missing, or out-of-bound values are not copied.

The Session reports archive size through an advisory threshold of 100 MiB by
default.  This warning is operational only: it does not change the calculation,
prevent writing, or define a schema limit.

Fit success versus scientific adequacy
--------------------------------------

Numerical convergence means that the selected backend met its termination
criteria.  It does not establish that:

- the chosen EOS family is appropriate;
- the data cover enough compression or temperature range;
- a high-order parameter is resolved;
- residuals are structureless;
- uncertainties are calibrated;
- extrapolation is safe;
- parameter correlations are acceptable.

Quantas therefore preserves and reports several diagnostic layers.

Parameter diagnostics
   Standard errors, bounds, bound hits, fixed/implied state, covariance, and
   correlation.

Residual diagnostics
   Physical residuals, standardized residuals when uncertainties permit,
   RMSE, weighted objective, chi-square, and reduced chi-square.

Model-specific diagnostics
   Finite strain and normalized pressure for compatible P--V models, P--T
   coverage for P--V--T, and validity checks for derived states.

Solver diagnostics
   Backend, termination category, native status, evaluation counts, starting
   and last parameter states, numerical ranges, and debug trace metadata.

A successful record should be accepted only after the relevant diagnostic
plots and tables have been examined.

Failed fits and last iterates
-----------------------------

Quantas does not reinterpret a finite last iterate as a successful result.
Failures are categorized, for example, as:

- invalid input or model state;
- unavailable or invalid required runtime backend;
- evaluation or iteration limit;
- backend exception;
- backend-declared non-convergence;
- non-finite or out-of-domain model evaluation.

At standard verbosity the report emphasizes the immediate scientific cause.
Debug verbosity adds backend status, numerical ranges, parameter trajectories,
bounds, and the last valid state where available.

The last iterate may help diagnose or initialize a later fit, but it is never
accepted automatically and is never reported as a converged parameter set.

Post-fit diagnostics
--------------------

``quantas eos diagnose`` operates on the accepted record for a slot or on an
explicit record ID.  It can report and export:

- selected and excluded observations;
- predicted response and physical residual;
- standardized residual where defined;
- model-specific finite-strain quantities;
- normalized pressure where applicable;
- grouping and selection provenance.

An explicit record ID is valuable when comparing rejected or candidate fits.
When a slot is omitted, the archive must contain exactly one accepted result;
otherwise the choice would be ambiguous and Quantas requires the user to state
it.

Property calculation
--------------------

The post-fit calculator evaluates a stored model without refitting.  Depending
on the domain it can calculate:

- pressure or volume along P--V states;
- volume, length, or expansion quantities along V--T states;
- coupled P--V--T states, isotherms, or isobars;
- applicable derivatives and moduli exposed by the selected fitted model.

Calculation coordinates are validated against the model domain.  Numerical
root finding used for an inverse calculation is a post-fit evaluation detail;
it does not alter the accepted parameters or create a new fit record.

Calculated tables can be exported to CSV.  A calculated value outside the
sampled data range must still be interpreted as model extrapolation even when
the numerical evaluation succeeds.

Plotting implementation
-----------------------

EOS plotting is record-driven.  Quantas first builds frontend-neutral plot
specifications from the selected fit record; Matplotlib is only the renderer.
Available plots depend on the domain and model:

- observed data with fitted curve;
- physical residuals;
- standardized residuals;
- normalized-pressure diagnostic;
- P--T coverage;
- calculated isotherms;
- calculated isobars.

Excluded observations and data groups are preserved as separate plot series.
Error bars are shown only when corresponding uncertainties exist.

The default smooth-curve sampling is 300 points.  This value changes only the
rendered line smoothness; it does not change fitted parameters or diagnostic
residuals.  Increasing it is not a convergence test for the fit.

Defaults and rationale
----------------------

The direct command provides conservative, inspectable defaults:

.. list-table:: Principal direct-command defaults
   :header-rows: 1
   :widths: 30 20 50

   * - Control
     - Default
     - Rationale
   * - Domain
     - ``pv``
     - Most common standalone EOS use; must be changed explicitly for V--T or P--V--T.
   * - Target
     - ``volume``
     - Unambiguous volumetric baseline.
   * - P--V model
     - Birch--Murnaghan order 3
     - Widely applicable balance between flexibility and parameter count.
   * - V--T model
     - Berman
     - Simple reference thermal-expansion model; variant may be selected explicitly.
   * - P--V--T coupling
     - ``linear``
     - Minimal coupled baseline; more structured couplings remain explicit decisions.
   * - Solver
     - OLS
     - Never assumes uncertainty semantics merely because uncertainty columns are present.
   * - Failure policy
     - ``stop``
     - Prevents later jobs from obscuring the first failed required step.
   * - Covariance scaling for weighted fits
     - ``inflate-only``
     - EosFit-like policy that does not shrink covariance below supplied uncertainty scale.
   * - Accepted result
     - yes for a successful direct job
     - Makes a simple one-job archive immediately usable while retaining the immutable record.

Defaults are starting points, not automatic model recommendations for every
dataset.  The tutorial exercises demonstrate why the model and solver should be
compared when the data support such a comparison.

Performance and numerical controls
----------------------------------

EOS fitting is normally inexpensive compared with phonon or spherical-sampling
workflows.  Runtime is governed mainly by:

- number of selected observations;
- number of free parameters;
- model-evaluation cost;
- number of batch jobs;
- effective-variance outer cycles;
- ODR coordinate corrections;
- convergence limits and poor starting values;
- P--V--T thermal integrals and repeated derivative evaluation.

Useful performance practices are:

#. run ``--dry-run`` before a large specification;
#. start with OLS or WLS to validate the model and initial state;
#. add effective variance or ODR only when justified by uncertainty semantics;
#. compare a few candidate models before launching a large combinatorial plan;
#. avoid unnecessarily broad parameter bounds;
#. use explicit initial values for difficult coupled P--V--T models;
#. calculate dense curves only after an accepted fit exists.

``max_iterations``, ``ftol``, ``xtol``, and ``gtol`` control numerical
termination, not scientific accuracy.  Tightening them cannot repair an
underconstrained model or incorrect data normalization.

EOS has no batching control analogous to SEISMIC ``batch_size=512``.  The
number 512 is not an EOS convergence parameter.

Recommended staged workflow
---------------------------

For a new dataset, the following sequence is robust.

1. Validate the input
~~~~~~~~~~~~~~~~~~~~~

- inspect units, target columns, uncertainties, and groups;
- confirm that relative and absolute quantities are declared correctly;
- use a dry run to validate domain/model compatibility and selections.

2. Establish a baseline
~~~~~~~~~~~~~~~~~~~~~~~

- fit the simplest scientifically plausible formulation;
- use OLS when uncertainty semantics are not yet established;
- inspect residuals, parameter state, and data coverage.

3. Test statistical assumptions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- compare OLS with WLS or effective variance when uncertainties are available;
- use ODR only when orthogonal coordinate uncertainty is defensible;
- investigate large solver-dependent shifts rather than averaging them.

4. Test model sensitivity
~~~~~~~~~~~~~~~~~~~~~~~~~

- compare adjacent orders or alternative families;
- inspect whether additional parameters are resolved;
- compare interpolation inside the data range before considering
  extrapolation.

5. Accept deliberately
~~~~~~~~~~~~~~~~~~~~~~~

- accept one record per slot only after diagnostics are satisfactory;
- retain alternatives as candidate or rejected records with notes;
- never replace an accepted result implicitly.

6. Calculate and plot
~~~~~~~~~~~~~~~~~~~~~

- use the accepted result for ordinary downstream work;
- use explicit record IDs for scientific comparisons;
- label all extrapolated calculations clearly.

Decision guide
--------------

.. list-table:: EOS implementation decisions
   :header-rows: 1
   :widths: 28 27 45

   * - Situation
     - Recommended action
     - Reason
   * - One simple fit
     - Direct ``eos run`` options
     - Minimal reproducible archive with one accepted slot.
   * - Several models or solvers
     - Strict batch specification
     - Ordered plan, inherited defaults, explicit acceptance, and full provenance.
   * - Need interactive scientific decisions
     - Persistent Session API
     - Acceptance, rejection, notes, candidates, and guess reuse without CLI prompts.
   * - Multiple axes in one file
     - Separate target jobs
     - Each target has independent parameters and accepted slot.
   * - Suspected outlier
     - Fit explicit alternate selections
     - Preserve the observation and compare sensitivity; do not delete silently.
   * - Reliable response errors only
     - WLS
     - Uses the uncertainty model actually supported by the data.
   * - Errors in both coordinate and response
     - Effective variance, then optionally ODR
     - Compare derivative-projected and orthogonal interpretations.
   * - High-order parameter hits a bound
     - Reconsider order, range, and starting state
     - A bound hit is evidence, not a successful physical determination.
   * - P--V--T MGD fit
     - Declare normalization explicitly
     - Thermal pressure depends on the number of vibrational degrees of freedom.
   * - Several accepted slots exist
     - Pass ``--slot`` downstream
     - Avoid ambiguous record selection.
   * - Compare a rejected fit
     - Pass ``--record-id``
     - Immutable history remains available independently of acceptance.
   * - Need smoother plot lines
     - Increase curve points
     - Rendering-only control; does not change the fit.

Common interpretation errors
----------------------------

Avoid the following mistakes:

- treating the last fit attempt as the accepted result;
- assuming a successful solver termination proves the EOS is appropriate;
- selecting a weighted method solely because uncertainty columns exist;
- interpreting fixed parameters as measured by the fitted dataset;
- treating a bound hit as a well-determined estimate;
- deleting excluded observations from the source file and losing provenance;
- comparing RMSE values from different uncertainty models as though they were
  the same objective;
- interpreting normalized-pressure plots as a second independent fit;
- using dense calculated curves to imply dense experimental support;
- applying cell-volume MGD without an explicit vibrational normalization;
- treating an extrapolated post-fit calculation as validated because the root
  finder converged;
- changing plotting resolution or solver tolerance as a substitute for model
  sensitivity analysis.

Interoperability and outputs
----------------------------

EOS accepted records are consumed by EOS post-fit diagnostics, calculators,
and renderers.  Shared energy EOS formulations are also used internally by QHA
and Thermoelasticity, but standalone EOS fit archives are not silently
substituted into those workflows.  Cross-module use must remain explicit and
normalization-compatible.

The public Python boundary is :mod:`quantas.api.eos`.  It exposes passive
requests and options, fitting and batch services, archives and sessions,
post-fit calculation and diagnostics, and neutral plotting contracts.  Click
owns only argument conversion, validation, dispatch, and terminal presentation.

Further reading
---------------

- :doc:`../theory/eos` -- equations, assumptions, strengths, and limitations;
- :doc:`../tutorials/eos` -- complete P--V, V--T, and P--V--T worked examples;
- :doc:`../formats/eos_input` -- keyword-directed dataset format;
- :doc:`../formats/eos_spec` -- strict batch specification;
- :doc:`../formats/eos_hdf5` -- archive organization and history;
- :doc:`../cli/eos` -- exact command syntax;
- :doc:`../api/eos` -- public Python contracts.
