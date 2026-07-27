Thermoelasticity: implementation and workflow
==============================================

Purpose and scope
-----------------

The Thermoelasticity workflow combines two independently calculated datasets:

- a static series of hydrostatically pre-stressed elastic tensors sampled at
  different volumes;
- a QHA result providing the equilibrium volume and, when available, the
  thermal fields required for an isothermal-to-adiabatic conversion.

Quantas applies the quasi-static approximation (QSA):

.. math::

   C^T_{IJ}(P,T)
   \simeq
   C^{\mathrm{cold}}_{IJ}\!\left[V_{\mathrm{QHA}}(P,T)\right].

Temperature therefore enters the isothermal elastic tensor through the QHA
volume path.  The explicit strain derivatives of the phonon frequencies that
appear in a complete quasi-harmonic theory of elasticity are not calculated.
The physical derivation and the distinction between full quasi-harmonic
elasticity and the QSA are discussed in :doc:`../theory/thermoelasticity`.

This page explains how Quantas implements that approximation, why the workflow
is divided into calibration and analysis stages, how its scientific controls
should be selected, and which diagnostics must be checked before interpreting
a pressure-temperature or depth-dependent tensor field.

Questions addressed by the workflow include:

- how should independently calculated elastic tensors be placed in one common
  Cartesian frame?;
- why is the static energy EOS fitted separately from the elastic components?;
- when is a second- or third-order Eulerian finite-strain description
  justified?;
- how does Quantas distinguish numerical convergence from adequate scientific
  support?;
- which uncertainties are propagated, and which are not?;
- how are arbitrary P--T states and geological profiles reconstructed without
  repeating the fits?;
- when can an adiabatic tensor be produced?;
- what makes a requested state an interpolation, a coordinate extrapolation,
  or an elastic-volume extrapolation?;

Computational pipeline
----------------------

The complete workflow is intentionally staged:

.. code-block:: text

   CRYSTAL elastic outputs
       │
       ├─ require PRESSURE-corrected Wallace coefficients
       ├─ validate phase, symmetry, atom order, volumes, and pressures
       ├─ co-rotate all tensors into one reference Cartesian frame
       └─ write a normalized thermoelastic YAML

   QHA input or QHA HDF5
       │
       ├─ validate primitive-cell normalization
       ├─ require sampled static E(V) and equilibrium V(P,T)
       └─ optionally collect C_V and the Cartesian expansion tensor

   calibration: thermoelasticity run
       │
       ├─ fit one cold static reference EOS
       ├─ identify and symmetry-average independent C_IJ observations
       ├─ fit each active component as C_IJ(V)
       ├─ assess scientific support and propagate calibration covariance
       └─ write a reusable fit HDF5 without reconstructing a P--T tensor grid

   post-fit analysis
       │
       ├─ interpolate archived QHA fields to requested states
       ├─ evaluate independent cold finite-strain components at V(P,T)
       ├─ reconstruct complete isothermal stiffness matrices
       ├─ evaluate mechanical stability
       ├─ optionally convert C^T to C^S
       └─ produce point, rectangular-grid, or geological-profile results

   presentation and exchange
       │
       └─ inspect / report / table / export / plot / Elasticity / SEISMIC

The separation is scientifically useful.  The reference EOS and
:math:`C_{IJ}(V)` models are calibrated once.  Changing the requested P--T
grid, evaluating a new depth profile, exporting another tensor condition, or
creating additional figures does not repeat those fits.

What the workflow calculates
----------------------------

The calibrated model predicts:

- symmetry-independent isothermal components;
- complete symmetric :math:`6\times6` isothermal Wallace stiffness matrices;
- propagated one-standard-deviation uncertainties and local component
  covariance matrices;
- density from the QHA equilibrium volume and the normalized cell mass;
- generic positive-definiteness diagnostics;
- optional adiabatic stiffness matrices, thermal-stress vectors, and
  adiabatic corrections;
- independent masks for QHA-coordinate and elastic-volume extrapolation.

It does not calculate:

- phonon frequencies under general anisotropic strain;
- explicit vibrational contributions to individual elastic components;
- intrinsic anharmonic elasticity at fixed volume;
- electronic, magnetic, defect, or configurational contributions unless they
  are already represented by the source calculations;
- phase stability or phase transformations along a profile;
- systematic electronic-structure uncertainty.

Input normalization
-------------------

Why CRYSTAL ``PRESSURE`` is required
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The elastic tensors must be the stress--strain coefficients appropriate to a
hydrostatically pre-stressed solid.  Quantas therefore requires CRYSTAL
elastic calculations performed with the ``PRESSURE`` keyword.  The input
creator checks that:

- the keyword was present;
- its value agrees with the pressure reported for the corrected elastic
  properties within 0.05 GPa by default;
- every sampled point contains a finite symmetric stiffness matrix.

The pressure-dependent term in the cold finite-strain formula is part of the
constitutive expansion.  It is **not** a second pressure correction applied to
input data.  Quantas fits the already corrected Wallace coefficients directly.

Passing uncorrected energy Hessians can produce smooth fits and apparently
reasonable values while remaining inconsistent with finite-pressure wave
propagation and mechanical-stability analysis.

Consistency of the elastic-volume series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The input creator sorts all points by increasing primitive-cell volume and
requires:

- unique resolved volumes;
- identical primitive atomic species and atom ordering;
- the same space-group number;
- compatible Hall number and crystallographic setting when available;
- the same detected elastic symmetry;
- a common structural branch rather than unrelated polymorphs.

The default maximum ordered-atom displacement along the path is 0.5 Å.  The
displacement is evaluated after wrapping fractional differences into the
nearest periodic image and expressing them in the reference lattice.  This
criterion allows ordinary internal relaxation while rejecting a structurally
unrelated path.

A phase transformation, symmetry change, change in atom order, or discontinuous
reconstruction must be treated as a separate calibration domain.  A single
finite-strain fit is not intended to bridge such a discontinuity.

The input-generation defaults are:

.. list-table:: Input-normalization tolerances
   :header-rows: 1
   :widths: 38 22 40

   * - Control
     - Default
     - Role
   * - Structural symmetry ``symprec``
     - :math:`10^{-5}` Å
     - Cartesian tolerance passed to structural symmetry analysis
   * - Angle tolerance
     - -1 degree
     - Request the symmetry library's automatic angle treatment
   * - Elastic-symmetry tolerance
     - :math:`10^{-3}` GPa
     - Detect one common elastic crystal-system pattern
   * - PRESSURE consistency tolerance
     - 0.05 GPa
     - Compare the keyword and reported corrected-tensor pressures
   * - Ordered-atom path tolerance
     - 0.5 Å
     - Reject a structurally unrelated volume path

These values are validation tolerances, not uncertainty estimates.  Loosening
them can permit incompatible source calculations to enter the same fit.

Reference point and frame normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When no reference point is specified, Quantas selects the volume-sorted elastic
point with the smallest absolute applied pressure.  This point defines:

- the reference Cartesian frame;
- the compact reference structure;
- the atom-order comparison;
- the crystallographic labels attached to later directions.

It does **not** define the fitted static :math:`V_0`.  The thermodynamic
reference volume comes from the QHA static energy EOS described below.

Different electronic-structure calculations can print equivalent structures
in Cartesian frames related by a rigid rotation.  Component-wise fitting would
be meaningless if :math:`C_{11}` at one volume referred to a different axis
than :math:`C_{11}` at another volume.  Quantas therefore constructs the
lattice deformation gradient relative to the selected reference and performs
a right polar decomposition:

.. math::

   \mathbf F = \mathbf R\mathbf U.

The proper rotation :math:`\mathbf R` is removed from the lattice and stiffness
components, while the symmetric stretch :math:`\mathbf U` is retained.  The
normalization verifies that:

- :math:`\det\mathbf R=1` within numerical tolerance;
- the stretch is positive definite;
- the sampled cell volume is preserved;
- every co-rotated tensor remains expressed in the same reference frame.

The removed rotation, principal logarithmic strains, source lattice, and
rotation matrix are stored in the YAML and HDF5 provenance.

What is and is not fitted from the elastic outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The elastic-component regressions use **volume** as the independent variable.
The CRYSTAL pressures are retained for validation, provenance, and physical
inspection, but are not the regression coordinate.

Likewise, the static energies printed in the individual SOEC outputs are
retained in the normalized input but do not define the authoritative reference
EOS.  The EOS is fitted to the sampled static energy-volume field stored by
the QHA workflow.  This keeps the cold reference state, cell normalization,
and QHA equilibrium-volume surface internally consistent.

Coupling to QHA
---------------

Required cold-QSA fields
~~~~~~~~~~~~~~~~~~~~~~~~

A reusable calibration requires the QHA payload to contain:

- temperature coordinates;
- pressure coordinates;
- sampled static volumes;
- sampled static energies;
- equilibrium volume :math:`V(P,T)`.

The sampled static energy-volume arrays define the cold reference EOS.  The
equilibrium-volume field provides the thermodynamic path on which the static
elastic model is evaluated.

The QHA temperature and pressure coordinates must be finite and strictly
increasing.  The equilibrium-volume array must have shape
``(n_temperature, n_pressure)``.

Optional adiabatic fields
~~~~~~~~~~~~~~~~~~~~~~~~~

Adiabatic conversion additionally requires:

- isochoric heat capacity :math:`C_V(P,T)`;
- a symmetric Cartesian thermal-expansion tensor
  :math:`\boldsymbol\alpha(P,T)`;
- the same cell normalization as the QHA equilibrium volume.

Quantas converts the QHA heat capacity to J K\ :sup:`-1` per normalized cell
and stores the expansion tensor in K\ :sup:`-1`.  Missing or shape-incompatible
fields do not prevent an isothermal calibration unless ``adiabatic=require``
is selected.

Atomic and cell normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the QHA archive contains primitive atomic numbers, Quantas requires exact
agreement with the thermoelastic reference structure.  Density is not
interpolated from the CRYSTAL values.  Instead, Quantas infers the normalized
cell mass from every CRYSTAL density-volume pair:

.. math::

   m_i = \rho_i V_i,

uses the median mass, and reconstructs density from:

.. math::

   \rho(P,T)=\frac{m}{V_{\mathrm{QHA}}(P,T)}.

A relative mass spread above :math:`10^{-5}` generates a warning.  This check
is useful for detecting inconsistent primitive/conventional cells, formula-unit
normalizations, or density metadata.

QHA HDF5 versus QHA YAML
~~~~~~~~~~~~~~~~~~~~~~~~

``thermoelasticity run`` can consume either:

- an existing native QHA HDF5 archive;
- a QHA or phonon YAML that Quantas must calculate first.

Using an existing QHA HDF5 is normally preferable for production work because
it:

- avoids repeating QHA minimization;
- preserves an already inspected thermodynamic surface;
- makes sensitivity tests of the elastic calibration inexpensive;
- ensures that several thermoelastic runs use exactly the same QHA result.

When a YAML source is supplied, the ``--qha-*`` controls belong to the embedded
QHA calculation, not to the cold finite-strain component fit.  Their scientific
meaning is described in :doc:`qha`.

Volume support at coupling time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Quantas compares every archived QHA equilibrium volume with the sampled
elastic-volume interval.  States outside that interval are recorded in the
calibration archive, but calibration itself is still possible.  The selected
extrapolation policy is enforced later, when a point, grid, or profile is
actually evaluated.

This distinction matters because the QHA coordinate domain and the elastic
volume domain are different concepts:

- a state can lie inside the archived P--T rectangle while its equilibrium
  volume lies outside the elastic series;
- a state can lie outside the archived P--T rectangle yet extrapolate to a
  volume that remains inside the elastic interval;
- a state can violate both domains.

Calibration stage
-----------------

The calibration archive contains models and source fields.  It intentionally
does not contain a reconstructed P--T stiffness grid.

Cold reference EOS
~~~~~~~~~~~~~~~~~~

Quantas first fits an integrated energy EOS to the QHA sampled static
energy-volume data.  The resulting reference parameters are:

- :math:`V_0`, the static zero-pressure volume;
- :math:`K_0`, the static isothermal bulk modulus;
- :math:`K'_0`, its first pressure derivative.

The fit excludes:

- zero-point energy;
- thermal vibrational energy;
- the QHA equilibrium free-energy minimum.

Consequently the reference is a cold static 0 K, 0 GPa state.  Thermal effects
enter later through :math:`V_{\mathrm{QHA}}(P,T)` and, for :math:`C^S`, through
:math:`C_V` and :math:`\boldsymbol\alpha`.

Choosing BM2, BM3, or BM4
~~~~~~~~~~~~~~~~~~~~~~~~~

The CLI exposes three reference EOS choices.

``BM2``
   Fixes :math:`K'_0=4`.  It has the fewest free elastic-curvature parameters
   and can be useful when the static energy range is narrow or cannot support
   an independent :math:`K'_0`.  Its limitation is structural rigidity: a poor
   BM2 reference can bias every component fit through the shared
   :math:`V_0`, :math:`K_0`, and :math:`K'_0`.

``BM3`` — default
   Fits :math:`K'_0` and is the usual balance between flexibility and support
   for a well-sampled static energy curve.  It is the recommended first
   production model.

``BM4``
   Also fits :math:`K''_0`.  It requires a sufficiently broad and precise
   static compression range.  The cold component formula uses
   :math:`V_0`, :math:`K_0`, and :math:`K'_0`; therefore BM4 affects the
   thermoelastic calibration mainly by changing those estimates and their
   covariance.  It does not add a separate :math:`K''_0` term to the current
   :math:`C_{IJ}(V)` expression.

A lower EOS residual alone is not sufficient reason to select BM4.  Inspect:

- parameter uncertainties;
- correlation and covariance;
- physical plausibility of :math:`K'_0` and :math:`K''_0`;
- stability of :math:`V_0` and :math:`K_0` under model changes;
- consequences for the reconstructed elastic field.

Independent components and symmetry averaging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Quantas detects the elastic crystal system during input generation.  At each
sampled volume it then:

1. extracts the symmetry-independent component set;
2. applies the sign conventions associated with symmetry-equivalent entries;
3. averages equivalent entries;
4. records their spread as a diagnostic.

Derived relations such as

.. math::

   C_{66}=\frac{C_{11}-C_{12}}{2}

are reconstructed after fitting from the independent component covariance.
They are not fitted as unrelated observations.

Exact-zero components
~~~~~~~~~~~~~~~~~~~~~

If the maximum absolute symmetry-averaged value of a component does not exceed
``zero_tolerance``, the component is retained as an exact zero.  The default is

.. code-block:: text

   zero_tolerance = 1.0e-10 GPa

No optimizer is called and no artificial uncertainty is assigned.  This is
intended for symmetry-forbidden entries, not for suppressing a physically small
but non-zero component.  Increasing the threshold changes the scientific
model and should not be used as a numerical cleanup tool.

Cold finite-strain component model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each active independent component is represented by two fitted parameters:

- :math:`C^0_{IJ}`, its value at the cold reference state;
- :math:`C'^0_{IJ}=\partial C_{IJ}/\partial P`, its reference pressure
  derivative.

The fixed EOS parameters and the component-specific Wallace coefficient enter
the Eulerian finite-strain expression described in
:doc:`../theory/thermoelasticity`.

The component model is linear in :math:`C^0_{IJ}` and
:math:`C'^0_{IJ}` once :math:`V_0`, :math:`K_0`, and :math:`K'_0` are fixed.
Quantas therefore obtains exact conditional linear least-squares estimates for
support diagnostics and initialisation.  The production fit uses ordinary
least squares.

Why only OLS is currently available
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current CRYSTAL SOEC interface does not provide observation uncertainties
for each elastic component.  Quantas therefore does not pretend that a
statistically meaningful WLS or ODR problem exists.  Every sampled volume has
equal regression weight.

The reported residual sum of squares and reduced unweighted chi-square are
therefore descriptive measures in GPa\ :sup:`2`; they are not probability-based
chi-square statistics.  A small residual does not account for systematic DFT
error, incomplete strain convergence, or model-form discrepancy.

Second- and third-order finite strain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``finite_strain_order=2`` retains terms through first order in the Eulerian
strain inside the constitutive polynomial.  It is more constrained and can be
useful when:

- only a modest strain interval is sampled;
- the data do not resolve order sensitivity;
- a third-order description produces unstable parameters.

``finite_strain_order=3`` — default — includes the quadratic strain term
required by the thermodynamically consistent third-order cold finite-strain
formulation.  It is generally preferred when the volume range is sufficiently
broad and the response shows resolved curvature.

Both orders still fit only :math:`C^0_{IJ}` and :math:`C'^0_{IJ}`.  The
quadratic coefficient of the third-order expression is constrained by these
component parameters, the reference EOS, and the Wallace hydrostatic term; it
is not an additional freely fitted coefficient.

A third-order fit should not be accepted merely because it is available.
Compare:

- leave-one-out sensitivity;
- second-versus-third-order parameter changes;
- residual structure;
- reference-volume bracketing;
- behavior near the ends of the sampled volume interval.

Scientific support diagnostics
------------------------------

Numerical convergence and scientific support are evaluated separately.  Every
active component receives a quality classification:

``supported``
   No configured support criterion is violated.

``caution``
   The fit is numerically usable but one or more heuristic support criteria are
   weak.

``unsupported``
   The number or rank of observations cannot support the two-parameter model,
   or the design matrix is singular.

``not_applicable``
   The component was retained as an exact zero.

The assessment does not alter observations, fitted parameters, or predictions.
It records why a result may be fragile.

Default support criteria
~~~~~~~~~~~~~~~~~~~~~~~~

The standard defaults are:

.. list-table:: Standard component-fit support thresholds
   :header-rows: 1
   :widths: 48 22 30

   * - Criterion
     - Default
     - Interpretation
   * - Preferred minimum observations
     - 4
     - More observations than the two fitted parameters
   * - Minimum Eulerian strain span
     - 0.005
     - Avoid a nearly singular local-volume calibration
   * - Maximum normalized design condition number
     - :math:`10^6`
     - Detect poorly separated parameter directions
   * - Maximum relative symmetry-equivalent spread
     - 0.01
     - Detect disagreement among entries expected to be equivalent
   * - Maximum scale-aware leave-one-out parameter change
     - 0.5
     - Detect domination by one sampled point
   * - Maximum second/third-order parameter change
     - 0.5
     - Detect strong finite-strain-order sensitivity

The relative-residual diagnostic uses a denominator floor of
:math:`10^{-8}` GPa by default.  This prevents division by an almost-zero
component from creating meaningless enormous ratios; it does not change the
OLS objective or fitted parameters.

Additional recorded quantities include:

- degrees of freedom;
- minimum, maximum, and span of Eulerian strain;
- whether :math:`V_0` is bracketed by the elastic volumes;
- distance of an unbracketed :math:`V_0` from the interval;
- design rank and leverage of each observation;
- maximum relative residual;
- maximum absolute symmetry spread.

Reference-volume bracketing is important.  If :math:`V_0` lies outside the
elastic interval, :math:`C^0_{IJ}` is already an extrapolated parameter even
when every later QHA state lies inside the sampled elastic volumes.

Validation presets
~~~~~~~~~~~~~~~~~~

The CLI groups the support and policy controls into three named presets.

.. list-table:: Validation presets
   :header-rows: 1
   :widths: 20 27 27 26

   * - Control
     - ``standard``
     - ``strict``
     - ``exploratory``
   * - Extrapolation
     - warn
     - fail
     - allow
   * - Unsupported fit
     - warn
     - fail
     - allow
   * - Unstable tensor
     - warn
     - fail
     - allow
   * - Minimum observations
     - 4
     - 5
     - 3
   * - Minimum strain span
     - 0.005
     - 0.0075
     - 0
   * - Maximum design condition
     - :math:`10^6`
     - :math:`10^5`
     - :math:`10^{12}`
   * - Maximum symmetry spread
     - 0.01
     - 0.005
     - 1
   * - Maximum leave-one-out change
     - 0.5
     - 0.25
     - 10
   * - Maximum order change
     - 0.5
     - 0.25
     - 10

``standard``
   Balanced production behavior.  Weak support, extrapolation, and instability
   are retained with warnings and complete diagnostics.

``strict``
   Appropriate when an automated production pipeline must reject unsupported
   fits, any requested extrapolation, or unstable reconstructed states.  A
   ``caution`` classification remains a diagnostic; only ``unsupported`` is a
   fit-quality failure.

``exploratory``
   Retains finite numerical results for method development and sensitivity
   analysis.  It does not make extrapolated or weakly supported predictions
   more reliable.

Fit failure and fit quality are different
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A component fit can fail numerically, or it can converge but be scientifically
unsupported.

- ``fit_failure_policy`` controls failed active fits;
- ``quality_policy`` controls converged fits classified as unsupported.

All active components are attempted.  Any failed active component makes the
calibration payload incomplete and prevents post-fit analysis.  ``raise``
aborts with an exception; the non-raising policies preserve diagnostic records
for inspection.  They do not create a scientifically complete tensor by
silently replacing a failed component.

Uncertainty propagation
-----------------------

Quantas propagates several available first-order uncertainty sources.  The
result should be interpreted as the uncertainty of the fitted numerical model
under its stated assumptions, not as a complete uncertainty budget for the
material.

Component-fit covariance
~~~~~~~~~~~~~~~~~~~~~~~~

Each active OLS fit provides covariance for
:math:`(C^0_{IJ},C'^0_{IJ})`.  Analytical Jacobians propagate that covariance
to every evaluated volume.

Shared reference-EOS covariance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All elastic components share the same fitted :math:`V_0`, :math:`K_0`, and
:math:`K'_0`.  Their uncertainties therefore create correlations between
otherwise separate component predictions.

Quantas includes two contributions:

1. the direct analytical sensitivity of the component formula to the EOS
   parameters;
2. the indirect sensitivity caused because the best-fitting
   :math:`C^0_{IJ}` and :math:`C'^0_{IJ}` change when the fixed EOS changes.

The second contribution is evaluated by symmetric finite differences of exact
conditional linear least-squares refits.  The resulting shared covariance is
retained in the local independent-component covariance matrix.

Disable this contribution only for a deliberate sensitivity test with
``--no-eos-covariance``.  Removing it usually underestimates correlations and
uncertainties rather than improving the fit.

QHA equilibrium-volume uncertainty
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the QHA archive provides a valid uncertainty field matching
:math:`V(P,T)`, Quantas propagates it through the analytical derivative
:math:`\partial C_{IJ}/\partial V`.  This contribution is optional through
``--volume-uncertainty/--no-volume-uncertainty``.

Interpolated uncertainty values are constrained to remain non-negative.  An
absent, invalid, or shape-incompatible QHA volume uncertainty is excluded with
a warning rather than replaced by an invented value.

Derived components and tensor covariance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Full stiffness matrices are reconstructed through linear symmetry relations.
The complete independent-component covariance is used, so correlations are
retained in derived entries.  For example, the uncertainty of
:math:`C_{66}=(C_{11}-C_{12})/2` includes the covariance between
:math:`C_{11}` and :math:`C_{12}`.

Adiabatic uncertainty
~~~~~~~~~~~~~~~~~~~~~

When enabled, Quantas propagates available uncertainties in:

- the isothermal stiffness tensor;
- volume;
- :math:`C_V`;
- the Cartesian expansion tensor.

The current method is a first-order delta approximation.  Cross-covariances
between QHA heat capacity, volume, and expansion are not available in the
archive and are assumed zero.

Uncertainty not represented
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The propagated standard deviations do not include, unless already encoded by
the source data:

- exchange-correlation functional error;
- basis-set or k-point incompleteness;
- uncertainty in imposed strains or stress corrections;
- systematic phonon error;
- QSA model-form error;
- uncertainty caused by a possible phase transition;
- correlation between static elastic and QHA calculations.

A narrow plotted band must therefore not be interpreted as proof of physical
accuracy.

The reusable calibration archive
----------------------------------

A successful ``thermoelasticity run`` writes:

- the cold reference EOS and covariance;
- every independent-component fit and diagnostic;
- the common tensor-frame provenance;
- the original elastic-volume observations;
- the archived QHA P--T coordinates and equilibrium volumes;
- optional QHA volume uncertainties;
- optional per-cell :math:`C_V` and Cartesian expansion tensors;
- density and normalized-cell mass metadata;
- validation policies, warnings, and meaningful events.

It does **not** write :math:`C_{IJ}(P,T)` arrays.  This is intentional:

- the archive remains compact;
- alternative grids and profiles reuse the same calibration;
- calibration quality is not confused with an arbitrary presentation grid;
- Quantas GUI can evaluate only the states currently requested.

Use:

.. code-block:: console

   quantas thermoelasticity inspect FIT.hdf5

before post-fit analysis.  The inspector reports the workflow stage, fit
support, tensor frame, QHA and elastic support, available tensor conditions,
and suitable next actions.

Post-fit reconstruction
-----------------------

Archived-field interpolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The post-fit engine interpolates archived QHA fields on their rectilinear
:math:`(T,P)` coordinates using piecewise-linear interpolation.

For a rectangular target grid, Quantas evaluates the Cartesian product of the
requested temperature and pressure axes.  For a geological profile, it
evaluates aligned pairs :math:`(T_i,P_i)` directly.  A profile therefore does
not construct a hidden rectangular grid.

Outside the source coordinate interval, the numerical interpolator uses the
slope of the nearest endpoint interval and marks the state as extrapolated.  A
singleton source axis is treated as constant.  The selected extrapolation
policy determines whether the marked state is rejected, warned about, or
retained.

Piecewise-linear interpolation has useful properties:

- it reproduces archived nodes exactly;
- it introduces no high-order oscillation between nodes;
- it is transparent and inexpensive;
- it preserves arbitrary trailing tensor dimensions.

Its limitations are equally important:

- derivatives are discontinuous at source-grid nodes;
- a sparse QHA grid can hide curvature between nodes;
- linear endpoint extrapolation is not a physical high-P or high-T model;
- refining only the target grid cannot recover information absent from the
  source QHA grid.

Evaluation at the QHA volume
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For every requested state, Quantas performs the following sequence:

1. interpolate :math:`V(P,T)` and optional :math:`\sigma_V`;
2. evaluate all independent cold finite-strain component models at that
   volume;
3. propagate component, EOS, and volume covariance;
4. reconstruct the complete symmetric isothermal tensor;
5. compute density from normalized mass and volume;
6. evaluate mechanical stability;
7. optionally interpolate adiabatic inputs and calculate :math:`C^S`;
8. store QHA-coordinate and elastic-volume masks independently.

The static component model depends on the requested state only through volume.
Two different :math:`(P,T)` pairs with exactly the same equilibrium volume
therefore produce the same QSA isothermal tensor.  They can still produce
different adiabatic corrections because :math:`T`, :math:`C_V`, and
:math:`\boldsymbol\alpha` differ.

Point, grid, and profile analyses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``analysis point``
   Evaluates one P--T state.  It is the best choice for checking a proposed
   condition, creating an Elasticity/SEISMIC input, or comparing methods before
   calculating a large grid.  No one-point HDF5 is required; an optional text
   file contains the complete tensor and QHA-consistent density.

``analysis grid``
   Evaluates a rectangular P--T field and writes a native HDF5 archive.  Use it
   for contour maps, fixed-pressure/fixed-temperature comparisons, or
   systematic export.  The state count is
   :math:`N_TN_P`.

``analysis profile``
   Evaluates only the ordered states of one depth-dependent path.  It is the
   preferred choice for geological applications because its cost and storage
   scale with the number of depth samples rather than with an enclosing P--T
   rectangle.

Profiles can be read from a complete depth/P/T table, generated from a YAML
composition of pressure and temperature models, selected from a built-in
preset, or created as an explicitly requested linear test path.  The underlying
Earth models and their limitations are described in
:doc:`../theory/earth_profiles`.

Two independent extrapolation masks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every evaluated state carries:

``qha_extrapolation_mask``
   The requested temperature or pressure lies outside the archived QHA
   coordinate rectangle.

``elastic_extrapolation_mask``
   The interpolated equilibrium volume lies outside the sampled elastic-volume
   interval, or is non-finite.

The masks must not be merged when interpreting a result.  They diagnose two
different limitations:

- QHA coordinate extrapolation tests the thermodynamic source surface;
- elastic-volume extrapolation tests the cold finite-strain calibration.

The default ``warn`` policy retains the numerical value and records the
limitation.  ``fail`` rejects any marked requested state.  ``allow`` retains it
without a warning, but the masks remain stored.

Mechanical stability
--------------------

Quantas evaluates the eigenvalues of each reconstructed **isothermal** Wallace
stiffness matrix.  The generic positive-definiteness criterion is:

.. math::

   \lambda_{\min} > \lambda_{\mathrm{tol}},

with a default tolerance of 0 GPa.

Each isothermal state is classified as:

- stable;
- unstable;
- indeterminate because the tensor contains invalid values.

Quantas never repairs an unstable tensor by clipping eigenvalues or replacing
components.  The current result stores no separate adiabatic stability field.
For a valid state the implemented adiabatic correction is positive
semidefinite, so it cannot destabilize an already positive-definite
:math:`C^T`.  The selected policy is:

``warn`` — default
   Preserve the state and report the count.

``fail``
   Reject a point, grid, or profile containing any unstable or indeterminate
   state.

``allow``
   Preserve the result without a warning while retaining the masks and minimum
   eigenvalues.

A loss of stability can indicate:

- a real elastic instability;
- extrapolation beyond the calibrated phase branch;
- a poor component fit;
- an inconsistent stress convention;
- an incorrect tensor frame;
- numerical noise near a stability boundary.

The stability diagnostic cannot determine which explanation is correct.

Isothermal-to-adiabatic conversion
----------------------------------

Implemented identity
~~~~~~~~~~~~~~~~~~~~

Quantas converts the isothermal field using:

.. math::

   C^S_{IJ}
   =
   C^T_{IJ}
   +
   \frac{TV}{C_V}\lambda_I\lambda_J,
   \qquad
   \lambda_I=\sum_J C^T_{IJ}\alpha_J.

The expansion tensor is converted to engineering-Voigt strain order
``11, 22, 33, 23, 13, 12``.  The shear entries are doubled:

.. math::

   \alpha_4=2\alpha_{23},\quad
   \alpha_5=2\alpha_{13},\quad
   \alpha_6=2\alpha_{12}.

For valid positive-temperature states the correction is a positive-semidefinite
rank-one update.  It acts through the thermal-stress vector and is not an
independent fit.

Availability modes
~~~~~~~~~~~~~~~~~~

``auto`` — default
   Calculate :math:`C^S` when complete QHA :math:`C_V` and expansion-tensor
   fields are available.  Otherwise preserve an isothermal-only calibration.

``off``
   Do not calculate or archive adiabatic tensors, even if the QHA fields are
   available.

``require``
   Reject calibration if the mandatory QHA fields are absent and reject an
   analysis if any nonzero-temperature state has invalid conversion inputs.

When an adiabatic field can be constructed, exactly zero kelvin is treated as
the thermodynamic limit:

.. math::

   C^S=C^T.

If the complete QHA adiabatic fields are absent and ``auto`` is selected, no
adiabatic field is created, including at zero kelvin.  When the fields exist, a
nonzero-temperature state with non-finite or non-positive :math:`C_V`, an
invalid volume, or a non-finite expansion tensor is represented by NaN in
``auto`` mode and rejected in ``require`` mode.  Quantas does not silently
substitute the isothermal tensor.

Interpretation of :math:`C^S`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The adiabatic conversion is thermodynamically consistent with the available
QHA heat capacity and expansion tensor, but it does not restore the explicit
vibrational strain derivatives omitted by the QSA.  Therefore:

- :math:`C^T` is the cold QSA tensor evaluated at the thermally expanded volume;
- :math:`C^S-C^T` is the standard thermodynamic conversion of that tensor;
- neither quantity is a full anisotropic phonon-QHA elastic calculation.

Choosing the main scientific options
------------------------------------

.. list-table:: Initial option choices
   :header-rows: 1
   :widths: 31 30 39

   * - Situation
     - Recommended starting choice
     - Reason
   * - Ordinary well-sampled calibration
     - BM3, order 3, standard validation
     - Balanced reference EOS and thermodynamically consistent curvature
   * - Narrow static energy range
     - Compare BM2 and BM3
     - Determine whether :math:`K'_0` is genuinely resolved
   * - Broad high-quality compression range
     - Compare BM3 and BM4
     - Test whether the extra EOS curvature is supported
   * - Few elastic volumes or narrow strain span
     - Compare order 2 and 3; add data if possible
     - A higher-order expression cannot replace missing support
   * - Production requires seismic velocities
     - ``adiabatic=require``
     - Prevent accidental use of an isothermal tensor
   * - Only static QSA trends are required
     - ``adiabatic=off``
     - Avoid storing or interpreting an unnecessary conversion
   * - Automated conservative pipeline
     - strict validation and extrapolation ``fail``
     - Stop on unsupported fits or unsupported requested states
   * - Method development
     - exploratory, followed by explicit mask inspection
     - Retain finite results without confusing them with validated predictions

Recommended sensitivity study
-----------------------------

A defensible production analysis should vary the assumptions that can
materially change the result.  A compact study can be organized as follows.

Reference EOS and finite-strain order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Evaluate a few representative states with:

.. list-table:: Calibration sensitivity table to complete
   :header-rows: 1

   * - Reference EOS
     - Strain order
     - :math:`V_0`
     - :math:`K_0`
     - :math:`K'_0`
     - :math:`C_{11}`
     - :math:`C_{33}`
     - Minimum eigenvalue
   * - BM2
     - 2
     - ...
     - ...
     - 4 fixed
     - ...
     - ...
     - ...
   * - BM2
     - 3
     - ...
     - ...
     - 4 fixed
     - ...
     - ...
     - ...
   * - BM3
     - 2
     - ...
     - ...
     - ...
     - ...
     - ...
     - ...
   * - BM3
     - 3
     - ...
     - ...
     - ...
     - ...
     - ...
     - ...
   * - BM4
     - 3
     - ...
     - ...
     - ...
     - ...
     - ...
     - ...

Compare both central states and states near the edge of elastic support.
Differences often appear first in pressure derivatives and boundary behavior,
not at the reference state.

Observation support
~~~~~~~~~~~~~~~~~~~

Inspect for every component:

- residual pattern versus volume;
- leverage;
- maximum symmetry-equivalent spread;
- leave-one-out changes;
- order sensitivity;
- whether :math:`V_0` is bracketed.

A useful exercise is to remove one endpoint, recalibrate, and compare a few
states.  This is more informative than relying only on the full-data residual.
Do not remove points solely because they worsen a preferred model; first check
whether they reveal real curvature, structural change, or a source-calculation
problem.

QHA sensitivity
~~~~~~~~~~~~~~~

The thermoelastic result inherits every important assumption of the QHA volume
surface.  Compare, where scientifically possible:

- ``freq`` and ``td`` interpolation;
- polynomial and EOS minimization;
- alternative QHA polynomial degrees;
- thermal-expansion methods when adiabatic tensors are required.

If two QHA routes produce similar :math:`V(P,T)` but different expansion
tensors, the QSA isothermal stiffness may agree while the adiabatic correction
differs.

Interpolation resolution
~~~~~~~~~~~~~~~~~~~~~~~~

Refining the **target** grid tests only the presentation and localization of
features.  To test the interpolation itself, refine the **source QHA** grid and
repeat the comparison at identical physical states.

A practical table is:

.. list-table:: QHA-grid sensitivity table to complete
   :header-rows: 1

   * - Source QHA spacing
     - State
     - Volume
     - :math:`C_{11}^T`
     - :math:`C_{11}^S`
     - :math:`C_{44}^T`
   * - coarse
     - 0 GPa, 300 K
     - ...
     - ...
     - ...
     - ...
   * - coarse
     - 5 GPa, 1000 K
     - ...
     - ...
     - ...
     - ...
   * - refined
     - 0 GPa, 300 K
     - ...
     - ...
     - ...
     - ...
   * - refined
     - 5 GPa, 1000 K
     - ...
     - ...
     - ...
     - ...

Depth-profile sensitivity
~~~~~~~~~~~~~~~~~~~~~~~~~

A geological profile is an external path through the same calibrated field.
Compare plausible pressure and temperature models separately.  A change in the
profile is not a change in the mineral model, and a stable mineral tensor does
not prove that the phase is thermodynamically stable at that depth.

Performance and acceleration
----------------------------

Calibration cost
~~~~~~~~~~~~~~~~

Calibration includes:

- one nonlinear static energy EOS fit;
- one two-parameter OLS fit for each active independent component;
- leave-one-out exact conditional fits;
- one alternate-order exact conditional fit;
- reference-EOS sensitivity refits for uncertainty propagation.

For ordinary crystal symmetries and tens of elastic volumes, calibration is
usually inexpensive compared with generating the electronic-structure and QHA
data.  Reuse the fit archive rather than recalibrating for each plot or
profile.

``max_iterations`` is only a ceiling on model evaluations for the reference
EOS and component fits.  Increasing it can help a fit that stops prematurely;
it cannot cure a singular design, inadequate strain span, inconsistent data,
or an unsuitable model.

Post-fit cost and memory
~~~~~~~~~~~~~~~~~~~~~~~~

Post-fit evaluation is vectorized.  Its work scales approximately with:

.. math::

   N_{\mathrm{state}}N_{\mathrm{component}},

plus reconstruction and eigensolution of one :math:`6\times6` matrix per
state.  A rectangular grid has
:math:`N_{\mathrm{state}}=N_TN_P`; a profile has only the number of depth
samples.

Memory can be more important than CPU time because a grid archive may contain:

- independent component values and covariance;
- full isothermal tensors and uncertainties;
- full adiabatic tensors and uncertainties;
- correction tensors and thermal-stress vectors;
- stability and extrapolation masks.

For large studies:

1. inspect and calibrate once;
2. test a few points first;
3. use a profile rather than an enclosing grid when only one geological path is
   needed;
4. use a coarse exploratory grid before a publication grid;
5. disable adiabatic uncertainty propagation only when its absence is
   scientifically acceptable;
6. select ``adiabatic=off`` at calibration when only isothermal tensors are
   scientifically required; table and export selectors change presentation, not
   the tensor fields stored by an existing analysis archive;
7. generate fit and analysis plots after validating the numerical archive.

No ``512`` convergence parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Thermoelasticity has no scientific or numerical control with a default value of
512.  The post-fit engine evaluates requested state arrays vectorially and does
not expose a chunk-size control in the CLI.

The default value 512 discussed in :doc:`seismic` is the number of propagation
directions processed per vectorized SEISMIC batch.  It does not belong to the
thermoelastic calculation and must not be transferred here.

The controls that actually change thermoelastic resolution or support are:

- the elastic-volume sampling;
- the QHA source P--T grid;
- the reference EOS family;
- the finite-strain order;
- the requested point, grid, or profile coordinates.

Controls such as plot DPI, contour levels, line width, marker style, or CSV
format do not alter scientific results.

Diagnostics and reporting
-------------------------

``standard`` reporting
   Summarizes the reference EOS, component parameters, support classification,
   support bounds, and final warnings.

``extended`` reporting
   Adds detailed residual, uncertainty, frame, and stability information useful
   for scientific review.

``debug`` reporting
   Retains bounded solver traces and detailed policy decisions.  It is intended
   for diagnosing a problematic calibration, not as the default publication
   report.

Operational progress events are emitted while active components are processed.
Progress events are observer-only and are not persisted as scientific history.
Warnings, meaningful messages, structured fit results, and provenance are
stored in the HDF5 envelope.

Recommended production workflow
-------------------------------

A robust sequence is:

1. calculate a consistent hydrostatic elastic-volume series with CRYSTAL
   ``PRESSURE``;
2. generate and inspect the normalized thermoelastic YAML;
3. calculate and validate QHA independently;
4. confirm that the elastic volumes bracket the relevant QHA volume range;
5. calibrate with BM3, third-order strain, and standard validation;
6. inspect reference-EOS and component support diagnostics;
7. compare at least one alternative EOS or finite-strain order;
8. evaluate a few representative points with extrapolation ``fail``;
9. check isothermal/adiabatic availability and stability;
10. create the final grid or geological profile;
11. export selected states to Elasticity or SEISMIC;
12. retain the fit and analysis HDF5 archives as the authoritative results.

Common interpretation errors
----------------------------

**Treating CRYSTAL pressure as the fitted coordinate**
   Quantas fits elastic components against volume.  Pressure verifies the
   pre-stressed source state and remains provenance.

**Assuming the reference elastic point is the EOS reference volume**
   The selected point defines the common frame.  :math:`V_0` comes from the QHA
   static energy EOS.

**Applying a second pressure correction**
   Input tensors must already be PRESSURE-corrected Wallace coefficients.  The
   Wallace term in the finite-strain equation is not another correction of the
   observations.

**Interpreting a calibration archive as a completed P--T grid**
   Calibration stores models and QHA source fields.  Point, grid, or profile
   analysis must still be requested.

**Confusing coordinate and volume extrapolation**
   Always inspect both masks.  They diagnose different source limitations.

**Using a denser target grid to repair a sparse QHA source**
   Target-grid refinement interpolates the same information.  It cannot add
   missing thermodynamic curvature.

**Interpreting :math:`C^S` as full phonon-QHA elasticity**
   The adiabatic correction converts the QSA isothermal tensor; it does not add
   explicit phonon-strain elastic terms.

**Treating narrow propagated bands as total material uncertainty**
   Current uncertainties describe fitted numerical inputs under first-order
   assumptions, not all systematic errors.

**Smoothing away a stability loss or profile discontinuity**
   Quantas stores physical masks and layer transitions.  Plot smoothing should
   not be used to conceal them.

Results and interoperability
----------------------------

A point analysis can write a shared text input containing:

- the complete symmetric stiffness matrix;
- the QHA-consistent density;
- the requested P and T in the title;
- the selected isothermal or adiabatic condition.

That file can be passed directly to:

.. code-block:: console

   quantas elasticity run STATE.dat
   quantas seismic run STATE.dat

The public Python API provides the same transformation without a temporary
file:

.. code-block:: python

   from quantas.api import interop, thermoelasticity

   fit = thermoelasticity.read_result("dolomite_fit.hdf5")
   seismic_input = interop.thermoelastic_to_seismic(
       fit,
       pressure=5.0,
       temperature=800.0,
       tensor_condition="adiabatic",
       extrapolation_policy="fail",
   )

The exact CLI syntax is documented in :doc:`../cli/thermoelasticity`.  A
complete guided calculation is provided in
:doc:`../tutorials/thermoelasticity`, and the native schemas are described in
:doc:`../formats/thermoelastic_input` and
:doc:`../formats/thermoelastic_hdf5`.
