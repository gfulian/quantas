Quasi-Harmonic Approximation: implementation and workflow
=========================================================

Purpose and scope
-----------------

The QHA workflow converts a multi-volume phonon dataset into equilibrium and
thermodynamic properties on a pressure-temperature grid.  The physical model is
described in :doc:`../theory/qha`; this page explains how Quantas represents the
free-energy surface, how the available choices differ, and how to assess their
numerical reliability.

QHA contains several legitimate routes because no single representation is
optimal for every dataset.  The main decisions are:

- interpolate **mode frequencies** or integrated **thermodynamic properties**;
- minimize a **polynomial** free-energy representation or an integrated
  **equation of state**;
- derive :math:`\alpha_V` from a mixed derivative, from mode Grüneisen
  parameters, or from numerical differentiation of :math:`V(P,T)`;
- evaluate polynomial thermoelastic derivatives analytically or from a local
  volume grid.

A production calculation should treat these alternatives as sensitivity tests,
not merely as interchangeable command-line settings.

Complete computational pipeline
-------------------------------

.. code-block:: text

   multi-volume phonon YAML
       │
       ├─ validate V, U0, frequencies, weights, structure, and mode continuity
       ├─ inspect static E(V) and estimate the sampled pressure interval
       ├─ calculate HA properties on the sampled T × V grid
       ├─ fit one F(V,T) representation for every temperature
       ├─ evaluate the minimum of F(V,T) + P V for every requested pressure
       ├─ reconstruct thermodynamic properties at V(P,T)
       ├─ select alpha_V and calculate Cp, Ks, and macroscopic gamma
       ├─ reconstruct cell parameters and anisotropic expansion when available
       └─ persist arrays, diagnostics, warnings, and provenance in HDF5

The pressure-temperature result arrays use temperature as the first axis and
pressure as the second axis.

Preflight inspection
--------------------

Before a full run, use

.. code-block:: console

   quantas qha inspect input.yaml --eos BM3

The inspector fits the **static** energy-volume data using both a polynomial
and, when requested, an integrated energy EOS.  It reports the pressure
associated with every sampled volume and the pressure range implied by each
representation.

This preview is valuable because QHA reliability is controlled first by the
volume support of the input.  A dense pressure grid does not compensate for an
equilibrium volume outside the sampled interval.

Inspect at least:

- whether polynomial and EOS pressure estimates have the same trend;
- whether either fit is classified as poor;
- whether the intended pressure range is bracketed by the sampled volumes;
- whether the reference volume lies inside a well-resolved energy basin;
- whether the outermost volumes are sufficiently far from the expected minima.

Input requirements and mode continuity
---------------------------------------

QHA uses the same normalized phonon arrays as HA, but requires a genuine
multi-volume series for interpolation and minimization.  The input may also
store a volume-constrained structural path used to reconstruct equilibrium
lattice parameters and anisotropic expansion.

The scientific construction of this normalized dataset is described in
:doc:`phonon_input_generation`.  In particular, independent CRYSTAL phonon
outputs can be checked by a backend-neutral eigenvector tracker before their
frequency arrays are assembled.

The frequency scheme depends on the public ``mode_continuity`` status:

``verified``
   Branch correspondence was established by a documented procedure.  Inspect
   ``mode_continuity_metadata`` to distinguish Quantas eigenvector tracking from
   source-managed continuity such as a native CRYSTAL QHA calculation.

``assumed``
   The stored array order is accepted as continuous, but no explicit
   verification is recorded.  Quantas permits frequency QHA and emits a
   warning.

``unknown``
   Continuity was not established, commonly because the required eigenvectors
   were unavailable.  Frequency QHA is rejected.

``unreliable``
   At least one mode assignment or degenerate subspace remains unresolved.
   Frequency QHA is rejected.

The thermodynamic ``td`` scheme does not require mode-by-mode continuity
because it interpolates harmonic quantities after the mode summation.

Kieffer acoustic enrichment
---------------------------

The Python API can add a volume-resolved ``KiefferVolumeSeries`` to primitive,
Gamma-only QHA data.  Every sampled QHA volume must match exactly one direct
cutoff state under the documented volume-tolerance policy.  File order is not
used as a scientific association: Quantas records and applies an explicit
one-to-one volume mapping.

The three Kieffer branches are added to the calculated Gamma phonons.  No
Gamma mode is removed or replaced.  The acoustic contribution is retained
separately on the sampled temperature-volume grid and persisted in the native
HDF5 result together with cutoff frequencies, effective velocities, and
matching diagnostics.

This integration currently applies to ``scheme=td``.  In that scheme the
combined harmonic-plus-acoustic properties are fitted and interpolated through
the existing thermodynamic QHA path, so the Kieffer contribution enters both
the free-energy minimization and the final pressure-temperature properties.
``scheme=freq`` is rejected explicitly: supporting it requires a documented
cutoff-volume evaluator used consistently by local minimization and by every
equilibrium thermodynamic property.

.. warning::

   Do not change ``unknown`` or ``unreliable`` to ``assumed`` merely to make a
   ``scheme=freq`` run start.  That edits the scientific claim made by the input
   without adding evidence.  Either establish continuity from the source data,
   extend or improve the phonon sampling, or use ``scheme=td`` when
   mode-resolved information is not required.

.. important::

   ``mode_continuity: verified`` addresses branch correspondence only.  It does
   not establish dynamical stability, force-constant convergence, adequacy of
   the volume range, or validity of the quasi-harmonic approximation.

The sampled harmonic stage
--------------------------

For every requested temperature, Quantas first computes HA properties at all
sampled volumes.  The resulting total Helmholtz free energy

.. math::

   F(V_i,T)=U_0(V_i)+F_{\mathrm{vib}}(V_i,T)

is the common starting surface for both QHA schemes.

If harmonic thermodynamics cannot be evaluated, Quantas can continue with the
static energy repeated at every temperature.  This behavior preserves a useful
static pressure-volume minimization and records a warning, but it is a
**degraded static-only calculation**, not a complete QHA result.  Temperature-
dependent properties from such a run must not be interpreted as physical QHA
predictions.

Frequency and thermodynamic schemes
-----------------------------------

``freq``: mode-resolved frequency interpolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For every q-point and mode, Quantas fits

.. math::

   \nu_{qj}(V)

with a polynomial of degree ``frequency_degree``.  At each equilibrium volume,
the frequencies are re-evaluated and all harmonic thermodynamic functions are
recalculated from the oscillator expressions.

Advantages
^^^^^^^^^^

- preserves mode-resolved information;
- enables mode Grüneisen parameters;
- allows thermodynamic properties to be recalculated from one internally
  consistent interpolated spectrum;
- for local polynomial derivatives, permits free energies to be regenerated
  around the equilibrium volume rather than read only from a global fit.

Limitations
^^^^^^^^^^^

- requires reliable branch correspondence across volumes;
- a mode crossing or accidental reordering can corrupt an otherwise smooth
  polynomial fit;
- one fit is required for every q-point and mode;
- fitted frequencies may become non-positive in extrapolation;
- mode-resolved fitting and Grüneisen analysis increase time and memory use.

Use ``freq`` when the phonon branches are continuous and mode information is
scientifically important.

``td``: thermodynamic-property interpolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Quantas first calculates integrated HA properties at every sampled volume and
then fits each property as a function of volume independently at every
temperature.  The fitted property is evaluated at the QHA equilibrium volume.

Advantages
^^^^^^^^^^

- does not require mode-by-mode branch continuity;
- is insensitive to harmless permutations of phonon modes;
- usually provides a robust route to integrated thermodynamic quantities;
- avoids storing and differentiating mode-resolved fits.

Limitations
^^^^^^^^^^^

- loses mode-resolved interpretation;
- cannot provide mode Grüneisen parameters;
- different thermodynamic properties are fitted separately, so consistency
  depends on the quality of all property fits;
- high-order volume derivatives remain sensitive to interpolation choices.

Use ``td`` when branch tracking is unavailable or the scientific objective is
limited to integrated thermodynamic properties.

Comparison
~~~~~~~~~~

.. list-table:: Frequency and thermodynamic QHA schemes
   :header-rows: 1
   :widths: 28 36 36

   * - Aspect
     - ``freq``
     - ``td``
   * - Interpolated object
     - Every mode frequency
     - Integrated HA properties
   * - Mode continuity
     - Required
     - Not required
   * - Mode Grüneisen parameters
     - Available
     - Unavailable
   * - Main sensitivity
     - Branch ordering and frequency fits
     - Property-volume fits
   * - Typical computational cost
     - Higher
     - Lower
   * - Preferred use
     - Mode-resolved analysis
     - Robust integrated thermodynamics

Agreement between the two schemes for one material is evidence of robustness,
not a universal guarantee.  The method-comparison exercise in
:doc:`../tutorials/qha` demonstrates the recommended check.

Representing the free-energy curve
----------------------------------

At each temperature Quantas fits one free-energy model over the sampled
volumes.  The fitted model is reused for all requested pressures at that
temperature.  Increasing the number of pressure points therefore does not
require refitting :math:`F(V,T)`, although later property evaluation still
scales with the number of pressure states.

Polynomial minimization
~~~~~~~~~~~~~~~~~~~~~~~

The polynomial path fits :math:`F(V,T)` in a centered and scaled volume
coordinate.  Scaling improves numerical conditioning without changing the
physical volume derivatives.

At pressure :math:`P`, Quantas finds positive-curvature roots of

.. math::

   \frac{\partial F}{\partial V}+P=0.

If several local minima exist, the implementation selects the one connected to
the sampled free-energy basin: it first identifies the sampled volume with the
lowest :math:`F+PV`, then chooses the closest stationary minimum.  Objective
values break an exact distance tie.

Advantages
^^^^^^^^^^

- flexible local representation near the sampled basin;
- no commitment to a specific EOS family;
- efficient when the relevant minimum is well bracketed;
- convenient for comparing alternative interpolation degrees.

Limitations
^^^^^^^^^^^

- behavior outside the sampled interval is uncontrolled;
- high polynomial degree can produce extra extrema or edge oscillations;
- :math:`K'_T` depends on a third volume derivative and is more sensitive than
  the equilibrium volume;
- polynomial coefficients are numerical parameters, not material constants.

The default degree is three.  The fit requires more volume points than fitted
coefficients, and a scientifically useful fit normally needs additional
redundancy beyond the mathematical minimum.

EOS minimization
~~~~~~~~~~~~~~~~

The EOS path fits an integrated energy EOS to :math:`F(V,T)` at every
temperature and evaluates that fitted model at each pressure.

Advantages
^^^^^^^^^^

- imposes a physically structured compression curve;
- gives interpretable parameters such as a reference volume and bulk modulus;
- usually behaves more regularly over a wider compression interval;
- supports covariance-based uncertainty propagation when the fit covariance is
  usable.

Limitations
^^^^^^^^^^^

- results depend on EOS family and order;
- parameters may be strongly correlated for a narrow volume interval;
- a good residual norm does not guarantee reliable high-order derivatives;
- EOS extrapolation remains extrapolation and must be flagged as such.

The default EOS tag is ``BM3``.  Alternative energy EOS forms are described in
:doc:`../theory/eos` and can be explored through the CLI reference.

Choosing polynomial or EOS minimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Initial minimization choice
   :header-rows: 1
   :widths: 42 24 34

   * - Dataset or objective
     - Start with
     - Why
   * - Dense, symmetric sampling around the minimum
     - polynomial
     - Provides a flexible local representation.
   * - Broad compression interval
     - EOS
     - A structured compression law is usually more stable.
   * - Few sampled volumes
     - low-order model
     - Extra parameters are not independently resolved.
   * - Strong disagreement between polynomial and EOS preview
     - neither blindly
     - Reinspect input range, convergence, and outliers.
   * - Properties requested outside sampled volume support
     - extend the input
     - Changing the fit form cannot create missing physical information.

Polynomial degrees
------------------

The CLI exposes two principal degree controls.

``--energy-degree``
   Sets the degree used for static energy, free energy, and thermodynamic-
   property volume fits.  The default is three.

``--frequency-degree``
   Sets the degree of every mode-resolved frequency-volume fit in the ``freq``
   scheme.  The default is three.

The public API also has a structural-path degree, default three, for the
volume-dependent cell-shape reconstruction.

Increasing a degree may reduce residuals while worsening extrapolation and
higher derivatives.  Use the smallest degree that resolves reproducible
curvature.  Inspect residuals mode-by-mode or temperature-by-temperature, and
compare the final physical properties rather than selecting a model from
:math:`R^2` alone.

Polynomial thermoelastic derivatives
------------------------------------

After polynomial minimization, Quantas must calculate

.. math::

   K_T=V\frac{\partial^2F}{\partial V^2},

and

.. math::

   K'_T=-1-V
   \frac{\partial^3F/\partial V^3}{\partial^2F/\partial V^2}.

Two implementations are available.

``local_grid`` — default
~~~~~~~~~~~~~~~~~~~~~~~~

A small symmetric volume grid is built around every equilibrium volume.  The
default contains **5 points**, with adjacent points separated by **0.05%** of
the central volume.  The local free energies are fitted again to obtain the
bulk properties.

With ``scheme=freq``, local free energies are regenerated from the fitted
phonon frequencies and static-energy model.  With ``scheme=td``, the local
values come from the fitted global polynomial.

Advantages
^^^^^^^^^^

- evaluates curvature in the immediate neighborhood of the state;
- provides a useful cross-check against derivatives of the global fit;
- in the frequency scheme, uses the locally reconstructed oscillator free
  energy.

Limitations
^^^^^^^^^^^

- adds local evaluations at every pressure-temperature state;
- too small a separation can amplify floating-point or fit noise;
- too large a separation no longer represents a local derivative;
- :math:`K'_T` remains intrinsically sensitive to the local polynomial degree.

``analytic``
~~~~~~~~~~~~

The second and third derivatives are evaluated directly from the global
free-energy polynomial.

Advantages
^^^^^^^^^^

- fastest route;
- exact derivative of the selected global polynomial;
- no additional local grid settings.

Limitations
^^^^^^^^^^^

- inherits every high-order feature of the global polynomial;
- can make :math:`K'_T` especially sensitive to edge behavior or excessive
  degree.

Convergence test for local derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At a small set of representative states, compare for example:

.. list-table:: Suggested derivative convergence table
   :header-rows: 1
   :widths: 18 22 20 20 20

   * - Points
     - Separation (%)
     - :math:`V`
     - :math:`K_T`
     - :math:`K'_T`
   * - 3
     - 0.05
     - …
     - …
     - …
   * - 5
     - 0.05
     - …
     - …
     - …
   * - 7
     - 0.05
     - …
     - …
     - …
   * - 5
     - 0.025
     - …
     - …
     - …
   * - 5
     - 0.10
     - …
     - …
     - …
   * - analytic
     - —
     - …
     - …
     - …

The equilibrium volume should not change when only the derivative method is
changed.  Material changes in :math:`K_T` or :math:`K'_T` indicate that the
curvature is not robustly resolved.

.. note::

   There is no HA/QHA numerical control with a default value of ``512`` in the
   current source.  The relevant QHA local derivative defaults are 5 points and
   0.05% separation.  The main global resolutions are the user-defined
   temperature and pressure grids.  Do not tune an unrelated ``512`` value as
   though it controlled QHA convergence.

Three routes to volumetric thermal expansion
--------------------------------------------

Quantas stores the available estimates separately and selects one requested
method as the authoritative :math:`\alpha_V`.  If the requested mixed or
mode-Grüneisen value is unresolved at an individual state, the implementation
uses the numerical volume derivative at that state when available.  A source
code array records which method supplied every value.

Mixed derivative — default
~~~~~~~~~~~~~~~~~~~~~~~~~~

The implementation uses the Maxwell relation

.. math::

   \alpha_V
   =\frac{1}{K_T}\left(\frac{\partial S}{\partial V}\right)_T.

At every temperature, entropy on the sampled volume grid is fitted by a
polynomial.  Its analytical volume derivative is evaluated at every equilibrium
volume and divided by :math:`K_T`.

Why it is the default
^^^^^^^^^^^^^^^^^^^^^

- available for both ``freq`` and ``td``;
- based directly on the fitted free-energy thermodynamics;
- avoids differentiating the final :math:`V(T)` curve;
- provides a smooth pressure-dependent estimate when the entropy-volume surface
  is well resolved.

Limitations
^^^^^^^^^^^

- depends on the entropy-volume polynomial and its derivative;
- inherits uncertainty in :math:`K_T`;
- can become unreliable near volume boundaries or with too few volumes.

Mode-Grüneisen route
~~~~~~~~~~~~~~~~~~~~

For ``scheme=freq``, Quantas differentiates each fitted frequency polynomial:

.. math::

   \gamma_{qj}(V)
   =-\frac{V}{\nu_{qj}}
   \frac{\partial\nu_{qj}}{\partial V}.

The modes are weighted by their harmonic contribution to :math:`C_V`, and the
thermal expansion is obtained from

.. math::

   \alpha_V=\frac{\bar\gamma C_V}{K_TV}.

Advantages
^^^^^^^^^^

- supplies mode-resolved interpretation;
- identifies modes driving positive or negative expansion;
- provides an independent thermodynamic consistency check.

Limitations
^^^^^^^^^^^

- requires mode continuity and positive fitted frequencies;
- non-positive or failed modes are excluded by default and recorded;
- sensitive to derivatives of every frequency-volume fit;
- more expensive than an integrated-property route.

Numerical volume derivative
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The numerical route evaluates

.. math::

   \alpha_V=\frac{1}{V}
   \left(\frac{\partial V}{\partial T}\right)_P

from the final equilibrium-volume columns.

Advantages
^^^^^^^^^^

- intuitive and independent of entropy or mode derivatives;
- useful as a cross-check;
- remains available as the pointwise fallback for unresolved primary methods.

Limitations
^^^^^^^^^^^

- requires at least two temperatures;
- depends on temperature spacing and endpoint differentiation;
- can amplify small irregularities in the minimized volume curve;
- cannot supply a value for a single-temperature calculation.

Choosing an expansion method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Thermal-expansion decision guide
   :header-rows: 1
   :widths: 35 25 40

   * - Situation
     - Primary choice
     - Recommended check
   * - General QHA calculation
     - mixed derivative
     - Compare with numerical differentiation.
   * - Verified continuous modes
     - mode Grüneisen
     - Compare weighted and macroscopic :math:`\gamma`.
   * - No trustworthy mode ordering
     - mixed derivative
     - Use ``td`` and compare with numerical values.
   * - Coarse temperature grid
     - mixed derivative
     - Numerical :math:`\alpha_V` may be too step-dependent.
   * - Disagreement among methods
     - no automatic winner
     - Inspect volume support, fit degrees, mode continuity, and derivatives.

Derived thermodynamic quantities
--------------------------------

After selecting :math:`\alpha_V`, Quantas calculates

.. math::

   C_P-C_V=\alpha_V^2K_TVT,

.. math::

   C_P=C_V+(C_P-C_V),

and

.. math::

   K_S=K_T\frac{C_P}{C_V}.

At :math:`T=0`, and wherever :math:`C_V` is too small for a stable ratio,
Quantas sets :math:`K_S=K_T`.

The optional macroscopic Grüneisen parameter is

.. math::

   \gamma=\frac{\alpha_VK_TV}{C_V}.

The default low-heat-capacity threshold is 1% of the Dulong–Petit value
:math:`3NR`.  Below this threshold the ratio is reported as unresolved rather
than allowing a numerically unstable value.  At exactly zero temperature the
stored macroscopic value is set to zero by convention.

Structural properties
---------------------

When the input contains a volume-constrained structural series, Quantas builds
a one-dimensional structural path in volume and evaluates it at
:math:`V(P,T)`.  Cell parameters and the deformation gradient are therefore
reconstructed from static structures, while temperature enters through the QHA
equilibrium volume.

Axial expansion is evaluated by the chain rule, for example

.. math::

   \alpha_a=
   \frac{\partial\ln a}{\partial\ln V}\,\alpha_V.

The tensor trace reproduces the selected volumetric expansion by construction.
This is not a full anisotropic phonon QHA: it assumes that the cell-shape path
is governed by the sampled static volume-constrained structures.

Volume support, extrapolation, and fit quality
----------------------------------------------

Every equilibrium volume is classified as:

``inside``
   More than 5% of the sampled interval width from either boundary.

``near_boundary``
   Inside the interval but within 5% of one boundary.

``outside``
   Outside the sampled volume interval.

The default API policy is to warn on extrapolation.  A strict application can
request failure instead.  Near-boundary points deserve scrutiny because
curvature and structural derivatives are already weakly constrained there.

Polynomial fits are classified as poor when, among other diagnostics, the
design matrix is rank deficient, its condition number exceeds
:math:`10^{12}`, or :math:`R^2<0.95`.  EOS fits use the common Quantas fit
quality contract.  A converged fit with warnings is retained by the default
``warn`` quality policy; a stricter API workflow can stop on poor quality.

Failure policies and partial results
------------------------------------

A local state can fail because:

- the free-energy fit failed;
- no physical minimum was found;
- thermoelastic derivatives were invalid;
- the selected EOS could not be evaluated at the requested pressure;
- an equilibrium volume violated a strict extrapolation policy.

The CLI exposes three fit-failure policies:

``continue``
   Record the failed state and proceed through the grid.

``stop`` — default
   Stop after the configured number of consecutive failures, default five.
   Previously completed states remain in the result.

``raise``
   Raise immediately for applications that require all states to succeed.

A partially completed HDF5 result is useful diagnostic evidence, but failed or
masked points must not be silently interpreted as a complete thermodynamic
surface.

Uncertainty behavior
--------------------

The default public options request covariance propagation where supported.  In
the EOS minimization path, the covariance of the temperature-specific EOS fit
can be propagated to :math:`V`, :math:`K_T`, and :math:`K'_T`.

Advanced API users may request Monte Carlo propagation and control sample count,
seed, confidence level, and the minimum accepted fraction of physical samples.
Bootstrap propagation is not currently available for EOS pressure states; the
workflow records a warning and omits those uncertainties rather than silently
substituting another method.

Polynomial minimization does not provide an equivalent full parameter-
covariance propagation for all final state quantities.  Do not compare missing
polynomial uncertainties with EOS uncertainties as though they represented the
same statistical model.

Performance and practical acceleration
--------------------------------------

The total cost has several components.

Sampled harmonic thermodynamics
   Scales approximately as
   :math:`N_TN_qN_{\mathrm{mode}}N_V`.

Frequency fitting
   In ``freq``, one volume polynomial is fitted for every q-point and mode.

Free-energy fitting
   One polynomial or EOS is fitted for every temperature, not every pressure.

State evaluation
   Scales with :math:`N_TN_P`; frequency thermodynamics and mode Grüneisen
   analysis add mode-resolved work at these states.

Structural reconstruction
   Scales with the final pressure-temperature grid and is usually secondary.

A disciplined acceleration strategy is:

1. run ``qha inspect`` before any large grid;
2. test one or a few temperatures and pressures;
3. disable mode Grüneisen analysis when it is not required;
4. use ``td`` when mode-resolved information is unnecessary or unreliable;
5. compare polynomial and EOS minimization on a reduced grid;
6. validate local derivative settings at representative states;
7. only then expand the final pressure-temperature domain.

Adding pressure points is often cheaper than adding temperature points because
the fitted free-energy representation is shared across all pressures at one
temperature.  Nevertheless, every final state still requires property
reconstruction, so very dense pressure grids are not free.

The following changes accelerate a run without changing the source phonon
dataset:

- coarser preliminary temperature and pressure steps;
- ``--no-mode-gruneisen`` when mode analysis is not needed;
- ``--no-gruneisen`` when the macroscopic ratio is not needed;
- analytic polynomial derivatives for a controlled sensitivity test;
- the ``td`` scheme for integrated properties.

Do not accelerate a production calculation by deleting volumes until the
impact on fitted minima and derivatives has been quantified.

Recommended staged workflow
---------------------------

.. list-table:: Recommended QHA workflow
   :header-rows: 1
   :widths: 12 34 54

   * - Stage
     - Action
     - Acceptance question
   * - 1
     - Validate the normalized input and mode-continuity status.
     - Are arrays, units, q-point weights, and structures consistent?
   * - 2
     - Run ``qha inspect`` with polynomial and EOS previews.
     - Is the intended pressure range bracketed by sampled volumes?
   * - 3
     - Run a small ``T × P`` grid with default settings.
     - Are all minima physical and inside the volume interval?
   * - 4
     - Compare ``freq`` and ``td`` where both are meaningful.
     - Are integrated properties insensitive to branch representation?
   * - 5
     - Compare polynomial and EOS minimization.
     - Are :math:`V`, :math:`K_T`, and :math:`\alpha_V` robust?
   * - 6
     - Compare thermal-expansion methods.
     - Is disagreement traceable to fits, modes, or grid resolution?
   * - 7
     - Test derivative and grid convergence.
     - Are high-order properties stable within the scientific tolerance?
   * - 8
     - Run the production domain and validate the result.
     - Are warnings, masks, uncertainty coverage, and provenance acceptable?

Programmatic validation
-----------------------

The public QHA API includes result-validation and comparison helpers.  They can
be used to compare two calculations on the same grid, for example ``freq``
against ``td`` or polynomial against BM3 minimization, without writing custom
array-alignment code.

A method comparison should focus on physical quantities, not only fit metrics:

- equilibrium volume;
- :math:`K_T` and :math:`K'_T`;
- :math:`\alpha_V` and its source code;
- :math:`C_P-C_V`;
- :math:`K_S`;
- macroscopic and mode-weighted Grüneisen parameters;
- structural extrapolation masks.

Related documentation
---------------------

- Input generation and mode continuity: :doc:`phonon_input_generation`
- Scientific theory: :doc:`../theory/qha`
- Complete worked example and method exercise: :doc:`../tutorials/qha`
- Input specification: :doc:`../formats/phonon_yaml`
- CLI syntax: :doc:`../cli/qha`
- Public Python surface: :doc:`../api/qha`
- EOS background: :doc:`../theory/eos`
