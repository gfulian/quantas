Quasi-Harmonic Approximation: implementation and workflow
=========================================================

Purpose and scope
-----------------

The QHA workflow turns a multi-volume phonon dataset into equilibrium and
thermodynamic properties on a pressure-temperature grid.  The physical model is
introduced in :doc:`../theory/qha`; here the focus is practical: how the
free-energy surface is represented, when the available numerical routes differ,
and which diagnostics deserve attention.

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

The CLI activates the embedded series explicitly with
``quantas qha run INPUT --kieffer``.  Omitting the flag ignores the optional
block and runs the normal phonon-only model.  This opt-in boundary separates
input enrichment from scientific model selection and permits direct paired
calculations from one immutable YAML file.

Raw elastic tensors and hydrostatic pre-stress
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Elastic constants obtained from a CRYSTAL energy--strain calculation without
an explicit finite-pressure correction cannot be passed directly to the
Christoffel solver.  Once a hydrostatic pressure has been assigned to every
state, the CRYSTAL interface converts the raw coefficients with the same
finite-pressure relation used by CRYSTAL ``PRESSURE``/``PRESSEOS``, as
documented for CRYSTAL finite-pressure elasticity
[#erba_mahmoud_belmonte_dovesi_2014]_:

.. math::

   B_{ijkl}(V_i)=C^{\mathrm{raw}}_{ijkl}(V_i)
   +\frac{P_i}{2}\left(2\delta_{ij}\delta_{kl}
   -\delta_{il}\delta_{jk}-\delta_{ik}\delta_{jl}\right).

VASP follows a separate ingestion rule. ``TOTAL ELASTIC MODULI`` are retained
as ``raw_stress_strain`` coefficients together with the unstrained reference
stress. After the selected hydrostatic pressure has been attached, the VASP
interface applies its own pressure adjustment and only the resulting
``wallace_hydrostatic`` state is admitted to Christoffel/Kieffer acoustics.
Quantas does not reuse the CRYSTAL/Erba transformation for VASP.

Pressure is positive in compression.  This CRYSTAL adapter rule is distinct
from the finite-strain ``wallace_delta`` term used by the QSA equations.  The
source of every :math:`P_i` is part of the data contract: it may come from the
output stress, a manually supplied value, an integrated energy EOS, or a
polynomial derivative of the QHA input's static :math:`E(V)` series. For the
latter two routes, the importer first treats the tensors as raw, matches elastic
and phonon volumes explicitly, evaluates :math:`P(V)=-dE/dV`, and only then
applies the selected backend correction. The generated input records the selected EOS
tag or polynomial degree, fit diagnostics, units,
evaluated pressures, and volume associations. The correction produces a new
elastic state and records the source and target tensor kinds, method, pressure
source, and software applying it. An already incremental tensor is rejected,
which prevents accidental double correction.

The backend-neutral pressure-assignment service is
``assign_hydrostatic_pressures()``.  Quantas also exposes the internal
Eulerian operators ``eulerian_hydrostatic_incremental_stiffness()``,
``convert_eulerian_hydrostatic_elastic_state()``, and
``convert_eulerian_hydrostatic_elastic_series()`` for tensors whose derivative
definition is explicitly compatible with that convention.  These operators
are **not** used as a substitute for external-code adapters: raw CRYSTAL tensors
are converted by ``crystal_hydrostatic_stiffness()`` in
:mod:`quantas.interfaces.crystal`, while raw VASP stress--strain tensors are
converted by ``vasp_hydrostatic_incremental_stiffness()`` in
:mod:`quantas.interfaces.vasp`.  Historical ``hydrostatic_wallace_*`` and
``correct_hydrostatic_*`` names remain compatibility aliases only.  A
correctly converted series can be passed directly to
``build_kieffer_volume_series()``.

Both QHA schemes include the acoustic contribution consistently.  With
``scheme=td``, the combined harmonic-plus-acoustic properties are fitted and
interpolated through the thermodynamic QHA path.  With ``scheme=freq``, Quantas
fits each of the three cutoff frequencies against volume using
``frequency_degree``.  Those fitted cutoffs are evaluated both during local
free-energy minimization and at the final equilibrium volumes before the
thermodynamic properties are recalculated.

Kieffer-enriched frequency QHA does not currently support the
``mode_gruneisen`` thermal-expansion route or the optional mode-Gruneisen
analysis.  A phonon-only weighted average would omit the acoustic branches and
is therefore rejected.  At the CLI, ``--kieffer`` automatically
turns off the otherwise enabled-by-default mode-Gruneisen output and records
that resolution in the options; an explicit ``--mode-gruneisen`` request is an
error.  The default ``mixed_derivative`` route and the numerical
thermal-expansion route remain available.

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

For every requested temperature, the workflow first computes HA properties at
all sampled volumes.  The resulting total Helmholtz free energy

.. math::

   F(V_i,T)=U_0(V_i)+F_{\mathrm{vib}}(V_i,T)

is the common starting surface for both QHA schemes.

If harmonic thermodynamics cannot be evaluated, the workflow can continue with the
static energy repeated at every temperature.  This behavior preserves a useful
static pressure-volume minimization and records a warning, but it is a
**degraded static-only calculation**, not a complete QHA result.  Temperature-
dependent properties from such a run must not be interpreted as physical QHA
predictions.

Frequency and thermodynamic schemes
-----------------------------------

The two QHA schemes differ in **what is interpolated across volume**. The
choice should follow the scientific question and the quality of branch
continuity rather than a preference for one numerical route.

.. list-table:: Frequency and thermodynamic QHA schemes
   :header-rows: 1
   :widths: 22 39 39

   * - Aspect
     - ``freq``
     - ``td``
   * - Interpolated object
     - Every :math:`\nu_{qj}(V)` branch
     - Integrated HA properties after mode summation
   * - Mode continuity
     - Required
     - Not required
   * - Main strength
     - Preserves mode-resolved information and permits mode Gruneisen analysis
     - Robust to harmless mode permutations and cheaper for bulk thermodynamics
   * - Main sensitivity
     - Branch identity, polynomial frequency fits, and extrapolated frequencies
     - Independent property-volume fits and their derivatives
   * - Prefer when
     - Continuous branches and mode-resolved interpretation are part of the study
     - Integrated properties are the goal or reliable branch tracking is absent

For ``freq``, each fitted branch is evaluated at the equilibrium volume and the
harmonic thermodynamics are recalculated from that interpolated spectrum. With
Kieffer enrichment, the three acoustic cutoffs are fitted in the same way and
must remain finite and positive over every volume reached by minimization.

With ``td``, the harmonic sum is first evaluated at each sampled volume and
the resulting thermodynamic quantities are then interpolated. This sacrifices
mode-resolved interpretation but avoids making the QHA result depend on a
branch assignment that the input cannot justify.

Agreement between ``freq`` and ``td`` is a useful robustness check when both
are scientifically meaningful. The worked comparison belongs to
:doc:`../tutorials/qha`, not to this implementation chapter.

Representing the free-energy curve
----------------------------------

At each temperature, one free-energy model is fitted over the sampled
volumes.  The fitted model is reused for all requested pressures at that
temperature.  Increasing the number of pressure points therefore does not
require refitting :math:`F(V,T)`, although later property evaluation still
scales with the number of pressure states.

Polynomial and EOS minimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At each temperature, the fitted free-energy representation is reused for all
requested pressures. The two available representations answer the same
minimization problem but impose different structure on :math:`F(V,T)`.

.. list-table:: Free-energy minimization choices
   :header-rows: 1
   :widths: 24 38 38

   * - Choice
     - Polynomial
     - Integrated EOS
   * - Representation
     - Centered and scaled local polynomial
     - Selected physical Energy EOS
   * - Main strength
     - Flexible description of a well-bracketed local basin
     - Structured compression curve with interpretable parameters
   * - Main risk
     - Edge oscillations, extra extrema, and sensitive high derivatives
     - Model/order correlation and apparently good fits outside useful support
   * - Prefer when
     - Sampling is dense and centered around the relevant minimum
     - A broader compression interval supports a physical EOS description

For the polynomial route, the minimum follows from

.. math::

   \frac{\partial F}{\partial V}+P=0

and retains a positive-curvature stationary point connected to the sampled
free-energy basin. Polynomial coefficients are numerical interpolation
parameters rather than material constants.

For the EOS route, an integrated Energy EOS is fitted at every temperature
and evaluates it at the requested pressures. Covariance-based uncertainty
propagation is available when the fitted covariance is usable. The default
model is ``BM3``; alternative Energy EOS families are described in
:doc:`../theory/eos`.

A disagreement between polynomial and EOS minimization is diagnostic evidence.
Changing the fit form cannot replace missing volume support; inspect the sampled
range, residuals, and location of the equilibrium state before selecting a
production model.

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

After polynomial minimization, the thermoelastic reconstruction needs the second and third volume
derivatives of free energy to obtain :math:`K_T` and :math:`K'_T`. Two routes
are available:

.. list-table:: Polynomial derivative methods
   :header-rows: 1
   :widths: 22 39 39

   * - Method
     - ``local_grid`` (default)
     - ``analytic``
   * - Evaluation
     - Reconstruct and refit a small volume neighborhood around each equilibrium state
     - Differentiate the global free-energy polynomial directly
   * - Strength
     - Probes curvature close to the requested state and provides an independent sensitivity test
     - Fast and exactly consistent with the selected global polynomial
   * - Sensitivity
     - Local spacing, local degree, and numerical noise
     - High-order behavior of the global polynomial, especially near edges

The default local grid contains five points separated by 0.05% of the central
volume. With ``scheme=freq`` the local free energies are regenerated from the
fitted spectrum; with ``scheme=td`` they are evaluated from the global property
fits.

Convergence should be checked at representative states by changing the number
and spacing of local points and comparing :math:`K_T` and :math:`K'_T` with the
analytic route. The equilibrium volume itself should not change when only the
derivative method changes.

Three routes to volumetric thermal expansion
--------------------------------------------

The available estimates are stored separately, together with the method that
supplies the authoritative :math:`\alpha_V` at each state.

.. list-table:: Thermal-expansion routes
   :header-rows: 1
   :widths: 25 36 39

   * - Method
     - Definition and availability
     - Main use and sensitivity
   * - ``mixed_derivative`` (default)
     - :math:`\alpha_V=K_T^{-1}(\partial S/\partial V)_T`; available for both schemes
     - Smooth thermodynamic route; depends on the entropy-volume derivative and :math:`K_T`
   * - ``mode_gruneisen``
     - Heat-capacity-weighted :math:`\gamma_{qj}`; available only for ``freq`` without Kieffer enrichment
     - Provides mode interpretation; sensitive to every frequency-volume derivative
   * - numerical volume derivative
     - :math:`V^{-1}(\partial V/\partial T)_P` from the final equilibrium-volume columns
     - Independent cross-check; sensitive to temperature spacing and endpoints

If the selected mixed or mode-Gruneisen value is unresolved at one state, the
numerical volume derivative is used there when it is available and the source
code records the fallback. A disagreement among methods is a reason to inspect
volume support, interpolation degree, branch continuity, and derivative
resolution rather than to select an automatic winner.

Derived thermodynamic quantities
--------------------------------

After selecting :math:`\alpha_V`, the workflow calculates

.. math::

   C_P-C_V=\alpha_V^2K_TVT,

.. math::

   C_P=C_V+(C_P-C_V),

and

.. math::

   K_S=K_T\frac{C_P}{C_V}.

At :math:`T=0`, and wherever :math:`C_V` is too small for a stable ratio,
The result then uses :math:`K_S=K_T`.

The optional macroscopic Grüneisen parameter is

.. math::

   \gamma=\frac{\alpha_VK_TV}{C_V}.

The default low-heat-capacity threshold is 1% of the Dulong–Petit value
:math:`3NR`.  Below this threshold the ratio is reported as unresolved rather
than allowing a numerically unstable value.  At exactly zero temperature the
stored macroscopic value is set to zero by convention.

Structural properties
---------------------

When the input contains a volume-constrained structural series, the workflow builds
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

Performance notes
-----------------

The dominant costs are the sampled harmonic calculation, mode fitting in
``freq``, one free-energy fit per temperature, and reconstruction of the final
P--T states. Adding pressure points is usually cheaper than adding temperature
points because all pressures at one temperature share the same fitted
free-energy representation.

For exploratory work, use ``qha inspect`` first, test a small P--T grid, disable
mode-Gruneisen analysis when it is not required, and prefer ``td`` when the
scientific objective is limited to integrated thermodynamics. Expand the
production domain only after the volume support and derivative choices have
been checked. Do not reduce the sampled volume set merely to accelerate the
calculation unless the effect on minima and derivatives has been quantified.

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

References
----------

.. include:: ../_generated/references/workflows_qha.inc
