HA and QHA validation
=====================

Scope
-----

HA and QHA validation is divided into two complementary questions:

#. do the normalized phonon inputs preserve the quantities and branch
   relationships required by the scientific model?;
#. once the input is accepted, do the harmonic and quasi-harmonic numerical
   workflows preserve the expected formulas, limits, shapes, and cross-method
   consistency?

This page concentrates on the first question where the recent CRYSTAL interface
and phonon-mode continuity work required new scientific validation.  The HA/QHA
thermodynamic equations and workflow choices are documented in
:doc:`../theory/ha`, :doc:`../theory/qha`, :doc:`../workflows/ha`, and
:doc:`../workflows/qha`.

The authoritative implementation of input generation is described in
:doc:`../workflows/phonon_input_generation`.

Kieffer acoustic thermodynamics
-------------------------------

The Kieffer sine-wave model is validated as a statistical-thermodynamics core
and is connected to the single-volume HA Python API and both multi-volume QHA
schemes.  End-to-end CLI tests now verify that ``run --kieffer`` reads the
embedded YAML series, evaluates the acoustic contribution, and preserves it in
native HA and QHA HDF5 results.

The validation uses ordinary cutoff frequencies in hertz and the nonsingular
integration variable

.. math::

   \theta=\arcsin(\nu/\nu_{\max}), \qquad 0\leq\theta\leq\frac{\pi}{2}.

It verifies the following properties:

* the historical Helmholtz-above-0-K and heat-capacity results are reproduced;
* the zero-point energy is checked against the analytical mean frequency

  .. math::

     \langle\nu\rangle = \nu_{\max}\frac{24(\pi-2)}{\pi^3};

* :math:`S=-\partial F/\partial T`, :math:`C_V=T\partial S/\partial T`, and
  :math:`U_{\mathrm{th}}=F_{\mathrm{th}}+TS` hold numerically;
* entropy and all thermal contributions vanish at zero temperature;
* three acoustic branches approach :math:`3R` in heat capacity at high
  temperature;
* multi-volume inputs preserve ``float64`` values and the ``(T,V)`` result
  shape.

Single-volume HA composition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first workflow integration accepts an explicit ``KiefferVolumeSeries``
through ``quantas.api.ha.run(..., kieffer_cutoffs=...)``.  It requires exactly
one direct cutoff state and validates all of the following before evaluating
the acoustic thermodynamics:

* one sampled volume and one q-point;
* explicit Gamma coordinates, modulo reciprocal-lattice vectors;
* the identity phonon-supercell matrix;
* primitive, single-repetition structural normalization when structural
  metadata are present;
* an explicit unique match between the HA and cutoff primitive-cell volumes.

The three sine-wave branches are **additional** acoustic contributions.  No
calculated Gamma frequency is selected, removed, or replaced.  Tests compare a
normal HA result with an enriched result and verify that their difference is
exactly the independently evaluated Kieffer contribution, even when the three
lowest Gamma frequencies are small and positive.

The total HA arrays contain the composed thermodynamic properties, while the
acoustic zero-point energy, thermal energy, entropy, heat capacity, and
Helmholtz energy are retained separately under ``kieffer_contribution``.  The
same separation is preserved in the native HDF5 payload together with cutoff
frequencies, effective velocities, composition policy, and volume-match
diagnostics.

Multi-volume thermodynamic QHA composition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The QHA validation uses input volumes deliberately ordered differently from
the increasing cutoff series.  It verifies that explicit one-to-one volume
matches, rather than array indices, determine the acoustic association.  The
enriched-minus-harmonic sampled Helmholtz surface must equal an independent
Kieffer evaluation at every temperature and volume, while the original Gamma
frequency array remains unchanged.

Negative tests cover missing or non-Gamma coordinates, non-identity phonon
supercells, and incomplete or mismatched cutoff volume sets.  Frequency-scheme
tests fit all three cutoffs against volume, compare their independently
evaluated thermodynamic contribution at arbitrary volumes, and verify that the
same acoustic surface affects both minimization and final equilibrium
properties.  Mode-Gruneisen analysis is rejected until the acoustic branches
can be included in its heat-capacity-weighted average.  A public API and HDF5
round-trip test confirms that the sampled acoustic component and its provenance
survive the complete QHA lifecycle.  CLI tests additionally verify opt-in
forwarding, automatic resolution of the frequency-scheme modal default,
rejection of an explicit incompatible request, Kieffer citation output, and
complete enriched-input-to-HDF5 execution for both HA and QHA.

Historical entropy defect
~~~~~~~~~~~~~~~~~~~~~~~~~

The original Quantas ``Kieffer.entropy`` routine squared the Bose occupation
denominator in the first entropy term.  The published equation contains

.. math::

   \frac{x}{e^x-1}

rather than :math:`x/(e^x-1)^2`.  The historical value is retained as a frozen
characterization datum, but the new core intentionally uses the published
formula.  The corrected result is independently constrained by
:math:`S=-\partial F/\partial T`; no compatibility switch preserves the defect.

Acoustic velocity averages
~~~~~~~~~~~~~~~~~~~~~~~~~~

The acoustic cutoff path reuses the SEISMIC Christoffel solver.  Phase
velocities are integrated with Gauss--Legendre quadrature in
:math:`\mu=\cos\theta` and a periodic uniform quadrature in :math:`\phi`:

.. math::

   u_i = \left[\frac{1}{4\pi}\int_{4\pi}v_i^{-3}\,d\Omega\right]^{-1/3}.

The modes follow the existing SEISMIC convention of local ascending phase
speed: slow quasi-shear, fast quasi-shear, and quasi-longitudinal.  Exact shear
degeneracy is therefore harmless for the integral because the two coincident
speeds make the local labelling immaterial.  Degenerate directions and clamped
eigenvalues remain explicit diagnostics.

An isotropic analytical test recovers both shear velocities and the
longitudinal velocity independently of quadrature order.  An anisotropic
hydroxylapatite test checks convergence under simultaneous refinement of both
angular orders.

For primitive-cell volume :math:`V` the cutoff validation uses

.. math::

   K_{\max}=\left(\frac{6\pi^2}{V}\right)^{1/3},\qquad
   \omega_{i,\max}=\frac{2}{\pi}u_iK_{\max},\qquad
   \nu_{i,\max}=\frac{\omega_{i,\max}}{2\pi}.

Tests independently verify the conversions from cubic angstrom to cubic metre,
from km/s to m/s, and from hertz to inverse centimetre.

Elastic-state and cutoff provenance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The shared model layer distinguishes raw energy--strain stiffness matrices
from hydrostatic Wallace and full-stress incremental tensors.  A correction
record contains the source tensor convention, pressure value and origin,
correction method, and the component that applied it.  Validation rejects a
correction whose source is already incremental, preventing silent double
application before the data reach Christoffel acoustics.

Pressure provenance distinguishes applied pre-stress, parsed output stress,
manual pressure, energy-EoS pressure, and energy-polynomial pressure.  Raw
tensors remain representable because a later enrichment stage may correct
them, but the acoustic eligibility check rejects them until that operation has
been completed explicitly.

Volume-resolved cutoff states retain the source elastic-state indices and mark
each value as direct or interpolated.  Interpolated series require a named
interpolation method.  Exact matching between independently printed QHA and
elastic volumes uses a stored relative and absolute tolerance and reports both
differences for every association; missing or ambiguous matches are errors.

Direct elastic-to-cutoff composition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An incremental elastic series can now be transformed directly into a Kieffer
cutoff series.  Every volume is processed independently through the shared
Christoffel solver, inverse-cube spherical average, equal-volume Brillouin
sphere, and sine-wave cutoff equations.  The result preserves the source-state
index, tensor convention, pressure origin, quadrature orders, refinement
change, degeneracy count, and clamped-eigenvalue count.

The workflow validates the complete tensor series before processing its first
volume.  A raw or unknown tensor therefore fails without producing a partial
cutoff series.  An isotropic two-volume test independently verifies the
expected density dependence :math:`u\propto\rho^{-1/2}` and the combined cutoff
scaling :math:`\nu_{\max}\propto uV^{-1/3}`.

Hydrostatic correction of raw elastic states
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The raw-to-incremental path is tested separately from parsing and acoustics.
Analytical component tests verify the normal, normal-coupling, and shear terms
of ``C_raw - P * Delta`` with pressure positive in compression.  State and
series tests require complete pressure provenance, preserve structural and
energy metadata, and confirm that applying the operation to an incremental
tensor is rejected.  The corrected synthetic series is then passed through the
complete Christoffel averaging and Kieffer cutoff construction to verify that
the result is immediately acoustic-ready.

The comparison against output frozen from the historical implementation uses
a cross-platform tolerance for adaptive quadrature.  This tolerance is not
used by the analytical zero-point, high-temperature, or thermodynamic-identity
tests.

OHAp multi-volume characterization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The complete hydroxylapatite (OHAp) characterization starts from ten
independent CRYSTAL Gamma-frequency calculations and ten raw ELAPIEZO elastic
calculations.  Every volume has both data types.  The normalized QHA input
contains 44 atoms, 132 Gamma modes, two formula units, and volumes from
482.2593 to 566.1216 angstrom cubed.  The first three frequencies are exactly
zero at every volume; they remain in the input, while the three Kieffer
branches are added as a separate continuous acoustic contribution.

The source archive used for this characterization has SHA-256
``b27ee27f1bb3832816edfe693513c16c278b5a9b0dab0c021b414bee024f8563``.
The aggregate digest over the two portable list files and their 20 referenced
outputs is
``2c1c1483ac67afe88b11eabce03e147366ae4d19f8cae9611c0fc14e154e428b``.
The compact numerical subset retained by the automated regression is
``tests/modules/qha/data/ohap_kieffer_reference.json``; it records both
digests so that results from a different archive are not silently compared.

Reproducible validation driver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The developer utility ``tools/validate_kieffer_ohap.py`` reproduces the full
chain without modifying the source outputs:

.. code-block:: console

   python tools/validate_kieffer_ohap.py /path/to/ohap_complete_kieffer \
      --source-archive /path/to/ohap_complete_kieffer.zip \
      --output-dir ohap_kieffer_validation

It regenerates the base QHA YAML, creates all pressure-source variants,
performs the acoustic quadrature study, runs phonon-only and Kieffer-enriched
QHA, round-trips the primary pair through HDF5, exercises expected failures,
and writes Markdown, JSON, and long-form CSV evidence.  Existing results are
protected unless ``--force`` is stated explicitly.

Pressure routes and acoustic convergence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``auto``, ``output-stress``, and manual pressures reproduce exactly the same
cutoffs; the result is also invariant to reversing the elastic-file list.
``auto`` chooses output stress here because all ten CRYSTAL outputs contain a
complete unstrained-stress pressure.  The two energy-derived alternatives
exercise the path needed when those values are absent:

.. list-table:: OHAp pressure-source characterization
   :header-rows: 1
   :widths: 30 22 22 26

   * - Pressure source
     - Pressure range (GPa)
     - Maximum difference from output stress (GPa)
     - Maximum cutoff difference from output stress (%)
   * - CRYSTAL output stress
     - 11.5300 to -7.4420
     - 0
     - 0
   * - Energy BM3
     - 11.0863 to -7.6862
     - 0.4437
     - 0.4745
   * - Degree-three energy polynomial
     - 10.9800 to -7.5801
     - 0.5500
     - 0.3867

The integrated BM3 fit is classified ``good`` with
:math:`R^2=0.9998476` and an energy RMSE of
:math:`2.0103\times10^{-4}` Ha.  The centred and scaled cubic polynomial is
also classified ``good`` with :math:`R^2=0.9998586` and an RMSE of
:math:`1.9364\times10^{-4}` Ha.  These numbers compare independent pressure
estimators; neither fitted curve is forced to reproduce the output stress.

The default 12 by 24 quadrature, refined to 24 by 48, changes by at most
:math:`1.0501\times10^{-3}` internally.  Comparing its final cutoffs with a
24 by 48 calculation refined to 48 by 96 gives a maximum relative change of
:math:`7.1339\times10^{-4}`.  A 48 by 96 calculation refined to 96 by 192
reduces both its reported refinement change and the dense-to-fine difference
to :math:`1.2787\times10^{-4}`.  This fine series is used for the reference
QHA comparison, rather than treating the default as converged by assumption.

At the static-energy minimum, :math:`V=527.75215929` angstrom cubed, the fine
effective velocities are 3.75295, 4.14359, and 7.57286 km/s.  The corresponding
ordinary cutoff frequencies are 1.83406, 2.02497, and 3.70084 THz.

Acoustic and optical temperature dependence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cleanest composition diagnostic evaluates both contributions at the same
sampled volume.  For each positive additive property the reported percentage
is

.. math::

   100\frac{X_{\mathrm{acoustic}}}
   {X_{\mathrm{optical}}+X_{\mathrm{acoustic}}}.

.. list-table:: Direct Kieffer share at the OHAp static minimum
   :header-rows: 1
   :widths: 16 28 28 28

   * - Temperature (K)
     - Thermal energy (%)
     - Entropy (%)
     - :math:`C_V` (%)
   * - 1
     - 100.0000
     - 100.0000
     - 100.0000
   * - 5
     - 85.4191
     - 87.7057
     - 68.2755
   * - 10
     - 30.0217
     - 32.0743
     - 25.7445
   * - 20
     - 29.9219
     - 29.8223
     - 31.8078
   * - 50
     - 23.2324
     - 24.6772
     - 16.9299
   * - 100
     - 12.0008
     - 14.4319
     - 7.4351
   * - 300
     - 5.0527
     - 6.9928
     - 3.3478
   * - 1000
     - 3.0269
     - 4.4616
     - 2.4228
   * - 1500
     - 2.7667
     - 4.0734
     - 2.3483

The expected broad trend is present: the acoustic fraction dominates at low
temperature and decreases as the 129 positive Gamma modes become populated.
It must not, however, be encoded as strict point-by-point monotonicity.  OHAp
has optical modes beginning near 39 inverse centimetres, and their activation
in the same low-temperature window produces a genuine 10--20 K crossover in
the :math:`C_V` fraction.  At the classical limit the acoustic fraction does
not vanish; it approaches the branch-count ratio

.. math::

   \frac{3}{129+3}\,100 = 2.272727\%.

The numerical result at :math:`10^7` K is 2.2727273%.  Equivalently, the
optical share grows from approximately zero at 1 K toward 97.727273%.

Net effect on equilibrium QHA properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All four supported combinations--frequency or thermodynamic-property scheme,
with polynomial or EOS minimization--complete from 0 to 1500 K at zero
pressure with finite properties, no failed states, and equilibrium volumes
inside the sampled interval.  The primary frequency-polynomial comparison is:

.. list-table:: Kieffer-enriched minus phonon-only QHA at zero pressure
   :header-rows: 1
   :widths: 13 18 18 17 17 17

   * - Temperature (K)
     - :math:`\Delta V` (angstrom cubed)
     - :math:`\Delta V/V` (%)
     - :math:`\Delta C_V/C_V` (%)
     - :math:`\Delta S/S` (%)
     - :math:`\Delta K_T/K_T` (%)
   * - 0
     - 0.0552866
     - 0.010470
     - --
     - --
     - -0.050268
   * - 300
     - 0.3415710
     - 0.064421
     - 3.483974
     - 7.608327
     - -0.483698
   * - 1000
     - 1.3872771
     - 0.256172
     - 2.493718
     - 4.876451
     - -1.910547
   * - 1500
     - 2.2929304
     - 0.416414
     - 2.412688
     - 4.539445
     - -3.130525

These percentages are **not** the direct acoustic fractions in the preceding
table.  They compare two independently minimized QHA states, so they also
contain the volume shift caused by the acoustic free energy.  Keeping the two
definitions separate avoids attributing an equilibrium-volume feedback to a
fixed-volume branch contribution.  The enriched :math:`C_V` moves closer to
the complete 132-branch Dulong--Petit limit: at 1500 K the phonon-only value is
94.5675% of that limit and the enriched value is 96.8491%.

Stress and failure characterization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A second run over -2 to 5 GPa in 1 GPa steps and 0 to 1500 K in 300 K steps
passes for all four QHA combinations.  Its equilibrium volumes remain between
506.7 and 564.2 angstrom cubed, within the source interval.  It can be
reproduced with:

.. code-block:: console

   python tools/validate_kieffer_ohap.py /path/to/ohap_complete_kieffer \
      --temperature 0 1500 300 --pressure -2 5 1 \
      --output-dir ohap_kieffer_pressure_validation

The same driver deliberately confirms rejection of a base input without an
embedded Kieffer block, a non-Gamma q-point, a phonon supercell, an incomplete
or mismatched cutoff series, a negative cutoff, an incompatible mode-Gruneisen
request, duplicate elastic sources, in-place enrichment, silent replacement
of an existing Kieffer block, and an incomplete elastic-volume series.  The
primary phonon-only and enriched HDF5 results reproduce every public property
exactly after a write/read round trip.

This characterization establishes numerical consistency, provenance, and
workflow behavior for a realistic material.  It is not an experimental
validation of the chosen electronic-structure method, the quasi-harmonic
approximation, or the hydrostatic treatment of a state with appreciable
deviatoric stress.

Validation hierarchy
--------------------

The phonon input path is tested at four levels:

``parser characterization``
   Real CRYSTAL outputs establish the syntax, block structure, complex phase
   reconstruction, atom ordering, and normalization actually produced by the
   external code.

``synthetic invariance tests``
   Controlled eigenvectors test properties that should be exact: permutation
   recovery, invariance to complex phase, basis rotation inside degenerate
   subspaces, and reference-independent local assignments.

``adversarial continuity tests``
   Deliberately weak or inconsistent mode links verify that the tracker does
   not silently promote every difficult assignment to ``verified``.

``real multi-volume regression``
   A seven-volume dolomite dataset exercises the complete parser and tracking
   path on a realistic dispersion calculation with many branches and q-points.

.. important::

   Validation is designed to preserve scientific failure modes.  A useful mode
   tracker must be able to return ``unreliable`` when the data do not justify a
   branch assignment.  A test suite in which every difficult case is forced to
   pass would validate convenience, not scientific reliability.

CRYSTAL phonon eigenvector parser
---------------------------------

Block structure and dimensions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CRYSTAL parser is characterized against the distributed dolomite phonon
output.  It reads

.. code-block:: text

   27 q-points
   30 modes per q-point
   10 atoms per eigenvector

and returns arrays with shapes

.. math::

   \nu:\;(27,30),

and

.. math::

   e:\;(27,30,10,3).

A separate synthetic block verifies that the final printed eigenvector block
may contain fewer than six modes.  This protects against an implementation
that assumes all CRYSTAL blocks have equal width.

Real and complex eigenvectors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At Gamma, parsed eigenvectors are real.  At non-zero q-points, the
characterization output contains non-zero imaginary components reconstructed
from CRYSTAL's in-phase and anti-phase sections when the q-point is complex;
real q-points may provide only the in-phase section.

The parser rejects an inconsistent real/imaginary decomposition if the two
sections use different frequencies or atom orderings.

Mass-weighted normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every parsed mode is required to satisfy

.. math::

   \sum_{\kappa\alpha}|e_{qj,\kappa\alpha}|^2 = 1

within numerical tolerance.  The real dolomite fixture verifies unit norms for
all :math:`27\times30` modes.

The native MgO QHA output provides an additional physical check.  The relative
translation amplitudes of Mg and O in a mass-weighted acoustic mode reproduce
the expected factor

.. math::

   \left|\frac{e_{\mathrm{Mg}}}{e_{\mathrm O}}\right|
   = \sqrt{\frac{m_{\mathrm{Mg}}}{m_{\mathrm O}}}.

This confirms that the transformation from CRYSTAL classical-amplitude
coordinates to the Quantas mass-weighted unit representation is not merely
norm preserving; it has the expected mass dependence.

Synthetic mode-tracking invariants
----------------------------------

Permutation recovery
~~~~~~~~~~~~~~~~~~~~

A deliberately permuted orthonormal basis is recovered exactly by the global
one-to-one overlap assignment.  This verifies that branch restoration is based
on the overlap matrix rather than on the raw printed mode indices.

Complex phase invariance
~~~~~~~~~~~~~~~~~~~~~~~~

If one eigenvector is multiplied by an arbitrary phase

.. math::

   e_j \rightarrow e^{i\phi_j}e_j,

its physical mode identity must not change.  Because Quantas uses

.. math::

   |\langle e_i|e_j\rangle|,

synthetic complex-phase rotations leave the selected overlaps equal to one.

Degenerate-subspace invariance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A two-dimensional degenerate eigenspace is rotated by an arbitrary unitary
basis transformation between adjacent states.  Individual vectors no longer
match one-to-one, but the subspace singular values remain unity and the dataset
is verified.

This test protects the physically essential distinction between an eigenvector
and the eigenspace of a degenerate eigenvalue.

.. warning::

   A branch tracker that validates degenerate modes through individual scalar
   products will produce basis-dependent results.  Such an implementation can
   fail even when the physical eigenspace is exactly unchanged.

Reference independence
~~~~~~~~~~~~~~~~~~~~~~

The same multi-state synthetic dataset is tracked with different
``reference_index`` values.  Every adjacent-volume local permutation and
selected overlap remains unchanged.  Only the final branch labels differ.

This explicitly characterizes the design rule

.. math::

   \text{local matching} \perp \text{reference choice}.

The tracker also accepts input files in non-monotonic volume order and confirms
that the actual comparison path follows sorted adjacent volumes.

Ambiguity and failure semantics
-------------------------------

Ambiguity is a caution
~~~~~~~~~~~~~~~~~~~~~~

A synthetic mixed pair with strong but insufficiently separated scalar
products is retained as ``verified`` with cautions.  This verifies that the
configured ambiguity margin is a diagnostic of non-uniqueness rather than an
automatic rejection criterion.

Weak overlap without independent evidence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With only two sampled states, a deliberately weak overlap cannot be validated
from an independent frequency path.  The tracker therefore returns
``unreliable``.

This test is important because any two points can be joined exactly by a line;
such a fit contains no independent continuity information.

Leave-one-out rescue of a smooth branch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A five-volume synthetic series contains a weak final eigenvector overlap but a
smooth frequency path.  The endpoint is removed from the fit, predicted from
the remaining states, and accepted only when the independent residual satisfies
the configured threshold.  The link remains a ``caution`` and the full dataset
is ``verified``.

Adversarial leave-one-out rejection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A companion test modifies the weak endpoint so that the global polynomial fit
is still classified as smooth by the descriptive thresholds, while the
leave-one-out prediction of the omitted point fails.

The expected result is

.. code-block:: text

   global branch diagnostic: supported
   leave-one-out endpoint:    not supported
   continuity status:         unreliable

.. important::

   This adversarial case demonstrates that the global :math:`\nu(V)` fit cannot
   rescue its own outlier.  The distinction is methodological, not cosmetic:
   low-overlap acceptance requires evidence from a model fitted **without** the
   point being tested.

Native CRYSTAL QHA characterization: MgO
----------------------------------------

The distributed CRYSTAL17 MgO QHA output contains eleven volume-dependent
Gamma-mode sets for a 64-atom phonon supercell.  The eigenvector parser returns

.. math::

   (11,192,64,3)

complex mode arrays and unit norms for every mode.

The native output is independently checked for CRYSTAL's explicit statement
that frequency continuity with volume was found.  The normalized YAML therefore
records

.. code-block:: yaml

   mode_continuity: verified
   mode_continuity_metadata:
     method: crystal-qha
     source: crystal

As a characterization cross-check, the parsed mode sets are also passed through
the Quantas overlap tracker.  The result contains no ambiguous or low-overlap
assignments.  In this highly degenerate Gamma dataset, individual minimum
overlap is not a meaningful scalar diagnostic; the minimum matched-subspace
singular value remains greater than 0.9.

.. note::

   The generated ``crystal-qha`` YAML deliberately retains ``method:
   crystal-qha``.  Running the Quantas tracker as a characterization test does
   not change the provenance of the scientific result that the source workflow
   itself established.

Real multi-volume regression: dolomite
--------------------------------------

Dataset
~~~~~~~

The real regression series contains seven independent CRYSTAL23 phonon
dispersion calculations for dolomite.  Each normalized state contains

.. code-block:: text

   primitive atoms:       10
   phonon q-points:       27
   modes per q-point:     30
   sampled volumes:        7

The input file order is not monotonic in volume.  Quantas sorts the states to

.. code-block:: text

   source indices: 6, 4, 1, 3, 0, 2, 5

before constructing the six adjacent-volume matching problems.

Tracking result
~~~~~~~~~~~~~~~

The validated regression result is

.. list-table:: Dolomite mode-continuity regression
   :header-rows: 1
   :widths: 42 22 36

   * - Quantity
     - Value
     - Interpretation
   * - Adjacent-volume assignments
     - 4860
     - :math:`6\times27\times30` local mode links.
   * - Ambiguous assignments
     - 274
     - Retained for inspection.
   * - Low-overlap assignments
     - 6
     - All occur in the widest final volume interval.
   * - Caution assignments
     - 274
     - Difficult but accepted local links.
   * - Unresolved assignments
     - 0
     - No continuity failure remains.
   * - Final reordered assignments
     - 420
     - Final branch labels differ from raw mode indices.
   * - Local reordered assignments
     - 261
     - Non-identity adjacent-volume assignments.
   * - Minimum non-degenerate overlap
     - 0.4209753995
     - Weakest selected scalar product.
   * - Minimum subspace singular value
     - 0.8821327956
     - Weakest matched degenerate eigenspace.

The interval summary shows a physically sensible degradation with increasing
volume separation.  The first five intervals contain no overlap below 0.5.  The
six low-overlap links are confined to the final interval

.. math::

   110.788757\rightarrow117.488206\ \mathrm{angstrom}^3,

which is deliberately wider than the preceding sampling steps.

This localization matters.  A parser or atom-ordering failure would normally
produce widespread loss of overlap rather than a controlled deterioration at
the largest structural separation.

Global branch diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~

For the seven-volume series, the global post-tracking diagnostic uses cubic
frequency paths with three residual degrees of freedom.  Across all
:math:`27\times30=810` q-point/branch paths, the current regression metadata
reports

.. code-block:: text

   maximum RMSE:          1.058647 cm^-1
   maximum abs residual:  1.882649 cm^-1
   diagnostic support:    810 / 810 branches

The minimum :math:`R^2` is approximately 0.118.  This apparently poor value
belongs to a nearly flat frequency path and is accompanied by a small absolute
error.  The case demonstrates why Quantas records :math:`R^2`, RMSE, and the
maximum residual rather than using :math:`R^2` alone.

Leave-one-out validation
~~~~~~~~~~~~~~~~~~~~~~~~

The six low-overlap links trigger the symmetric endpoint validation.  With
seven states, the predictive fit is cubic and each leave-one-out training fit
retains two residual degrees of freedom.  All six weak assignments satisfy the
independent prediction criterion at both endpoints:

.. code-block:: text

   supported low-overlap assignments: 6
   unresolved low-overlap assignments: 0

The dataset is therefore classified as

.. code-block:: yaml

   mode_continuity: verified

with explicit cautions retained in the metadata and debug report.

.. important::

   The dolomite result is **verified with cautions**, not "perfectly matched".
   The six low overlaps are scientifically visible and remain part of the
   provenance.  Verification means that no assignment is unresolved after the
   independent checks; it does not erase evidence of strong mode mixing.

Structural reconstruction regression
------------------------------------

The same dolomite input exercises the supercell-to-primitive reconstruction.
The CRYSTAL phonons were evaluated in a :math:`3\times3\times3` supercell,
while the normalized QHA structure contains ten primitive atoms.  Every source
state reconstructs exactly 27 translational copies per primitive atom.

The recorded maximum and RMS translational residuals are of order
:math:`10^{-12}` angstrom, far below the ``exact`` reconstruction threshold.
The normalized structural volumes and top-level thermodynamic volumes are
therefore consistent representations of the same sampled states.

Scientific limits of this validation
------------------------------------

The tests above establish that Quantas:

- parses the characterized CRYSTAL eigenvector representations consistently;
- preserves complex phase invariance;
- restores mass-weighted unit normalization;
- compares independent volumes on compatible q meshes;
- tracks non-degenerate modes by one-to-one overlap assignment;
- treats numerical degeneracies as eigenspaces;
- keeps local matching independent of reference labels;
- distinguishes cautions from unresolved assignments;
- requires independent leave-one-out evidence before rescuing weak overlaps;
- preserves real multi-volume continuity for the validated dolomite dataset.

They do **not** establish that every upstream phonon calculation is converged or
that QHA is physically adequate for every material.

.. warning::

   Mode-continuity validation must not be used as a substitute for phonon
   convergence and stability tests.  Quantas can track an unstable imaginary
   branch perfectly.  It can also track frequencies generated from an
   insufficient k mesh, basis set, supercell, or force-constant threshold.  The
   upstream electronic-structure and lattice-dynamical calculation remains part
   of the validation chain.

.. warning::

   The validated thresholds are operational criteria for this implementation,
   not universal constants of lattice dynamics.  A future change to the
   ambiguity margin, overlap threshold, degeneracy tolerance, subspace
   criterion, or leave-one-out tolerance is a scientific numerical change.  It
   requires new characterization tests and real-data comparison before being
   accepted.

Traceability to tests
---------------------

The principal regression and characterization coverage is located in:

``tests/interfaces/test_crystal_phonon_modes.py``
   CRYSTAL real/complex eigenvectors, normalization, partial blocks, native MgO
   QHA continuity, and source-mode characterization.

``tests/modules/qha/test_mode_tracking.py``
   Permutations, phase invariance, degenerate rotations, ambiguity semantics,
   reference independence, leave-one-out support/rejection, and the real
   dolomite regression.

``tests/modules/ha/test_input_generation.py``
   Multi-file generator integration, continuity metadata, neutral diagnostic
   tables, unit labels, YAML presentation, and serialization equivalence.

The staged project runner remains the authoritative whole-package regression
check before merging this work.

Related pages
-------------

- :doc:`../workflows/phonon_input_generation`
- :doc:`../workflows/ha`
- :doc:`../workflows/qha`
- :doc:`../formats/phonon_yaml`
- :doc:`strategy`
