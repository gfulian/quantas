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

The most extensive recent validation concerns the first question, because the
CRYSTAL interface and phonon-mode continuity work introduced new scientific
contracts that needed direct evidence.  The HA/QHA
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

One QHA validation case deliberately scrambles the input-volume order relative
to the increasing cutoff series.  The test then verifies that explicit
one-to-one volume matches, rather than array indices, determine the acoustic
association.  The
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
characterization datum, while the new core uses the published formula.  The corrected result is independently constrained by
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

The same driver also checks that the workflow rejects a base input without an
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

Phonon parsing and mode-continuity validation
---------------------------------------------

The CRYSTAL phonon input path is validated at four complementary levels. The
algorithm and its user-facing diagnostics are described in
:doc:`../workflows/phonon_input_generation`; this page records the evidence that
the implementation satisfies those contracts.

.. list-table:: Mode-continuity validation layers
   :header-rows: 1
   :widths: 25 36 39

   * - Layer
     - Evidence
     - What it establishes
   * - Parser characterization
     - Real CRYSTAL Gamma and non-Gamma outputs, including partial final blocks
     - Array dimensions, atom ordering, complex phase reconstruction, and mass-weighted normalization
   * - Synthetic invariants
     - Permutations, arbitrary complex phases, rotated degenerate subspaces, and changed reference indices
     - Assignment invariance to representation choices that must not change physical mode identity
   * - Adversarial continuity cases
     - Weak overlaps, ambiguous links, leave-one-out support and deliberately unsupported endpoints
     - Difficult assignments remain cautions or become ``unreliable`` instead of being forced to pass
   * - Real multi-volume regression
     - Native MgO QHA output and seven-volume dolomite phonon dispersions
     - Complete parser/tracker behavior on realistic degeneracies, many branches, and non-monotonic source ordering

Parser characterization
~~~~~~~~~~~~~~~~~~~~~~~

The distributed dolomite fixture contains 27 q-points, 30 modes per q-point,
and 10 atoms. It produces frequency arrays of shape ``(27, 30)`` and complex
eigenvector arrays of shape ``(27, 30, 10, 3)``. Gamma eigenvectors are real;
non-zero q-points reconstruct complex components when CRYSTAL provides
in-phase and anti-phase sections. Inconsistent real/imaginary blocks are
rejected.

Every parsed mode is normalized to unit mass-weighted norm. A native MgO QHA
fixture provides an independent physical check: the relative acoustic
translation amplitudes reproduce the expected square-root mass ratio for Mg and
O.

Synthetic invariants and failure semantics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Controlled tests verify exact permutation recovery, invariance to arbitrary
complex phase, and basis invariance inside a degenerate eigenspace. Local
matching is independent of the chosen final reference labels and follows sorted
adjacent volumes even when the source files are supplied out of order.

Ambiguity is retained as a caution when the assignment remains usable. A weak
overlap with no independent evidence is ``unreliable``. When enough volumes are
available, a weak link may be accepted only if a leave-one-out frequency model
fitted **without the tested endpoint** predicts that endpoint within tolerance.
An adversarial fixture confirms that a global fit cannot rescue its own
outlier.

Native MgO and dolomite regressions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The native CRYSTAL17 MgO QHA output contains eleven Gamma-mode sets for a
64-atom phonon supercell. CRYSTAL explicitly reports frequency continuity with
volume; the normalized YAML therefore records source-managed ``verified``
continuity. Running the Quantas overlap tracker as an independent
characterization produces no unresolved assignments, without replacing that
source provenance.

The real dolomite regression contains seven independent CRYSTAL23 dispersion
calculations with 27 q-points and 30 modes at each volume. The source files are non-monotonic in volume, which prevents the regression
from relying accidentally on file order. Across the six adjacent-volume matching
steps the regression contains 4860 local links, including 274 cautions and six
low-overlap assignments. All six weak links occur in the widest volume
interval and are independently supported by leave-one-out checks; none remains
unresolved. The dataset is therefore classified ``verified`` **with cautions**,
not as a perfect one-to-one scalar match.

The same regression also validates supercell-to-primitive structural
reconstruction. The 3x3x3 phonon supercell maps to ten primitive atoms with 27
translational copies per atom, and translational residuals remain at numerical
round-off.

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
