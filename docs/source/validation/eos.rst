EOS validation
==============

The EOS validation strategy combines analytical identities, synthetic parameter
recovery, real-data end-to-end workflows, and comparison with established EOS
implementations where suitable.  The complete release-candidate matrix for the
experimental P--V, V--T, and P--V--T domains is still being consolidated.  The
``2.0.0b11`` Energy EOS tranche already has a dedicated validation set because
it introduces the public ``ev/energy`` workflow.

Energy EOS analytical checks
----------------------------

Every integrated Energy EOS is tied to its pressure representation through

.. math::

   P(V)=-\frac{dE}{dV}.

The test suite checks this relation analytically or numerically for the shared
integrated models.  Modified Tait is additionally compared with direct numerical
integration of its independently evaluated pressure equation, including the
removable logarithmic limit.  SJEOS tests cover its physical equilibrium
parameterization, analytical pressure and bulk-modulus derivatives, and the
implied second pressure derivative.

Synthetic public-workflow recovery
----------------------------------

The standalone E--V adapter is tested on DFT-scale energies rather than on
small values centred artificially around zero.  Synthetic BM3, T3, and SJEOS
datasets use absolute energies near -275 Ha while the physically relevant
energy differences are many orders of magnitude smaller.  Public fitting must
recover ``E0``, ``V0``, ``K0``, ``KP``, and ``KPP`` in the documented units and
preserve the derived arrays for energy, pressure, bulk modulus, and its pressure
derivatives.

The tests also characterize:

- OLS fitting without invented statistical energy errors;
- WLS when a genuine ``sigma_energy`` is supplied;
- ``ev/energy`` HDF5 round trips;
- diagnostics and optional comparison with an independently supplied pressure;
- forward ``V -> E,P,K,K',K''`` and inverse ``P -> V`` calculator paths;
- fit, pressure, and residual plot inventories;
- direct CLI execution and ``[defaults.ev]`` specification resolution.

Real MgO end-to-end regression
------------------------------

A seven-volume CRYSTAL/PBE MgO series is used as the real-data regression for
the complete public Energy EOS path.  The normalized dataset spans approximately
17.11--20.50 angstrom cubed and is processed through the same text reader used
by ordinary users.

A representative BM3 OLS fit in the reference environment gives:

.. list-table:: MgO Energy EOS regression target
   :header-rows: 1
   :widths: 30 30 40

   * - Quantity
     - Reference value
     - Unit
   * - ``E0``
     - -275.173937178
     - Ha
   * - ``V0``
     - 18.817428245
     - angstrom cubed
   * - ``K0``
     - 178.761458
     - GPa
   * - ``KP``
     - 3.815504
     - 1
   * - ``KPP``
     - -0.02091296
     - GPa :math:`^{-1}`
   * - energy RMSE
     - 3.28e-6
     - Ha

The same MgO static surface had already been used independently while validating
energy-derived pressure for thermoelastic input generation.  The cross-platform
regression gate therefore checks the numerically robust ``E0``, ``V0``, ``K0``,
and energy-RMSE observables.  ``KP`` and the BM3-implied ``KPP`` are reported
above as representative reference-environment values, but they are not used as
tight real-data CI gates because higher pressure derivatives are more sensitive
to nonlinear-solver termination across supported SciPy/BLAS combinations.
Their formulas and parameter recovery are characterized independently by the
synthetic Energy EOS and core-physics tests.

Theoretical structural-response regression
------------------------------------------

The Energy EOS structural extension is characterized independently of the
pressure-only and thermal EOS paths.  Synthetic cubic data verify the exact
identity :math:`d\ln a/d\ln V=1/3` and therefore :math:`M_a=3K`; a synthetic
tetragonal path verifies independent ``a`` and ``c`` logarithmic responses and
axial moduli.  The same tests exercise SJEOS as the primary Energy EOS, proving
that the structural response depends on the analytical E(V) derivatives rather
than on availability of a direct experimental P--V fitting surface.

The curated MgO series is stored in both primitive and crystallographic
normalizations.  The input data verify the extensive normalization exactly:
primitive-cell energy and volume are multiplied by four while the FCC
primitive lattice is transformed to the conventional cubic cell.  Independent
nonlinear fits preserve the corresponding ``E0``/``V0`` scaling and the
intensive ``K0`` response within a tight cross-platform solver tolerance.
This distinction avoids treating platform-dependent least-squares termination
at the last few digits as a scientific regression.  The crystallographic
response remains approximately ``a0 = 4.22221`` angstrom and
``M_a = 536.28`` GPa, with ``M_a = 3 K0`` to numerical precision.

Optional secondary pressure-form axial EOS fits are characterized separately.
A primary SJEOS E(V) fit may feed a BM3 :math:`P(a^3)` fit; the full covariance
of the derived pressure vector survives the HDF5 round trip, while the current WLS
fit is explicitly labelled as a diagonal marginal-uncertainty approximation.
The primary ``a0``, ``eta_a``, and ``M_a`` results do not depend on requesting
this secondary parameterization.

CRYSTAL input-generation characterization
-----------------------------------------

The CRYSTAL Energy EOS interface is tested separately from the fit.  It must
normalize three source patterns to the same backend-neutral structure--energy
contract:

- one static calculation -> one state;
- one completed geometry optimization -> one final state;
- one native CRYSTAL ``EOS`` calculation -> multiple volume states.

Real urea calculations containing DFT-D3 and DFT-D3+gCP corrections are used to
characterize energy semantics.  The final CRYSTAL volume--energy summary defines
which states belong to the native EOS curve, while each state is independently
matched to its final optimized geometry and to the authoritative high-precision
corrected total energy.  ``OPT END`` or an intermediate SCF energy is not used
when a later corrected total is available.

Mixed atom counts, chemical compositions, total-energy correction signatures,
and duplicate volumes are rejected.  Nearby but genuinely distinct volumes are
preserved.  This protects the fitted E(V) surface from combining calculations
that are individually parseable but scientifically incompatible.

Experimental P--V reference regression
-------------------------------------------

The pressure--volume path has an external-reference regression for third-order
Birch--Murnaghan fits. Quartz and topaz datasets are compared with frozen
EosFit7-compatible results using both ordinary least squares and the effective
variance treatment used for uncertainties in both pressure and volume.

Representative reference parameters are:

.. list-table:: BM3 P--V external-reference checkpoints
   :header-rows: 1
   :widths: 18 18 18 18 18 10

   * - Dataset / solver
     - ``V0``
     - ``K0`` (GPa)
     - ``KP``
     - ``KPP`` (GPa :math:`^{-1}`)
     - reduced :math:`\chi^2`
   * - Quartz / OLS
     - 112.96752
     - 37.28543
     - 5.93351
     - -0.25642
     - --
   * - Quartz / effective variance
     - 112.98088
     - 37.12600
     - 5.98823
     - -0.26478
     - 0.95
   * - Topaz / OLS
     - 346.97214
     - 135.74806
     - 4.38050
     - -0.03252
     - --
   * - Topaz / effective variance
     - 345.50726
     - 161.99034
     - 2.97647
     - -0.02416
     - 18.55

The regression also checks parameter standard errors, the largest residual, and
the inflate-only covariance policy. Parameter tolerances are typically a few
times :math:`10^{-4}` in the reported units; the more nonlinear topaz effective-
variance ``V0`` comparison allows an absolute tolerance of 0.002. Those values
are regression tolerances for reproducing the external reference calculation,
not statements about experimental accuracy.

Validation still in progress
----------------------------

.. admonition:: Work in progress

   The V--T and P--V--T domains already have analytical, synthetic, workflow,
   and tutorial regressions, but their consolidated external/reference matrix
   has not yet been assembled to the release-candidate standard defined in
   :doc:`strategy`. They therefore remain work in progress even though the
   implementations are extensively tested.

Traceability
------------

The principal EOS validation coverage is located in:

``tests/physics/eos/test_energy_models.py`` and
``tests/physics/eos/test_energy_pressure.py``
   Integrated E(V) models, analytical pressure derivatives, and energy/pressure
   consistency.

``tests/modules/eos/test_energy_volume_workflow.py`` and
``tests/modules/eos/test_energy_structural_response.py``
   Public E--V fitting, structural normalization, derived properties, and
   secondary axial response.

``tests/interfaces/test_crystal_energy_volume.py``
   CRYSTAL state extraction and authoritative energy semantics.

``tests/examples/test_curated_examples.py``
   Curated MgO end-to-end Energy EOS regression.

``tests/modules/eos/test_eosfit_reference.py``
   Quartz/topaz BM3 OLS and effective-variance regression against the frozen
   EosFit7-compatible reference results listed above.

The remaining V--T and P--V--T tests protect implementation behaviour while
their public validation record is completed.
