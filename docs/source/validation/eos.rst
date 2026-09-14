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

A BM3 OLS fit gives:

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
energy-derived pressure for thermoelastic input generation.  Reproducing the
same equilibrium parameters through ``quantas eos run --domain ev`` therefore
checks that the public workflow has not changed the underlying numerical Energy
EOS service.

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
normalizations.  Multiplying primitive-cell energy and volume by four while
transforming the FCC primitive lattice to the conventional cubic cell changes
``E0`` and ``V0`` by exactly the same factor but leaves ``K0``, ``KP``, and the
axial response invariant.  The crystallographic regression target is
approximately ``a0 = 4.22221`` angstrom and ``M_a = 536.28`` GPa, with
``M_a = 3 K0`` to numerical precision.

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

Release-candidate validation still to close
-------------------------------------------

Before ``2.0.0rc1`` the broader EOS validation matrix will also collect the
existing P--V, V--T, and P--V--T analytical and real-data comparisons in one
place, including EosFit7 reference cases where appropriate.  That documentation
work does not change the already characterized Energy EOS formulas or public
``ev/energy`` contract.
