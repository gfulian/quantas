Citing Quantas
==============

For the base use of Quantas, cite:

   G. Ulian and G. Valdrè, “QUANTAS, a Python software for the analysis of
   solids from ab initio quantum mechanical simulations and experimental
   data”, *Journal of Applied Crystallography* **55**, 386--396 (2022),
   doi:10.1107/S1600576722000085.

For phase and group acoustic velocities produced by SEISMIC, cite:

   G. Ulian and G. Valdrè, “SEISMIC, a Python-based code of the Quantas package
   to calculate the phase and group acoustic velocities in crystals”,
   *Computers & Geosciences* **188**, 105615 (2024),
   doi:10.1016/j.cageo.2024.105615.

Module reports may provide additional method-specific references. Users should
cite both Quantas and the scientific methods used in a calculation.


Harmonic and quasi-harmonic thermodynamics
-------------------------------------------

For the statistical thermodynamics of independent harmonic oscillators, cite:

   D. A. McQuarrie and J. D. Simon, *Physical Chemistry: A Molecular
   Approach*, University Science Books, Sausalito, California (1997).

For the general high-pressure thermodynamic framework used in QHA, cite:

   O. L. Anderson, *Equations of State of Solids for Geophysics and Ceramic
   Science*, Oxford Monographs on Geology and Geophysics **31**, Oxford
   University Press, New York (1995).

   O. L. Anderson, K. Masuda, and D. G. Isaak, “A new thermodynamic approach
   for high-pressure physics”, *Physics of the Earth and Planetary Interiors*
   **91**, 3--16 (1995), doi:10.1016/0031-9201(95)03044-W.

For the combination of pressure and temperature effects with standard
first-principles calculations, and for practical assessment of QHA behaviour,
cite:

   A. Erba, “On combining temperature and pressure effects on structural
   properties of crystals with standard ab initio techniques”, *Journal of
   Chemical Physics* **141**, 124115 (2014), doi:10.1063/1.4896228.

   A. Erba, M. Shahrokhi, R. Moradian, and R. Dovesi, “On how differently the
   quasi-harmonic approximation works for two isostructural crystals: Thermal
   properties of periclase and lime”, *Journal of Chemical Physics* **142**,
   044114 (2015), doi:10.1063/1.4906422.

Quantas reports also include the Quantas software citation and, for the QHA
workflow, the application reference by Ulian and Valdrè (2018) registered in
``quantas.references``.


Equation-of-state methods
-------------------------

For isothermal, thermal-expansion, P--V--T, linear-EOS conventions, and the
effective-variance weighting used by EOS, cite:

   R. J. Angel, J. Gonzalez-Platas, and M. Alvaro, “EosFit7c and a Fortran
   module (library) for equation of state calculations”, *Zeitschrift für
   Kristallographie* **229**, 405--419 (2014),
   doi:10.1515/zkri-2013-1711.

The effective-variance approach is attributed to:

   J. Orear, “Least squares when both variables have uncertainties”,
   *American Journal of Physics* **50**, 912--916 (1982).


For weighted orthogonal distance regression and ODRPACK, cite:

   P. T. Boggs, R. H. Byrd, J. E. Rogers, and R. B. Schnabel,
   *User's Reference Guide for ODRPACK Version 2.01: Software for Weighted
   Orthogonal Distance Regression*, NISTIR 4834 (1992),
   doi:10.6028/NIST.IR.4834.

For the bound-constrained ODRPACK95 implementation used by the Quantas
runtime, cite:

   J. W. Zwolak, P. T. Boggs, and L. T. Watson, “Algorithm 869: ODRPACK95:
   A weighted orthogonal distance regression code with bound constraints”,
   *ACM Transactions on Mathematical Software* **33** (4), Article 27 (2007),
   doi:10.1145/1268776.1268782.

Thermal-expansion models audited for EOS
----------------------------------------

The individual V--T formulations and their historical parameterizations are
attributed to the following sources, as reviewed and revalidated by Angel,
Gonzalez-Platas, and Alvaro (2014):

* R. G. Berman, *Journal of Petrology* **29**, 445--522 (1988).
* Y. Fei, “Thermal expansion”, in *Mineral Physics & Crystallography: A
  Handbook of Physical Constants*, Volume 2, 29--44 (1995).
* A. R. Pawley, S. A. T. Redfern, and T. J. B. Holland,
  *American Mineralogist* **81**, 335--340 (1996).
* E. K. H. Salje, B. Wruck, and H. Thomas,
  *Zeitschrift für Physik B* **82**, 399--404 (1991).
* T. J. B. Holland and R. Powell, *Journal of Metamorphic Geology* **29**,
  333--383 (2011).
* G. Hellfrich and J. A. D. Connolly, *American Mineralogist* **94**,
  1616--1620 (2009), for Anderson--Grüneisen coupling of thermal expansion
  and the zero-pressure bulk modulus.
* H. Kroll, A. Kirfel, R. Heinemann, and B. Barbier,
  *European Journal of Mineralogy* **24**, 935--956 (2012).

The Angel et al. (2014) formulations and the independent EOS reference
snapshot define the Quantas implementation target. The public validation
record is maintained in :doc:`../validation/eos`.

Terrestrial pressure--temperature profiles
-------------------------------------------

For pressure reconstructed from the Preliminary Reference Earth Model, cite:

   A. M. Dziewonski and D. L. Anderson, “Preliminary reference Earth model”,
   *Physics of the Earth and Planetary Interiors* **25** (4), 297--356 (1981),
   doi:10.1016/0031-9201(81)90046-7.

For the layered continental conductive framework and representative
lithospheric geotherms, cite:

   D. Hasterok and D. S. Chapman, “Heat production and geotherms for the
   continental lithosphere”, *Earth and Planetary Science Letters* **307**
   (1--2), 59--70 (2011), doi:10.1016/j.epsl.2011.04.034.

For oceanic half-space and finite-plate cooling references, cite:

   B. Parsons and J. G. Sclater, “An analysis of the variation of ocean floor
   bathymetry and heat flow with age”, *Journal of Geophysical Research* **82**
   (5), 803--827 (1977), doi:10.1029/JB082i005p00803.

For the dry-pyrolite mantle adiabat and its archived scripts, cite both:

   T. Katsura, “A revised adiabatic temperature profile for the mantle”,
   *Journal of Geophysical Research: Solid Earth* **127** (2), e2021JB023562
   (2022), doi:10.1029/2021JB023562.

   T. Katsura, *Matlab scripts of “A revised adiabatic temperature profile for
   the mantle”*, Version 1.1.0, Zenodo (2022),
   doi:10.5281/zenodo.5903286.

The Quantas Katsura implementation is a deterministic reconstruction from the
published temperature and gradient constraints. It is not a re-execution of
the complete Monte Carlo MATLAB workflow.

Thermoelastic and adiabatic elastic tensors
--------------------------------------------

For the general quasi-harmonic formulation of thermoelastic stiffness,
the quasi-static approximation, and the conversion from isothermal to
adiabatic elastic constants, cite:

   M. Destefanis, C. Ravoux, A. Cossard, and A. Erba,
   “Thermo-Elasticity of Materials from Quasi-Harmonic Calculations”,
   *Minerals* **9**, 16 (2019), doi:10.3390/min9010016.

For the thermodynamically self-consistent Eulerian finite-strain derivation of
the cold elastic tensor and its quasi-harmonic extension, cite:

   L. Stixrude and C. Lithgow-Bertelloni, “Thermodynamics of mantle
   minerals—I. Physical properties”, *Geophysical Journal International*
   **162**, 610--632 (2005), doi:10.1111/j.1365-246X.2005.02642.x.

For the foundational quasi-harmonic treatment of elastic moduli under
hydrostatic pre-stress, cite:

   G. F. Davies, “Effective elastic moduli under hydrostatic stress—I.
   Quasi-harmonic theory”, *Journal of Physics and Chemistry of Solids* **35**,
   1513--1520 (1974), doi:10.1016/S0022-3697(74)80279-9.

For the thermodynamic relation between anisotropic isothermal and adiabatic
elastic tensors, cite:

   M. J. Waters and A. W. Bielawski, “Isothermal and adiabatic elastic
   tensors”, arXiv:1605.06548 (2016), doi:10.48550/arXiv.1605.06548.

The tensor thermodynamics and notation also follow:

   D. C. Wallace, *Thermodynamics of Crystals*, John Wiley & Sons, New York
   (1972).

Thermoelastic validation systems
--------------------------------

For experimental and first-principles pressure-dependent elasticity of the
cubic MgO validation system, see:

   B. B. Karki, L. Stixrude, S. J. Clark, M. C. Warren, G. J. Ackland, and
   J. Crain, “Structure and elasticity of MgO at high pressure”, *American
   Mineralogist* **82**, 51--60 (1997), doi:10.2138/am-1997-1-207.

   S. V. Sinogeikin and J. D. Bass, “Single-crystal elasticity of MgO at high
   pressure”, *Physical Review B* **59**, R14141--R14144 (1999),
   doi:10.1103/PhysRevB.59.R14141.

For ambient single-crystal elasticity of the low-trigonal dolomite validation
system, see:

   F. Jiang, S. Speziale, and T. S. Duffy, “Elasticity of magnesite and
   dolomite from a genetic algorithm for inverting Brillouin spectroscopy
   measurements”, *Physics of the Earth and Planetary Interiors* **155**,
   1--20 (2006), doi:10.1016/j.pepi.2005.08.004.

Elasticity and seismic-wave analysis
-------------------------------------

For tensor notation and the transformation properties of elastic coefficients,
see:

   J. F. Nye, *Physical Properties of Crystals: Their Representation by
   Tensors and Matrices*, second edition, Oxford University Press, Oxford
   (1985).

For the Voigt--Reuss bounds and the Hill average of a crystalline aggregate,
cite:

   R. Hill, “The elastic behaviour of a crystalline aggregate”,
   *Proceedings of the Physical Society. Section A* **65**, 349--354 (1952),
   doi:10.1088/0370-1298/65/5/307.

For directional elastic-property surfaces and their interpretation, cite:

   R. Gaillac, P. Pullumbi, and F.-X. Coudert, “ELATE: an open-source online
   application for analysis and visualization of elastic tensors”,
   *Journal of Physics: Condensed Matter* **28**, 275201 (2016),
   doi:10.1088/0953-8984/28/27/275201.

For phase velocity, group velocity, polarization, and acoustic enhancement from
the Christoffel equation, cite:

   J. W. Jaeken and S. Cottenier, “Solving the Christoffel equation: Phase and
   group velocities”, *Computer Physics Communications* **207**, 445--451
   (2016), doi:10.1016/j.cpc.2016.06.014.
