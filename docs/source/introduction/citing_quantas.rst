Citing Quantas
==============

For the base use of Quantas, cite the main Quantas application paper
[#quantas_2022]_.  For phase and group acoustic velocities produced by SEISMIC,
also cite the dedicated SEISMIC paper [#seismic_ulian_valdre_2024]_.

Module reports may provide additional method-specific references. Users should
cite both Quantas and the scientific methods used in a calculation.


Harmonic and quasi-harmonic thermodynamics
-------------------------------------------

For the statistical thermodynamics of independent harmonic oscillators, cite
McQuarrie and Simon [#mcquarrie_simon_1997]_.

For the general high-pressure thermodynamic framework used in QHA, cite
Anderson [#anderson_1995]_ and Anderson, Masuda, and Isaak
[#anderson_masuda_isaak_1995]_.  For the combination of pressure and
temperature effects with standard first-principles calculations, and for
practical assessment of QHA behaviour, cite the corresponding studies by Erba
[#erba_2014]_ and Erba *et al.* [#erba_shahrokhi_moradian_dovesi_2015]_.

Quantas reports also include the Quantas software citation and, for the QHA
workflow, the Quantas QHA application paper [#qha_ulian_valdre_2018]_.


Equation-of-state methods
-------------------------

For isothermal, thermal-expansion, P--V--T, linear-EOS conventions, and the
EosFit7 framework used as a principal EOS reference, cite Angel, Gonzalez-Platas,
and Alvaro [#eosfit7_angel_gonzalez_platas_alvaro_2014]_.  The same reference
defines the modified Tait pressure parameterization used by Quantas.
Volume-integrated E--V forms retain the corresponding physical EOS parameters
and are constrained by :math:`P(V)=-\mathrm dE/\mathrm dV`; the integrated Tait
expression is obtained analytically from that registered pressure form.

For the stabilized-jellium energy EOS (SJEOS), also cite Alchagirov *et al.*
[#sjeos_alchagirov_perdew_boettger_albers_fiolhais_2001]_.

The effective-variance approach is attributed to Orear [#orear_1982]_.  For
weighted orthogonal distance regression and ODRPACK, cite the ODRPACK reference
guide [#boggs_byrd_rogers_schnabel_1992]_ and the ODRPACK95 implementation
paper [#zwolak_boggs_watson_2007]_.


Thermal-expansion models audited for EOS
----------------------------------------

The individual V--T formulations and their historical parameterizations are
attributed to Berman [#berman_1988]_, Fei [#fei_1995]_, Pawley *et al.*
[#pawley_redfern_holland_1996]_, Salje *et al.* [#salje_wruck_thomas_1991]_,
Holland and Powell [#holland_powell_2011]_, Helffrich and Connolly
[#helffrich_connolly_2009]_, and Kroll *et al.*
[#kroll_kirfel_heinemann_barbier_2012]_.  These formulations are reviewed and
revalidated in the EosFit7 reference above
[#eosfit7_angel_gonzalez_platas_alvaro_2014]_.  The public Quantas validation
record is maintained in :doc:`../validation/eos`.


Terrestrial pressure--temperature profiles
-------------------------------------------

For pressure reconstructed from the Preliminary Reference Earth Model, cite
Dziewonski and Anderson [#prem_dziewonski_anderson_1981]_.  For the layered
continental conductive framework and representative lithospheric geotherms,
cite Hasterok and Chapman [#hasterok_chapman_2011]_.  For oceanic half-space
and finite-plate cooling, cite Parsons and Sclater [#parsons_sclater_1977]_.

For the dry-pyrolite mantle adiabat and its archived scripts, cite both the
published profile [#katsura_2022]_ and the corresponding software archive
[#katsura_software_2022]_.  The Quantas Katsura implementation is a
deterministic reconstruction from the published temperature and gradient
constraints. It is not a re-execution of the complete Monte Carlo MATLAB
workflow.


Thermoelastic and adiabatic elastic tensors
--------------------------------------------

For the general quasi-harmonic formulation of thermoelastic stiffness, the
quasi-static approximation, and the conversion from isothermal to adiabatic
elastic constants, cite Destefanis *et al.*
[#destefanis_ravoux_cossard_erba_2019]_.

For CRYSTAL elastic tensors evaluated or reconstructed under hydrostatic
pre-stress, cite the finite-pressure formulation implemented by CRYSTAL
[#erba_mahmoud_belmonte_dovesi_2014]_.  For the thermodynamically
self-consistent Eulerian finite-strain derivation of the cold elastic tensor
and its quasi-harmonic extension, cite Stixrude and Lithgow-Bertelloni
[#stixrude_lithgow_bertelloni_2005]_.

For the foundational quasi-harmonic treatment of elastic moduli under
hydrostatic pre-stress, cite Davies [#davies_1974]_.  For the thermodynamic
relation between anisotropic isothermal and adiabatic elastic tensors, cite
Waters and Bielawski [#waters_bielawski_2016]_.  The tensor thermodynamics and
notation also follow Wallace [#wallace_1972]_.


Thermoelastic validation systems
--------------------------------

For experimental and first-principles pressure-dependent elasticity of the
cubic MgO validation system, see Karki *et al.*
[#karki_stixrude_clark_warren_ackland_crain_1997]_ and Sinogeikin and Bass
[#sinogeikin_bass_1999]_.  For ambient single-crystal elasticity of the
low-trigonal dolomite validation system, see Jiang, Speziale, and Duffy
[#jiang_speziale_duffy_2006]_.


Elasticity and seismic-wave analysis
-------------------------------------

For tensor notation and the transformation properties of elastic coefficients,
see Nye [#nye_1985]_.  For the Voigt--Reuss bounds and the Hill average of a
crystalline aggregate, cite Hill [#hill_1952]_.  For directional
elastic-property surfaces and their interpretation, cite Gaillac, Pullumbi,
and Coudert [#elate_gaillac_pullumbi_coudert_2016]_.  For phase velocity, group
velocity, polarization, and acoustic enhancement from the Christoffel equation,
cite Jaeken and Cottenier [#jaeken_cottenier_2016]_.


.. include:: ../_generated/references/introduction_citing_quantas.inc
