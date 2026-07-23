Equations of state
==================

The scientific reference for the equation-of-state framework adopted by
Quantas is the review and revalidation by Angel, Gonzalez-Platas, and Alvaro,
*EosFit7c and a Fortran module (library) for equation of state calculations*
[#eosfit7_angel_gonzalez_platas_alvaro_2014]_.  That work reviews the
connection between linear elasticity and equations of state, derives the
implemented formulae, and discusses their ranges of validity.  The complete
citation guidance is collected in :doc:`../introduction/citing_quantas`.

An equation of state (EOS) describes how a structural quantity changes with
pressure, temperature, or both.  For a volume EOS, the fundamental variables
are pressure :math:`P`, volume :math:`V`, and temperature :math:`T`.  Quantas
also applies the same formalism to linear quantities, such as lattice
parameters, through the conventions described at the end of this chapter.

The presentation follows the order used by Angel *et al.* (2014):

#. pressure--volume equations at fixed temperature;
#. volume--temperature equations at fixed pressure;
#. pressure--volume--temperature equations obtained by coupling compression
   and thermal expansion.

The equations below define the scientific models.  Parameter estimation,
uncertainty treatment, diagnostics, and command-line options belong to the
:doc:`../workflows/eos` workflow rather than to this chapter.

Reference state and notation
----------------------------

At a selected reference temperature :math:`T_{\mathrm{ref}}` and zero
pressure, let

.. math::

   V_0 = V(P=0,T_{\mathrm{ref}}),

.. math::

   K_0 = -V\left(\frac{\partial P}{\partial V}\right)_{T,
   P=0,T=T_{\mathrm{ref}}},

with pressure derivatives

.. math::

   K'_0 = \left(\frac{\partial K}{\partial P}\right)_{P=0,T=T_{\mathrm{ref}}},
   \qquad
   K''_0 = \left(\frac{\partial^2 K}{\partial P^2}\right)_{P=0,T=T_{\mathrm{ref}}}.

Pressure is positive under compression.  :math:`K_0` has the same unit as
pressure, :math:`K'_0` is dimensionless, and :math:`K''_0` has units of
inverse pressure.

The order assigned to an EOS identifies which terms of its defining expansion
are retained.  A lower-order equation fixes or implies one or more derivatives
of the bulk modulus; a higher-order equation introduces additional freedom,
but that freedom is useful only when the data span a sufficient pressure or
volume range.

Pressure--volume equations of state
-----------------------------------

An isothermal P--V equation extends infinitesimal linear elasticity to finite
compression.  Linear elasticity defines the bulk modulus at one state, but it
does not uniquely determine how :math:`K` changes over a finite pressure
interval.  Different EOS families therefore embody different assumptions
about the evolution of strain energy or bulk modulus.

Murnaghan equation
~~~~~~~~~~~~~~~~~~

The Murnaghan EOS assumes a bulk modulus that varies linearly with pressure,

.. math::

   K(P)=K_0+K'_0P.

Integration gives

.. math::

   P(V)=\frac{K_0}{K'_0}
   \left[
   \left(\frac{V_0}{V}\right)^{K'_0}-1
   \right],

or, equivalently,

.. math::

   V(P)=V_0
   \left(1+\frac{K'_0P}{K_0}\right)^{-1/K'_0}.

**Advantages**

- It is algebraically simple and exactly invertible between :math:`P` and
  :math:`V`.
- Its parameters have an immediate physical interpretation.
- It is useful for modest compressions and for applications that repeatedly
  require both :math:`P(V)` and :math:`V(P)`.

**Limitations**

- The assumption implies :math:`K''=0` at every pressure.
- Most solids exhibit a small negative :math:`K''_0`, so the equation becomes
  progressively less realistic as compression increases.
- Angel *et al.* note that it is normally reliable only over relatively small
  compressions, approximately up to :math:`V/V_0\simeq0.9` for many solids.

Tait equation
~~~~~~~~~~~~~

The modified Tait EOS is an invertible generalization of Murnaghan that allows
non-zero curvature in :math:`K(P)`.  In the parameterization used by Angel
*et al.*, define

.. math::

   a=\frac{1+K'_0}{1+K'_0+K_0K''_0},

.. math::

   b=\frac{K'_0}{K_0}-\frac{K''_0}{1+K'_0},

.. math::

   c=\frac{1+K'_0+K_0K''_0}
   {(K'_0)^2+K'_0-K_0K''_0}.

The pressure form is

.. math::

   P(V)=\frac{1}{b}
   \left\{
   \left[
   \frac{V/V_0+a-1}{a}
   \right]^{-1/c}-1
   \right\},

and the inverse relation is

.. math::

   V(P)=V_0\left[1-a\left(1-(1+bP)^{-c}\right)\right].

Quantas supports second-, third-, and fourth-order Tait forms.  The
third-order form uses the common truncation
:math:`K''_0=-K'_0/K_0`; the fourth-order form allows :math:`K''_0` to remain
independent.

**Advantages**

- It retains the practical invertibility of Murnaghan.
- It represents the observed curvature of the bulk modulus more realistically.
- Angel *et al.* report that Tait and Birch--Murnaghan commonly yield mutually
  consistent parameters for ordinary solid P--V datasets.

**Limitations**

- It is primarily an empirical compressional representation rather than an
  expansion derived from a specific finite-strain energy.
- Its auxiliary parameters :math:`a`, :math:`b`, and :math:`c` can become
  singular for non-physical combinations of :math:`K_0`, :math:`K'_0`, and
  :math:`K''_0`.
- A fourth-order form is not automatically preferable: :math:`K''_0` is often
  weakly constrained and strongly correlated with the lower derivatives.

Birch--Murnaghan equation
~~~~~~~~~~~~~~~~~~~~~~~~~

The Birch--Murnaghan EOS is derived by expanding the strain energy in the
Eulerian finite strain

.. math::

   f_E=\frac{1}{2}
   \left[
   \left(\frac{V_0}{V}\right)^{2/3}-1
   \right].

To fourth order in the strain expansion, the pressure is

.. math::

   P=3K_0f_E(1+2f_E)^{5/2}
   \left\{
   1+\frac{3}{2}(K'_0-4)f_E
   +\frac{3}{2}
   \left[
   K_0K''_0+(K'_0-4)(K'_0-3)+\frac{35}{9}
   \right]f_E^2
   \right\}.

The normalized pressure

.. math::

   F_E=\frac{P}{3f_E(1+2f_E)^{5/2}}

is a polynomial in :math:`f_E`, which makes an :math:`F_E`--:math:`f_E` plot a
useful diagnostic of EOS order.

The supported truncations are:

- **second order:** :math:`K'_0=4`;
- **third order:** :math:`K'_0` is independent and
  :math:`K''_0` is implied by

  .. math::

     K''_0=-\frac{1}{K_0}
     \left[(3-K'_0)(4-K'_0)+\frac{35}{9}\right];

- **fourth order:** :math:`K''_0` is independent.

**Advantages**

- It is grounded in finite-strain theory and has a clear hierarchy of orders.
- It is widely used for solids at moderate compression.
- The normalized-pressure representation provides a direct visual test of
  whether higher-order curvature is supported.
- It is also naturally connected to integrated energy--volume equations used
  in first-principles calculations.

**Limitations**

- A finite-order Eulerian expansion is not expected to remain accurate under
  extreme compression; Angel *et al.* identify approximately
  :math:`V/V_0<0.6` as a regime where finite-strain EOS can become inadequate
  for many solids.
- Higher-order terms require a sufficiently broad and precise dataset.
- Apparent improvement from a fourth-order fit may reflect parameter
  correlation rather than physically resolved curvature.

Natural-strain or Poirier--Tarantola equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The natural-strain EOS uses the Hencky strain

.. math::

   f_N=\frac{1}{3}\ln\left(\frac{V_0}{V}\right),

chosen here to be positive on compression.  Its fourth-order pressure form is

.. math::

   P=3K_0\frac{V_0}{V}f_N
   \left\{
   1+\frac{3}{2}(K'_0-2)f_N
   +\frac{3}{2}
   \left[
   1+K_0K''_0+(K'_0-2)+(K'_0-2)^2
   \right]f_N^2
   \right\}.

The corresponding lower-order forms imply:

- **second order:** :math:`K'_0=2`;
- **third order:**

  .. math::

     K''_0=-\frac{1}{K_0}
     \left[1+(K'_0-2)+(K'_0-2)^2\right];

- **fourth order:** :math:`K''_0` is independent.

**Advantages**

- It provides a mathematically consistent alternative finite-strain measure.
- It can be useful when the natural strain produces a simpler normalized
  relation for a particular material or pressure interval.
- Like Birch--Murnaghan, it offers a systematic order hierarchy.

**Limitations**

- The implied :math:`K''_0` of the third-order form is commonly more negative
  in magnitude than the Birch--Murnaghan value.
- Angel *et al.* observe that this often gives a poorer description of ordinary
  solid P--V data.
- Agreement or disagreement with Birch--Murnaghan is dataset dependent and
  should be interpreted as model sensitivity, not merely as a statistical
  ranking.

Vinet equation
~~~~~~~~~~~~~~

The Vinet EOS was derived from a generalized interatomic potential rather than
from a finite-order strain-energy polynomial.  Define

.. math::

   f_V=1-\left(\frac{V}{V_0}\right)^{1/3},
   \qquad
   \eta=\frac{3}{2}(K'_0-1).

Then

.. math::

   P=3K_0\frac{f_V}{(1-f_V)^2}\exp(\eta f_V).

The standard form contains :math:`V_0`, :math:`K_0`, and :math:`K'_0` and is
conventionally called third order.  Quantas also provides a second-order
representation with the implied value :math:`K'_0=1`.

**Advantages**

- It is designed to behave more realistically than finite-strain polynomials
  at very large compression.
- It is compact and often effective over broad pressure ranges.
- Its normalized-pressure form can be used to assess deviations from the
  assumed exponential behavior.

**Limitations**

- There is no compelling theoretical truncation to a lower-order Vinet EOS;
  the second-order form is mainly a constrained representation.
- A model intended for high compression is not necessarily the best choice for
  a narrow low-pressure dataset.
- The ordinary three-parameter form fixes the implied :math:`K''_0`; resolving
  additional curvature requires a different extended model.

Choosing a P--V family
~~~~~~~~~~~~~~~~~~~~~~

No single P--V EOS is universally superior.  A scientifically defensible
choice considers the compression range, the physical origin of the data, and
whether the dataset genuinely resolves higher derivatives.

.. list-table:: P--V model summary
   :header-rows: 1
   :widths: 19 25 28 28

   * - Family
     - Main idea
     - Principal strength
     - Principal caution
   * - Murnaghan
     - Linear :math:`K(P)`
     - Simple and invertible
     - :math:`K''=0`; limited compression range
   * - Tait
     - Invertible nonlinear :math:`K(P)`
     - Flexible and practical
     - Empirical auxiliary parameters; high-order correlation
   * - Birch--Murnaghan
     - Eulerian finite-strain expansion
     - Standard solid-state EOS and clear order diagnostics
     - Finite-order expansion degrades at extreme compression
   * - Poirier--Tarantola
     - Natural/Hencky strain expansion
     - Alternative finite-strain description
     - Often implies stronger curvature and poorer ordinary P--V fits
   * - Vinet
     - Generalized interatomic-potential form
     - Good behavior at high compression
     - Lower-order truncation has limited theoretical basis

Volume--temperature equations of state
--------------------------------------

At fixed pressure, the volume thermal-expansion coefficient is

.. math::

   \alpha(T)=\frac{1}{V}
   \left(\frac{\partial V}{\partial T}\right)_P.

Integration gives

.. math::

   V(T)=V_{\mathrm{ref}}
   \exp\left[
   \int_{T_{\mathrm{ref}}}^{T}\alpha(t)\,\mathrm dt
   \right].

Thermodynamics requires :math:`\alpha\rightarrow0` and
:math:`\partial\alpha/\partial T\rightarrow0` as :math:`T\rightarrow0` for a
stable non-degenerate solid.  Simple empirical expressions often describe a
restricted high-temperature interval very well without satisfying this
low-temperature limit.  Conversely, forms designed for low-temperature
saturation may become unsuitable when extrapolated far above their intended
range.

Berman equation
~~~~~~~~~~~~~~~

Let :math:`\Delta T=T-T_{\mathrm{ref}}`.  The Berman form used by Quantas is

.. math::

   V(T)=V_{\mathrm{ref}}
   \left[1+\alpha_0\Delta T
   +\frac{1}{2}\alpha_1\Delta T^2\right].

The exact expansion coefficient of this truncated volume expression is

.. math::

   \alpha(T)=
   \frac{\alpha_0+\alpha_1\Delta T}
   {1+\alpha_0\Delta T+\frac{1}{2}\alpha_1\Delta T^2}.

**Advantages**

- It is simple, transparent, and effective for many moderate- or
  high-temperature datasets.
- The quadratic term captures smoothly varying non-linear expansion.
- :math:`V(T_{\mathrm{ref}})=V_{\mathrm{ref}}` and
  :math:`\alpha(T_{\mathrm{ref}})=\alpha_0` exactly.

**Limitations**

- It does not enforce the third-law low-temperature limit.
- :math:`\alpha_1` is a coefficient of the volume polynomial, not exactly
  :math:`(\partial\alpha/\partial T)_{T_{\mathrm{ref}}}`.
- The bracketed volume factor must remain positive; broad extrapolation can
  become non-physical.

Fei equation
~~~~~~~~~~~~

The Fei expansion coefficient is

.. math::

   \alpha(T)=\alpha_0+\alpha_1T+\alpha_2T^{-2}.

Integration yields

.. math::

   V(T)=V_{\mathrm{ref}}\exp\left[
   \alpha_0(T-T_{\mathrm{ref}})
   +\frac{1}{2}\alpha_1(T^2-T_{\mathrm{ref}}^2)
   -\alpha_2
   \left(\frac{1}{T}-\frac{1}{T_{\mathrm{ref}}}\right)
   \right].

The simplified form sets :math:`\alpha_2=0`.

**Advantages**

- The fitted coefficients describe the expansion function directly.
- They are defined against absolute temperature and are independent of the
  chosen reference temperature.
- The linear form :math:`\alpha_0+\alpha_1T` is useful for broad
  high-temperature trends.

**Limitations**

- The inverse-square term diverges as :math:`T\rightarrow0`.
- Even the simplified form does not force :math:`\alpha(0)=0`.
- It should therefore be regarded as a finite-temperature representation, not a
  universal low-temperature law.

Modified Holland--Powell equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The modified Pawley--Redfern--Holland expression adopted by EosFit7 and
Quantas is

.. math::

   V(T)=V_{\mathrm{ref}}
   \left\{
   1+\alpha_0(T-T_{\mathrm{ref}})
   -2(10\alpha_0+\alpha_1)
   \left(\sqrt{T}-\sqrt{T_{\mathrm{ref}}}\right)
   \right\}.

Setting the residual :math:`\alpha_1=0` gives the simplified
Holland--Powell form.

**Advantages**

- It is compact and can describe the approach toward approximately constant
  expansion at high temperature.
- The simplified form is useful for lower-resolution datasets with few
  independently resolvable thermal parameters.
- It is widely used in thermodynamic databases.

**Limitations**

- It is not a low-temperature EOS.  Below

  .. math::

     T_{\mathrm{limit}}=
     \left(\frac{10\alpha_0+\alpha_1}{\alpha_0}\right)^2,

  the predicted expansion becomes negative for the usual positive
  :math:`\alpha_0` case.  The simplified form has
  :math:`T_{\mathrm{limit}}=100\ \mathrm K`.
- The coefficient :math:`\alpha_0` is not generally the exact expansion at
  :math:`T_{\mathrm{ref}}`.

Salje equation
~~~~~~~~~~~~~~

The Salje form explicitly represents low-temperature saturation:

.. math::

   V(T)=\left[
   V_0^{1/3}
   +p_1\theta_{\mathrm{sat}}
   \left(
   \coth\frac{\theta_{\mathrm{sat}}}{T}-1
   \right)
   \right]^3.

It satisfies

.. math::

   V(0)=V_0,
   \qquad
   \alpha(0)=0.

**Advantages**

- It has the correct zero-temperature saturation of thermal expansion.
- It is particularly suitable for low-temperature datasets.
- The saturation scale :math:`\theta_{\mathrm{sat}}` has a clear qualitative
  role.

**Limitations**

- Above a few times :math:`\theta_{\mathrm{sat}}`, the predicted expansion
  becomes nearly temperature independent.
- That behavior is inconsistent with the approximately linear increase of
  :math:`\alpha(T)` observed at high temperature for many solids.
- It should not be used as an unrestricted high-temperature extrapolation.

Kroll form of Holland--Powell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Kroll--Holland--Powell form relates thermal expansion to an Einstein
oscillator while retaining parameters at a chosen reference temperature.  Let

.. math::

   \xi_0=
   \frac{(\theta_E/T_{\mathrm{ref}})^2
   \exp(\theta_E/T_{\mathrm{ref}})}
   {[\exp(\theta_E/T_{\mathrm{ref}})-1]^2},

.. math::

   A(T)=\alpha_{\mathrm{ref}}\frac{\theta_E}{\xi_0}
   \left[
   \frac{1}{\exp(\theta_E/T)-1}
   -\frac{1}{\exp(\theta_E/T_{\mathrm{ref}})-1}
   \right],

and, with :math:`k=K'_0`,

.. math::

   B=-\frac{1}{k(k+2)}.

The volume can be written

.. math::

   V(T)=V_{\mathrm{ref}}
   \left\{
   -k+(1+k)
   \left[
   1-\frac{k(k+2)}{k+1}A(T)
   \right]^B
   \right\}.

**Advantages**

- The Einstein function gives low-temperature saturation while retaining a
  realistic high-temperature trend.
- The reference identities
  :math:`V(T_{\mathrm{ref}})=V_{\mathrm{ref}}` and
  :math:`\alpha(T_{\mathrm{ref}})=\alpha_{\mathrm{ref}}` are exact.
- It offers a physically motivated bridge between purely empirical thermal
  expansion and lattice-vibrational models.

**Limitations**

- It is more complex than the simple polynomial forms.
- :math:`\theta_E` is frequently weakly constrained by ordinary V--T data.
- Its robustness depends on the availability of a credible compressional
  :math:`K'_0` and on remaining within a physically admissible temperature
  interval.

Choosing a V--T family
~~~~~~~~~~~~~~~~~~~~~~

The temperature interval is often the decisive criterion.

.. list-table:: V--T model summary
   :header-rows: 1
   :widths: 21 24 28 27

   * - Family
     - Intended behavior
     - Principal strength
     - Principal caution
   * - Berman
     - Smooth polynomial expansion
     - Simple and effective over restricted ranges
     - No low-temperature saturation
   * - Fei
     - Direct polynomial/inverse-square :math:`\alpha(T)`
     - Reference-independent expansion coefficients
     - Divergent or non-zero low-temperature limit
   * - Modified Holland--Powell
     - Approximately constant high-temperature expansion
     - Compact database-oriented form
     - Negative expansion below its limiting temperature
   * - Salje
     - Low-temperature saturation
     - Correct :math:`\alpha(0)=0`
     - Unrealistic unrestricted high-temperature behavior
   * - Kroll--Holland--Powell
     - Einstein-based low/high-temperature crossover
     - More physical global temperature dependence
     - Extra parameters may be weakly constrained

Pressure--volume--temperature equations of state
------------------------------------------------

A P--V--T EOS combines an isothermal compressional EOS with a thermal model.
The central thermodynamic issue is not only the increase of the zero-pressure
volume :math:`V_{0T}`, but also the change of the zero-pressure bulk modulus
:math:`K_{0T}` with temperature.  Angel *et al.* discuss three principal
coupling strategies: a linear temperature dependence of :math:`K_0`, an
Anderson--Gruneisen relation, and thermal pressure.

Linear variation of the zero-pressure bulk modulus
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest model is

.. math::

   K_{0T}=K_0+
   \left(\frac{\partial K_0}{\partial T}\right)
   (T-T_{\mathrm{ref}}),

combined with a selected V--T equation for :math:`V_{0T}` and a selected P--V
family evaluated with those temperature-dependent reference values.

**Advantages**

- It is transparent and introduces only one direct temperature derivative.
- In combination with non-linear thermal expansion and an isothermal EOS with
  non-zero :math:`K''_0`, it contains the principal second derivatives of
  volume with respect to pressure and temperature.
- It is often adequate over a restricted experimental temperature range.

**Limitations**

- A constant :math:`\partial K_0/\partial T` is a local approximation.
- Extrapolation can drive :math:`K_{0T}` toward zero or negative values.
- Hellfrich and Connolly showed that this coupling can predict non-physical
  negative thermal expansion at moderate pressure for many materials.

Anderson--Gruneisen coupling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Anderson--Gruneisen parameter :math:`\delta` relates the temperature
variation of the zero-pressure bulk modulus to thermal expansion
[#anderson_1995]_:

.. math::

   K_{0T}=K_0
   \exp\left[
   -\delta
   \int_{T_{\mathrm{ref}}}^{T}\alpha(t)\,\mathrm dt
   \right].

Using the zero-pressure thermal volume,

.. math::

   K_{0T}=K_0
   \left(\frac{V_0}{V_{0T}}\right)^{\delta}.

**Advantages**

- It couples thermal softening directly to expansion instead of imposing a
  constant temperature derivative.
- It generally behaves more smoothly under pressure and avoids some of the
  negative-expansion pathologies of linear coupling.
- The approximation :math:`\delta\simeq K'_0` can provide a first estimate
  when independent information is limited.

**Limitations**

- Treating :math:`\delta` as constant is itself an approximation.
- :math:`\delta\simeq K'_0` is not a thermodynamic identity.
- Broad extrapolation can fail when the pressure or temperature dependence of
  :math:`\delta` becomes important.

Holland--Powell Einstein thermal pressure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The thermal-pressure approach writes

.. math::

   P(V,T)=P(V,T_{\mathrm{ref}})+\Delta P_{\mathrm{th}}(T).

The reference term is any selected isothermal P--V EOS.  In the
Holland--Powell Einstein model,

.. math::

   \Delta P_{\mathrm{th}}(T)=
   \alpha_{\mathrm{ref}}K_0
   \frac{\theta_E}{\xi_0}
   \left[
   \frac{1}{\exp(\theta_E/T)-1}
   -\frac{1}{\exp(\theta_E/T_{\mathrm{ref}})-1}
   \right].

The thermal pressure vanishes on the reference isotherm.

**Advantages**

- It gives a direct physical interpretation: heating at fixed volume creates an
  additional pressure.
- The Einstein function produces low-temperature saturation and an
  approximately linear high-temperature trend.
- The formulation avoids several pathologies associated with prescribing a
  constant :math:`\partial K_0/\partial T`.

**Limitations**

- In this form, thermal pressure depends only on temperature, not explicitly on
  volume.
- The single Einstein frequency is an effective representation of the
  vibrational spectrum.
- :math:`\theta_E` may be poorly determined unless the temperature range is
  broad and precise.

Mie--Gruneisen--Debye thermal pressure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Quantas also includes a Mie--Gruneisen--Debye (MGD) coupling as an extension
beyond the model set reviewed in detail by Angel *et al.* (2014).  MGD is a
quasi-harmonic thermal-pressure EOS in which the vibrational internal energy is
represented by a Debye spectrum [#stixrude_lithgow_bertelloni_2005]_.  In
compact form,

.. math::

   P(V,T)=P(V,T_{\mathrm{ref}})
   +\frac{\gamma(V)}{V}
   \left[U_D(V,T)-U_D(V,T_{\mathrm{ref}})\right],

with the appropriate unit and atom-count normalization.

The full Quantas form uses

.. math::

   \gamma(V)=\gamma_0\left(\frac{V}{V_0}\right)^q,

and

.. math::

   \gamma(V)=-\frac{\partial\ln\theta_D}{\partial\ln V}.

For :math:`q\ne0`, integration gives

.. math::

   \theta_D(V)=\theta_{D0}
   \exp\left[
   \frac{\gamma_0-\gamma(V)}{q}
   \right].

The Debye thermal term uses

.. math::

   D_3(x)=\frac{3}{x^3}
   \int_0^x\frac{t^3}{\exp(t)-1}\,\mathrm dt.

Quantas also represents the explicitly named EosFit ``q-compromise``
approximation, in which :math:`\theta_D` and :math:`\gamma/V` are held
constant.

**Advantages**

- It introduces explicit volume dependence into the thermal pressure.
- It connects compression, the Gruneisen parameter, and lattice-vibrational
  energy within a quasi-harmonic framework.
- It is well suited to broad P--V--T applications when the material is
  reasonably represented by an effective Debye spectrum.

**Limitations**

- A single Debye temperature cannot reproduce all details of a real phonon
  density of states, especially at low temperature.
- :math:`\theta_{D0}`, :math:`\gamma_0`, and :math:`q` can be strongly
  correlated unless the dataset spans a broad P--T--V domain.
- Intrinsic anharmonicity at fixed volume, electronic excitations, magnetic
  effects, and phase transitions are outside the basic MGD model.
- The physical normalization must be consistent with the chosen cell or molar
  volume.

Choosing a P--V--T coupling
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: P--V--T coupling summary
   :header-rows: 1
   :widths: 22 24 27 27

   * - Coupling
     - Main assumption
     - Principal strength
     - Principal caution
   * - Linear :math:`K_0(T)`
     - Constant :math:`\partial K_0/\partial T`
     - Minimal and transparent
     - Local model; can become non-physical on extrapolation
   * - Anderson--Gruneisen
     - :math:`K_0` scales with thermal volume
     - Smooth pressure--temperature coupling
     - Constant :math:`\delta` is approximate
   * - Einstein thermal pressure
     - Temperature-only oscillator pressure
     - Physical low/high-temperature crossover
     - No explicit volume dependence of thermal pressure
   * - Mie--Gruneisen--Debye
     - Debye energy with volume-dependent :math:`\gamma` and :math:`\theta_D`
     - Broad quasi-harmonic P--V--T description
     - Effective-spectrum approximation and correlated parameters

Linear equations of state
-------------------------

The same EOS families can describe a lattice parameter or another linear
quantity :math:`x`, provided the finite-strain conventions remain consistent
with the volume EOS.  Following Angel *et al.*, Quantas uses the auxiliary
quantity

.. math::

   q=x^3

inside the volume-like equation.  The physical linear expansion and
compressibility are then

.. math::

   \alpha_x=\frac{1}{x}
   \left(\frac{\partial x}{\partial T}\right)_P,
   \qquad
   \beta_x=-\frac{1}{x}
   \left(\frac{\partial x}{\partial P}\right)_T.

For three orthogonal principal directions,

.. math::

   \alpha_1+\alpha_2+\alpha_3=\alpha_V,
   \qquad
   \beta_1+\beta_2+\beta_3=\frac{1}{K}.

A linear modulus is

.. math::

   M_x=\frac{1}{\beta_x}.

For an isotropic or cubic material, :math:`M_x=3K`.  The cubed-length
construction therefore preserves consistency between volume and linear
elasticity, but the numerical moduli reported for a linear EOS are physical
linear moduli, not volume-like auxiliary values.

Scientific interpretation and model selection
----------------------------------------------

An EOS is not only a curve through data.  Its derivatives determine bulk
modulus, expansion, and pressure--temperature response, so model choice should
be based on physical scope as well as residual size.

When comparing equations, consider:

- whether the sampled compression or temperature interval matches the intended
  domain of the model;
- whether an additional EOS order is supported by independent curvature rather
  than parameter correlation;
- whether the reference state lies inside or near the observed domain;
- whether extrapolated volumes, moduli, expansion coefficients, Debye
  temperatures, and thermal pressures remain physically admissible;
- whether conclusions are stable when another scientifically plausible EOS
  family is used.

A statistically converged representation can still be scientifically
inappropriate outside its calibration range.  Phase transitions, elastic or
vibrational instabilities, intrinsic anharmonicity, electronic transitions,
and changes in bonding cannot in general be repaired by increasing the order
of a single-phase EOS.

What Quantas provides
---------------------

Within the standalone EOS workflow, Quantas evaluates and analyzes:

- Murnaghan P--V;
- Birch--Murnaghan P--V, orders 2--4;
- natural-strain/Poirier--Tarantola P--V, orders 2--4;
- Vinet P--V, orders 2--3;
- Tait P--V, orders 2--4;
- Berman, Fei, modified Holland--Powell, Salje, and
  Kroll--Holland--Powell V--T models;
- linear, Anderson--Gruneisen, Einstein thermal-pressure, and MGD P--V--T
  couplings;
- volume and linear forms where the underlying scientific convention is
  applicable.

The workflow pages explain how these models are selected and assessed; the
present chapter defines the physical assumptions behind them.

.. include:: ../_generated/references/eos.inc
