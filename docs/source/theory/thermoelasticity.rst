Thermoelasticity
================

Thermoelasticity describes the coupled dependence of the elastic response of a
solid on pressure, temperature, and strain.  At a fixed thermodynamic state,
linear elasticity relates an infinitesimal stress increment to an infinitesimal
strain increment.  Thermoelasticity adds the question of how the corresponding
elastic coefficients evolve as the reference state itself changes.

A complete treatment must be thermodynamically consistent: the equation of
state, thermal expansion, heat capacities, elastic constants, and acoustic
velocities should be obtained from derivatives of the same thermodynamic
potential.  The formulation developed by L. Stixrude and C. Lithgow-Bertelloni
provides a particularly useful framework because it generalizes the fundamental
thermodynamic relation to anisotropic strain and combines Eulerian finite-strain
elasticity with quasi-harmonic lattice thermodynamics [#stixrude_lithgow_bertelloni_2005]_.

Quantas currently implements the **quasi-static approximation** (QSA) using
a **cold Eulerian finite-strain** model for the elastic coefficients.  The
complete thermoelastic framework is introduced first; the approximation used
by Quantas is then identified explicitly.

Thermodynamic foundation
------------------------

For hydrostatic thermodynamics, the Gibbs and Helmholtz free energies are
related by the Legendre transformation

.. math::

   G(P,T)=F(V,T)+PV.

A crystal can support deviatoric stress, so pressure and volume must be
replaced by tensorial stress and strain variables.  The stress tensor can be
written as

.. math::

   \sigma_{ij}=-P\delta_{ij}+\tau_{ij},

where :math:`P` is the hydrostatic pressure and :math:`\tau_{ij}` is the
deviatoric stress.

Thermodynamic properties follow from derivatives of a potential expressed in
its natural variables.  In particular, the isothermal compliance tensor is
obtained from the stress derivatives of the Gibbs free energy,

.. math::

   s^T_{ijkl}
   =-\frac{1}{V}
   \left(\frac{\partial^2 G}
   {\partial\sigma_{ij}\,\partial\sigma_{kl}}\right)_T,

and the stiffness tensor is its inverse,

.. math::

   \mathbf C^T=(\mathbf S^T)^{-1}.

This construction places elasticity on the same thermodynamic footing as
entropy, volume, heat capacity, bulk modulus, and thermal expansion.

Elastic coefficients under hydrostatic pre-stress
--------------------------------------------------

At finite pressure, the second derivative of the Helmholtz free energy is not,
by itself, the incremental stress--strain coefficient governing wave
propagation or mechanical response around the pre-stressed state.  For an
Eulerian strain increment :math:`S_{ij}`, the relevant Wallace stress--strain
coefficients are [#wallace_1972]_ [#barron_klein_1965]_

.. math::

   c^T_{ijkl}
   =\frac{1}{V}
   \left(\frac{\partial^2G}
   {\partial S_{ij}\,\partial S_{kl}}\right)_{P,T}
   =\frac{1}{V}
   \left(\frac{\partial^2F}
   {\partial S_{ij}\,\partial S_{kl}}\right)_T
   -P\,\delta^{ij}_{kl},

with

.. math::

   \delta^{ij}_{kl}
   =-\delta_{ij}\delta_{kl}
    -\delta_{il}\delta_{jk}
    -\delta_{jl}\delta_{ik}.

The pressure term is therefore part of the thermodynamic definition of the
finite-pressure stress--strain coefficients.  This distinction matters when
elastic tensors are used for stability analysis or acoustic-wave propagation.

Static and vibrational parts of the free energy
-----------------------------------------------

Stixrude and Lithgow-Bertelloni write the Helmholtz free energy as the sum of a
reference value, a cold compression term, and a quasi-harmonic vibrational
term,

.. math::

   F(E_{ij},T)
   =F_0+F_{\mathrm c}(E_{ij},T_0)
   +F_{\mathrm q}(E_{ij},T)-F_{\mathrm q}(E_{ij},T_0).

Here :math:`E_{ij}` is the Eulerian finite-strain tensor.  The cold part is
expanded as a power series in strain, while the quasi-harmonic part contains
the vibrational free energy with phonon frequencies that depend on strain.
This decomposition makes the origin of the elastic response explicit:

.. math::

   c^T_{ijkl}(P,T)
   =c^{\mathrm c}_{ijkl}(P,T)
   +\Delta c^{\mathrm q}_{ijkl}(P,T).

The first contribution describes the static response of the electronic ground
state to compression or expansion.  The second contains the explicit
vibrational response to strain.

Eulerian finite strain
----------------------

For an isotropically compressed reference state,

.. math::

   E_{ij}=-f\delta_{ij},

where the scalar Eulerian finite strain is

.. math::

   f=\frac{1}{2}
   \left[\left(\frac{V_0}{V}\right)^{2/3}-1\right].

This convention gives :math:`f=0` at the reference volume, positive strain on
compression, and negative strain on expansion.

Eulerian finite strain is especially useful at high pressure because the
resulting series generally converges more rapidly than a direct expansion in
volume or a Lagrangian-strain expansion.  In the benchmark examined by
Stixrude and Lithgow-Bertelloni for MgSiO\ :sub:`3` perovskite, the
third-order Eulerian expression reproduced first-principles bulk and shear
moduli much more accurately than the corresponding Lagrangian expression.  The
comparison also showed that the quadratic term in :math:`f` can be material at
large compression; omitting it changed the predicted shear modulus by several
percent in that example.

Cold finite-strain elastic coefficients
---------------------------------------

Let

- :math:`V_0` be the reference volume;
- :math:`K_0` be the reference isothermal bulk modulus;
- :math:`K'_0=(\partial K/\partial P)_0`;
- :math:`C^0_{IJ}` be one reference Wallace stiffness component;
- :math:`C'^0_{IJ}=(\partial C_{IJ}/\partial P)_0`;
- :math:`\Delta_{IJ}` be the corresponding Voigt component of
  :math:`\delta^{ij}_{kl}`.

The thermodynamically consistent third-order Eulerian finite-strain expression
for the cold stiffness component is

.. math::

   C^{\mathrm c}_{IJ}(f)
   =(1+2f)^{5/2}
   \left\{
   C^0_{IJ}
   +\left(3K_0C'^0_{IJ}-5C^0_{IJ}\right)f
   +\left[
   6K_0C'^0_{IJ}-14C^0_{IJ}
   -\frac{3}{2}K_0\Delta_{IJ}(3K'_0-16)
   \right]f^2
   \right\}.

The expression has several useful properties:

- it recovers :math:`C^0_{IJ}` exactly at :math:`f=0`;
- its initial pressure dependence is controlled by :math:`C'^0_{IJ}`;
- all tensor components share the same reference equation-of-state parameters;
- the hydrostatic pre-stress tensor enters the quadratic coefficient;
- truncating after the linear term gives the second-order form, while retaining
  :math:`f^2` gives the third-order form used for the thermodynamically
  consistent cold expansion.

The order refers to the truncation of the underlying free-energy expansion.
Because elastic constants are strain derivatives, a third-order finite-strain
model for the moduli contains a quadratic term in :math:`f`.

Explicit quasi-harmonic elastic contribution
--------------------------------------------

In the complete quasi-harmonic treatment, vibrational frequencies depend not
only on volume but on the full strain tensor.  Their strain sensitivity is
expressed through the generalized Grüneisen tensor

.. math::

   \gamma_{ij}
   =-\frac{1}{\nu_\lambda}
   \frac{\partial\nu_\lambda}{\partial S_{ij}},

and its strain derivative

.. math::

   \eta_{ijkl}
   =\frac{\partial\gamma_{ij}}{\partial S_{kl}}.

Under the Grüneisen approximations used in the finite-strain formulation, the
explicit vibrational correction to the elastic tensor has the form

.. math::

   \Delta c^{\mathrm q}_{ijkl}
   =\left[
   \gamma_{ij}\gamma_{kl}
   +\frac{1}{2}
   (\gamma_{ij}\delta_{kl}+\gamma_{kl}\delta_{ij})
   -\eta_{ijkl}
   \right]\rho\,\Delta U_{\mathrm q}
   -\gamma_{ij}\gamma_{kl}\rho\,\Delta(C_VT).

This equation shows why a full thermoelastic QHA calculation is more demanding
than an ordinary volume QHA calculation.  It requires derivatives of the
phonon frequencies with respect to independent strains, not only their
variation along a hydrostatic volume path.

The quasi-static approximation
------------------------------

The quasi-static approximation, as used in practical quasi-harmonic
thermoelastic workflows [#destefanis_ravoux_cossard_erba_2019]_, neglects
the explicit vibrational strain
contribution :math:`\Delta c^{\mathrm q}_{ijkl}`.  Temperature still changes
the elastic tensor, but only indirectly, because thermal vibrations alter the
QHA equilibrium volume.  The isothermal tensor is therefore approximated as

.. math::

   C^{T,\mathrm{QSA}}_{IJ}(P,T)
   =C^{\mathrm c}_{IJ}
   \left[V_{\mathrm{QHA}}(P,T)\right].

This approximation retains two distinct pieces of information:

1. a static finite-strain relation :math:`C^{\mathrm c}_{IJ}(V)`;
2. a quasi-harmonic equilibrium-volume surface :math:`V(P,T)`.

The QHA volume may include zero-point and thermal vibrational effects on the
equilibrium structure.  What is omitted is their **explicit anisotropic strain
derivative** in the elastic tensor.

The QSA is consequently more informative than applying an arbitrary linear
temperature correction to each elastic coefficient.  Pressure and temperature
remain coupled through one physical state variable, the equilibrium volume,
and the cold elastic coefficients retain the finite-strain structure required
by the reference equation of state.  However, the approximation cannot capture
materials for which explicit vibrational shear-strain coupling is large.

Isothermal and adiabatic elastic tensors
----------------------------------------

The cold finite-strain reconstruction gives an isothermal stiffness tensor.
Many acoustic and seismic measurements are closer to adiabatic conditions.  A
thermodynamic conversion [#davies_1974]_ [#waters_bielawski_2016]_ can be
written in Voigt notation as

.. math::

   C^S_{IJ}
   =C^T_{IJ}
   +\frac{TV}{C_V}\lambda_I\lambda_J,
   \qquad
   \lambda_I=\sum_K C^T_{IK}\alpha_K,

where :math:`\alpha_K` is the Voigt representation of the thermal-expansion
tensor.  The correction is positive semidefinite when :math:`T`, :math:`V`,
and :math:`C_V` are positive, so adiabatic stiffness is not smaller than the
corresponding isothermal stiffness in the associated quadratic form.

This conversion does not restore the explicit vibrational strain derivatives
omitted by the QSA.  It converts the thermodynamic boundary condition of an
already reconstructed tensor from isothermal to adiabatic.

Approximation implemented by Quantas
------------------------------------

The present Quantas thermoelastic workflow implements the following scientific
model:

#. determine a common cold reference :math:`V_0`, :math:`K_0`, and
   :math:`K'_0` from the static electronic energy--volume relation;
#. determine :math:`C^0_{IJ}` and :math:`C'^0_{IJ}` from static elastic tensors
   sampled at several hydrostatic volumes;
#. represent every independent stiffness component with the second- or
   third-order Eulerian cold finite-strain expression;
#. evaluate that relation at the QHA equilibrium volume
   :math:`V_{\mathrm{QHA}}(P,T)`;
#. optionally convert the reconstructed isothermal tensor to the adiabatic
   tensor using QHA heat capacity and thermal expansion.

Quantas does **not** presently evaluate the generalized Grüneisen tensor
:math:`\gamma_{ij}`, its strain derivative :math:`\eta_{ijkl}`, or the explicit
quasi-harmonic elastic contribution.  The result is therefore a quasi-static
thermoelastic tensor, not a complete strain-dependent QHA elastic tensor.

Scientific scope and limitations
--------------------------------

The QSA is most defensible when:

- the static elastic-volume relation is smooth and well constrained;
- the QHA equilibrium state remains within, or close to, the calibrated volume
  interval;
- thermal evolution is dominated by expansion along the hydrostatic volume
  path;
- explicit vibrational shear-strain coupling is modest;
- the phase remains mechanically and thermodynamically meaningful over the
  investigated domain.

The approximation does not describe intrinsic anharmonicity, phonon lifetimes,
magnetic or electronic thermal excitations, phase transformations, or explicit
strain derivatives of the phonon spectrum.  Agreement between a QSA result and
experiment is therefore material- and state-dependent and should be evaluated
rather than assumed.

The Eulerian finite-strain formulation is also a truncated series.  Its
third-order form is generally preferable over a broad compression interval,
but higher-order behavior may become important for some materials or at very
large strain.  Extrapolation outside the sampled elastic-volume and QHA domains
must be treated as a scientific prediction with correspondingly greater
uncertainty.

For operational details, input requirements, diagnostics, and supported
analyses, see :doc:`../workflows/thermoelasticity`.

.. include:: ../_generated/references/thermoelasticity.inc
