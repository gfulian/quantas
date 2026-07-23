Quasi-Harmonic Approximation
============================

The quasi-harmonic approximation (QHA) extends harmonic lattice dynamics by
allowing the phonon frequencies to depend on crystal volume.  At each fixed
volume the crystal is still treated harmonically, but the collection of
harmonic calculations defines a temperature-dependent free-energy surface from
which equilibrium properties can be obtained at finite pressure and
temperature [#anderson_1995]_ [#anderson_masuda_isaak_1995]_.

The treatment in this section describes the scientific model.  The numerical
representations, fit diagnostics, and execution choices available in Quantas
are documented in :doc:`../workflows/qha`.

Why the harmonic approximation is insufficient
-----------------------------------------------

In the harmonic approximation the vibrational frequencies are evaluated for a
fixed structure.  The equilibrium volume therefore does not change with
temperature, which implies

.. math::

   \alpha_V = 0,
   \qquad C_P = C_V.

Real crystals expand or contract because their vibrational spectrum changes
with volume.  The QHA introduces this effect while retaining a harmonic
phonon description at each sampled volume.

Quasi-harmonic free-energy surface
----------------------------------

Let :math:`U_0(V)` be the static electronic energy and let
:math:`F_{\mathrm{vib}}(V,T)` be the harmonic vibrational Helmholtz free energy
calculated from the volume-dependent phonon frequencies.  The total Helmholtz
free energy is

.. math::

   F(V,T) = U_0(V) + F_{\mathrm{vib}}(V,T).

The volume dependence of :math:`F_{\mathrm{vib}}` is the essential new element
of the QHA.  It describes the leading thermodynamic consequence of vibrational
anharmonicity without introducing explicit cubic or higher-order force
constants [#erba_2014]_.

Equilibrium at finite pressure
------------------------------

At an external hydrostatic pressure :math:`P`, the equilibrium state minimizes

.. math::

   G^*(V;P,T) = F(V,T) + PV.

The equilibrium volume :math:`V(P,T)` satisfies

.. math::

   \left(\frac{\partial G^*}{\partial V}\right)_{P,T} = 0,

or equivalently

.. math::

   P = -\left(\frac{\partial F}{\partial V}\right)_T.

Once the equilibrium volume has been determined, harmonic thermodynamic
properties are evaluated or interpolated at that volume to construct the
pressure--temperature state.

Thermodynamic potentials
------------------------

At the equilibrium volume, the enthalpy and Gibbs free energy are

.. math::

   H(P,T) = U(P,T) + PV,

and

.. math::

   G(P,T) = F(P,T) + PV.

The Gibbs free energy is the appropriate thermodynamic potential for comparing
phases at fixed pressure and temperature, provided all phases are represented
with mutually consistent electronic and vibrational free energies.

Isothermal bulk modulus
-----------------------

The isothermal bulk modulus is obtained from the curvature of the Helmholtz
free-energy surface,

.. math::

   K_T(V,T)
   = V\left(\frac{\partial^2 F}{\partial V^2}\right)_T
   = -V\left(\frac{\partial P}{\partial V}\right)_T.

At equilibrium it becomes a function of pressure and temperature,
:math:`K_T(P,T)`.  Its pressure derivative is

.. math::

   K'_T = \left(\frac{\partial K_T}{\partial P}\right)_T.

A positive :math:`K_T` is required for local mechanical stability against
hydrostatic volume fluctuations.

Thermal expansion
-----------------

The volumetric thermal-expansion coefficient is

.. math::

   \alpha_V(P,T)
   = \frac{1}{V}
     \left(\frac{\partial V}{\partial T}\right)_P.

Using

.. math::

   P(V,T) = -\left(\frac{\partial F}{\partial V}\right)_T,

and differentiating the constant-pressure condition gives the mixed-derivative
form

.. math::

   \alpha_V
   = -\frac{1}{K_T}
     \frac{\partial^2 F}{\partial V\,\partial T}.

Since

.. math::

   S = -\left(\frac{\partial F}{\partial T}\right)_V,

this can also be written as the Maxwell-relation form

.. math::

   \alpha_V
   = \frac{1}{K_T}
     \left(\frac{\partial S}{\partial V}\right)_T.

This relation is especially useful because it connects thermal expansion to the
volume dependence of the vibrational entropy.

Mode Grüneisen parameters
-------------------------

The volume sensitivity of a phonon mode is described by its mode Grüneisen
parameter,

.. math::

   \gamma_{qj}(V)
   = -\frac{\partial\ln\nu_{qj}(V)}
           {\partial\ln V}.

A positive mode Grüneisen parameter means that the mode softens upon expansion;
a negative value means that the mode stiffens upon expansion.

A heat-capacity-weighted average is

.. math::

   \bar{\gamma}(V,T)
   = \frac{
       \sum_{qj} w_q\,\gamma_{qj}(V) C_{V,qj}(V,T)
     }{
       \sum_{qj} w_q\,C_{V,qj}(V,T)
     }.

When the quasi-harmonic relations are internally consistent, the volumetric
thermal expansion can be expressed as

.. math::

   \alpha_V
   = \frac{\bar{\gamma} C_V}{K_T V},

with all quantities expressed on the same normalization basis.  Conversely,
the macroscopic thermodynamic Grüneisen parameter is

.. math::

   \gamma
   = \frac{\alpha_V K_T V}{C_V}.

The mode-resolved route requires physically continuous phonon branches across
the sampled volumes.  The macroscopic thermodynamic relation does not require
mode-by-mode tracking, but it remains sensitive to the quality of
:math:`\alpha_V`, :math:`K_T`, and :math:`C_V`.

Isochoric and isobaric heat capacities
--------------------------------------

The isochoric heat capacity is obtained from harmonic thermodynamics evaluated
at the equilibrium volume.  The standard thermodynamic relation between the
isobaric and isochoric heat capacities is

.. math::

   C_P-C_V = \alpha_V^2 K_T V T.

Therefore

.. math::

   C_P = C_V + \alpha_V^2 K_T V T.

The same normalization and unit basis must be used for :math:`C_V`,
:math:`K_T`, and :math:`V`.

Adiabatic bulk modulus
----------------------

The adiabatic and isothermal bulk moduli are related by

.. math::

   \frac{K_S}{K_T} = \frac{C_P}{C_V},

so that

.. math::

   K_S = K_T\frac{C_P}{C_V}.

At zero temperature, or whenever the thermal correction vanishes,
:math:`K_S` approaches :math:`K_T`.

Structural properties
---------------------

When a consistent set of relaxed crystal structures is available along the
sampled volume path, the equilibrium cell can be reconstructed at
:math:`V(P,T)`.  The linear expansion of a lattice parameter :math:`a` can be
written through the chain rule as

.. math::

   \alpha_a
   = \frac{1}{a}\left(\frac{\partial a}{\partial T}\right)_P
   = \frac{\partial\ln a}{\partial\ln V}\,\alpha_V.

Equivalent relations apply to the other lattice parameters and to the
symmetric Cartesian thermal-expansion tensor.  The trace of that tensor must be
consistent with the volumetric expansion coefficient,

.. math::

   \operatorname{tr}(\boldsymbol{\alpha}) = \alpha_V.

This structural reconstruction represents the cell-shape response along the
sampled volume-constrained structural path.  It is not a fully independent
anisotropic minimization of every strain degree of freedom at every
pressure--temperature point.

Two equivalent QHA viewpoints
-----------------------------

The QHA can be organized in two scientifically equivalent ways when the
underlying fits are accurate.

Frequency-based formulation
   The phonon frequencies are represented as functions of volume.  Harmonic
   thermodynamic functions are then recalculated at the equilibrium volume.
   This formulation retains mode-resolved information and permits direct
   calculation of mode Grüneisen parameters, but requires mode continuity.

Thermodynamic formulation
   Harmonic thermodynamic quantities such as :math:`F_{\mathrm{vib}}`,
   :math:`S`, and :math:`C_V` are first calculated at the sampled volumes and
   then represented as functions of volume.  This avoids mode-by-mode
   correspondence, but does not retain the same mode-resolved information.

Quantas supports both formulations.  Their numerical configuration belongs to
the workflow layer rather than to the scientific definition of the QHA.

Scientific scope and limitations
--------------------------------

The QHA captures the leading thermal effect arising from the volume dependence
of harmonic phonons.  It does **not** include explicit anharmonic interactions
at fixed volume.  Important limitations therefore include:

- phonon linewidths and lifetimes are not described;
- intrinsic temperature shifts of frequencies at fixed volume are neglected;
- strongly anharmonic crystals may not be represented accurately;
- dynamically unstable modes and soft-mode phase transitions require special
  treatment;
- results become unreliable when the equilibrium state lies far outside the
  sampled volume interval;
- electronic, magnetic, configurational, or defect contributions must be added
  separately when they are thermodynamically relevant.

The QHA is generally most reliable for mechanically stable crystalline phases
away from melting, strongly anharmonic regimes, and structural phase
transitions.  Its quality must be assessed from the phonon calculations, the
sampled volume range, the free-energy representation, and comparison with
experimental or independent theoretical data [#erba_shahrokhi_moradian_dovesi_2015]_.

Quantities derived by Quantas
-----------------------------

Depending on the supplied data and selected analysis, the QHA workflow can
provide:

- equilibrium volume :math:`V(P,T)`;
- zero-point, thermal, internal, Helmholtz, enthalpy, and Gibbs energies;
- entropy :math:`S(P,T)`;
- :math:`C_V(P,T)` and :math:`C_P(P,T)`;
- :math:`K_T(P,T)`, :math:`K_S(P,T)`, and :math:`K'_T(P,T)`;
- volumetric thermal expansion :math:`\alpha_V(P,T)`;
- thermodynamic and mode-weighted Grüneisen parameters;
- mode-resolved Grüneisen parameters when branch continuity is available;
- equilibrium lattice parameters, axial expansion coefficients, and the
  thermal-expansion tensor when a structural volume path is supplied.

The practical sequence for obtaining and validating these quantities is
described in :doc:`../workflows/qha`.

.. include:: ../_generated/references/qha.inc
