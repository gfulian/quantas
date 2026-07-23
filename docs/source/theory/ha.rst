Harmonic Approximation
======================

The harmonic approximation is the starting point for the lattice-dynamical
thermodynamics implemented in Quantas.  It describes the vibrations of a
crystal in terms of independent quantum harmonic oscillators whose frequencies
are obtained from the normal modes of the lattice [#mcquarrie_simon_1997]_.

The treatment in this section is purely scientific.  Input formats, numerical
options, and execution details are described separately in the
:doc:`../workflows/ha` workflow and in :doc:`../formats/phonon_yaml`.

Lattice vibrations in the harmonic approximation
-------------------------------------------------

Let :math:`U` be the potential energy of a crystal and let
:math:`u_{l\kappa\alpha}` denote the displacement of atom :math:`\kappa` in
cell :math:`l` along Cartesian direction :math:`\alpha`.  Expanding the
potential energy about an equilibrium structure gives

.. math::

   U = U_0
       + \sum_{l\kappa\alpha}
         \left(\frac{\partial U}{\partial u_{l\kappa\alpha}}\right)_0
         u_{l\kappa\alpha}
       + \frac{1}{2}
         \sum_{l\kappa\alpha,l'\kappa'\beta}
         \Phi_{l\kappa\alpha,l'\kappa'\beta}
         u_{l\kappa\alpha}u_{l'\kappa'\beta}
       + \cdots .

At equilibrium the first derivatives vanish.  The harmonic approximation
retains only the quadratic term and neglects cubic and higher-order force
constants.  Diagonalization of the mass-weighted dynamical matrix then
produces phonon frequencies :math:`\nu_{qj}`, where :math:`q` identifies a
wavevector in the Brillouin zone and :math:`j` identifies a phonon branch.

For a crystal containing :math:`N` atoms in the chosen cell, there are
:math:`3N` branches at each wavevector.  The three acoustic branches approach
zero frequency at the Brillouin-zone centre as required by translational
invariance; the remaining branches are optical when more than one atom is
present in the primitive cell.

Independent quantum oscillators
--------------------------------

For one harmonic mode of frequency :math:`\nu`, the energy levels are

.. math::

   \epsilon_n = h\nu\left(n+\frac{1}{2}\right),
   \qquad n=0,1,2,\ldots,

and the canonical partition function is

.. math::

   z(\nu,T)
   = \sum_{n=0}^{\infty} e^{-\beta\epsilon_n}
   = \frac{e^{-\beta h\nu/2}}
          {1-e^{-\beta h\nu}},

where

.. math::

   \beta = \frac{1}{k_{\mathrm B}T}.

Within the harmonic approximation the phonon modes are independent.  The
vibrational partition function of the crystal is therefore the product of the
single-mode partition functions,

.. math::

   Z_{\mathrm{vib}}
   = \prod_{qj} z(\nu_{qj},T)^{w_q},

where :math:`w_q` represents the normalized integration weight associated with
the wavevector :math:`q`.  Equivalently, the thermodynamic potentials are
weighted sums over wavevectors and phonon branches.

Zero-point energy
-----------------

Even at absolute zero the ground state of each harmonic oscillator has a
finite energy.  The total zero-point vibrational energy is

.. math::

   U_{\mathrm{zp}}(V)
   = \sum_{qj} w_q\,\frac{h\nu_{qj}(V)}{2}.

The zero-point contribution is independent of temperature but can depend on
volume because the phonon frequencies may change when the crystal is
compressed or expanded.

Thermal internal energy
-----------------------

The Bose--Einstein occupation number of one phonon mode is

.. math::

   \bar{n}_{qj}
   = \frac{1}
          {\exp\left(h\nu_{qj}/k_{\mathrm B}T\right)-1}.

The thermal vibrational contribution to the internal energy is therefore

.. math::

   U_{\mathrm{th}}(V,T)
   = \sum_{qj} w_q\,
     \frac{h\nu_{qj}(V)}
          {\exp\left[h\nu_{qj}(V)/(k_{\mathrm B}T)\right]-1}.

When a static electronic energy :math:`U_0(V)` is supplied, the total internal
energy used by Quantas is

.. math::

   U(V,T)
   = U_0(V) + U_{\mathrm{zp}}(V) + U_{\mathrm{th}}(V,T).

Helmholtz free energy
---------------------

The vibrational Helmholtz free energy follows directly from the partition
function:

.. math::

   F_{\mathrm{vib}}(V,T)
   = \sum_{qj} w_q\left[
       \frac{h\nu_{qj}(V)}{2}
       + k_{\mathrm B}T
         \ln\left(
           1-e^{-h\nu_{qj}(V)/(k_{\mathrm B}T)}
         \right)
     \right].

The total Helmholtz free energy is

.. math::

   F(V,T) = U_0(V) + F_{\mathrm{vib}}(V,T).

At :math:`T=0`, the thermal logarithmic term vanishes and the vibrational free
energy reduces to the zero-point energy.

Entropy
-------

Introducing the dimensionless variable

.. math::

   x_{qj} = \frac{h\nu_{qj}}{k_{\mathrm B}T},

one convenient expression for the vibrational entropy is

.. math::

   S(V,T)
   = k_{\mathrm B}\sum_{qj} w_q\left[
       \frac{x_{qj}}{e^{x_{qj}}-1}
       - \ln\left(1-e^{-x_{qj}}\right)
     \right].

This expression is consistent with

.. math::

   S = -\left(\frac{\partial F}{\partial T}\right)_V.

Isochoric heat capacity
-----------------------

The vibrational heat capacity at constant volume is

.. math::

   C_V(V,T)
   = k_{\mathrm B}\sum_{qj} w_q\,
     \frac{x_{qj}^2 e^{x_{qj}}}
          {\left(e^{x_{qj}}-1\right)^2}.

The heat capacity is a constant-volume quantity because the harmonic
frequencies are evaluated at a fixed structure and volume.

Low- and high-temperature limits
--------------------------------

For a dynamically stable crystal, the harmonic oscillator expressions satisfy
important limiting behaviours.

At low temperature
   High-frequency modes are only weakly populated.  The entropy and heat
   capacity tend to zero, while the internal and Helmholtz energies retain the
   zero-point contribution.

At high temperature
   Each fully excited harmonic mode approaches the classical equipartition
   limit.  For a cell containing :math:`N` atoms and a complete Brillouin-zone
   average, the molar heat capacity approaches the Dulong--Petit limit

   .. math::

      C_V \longrightarrow 3N R.

Scientific scope and limitations
--------------------------------

The harmonic approximation is appropriate when the vibrational potential is
well represented by its quadratic expansion around a mechanically stable
structure.  Its main assumptions are:

- phonon modes are independent;
- phonon frequencies do not change with temperature at fixed structure;
- explicit phonon--phonon interactions are neglected;
- thermal expansion is absent because the equilibrium structure is fixed;
- consequently, :math:`C_P=C_V` within a purely harmonic treatment.

Imaginary phonon frequencies indicate that the chosen structure is not a local
minimum on the harmonic potential-energy surface, or that the numerical phonon
calculation is insufficiently converged.  The standard harmonic thermodynamic
expressions are not physically defined for such modes and the underlying
structure or calculation should be assessed before interpreting the results.

Relation to the quasi-harmonic approximation
--------------------------------------------

The harmonic approximation can be applied independently at several volumes.
When the volume dependence of the phonon spectrum is used to construct a
free-energy surface :math:`F(V,T)`, the result is the quasi-harmonic
approximation described in :doc:`qha`.

Quantities derived by Quantas
-----------------------------

For each supplied structure or volume, the HA workflow can provide:

- zero-point energy :math:`U_{\mathrm{zp}}`;
- thermal vibrational energy :math:`U_{\mathrm{th}}`;
- total internal energy :math:`U`;
- vibrational and total Helmholtz free energies,
  :math:`F_{\mathrm{vib}}` and :math:`F`;
- entropy :math:`S`;
- isochoric heat capacity :math:`C_V`.

The corresponding execution sequence and output formats are documented in
:doc:`../workflows/ha`.

.. include:: ../_generated/references/ha.inc
