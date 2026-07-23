Elasticity
==========

Elasticity describes the reversible response of a solid to a sufficiently small
applied stress.  For a crystal this response is generally anisotropic: the
strain produced by a given load depends on the orientation of the load relative
to the crystallographic frame [#nye_1985]_.

The elastic tensor therefore contains more information than a single bulk or
Young's modulus.  It connects stress and strain, determines mechanical
stability, controls the response of polycrystalline aggregates, and provides
the starting point for the acoustic-wave analysis described in
:doc:`seismic`.

Stress, strain, and Hooke's law
-------------------------------

The infinitesimal strain tensor is

.. math::

   \varepsilon_{ij}
   = \frac{1}{2}\left(
   \frac{\partial u_i}{\partial x_j}
   + \frac{\partial u_j}{\partial x_i}
   \right),

where :math:`\mathbf u` is the displacement field.  In the linear elastic
regime, stress and strain are related by

.. math::

   \sigma_{ij}=C_{ijkl}\varepsilon_{kl},

or equivalently

.. math::

   \varepsilon_{ij}=S_{ijkl}\sigma_{kl},

where :math:`C_{ijkl}` is the stiffness tensor and
:math:`S_{ijkl}` is the compliance tensor.  In matrix form,
:math:`\mathbf S=\mathbf C^{-1}`.

The elastic strain-energy density is

.. math::

   U_{\mathrm{el}}
   =\frac{1}{2}C_{ijkl}\varepsilon_{ij}\varepsilon_{kl}
   =\frac{1}{2}S_{ijkl}\sigma_{ij}\sigma_{kl}.

The tensor symmetries

.. math::

   C_{ijkl}=C_{jikl}=C_{ijlk}=C_{klij}

reduce the general fourth-rank tensor to at most 21 independent coefficients.
Crystal symmetry reduces this number further.  A cubic crystal, for example,
has three independent second-order elastic constants, whereas a triclinic
crystal may retain all 21.

Voigt representation and shear convention
------------------------------------------

A symmetric second-rank tensor has six independent components.  Quantas uses
the conventional Voigt ordering

.. math::

   (11,22,33,23,13,12)\longleftrightarrow(1,2,3,4,5,6).

The strain vector uses engineering shear strains,

.. math::

   \boldsymbol\varepsilon^{\mathrm V}
   = (\varepsilon_{11},\varepsilon_{22},\varepsilon_{33},
      2\varepsilon_{23},2\varepsilon_{13},2\varepsilon_{12}),

while the stress vector is

.. math::

   \boldsymbol\sigma^{\mathrm V}
   = (\sigma_{11},\sigma_{22},\sigma_{33},
      \sigma_{23},\sigma_{13},\sigma_{12}).

Hooke's law can then be written as

.. math::

   \boldsymbol\sigma^{\mathrm V}
   =\mathbf C^{\mathrm V}\boldsymbol\varepsilon^{\mathrm V}.

The factors of two associated with engineering shear strain must be treated
explicitly when converting a compliance matrix to a full Cartesian tensor.
They are not optional formatting choices: omitting them changes directional
elastic properties.

Mechanical stability
--------------------

For an unstressed crystal in the linear elastic regime, a mechanically stable
elastic response requires the strain energy to be positive for every non-zero
strain.  In matrix form this means

.. math::

   \boldsymbol\varepsilon^{\mathrm V\,T}
   \mathbf C^{\mathrm V}
   \boldsymbol\varepsilon^{\mathrm V}>0
   \qquad
   \text{for all }\boldsymbol\varepsilon^{\mathrm V}\ne\mathbf 0.

Thus the stiffness matrix must be positive definite, or equivalently all its
eigenvalues must be positive [#mouhat_coudert_2014]_:

.. math::

   \lambda_i(\mathbf C)>0.

Symmetry-specific Born criteria provide algebraically equivalent conditions
when the crystal class and tensor convention are known.  The general
eigenvalue test is especially useful for a frontend-neutral analysis because it
does not require the user to choose a crystal system before assessing the
matrix.

At finite external stress, stability must be formulated using the appropriate
stress-corrected elastic coefficients.  A stiffness matrix reported under
pressure must therefore be interpreted together with the convention used by
the generating code.

Isotropic moduli from an anisotropic crystal
--------------------------------------------

A single crystal is generally anisotropic, but many applications require an
effective isotropic response, for example for a randomly oriented
polycrystalline aggregate.

Voigt approximation
^^^^^^^^^^^^^^^^^^^

The Voigt approximation assumes uniform strain in all grains.  The bulk and
shear moduli are

.. math::

   K_V=\frac{C_{11}+C_{22}+C_{33}
   +2(C_{12}+C_{13}+C_{23})}{9},

.. math::

   G_V=\frac{C_{11}+C_{22}+C_{33}
   -(C_{12}+C_{13}+C_{23})
   +3(C_{44}+C_{55}+C_{66})}{15}.

Reuss approximation
^^^^^^^^^^^^^^^^^^^^

The Reuss approximation assumes uniform stress and is written using the
compliance matrix:

.. math::

   K_R=\left[S_{11}+S_{22}+S_{33}
   +2(S_{12}+S_{13}+S_{23})\right]^{-1},

.. math::

   G_R=15\left[4(S_{11}+S_{22}+S_{33})
   -4(S_{12}+S_{13}+S_{23})
   +3(S_{44}+S_{55}+S_{66})\right]^{-1}.

Hill average
^^^^^^^^^^^^

For a macroscopically isotropic aggregate, the Voigt and Reuss values provide
upper and lower bounds under their idealized assumptions.  The Hill estimate
[#hill_1952]_ is the arithmetic mean,

.. math::

   K_H=\frac{K_V+K_R}{2},
   \qquad
   G_H=\frac{G_V+G_R}{2}.

Once :math:`K` and :math:`G` are known, the corresponding isotropic Young's
modulus and Poisson ratio are

.. math::

   E=\frac{9KG}{3K+G},

.. math::

   \nu=\frac{3K-2G}{2(3K+G)}.

These estimates summarize an aggregate response.  They do not replace the
full tensor when orientation-dependent behavior matters.

Directional elastic properties
------------------------------

Let :math:`\mathbf n` be a unit direction and :math:`\mathbf m` a unit vector
orthogonal to :math:`\mathbf n`.

Young's modulus
^^^^^^^^^^^^^^^

The directional Young's modulus is

.. math::

   E(\mathbf n)
   =\left(S_{ijkl}n_i n_j n_k n_l\right)^{-1}.

It measures the uniaxial stiffness along :math:`\mathbf n`.

Linear compressibility
^^^^^^^^^^^^^^^^^^^^^^

Under hydrostatic pressure, the linear compressibility along
:math:`\mathbf n` is

.. math::

   \beta(\mathbf n)
   =S_{ijkk}n_i n_j.

Unlike the bulk compressibility, directional linear compressibility can be
negative in anisotropic materials.  Negative values indicate expansion along a
particular direction during hydrostatic compression; they do not imply a
negative total volume compressibility.

Shear modulus
^^^^^^^^^^^^^

For shear involving longitudinal direction :math:`\mathbf n` and transverse
direction :math:`\mathbf m`, Quantas uses

.. math::

   G(\mathbf n,\mathbf m)
   =\left[4S_{ijkl}n_i m_j n_k m_l\right]^{-1}.

Poisson ratio
^^^^^^^^^^^^^

The Poisson ratio associated with axial direction :math:`\mathbf n` and
transverse measurement direction :math:`\mathbf m` is

.. math::

   \nu(\mathbf n,\mathbf m)
   =-\frac{S_{ijkl}n_i n_j m_k m_l}
   {S_{ijkl}n_i n_j n_k n_l}.

A negative value corresponds to auxetic response for that pair of directions.
The sign and magnitude can vary strongly over the unit sphere even when an
isotropic average appears ordinary.

Transverse extrema
^^^^^^^^^^^^^^^^^^

For a fixed :math:`\mathbf n`, both :math:`G` and :math:`\nu` vary as
:math:`\mathbf m` rotates in the plane normal to :math:`\mathbf n`.  Their
minimum and maximum are therefore properties of a two-dimensional transverse
quadratic form.  Diagonalizing that projected form yields exact transverse
extrema without relying on an angular search over :math:`\mathbf m`.

Global extrema and anisotropy
-----------------------------

Sampling all longitudinal directions produces global minima and maxima and the
associated Cartesian directions.  For a strictly positive property, a useful
anisotropy ratio is

.. math::

   A_X=\frac{X_{\max}}{X_{\min}}.

For signed properties such as linear compressibility or Poisson ratio, positive
and negative branches must be interpreted separately.  A single ratio cannot
summarize a surface that crosses zero.

Two- and three-dimensional representations
-------------------------------------------

A directional property can be represented radially [#elate_gaillac_pullumbi_coudert_2016]_:

.. math::

   \mathbf r(\mathbf n)=X(\mathbf n)\,\mathbf n.

For Young's modulus this produces one surface.  Shear modulus and Poisson ratio
produce minimum and maximum transverse branches.  Signed properties require
separate positive and negative surfaces to avoid confusing a negative value
with a reversed direction.

Two-dimensional sections are useful for identifying symmetry and comparing
specific planes.  Three-dimensional surfaces show the complete anisotropy, but
visual appearance depends on radial scale and viewing direction; numerical
extrema should always accompany the figure.

Elasticity at finite pressure and temperature
---------------------------------------------

The relations above describe the elastic response at one specified reference
state.  When pressure or temperature changes, the reference structure, density,
and elastic tensor also change.  A thermodynamically consistent treatment then
requires derivatives of a free energy with respect to strain and must account
for hydrostatic pre-stress [#stixrude_lithgow_bertelloni_2005]_.

Quantas treats this problem in the separate thermoelastic workflow.  Its
current approximation combines a QHA equilibrium-volume surface with the
Eulerian cold finite-strain evolution of the Wallace stiffness tensor.  See
:doc:`thermoelasticity` for the derivation and for the distinction between a
full strain-dependent QHA treatment and the quasi-static approximation.

Tensor rotations and reference frames
--------------------------------------

Elastic coefficients are components of a tensor in a selected Cartesian frame.
For a proper orthogonal transformation :math:`R`, the transformed stiffness is

.. math::

   C'_{ijkl}=R_{ia}R_{jb}R_{kc}R_{ld}C_{abcd}.

A rotation changes the numerical components and the coordinates of extrema,
but not the underlying physical tensor.  Comparisons between calculations are
meaningful only when the frames and Voigt conventions are consistent.

What Quantas derives from the elastic tensor
--------------------------------------------

The Elasticity workflow uses the concepts above to provide:

- stiffness and compliance matrices;
- positive-definiteness diagnostics;
- Voigt, Reuss, and Hill isotropic estimates;
- directional Young's modulus and linear compressibility;
- exact transverse extrema of shear modulus and Poisson ratio;
- global extrema, anisotropy measures, and extremal directions;
- optional tensor rotation;
- two-dimensional sections and three-dimensional directional surfaces.

The workflow describes static or state-specific elasticity.  Pressure- and
temperature-dependent tensors are treated by the Thermoelasticity workflow,
while acoustic propagation additionally requires density and is treated by
SEISMIC.

.. include:: ../_generated/references/elasticity.inc
