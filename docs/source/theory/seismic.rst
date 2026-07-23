Seismic-wave propagation
========================

The elastic tensor controls not only static deformation but also the long-
wavelength acoustic modes of a crystal.  In the continuum limit, combining the
stiffness tensor with the density determines the direction-dependent phase
velocities, polarizations, group velocities, and energy-flow geometry of the
three acoustic branches [#nye_1985]_.

SEISMIC is kept separate from the Elasticity workflow because the scientific
questions are different.  Elasticity characterizes deformation; SEISMIC
characterizes propagation.  The latter additionally requires density and must
distinguish wave-normal directions from ray or energy-flow directions.

Equation of motion and Christoffel matrix
-----------------------------------------

For a homogeneous elastic medium without body forces, the equation of motion is

.. math::

   \rho\,\ddot u_i
   = C_{ijkl}\frac{\partial^2u_k}{\partial x_j\partial x_l},

where :math:`\rho` is the density.  Consider a plane wave

.. math::

   u_i=A_i\exp\left[i\left(\mathbf k\cdot\mathbf x-\omega t\right)\right].

Writing :math:`\mathbf n=\mathbf k/|\mathbf k|` and
:math:`v=\omega/|\mathbf k|` gives the Christoffel eigenproblem

.. math::

   \Gamma_{ik}(\mathbf n)A_k=\rho v^2 A_i,

with

.. math::

   \Gamma_{ik}(\mathbf n)=C_{ijkl}n_jn_l.

Equivalently, if the density-normalized Christoffel matrix is used, its
eigenvalues are directly :math:`v^2`.

Acoustic modes and phase velocity
---------------------------------

For every wave-normal direction, the symmetric Christoffel matrix has three
real eigenvalues when the medium is mechanically stable.  Quantas orders the
corresponding speeds as

.. math::

   V_{S2}\le V_{S1}\le V_P,

where

- :math:`V_P` is the fastest, quasi-longitudinal branch;
- :math:`V_{S1}` is the faster quasi-shear branch;
- :math:`V_{S2}` is the slower quasi-shear branch.

The adjective *quasi* is important.  Away from symmetry directions, the
polarization of the fastest branch need not be exactly parallel to
:math:`\mathbf n`, and the shear polarizations need not be exactly
perpendicular to it.

The phase velocity describes propagation of a surface of constant phase.  Its
direction is the wave normal :math:`\mathbf n`, not necessarily the direction
of energy transport.

Polarization
------------

The normalized Christoffel eigenvectors define polarization axes.  Because an
eigenvector and its negative describe the same physical axis, the sign returned
by an eigensolver is arbitrary.  Moreover, when two eigenvalues approach each
other, their individual eigenvectors can rotate or exchange labels rapidly even
though the degenerate subspace remains well defined.

A direction-by-direction calculation therefore distinguishes two concepts:

- **local velocity ordering**, which always labels the slow shear, fast shear,
  and fastest branch at the current direction;
- **tracked polarization branches**, which choose signs and shear permutations
  to maximize continuity along a deterministic path.

Tracking changes the continuity labels of the polarization axes; it does not
change the locally ordered phase velocities.

Degeneracy and near-degeneracy
------------------------------

Two acoustic modes are degenerate when their Christoffel eigenvalues are equal.
At exact degeneracy, individual eigenvectors within the degenerate subspace are
not unique.  Close to degeneracy, they are numerically ill-conditioned and
small perturbations can produce large apparent rotations.

For adjacent modes :math:`a` and :math:`b`, a useful diagnostic is the absolute
gap

.. math::

   \Delta\lambda_{ab}=|\lambda_a-\lambda_b|,

and its value relative to the largest directional eigenvalue.  Numerical
degeneracy is identified using explicit absolute and relative tolerances.  A
degenerate solution may still have valid phase speeds, but mode-specific group
and polarization quantities may not be uniquely resolved.

Group velocity and ray direction
--------------------------------

The group velocity [#jaeken_cottenier_2016]_ is

.. math::

   \mathbf v_g=\nabla_{\mathbf k}\omega.

In a nondispersive elastic continuum, :math:`\omega=kv(\mathbf n)`, so the group
velocity can be obtained from the directional derivative of the Christoffel
eigenvalue.  If :math:`\lambda=v^2`, then

.. math::

   \mathbf v_g=\frac{1}{2v}\nabla_{\mathbf n}\lambda.

The group speed is :math:`|\mathbf v_g|`, while the unit ray direction is

.. math::

   \mathbf r=\frac{\mathbf v_g}{|\mathbf v_g|}.

In an isotropic material, :math:`\mathbf r=\mathbf n`.  In an anisotropic
crystal they generally differ.

Power-flow angle
----------------

The power-flow angle measures the deviation between wave-normal and ray
directions:

.. math::

   \psi=\cos^{-1}\left(\mathbf n\cdot\mathbf r\right).

A non-zero :math:`\psi` means that energy propagates in a different direction
from the normal to the phase front.  This distinction is essential when
interpreting ray paths, phonon focusing, and directional maps.

Phase and group velocity surfaces
---------------------------------

A phase-velocity surface is parameterized by the wave normal,

.. math::

   \mathbf x_{\mathrm{phase}}(\mathbf n)
   =v(\mathbf n)\mathbf n.

A group-velocity surface is parameterized by the ray vector,

.. math::

   \mathbf x_{\mathrm{group}}(\mathbf n)
   =\mathbf v_g(\mathbf n).

Although both surfaces can be plotted from the same sampled wave-normal grid,
they answer different questions.  The phase surface describes wavefront
kinematics; the group surface describes energy transport.

Shear-wave splitting
--------------------

An anisotropic crystal generally supports two distinct shear speeds.  Their
separation can be described by

.. math::

   \Delta V_S=V_{S1}-V_{S2},

or by a normalized splitting measure, for example

.. math::

   A_S=\frac{2(V_{S1}-V_{S2})}{V_{S1}+V_{S2}}.

The splitting vanishes along directions where the two shear modes are
degenerate.  The orientation and magnitude of shear splitting are often more
informative than a single isotropic shear velocity.

Acoustic enhancement and phonon focusing
-----------------------------------------

A uniform solid angle in wave-normal space does not generally map to a uniform
solid angle in ray-direction space.  The local Jacobian of the transformation
from :math:`\mathbf n` to :math:`\mathbf r` determines whether neighboring wave
normals converge or diverge in energy-flow space.

Let

.. math::

   J_{ij}=\frac{\partial r_i}{\partial n_j}.

The area factor used by SEISMIC is obtained from the cofactor matrix of
:math:`\mathbf J` acting on the wave normal.  The acoustic enhancement is
proportional to the inverse of this area factor.  Large enhancement indicates
focusing of acoustic energy; small enhancement indicates defocusing.

Because the inverse diverges when the area factor tends to zero, enhancement
must be interpreted together with numerical conditioning and degeneracy
information.

Caustic candidates
------------------

A caustic occurs where the mapping from wave-normal directions to ray
directions becomes locally singular.  Numerically, SEISMIC marks a **caustic
candidate** when the calculated area factor is sufficiently close to zero under
explicit absolute and relative tolerances.

The word *candidate* is deliberate.  A sampled point near a singular mapping is
a diagnostic that deserves refinement; it is not by itself proof of a complete
caustic curve or surface.  Confirmation requires examining neighboring
directions, sampling convergence, mode resolution, and the stability of the
area factor.

Isotropic reference velocities
------------------------------

For an isotropic medium with bulk modulus :math:`K`, shear modulus :math:`G`,
and density :math:`\rho`, the familiar reference speeds are

.. math::

   V_P=\sqrt{\frac{K+4G/3}{\rho}},

.. math::

   V_S=\sqrt{\frac{G}{\rho}}.

These are useful summaries when :math:`K` and :math:`G` are taken from a
Voigt--Reuss--Hill average, but they do not reproduce the complete directional
Christoffel solution of an anisotropic single crystal.

Physical validity and interpretation
------------------------------------

A meaningful acoustic calculation requires:

- a symmetric stiffness tensor in a known Cartesian and Voigt convention;
- positive density;
- a mechanically stable elastic response for the state being modeled;
- consistent units;
- sufficient angular resolution for directional conclusions.

Small negative eigenvalues may arise from floating-point noise and can be
clamped only within a documented tolerance.  Substantial negative eigenvalues
indicate a physically invalid direction or input tensor and cannot be repaired
by numerical presentation.

What Quantas derives in SEISMIC
-------------------------------

For sampled wave-normal directions, the SEISMIC workflow can provide:

- Christoffel eigenvalues and phase speeds;
- quasi-longitudinal and split quasi-shear branches;
- polarization axes and continuity-tracking diagnostics;
- absolute and relative eigenvalue gaps;
- degeneracy candidates;
- group-velocity vectors, group speeds, and ray directions;
- power-flow angles;
- shear-wave splitting and directional anisotropy;
- acoustic enhancement and area factors;
- caustic candidates;
- two-dimensional maps and three-dimensional phase/group surfaces.

The results describe the long-wavelength elastic continuum.  They do not
replace a full phonon-dispersion calculation at finite wave vector.  The
scientific scope of the Quantas SEISMIC module is described in
[#seismic_ulian_valdre_2024]_.

.. include:: ../_generated/references/seismic.inc
