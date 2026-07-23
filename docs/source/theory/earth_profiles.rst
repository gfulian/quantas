Terrestrial pressure--temperature depth profiles
================================================

Scientific purpose
------------------

Thermoelastic properties are functions of thermodynamic state.  In geological
and geophysical applications, however, the scientific question is often not the
response on an arbitrary rectangular pressure--temperature grid, but the
response along a physically motivated path through a crust, lithosphere, or
mantle.  Such a geothermobarometric depth profile supplies that path by
associating each depth with a
pressure and a temperature.

Such a profile can be used to examine, for example, how density, elastic
coefficients, isotropic moduli, or seismic velocities evolve along a selected
continental geotherm, an oceanic cooling model, or a mantle adiabat.  It is a
sampling path through a material-property surface, not an additional equation
of state and not a substitute for phase-equilibrium analysis.

Physical basis of terrestrial profiles
---------------------------------------

Pressure and temperature arise from different physical constraints and should
not be conflated into one universal depth law.

Pressure
   In a spherically symmetric Earth model, pressure follows hydrostatic
   equilibrium from density, gravity, and enclosed mass.  In shallow local
   models it may instead be represented by lithostatic integration through
   explicit density layers.

Continental temperature
   A first reference for stable continental lithosphere is steady conductive
   heat transfer through layers with specified conductivity and radiogenic heat
   production.  Advection, tectonic transients, and magmatic heat transport
   require more specialized models.

Oceanic temperature
   The thermal structure of oceanic lithosphere is commonly represented as a
   transient diffusion problem.  Half-space cooling describes early conductive
   cooling, whereas finite-plate models impose a basal temperature and avoid
   indefinite thickening of the thermal boundary layer.

Mantle temperature
   Away from major thermal boundary layers, mantle profiles are commonly
   approximated as adiabatic paths.  Phase transitions and changes in mineral
   assemblage may introduce slope changes or finite temperature offsets that
   must not be erased by indiscriminate smoothing.

Layering, discontinuities, and mineralogical interpretation
-----------------------------------------------------------

The Earth is radially stratified, and local geological columns may contain an
arbitrary number of scientifically defined layers.  Interfaces can represent
changes in density, heat production, conductivity, tectonic regime, or mineral
assemblage.  A useful profile representation must therefore preserve declared
layer boundaries and make every interpolation, offset, discontinuity, or blend
explicit.

A pressure--temperature path alone does not determine which phase is stable at
each depth.  Before evaluating a thermoelastic tensor along the path, the user
must verify that the phase remains meaningful and that the requested state lies
within the calibrated QHA and elastic-volume domains.  Discontinuities caused
by phase transformations belong to the mineralogical model, even when the
pressure and temperature fields themselves remain continuous.

Representation used by Quantas
-------------------------------

Quantas represents a geological path as three aligned float64 arrays,

.. math::

   z_i,\qquad P(z_i),\qquad T(z_i),

where depth is in km, pressure in GPa, and absolute temperature in K.  Pressure
and temperature are deliberately separate models.  They may therefore be
combined, for example, as PREM pressure plus a user-supplied geotherm, without
embedding a thermoelastic calculator or a frontend in the scientific core.

These paths are reference models, not phase-equilibrium calculations.  Before
a mineral tensor is interpreted along a path, the user must still assess phase
stability and whether the requested pressure, temperature, and volume lie
inside the calibrated thermoelastic domain.

PREM pressure
-------------

The built-in ``prem`` model reconstructs pressure from the density polynomials
of the Preliminary Reference Earth Model (PREM)
[#prem_dziewonski_anderson_1981]_.  Quantas integrates enclosed
mass from the centre,

.. math::

   \frac{dm}{dr}=4\pi r^2\rho(r),

calculates gravity,

.. math::

   g(r)=\frac{Gm(r)}{r^2},

and integrates hydrostatic pressure inward from the surface,

.. math::

   \frac{dP}{dr}=-\rho(r)g(r).

The complete radial density model is used internally so mantle gravity includes
the core mass.  Public mineral-profile presets stop at or above the
core--mantle boundary.  The default 0.25 km integration spacing gives a mass of
approximately :math:`5.973\times10^{24}` kg, surface gravity of approximately
9.82 m s\ :sup:`-2`, and pressures near 13.74, 23.45, and 135.85 GPa at 410,
660, and 2891 km, respectively.  These values are numerical regression anchors,
not replacements for local crustal density information.

PREM is a one-dimensional global-average Earth reference.  For a local shallow
crust, ``LayeredLithostaticPressureModel`` can instead integrate explicit
constant-density layers with a chosen gravitational acceleration.



Continental conductive geotherms
--------------------------------

``ContinentalConductiveGeotherm`` solves steady one-dimensional conduction in
piecewise-homogeneous layers [#hasterok_chapman_2011]_.  Within a layer with
constant conductivity
:math:`k` and heat production :math:`A`,

.. math::

   \frac{d}{dz}\left(k\frac{dT}{dz}\right)+A=0.

For a layer entered with downward conductive heat flux :math:`q_t`, the local
solution is

.. math::

   T(z)=T_t+\frac{q_t}{k}\Delta z-\frac{A}{2k}\Delta z^2,

and the heat flux entering the next layer is reduced by the radiogenic
production integrated over the layer.  Temperature and heat flux are therefore
propagated consistently through the stack.

The built-in ``continental-cratonic``, ``continental-reference``, and
``continental-active`` paths are explicit representative parameterizations of
cool, intermediate, and warm conductive lithospheres.  They are not digitized
site-specific curves from the cited article.  Their surface heat flow, layer
thicknesses, conductivities, and heat-production values are archived in profile
metadata and can be replaced through a profile specification.

The model neglects advection, transient tectonics, lateral variations, and
temperature-dependent conductivity.  These approximations must be reconsidered
for rifts, subduction zones, magmatic regions, or transiently cooling crust.



Oceanic cooling models
-----------------------

The half-space and finite-plate cooling formalisms used here follow the
classical oceanic-lithosphere treatment [#parsons_sclater_1977]_.  The
half-space model evaluates

.. math::

   T(z,t)=T_s+(T_m-T_s)\operatorname{erf}\left[
   \frac{z}{2\sqrt{\kappa t}}\right],

for a semi-infinite medium.  The finite plate model imposes fixed temperatures
at the surface and at a basal depth :math:`L`:

.. math::

   \frac{T-T_s}{T_m-T_s}=
   \frac{z}{L}+\sum_{n=1}^{N}
   \frac{2}{n\pi}\sin\left(\frac{n\pi z}{L}\right)
   \exp\left(-\frac{\kappa n^2\pi^2t}{L^2}\right).

At and below the plate base, Quantas returns the basal mantle temperature.
Built-in scientific paths are supplied at 10, 50, and 100 Ma using an explicit
125 km plate thickness.  Age, diffusivity, boundary temperatures, thickness,
and Fourier truncation are stored in provenance metadata.

These are one-dimensional conductive references.  They do not model ridge-axis
hydrothermal circulation, subduction, sediment blanketing, melt transport, or
lateral mantle-temperature variations.



Katsura (2022) mantle adiabat
-----------------------------

``Katsura2022MantleAdiabat`` is a deterministic reconstruction of the
most-probable dry-pyrolite profile reported by Katsura [#katsura_2022]_.
Quantas uses monotonic
cubic Hermite segments constrained by published temperatures and local
adiabatic gradients.  The reconstructed anchors are:

.. list-table:: Published constraints used by Quantas
   :header-rows: 1
   :widths: 30 35 35

   * - Depth (km)
     - Shallow-side temperature (K)
     - Deep-side temperature (K)
   * - 50
     - 1646
     - --
   * - 410
     - 1799
     - 1860
   * - 520
     - 1899
     - 1942
   * - 660
     - 1994
     - 1960
   * - 2400
     - 2490
     - --
   * - 2800
     - 2587
     - --

Paired depth samples are inserted immediately above and below 410, 520, and
660 km so abrupt phase-boundary changes are not erased by a coarse regular
grid.  At a depth exactly equal to a transition, the deeper assemblage is
selected.

This implementation does **not** re-run the full Monte Carlo MATLAB workflow.
It is an explicitly labelled published-constraint reconstruction.  The article
and archived software [#katsura_software_2022]_ are both recorded in every
generated profile.



Piecewise and tabulated models
------------------------------

A user profile can compose any pressure model and temperature model that share
a depth interval.  A layered specification may contain as many ordered layers
or segments as required by the geological model; each boundary and material
parameter remains explicit.  Tabulated fields support linear or shape-preserving
PCHIP interpolation; extrapolation is never implicit.  Piecewise temperature segments
support three explicit join policies:

``direct``
   Preserve both source models exactly, including a finite temperature jump.

``continuous-offset``
   Add a recorded constant offset to the entering segment so the join is
   continuous.

``blend``
   Linearly blend both models across a declared finite width.  Both source
   models must support the complete blend interval.

Every offset, blend width, table source, interpolation method, citation, and
nested model is retained in the profile metadata and therefore in the
thermoelastic HDF5 result.

Built-in scientific presets
---------------------------

The scientific presets currently exposed by the core and thermoelastic CLI are:

``continental-cratonic``
   Cool 0--200 km conductive continental reference, 40 mW m\ :sup:`-2`.

``continental-reference``
   Intermediate 0--120 km conductive continental reference, 55 mW m\ :sup:`-2`.

``continental-active``
   Warm 0--80 km conductive continental reference, 70 mW m\ :sup:`-2`.

``oceanic-10ma``, ``oceanic-50ma``, ``oceanic-100ma``
   Finite 125 km plate-cooling paths at the indicated ages.

``mantle-katsura-2022``
   PREM pressure plus the reconstructed dry-pyrolite mantle adiabat from 50 to
   2800 km.

Synthetic legacy presets are not exposed.  A deliberately linear path can
still be requested explicitly for numerical experiments, while editable
composed profiles should start from ``quantas thermoelasticity profile-template``.

.. include:: ../_generated/references/earth_profiles.inc
