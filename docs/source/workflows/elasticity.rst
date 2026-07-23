Elasticity analysis: implementation and workflow
================================================

Purpose and scope
-----------------

The Elasticity workflow analyzes one second-order elastic stiffness tensor at a
specified thermodynamic state.  It converts the supplied engineering-Voigt
matrix into a complete Cartesian tensor, evaluates isotropic aggregate
estimates and mechanical stability, and then determines directional mechanical
properties.

The workflow answers questions such as:

- is the supplied stiffness representation mechanically stable?;
- what are its Voigt, Reuss, and Hill aggregate moduli?;
- along which directions is the crystal stiffest, softest, most compressible,
  or auxetic?;
- how do the directional properties appear on principal-plane sections or
  complete three-dimensional surfaces?;
- how do the tensor components and directional results change when expressed
  in another Cartesian frame?

It does **not** calculate the pressure- or temperature-dependence of the tensor.
A tensor obtained at finite pressure and temperature may be analyzed here, but
its construction belongs to :doc:`thermoelasticity`.  Acoustic-wave
propagation additionally requires density and is handled by :doc:`seismic`.
The tensor definitions and physical equations are introduced in
:doc:`../theory/elasticity`.

Computational pipeline
----------------------

The workflow follows the sequence

.. code-block:: text

   stiffness input in engineering Voigt notation
       │
       ├─ validate shape, finiteness, and matrix symmetry
       ├─ optionally transform tensor components to another frame
       ├─ invert C to obtain compliance S
       ├─ detect the compatible elastic symmetry
       ├─ calculate Voigt, Reuss, and Hill estimates
       ├─ test positive definiteness
       ├─ if stable: calculate global directional extrema
       ├─ optionally sample the three principal Cartesian planes
       ├─ optionally sample complete 3D directional fields
       └─ build ResultData → HDF5 / report / plots / table export

The basic tensor, aggregate, stability, and global-extrema stages are performed
for every successful run.  Two- and three-dimensional fields are optional
because they create larger result arrays and are primarily needed for plotting
or tabular export.

Input validation and tensor convention
--------------------------------------

Quantas expects a finite ``6 × 6`` stiffness matrix in engineering Voigt
notation and in GPa.  The matrix must be symmetric within the historical
Elasticity tolerance

.. math::

   \lVert \mathbf C-\mathbf C^{\mathsf T}\rVert_F \leq 10^{-3}\ \mathrm{GPa}.

The tolerance permits harmless textual rounding, but it is not a license to
supply independently inconsistent upper and lower triangles.  Quantas validates
rather than silently symmetrizes the matrix.  A singular matrix cannot be
inverted and therefore fails before directional properties are evaluated.

Engineering shear factors are handled when the Voigt compliance is converted
to the Cartesian tensor :math:`S_{ijkl}`.  Users should therefore supply the
same engineering-Voigt convention described in :doc:`../formats/elasticity_input`
and must not manually alter shear rows or columns to mimic a Cartesian tensor.

Tensor rotation and analysis frame
----------------------------------

An optional proper orthogonal transformation is applied before any property is
calculated:

.. math::

   C'_{ijkl}=R_{ia}R_{jb}R_{kc}R_{ld}C_{abcd}.

This operation changes the **components and coordinate labels**, not the
physical material.  Scalar extrema and isotropic averages should remain
invariant apart from numerical noise, while the reported Cartesian directions
are expressed in the new analysis frame.

Use rotation when:

- the source code and the desired plotting frame use different axes;
- a tensor must be aligned with a laboratory, sample, or geological frame;
- results from different sources need a common component convention.

Do not use rotation to force a matrix into an ideal symmetry pattern.  If a
nominally symmetric tensor violates the expected equalities beyond ordinary
rounding, the source calculation should be inspected.

Symmetry detection and specialization
-------------------------------------

Quantas compares the stiffness matrix with the component patterns allowed for
cubic, hexagonal, trigonal, tetragonal, orthorhombic, monoclinic, and triclinic
systems.  Component equalities and required zeros are tested with a tolerance
of ``1e-3 GPa``.

The detected label is descriptive and also enables available specialized
expressions.  In the current implementation, orthorhombic tensors use optimized
closed expressions for some pointwise directional evaluations; general tensors
use full Cartesian contractions.  The scientific definitions are identical.

A detected lower symmetry does not necessarily mean that the material has
lower crystallographic symmetry.  It can also indicate:

- insufficient numerical convergence of the source tensor;
- a non-canonical coordinate frame;
- inconsistent rounding;
- finite-pressure coefficients supplied in an inappropriate convention.

Compliance and isotropic aggregate estimates
--------------------------------------------

The compliance matrix is obtained by direct inversion of the stiffness matrix.
Quantas then calculates the Voigt and Reuss bounds and their arithmetic Hill
mean.

These values have distinct meanings:

``Voigt``
   Uniform-strain estimate based on stiffness components.

``Reuss``
   Uniform-stress estimate based on compliance components.

``Hill``
   Arithmetic mean of the Voigt and Reuss values, commonly used as a practical
   estimate for an isotropic polycrystalline aggregate.

The reported :math:`K`, :math:`G`, :math:`E`, and :math:`\nu` are aggregate
properties.  They are not angular averages of the single-crystal directional
curves and they do not replace the directional analysis of an anisotropic
crystal.

A large Voigt--Reuss separation is itself useful information: it indicates
strong elastic anisotropy and warns that a single isotropic aggregate value may
hide substantial directional variability.

Mechanical stability
--------------------

The workflow tests positive definiteness of the symmetric ``6 × 6`` stiffness
representation by calculating its six eigenvalues.  The current acceptance
threshold is exactly zero:

.. math::

   \lambda_i>0\qquad i=1,\ldots,6.

For an unstressed crystal this is the generic Born stability condition.  At
finite hydrostatic pressure it is valid only when the matrix contains the
appropriate Wallace/Barron--Klein stress--strain coefficients.  An uncorrected
second derivative of energy under pre-stress must not be interpreted in the
same way.

If the tensor is not positive definite, Quantas still returns the matrix,
compliance when invertible, aggregate estimates, eigenvalues, and a warning.
Directional extrema and sampled fields are skipped because their compliance
denominators need not remain positive or physically meaningful.

The minimum eigenvalue is a stability margin, not a distance to every possible
structural transformation.  A small positive value can be physically
important and should be compared with numerical uncertainty and with the
thermodynamic state at which the tensor was obtained.

Global directional extrema
--------------------------

For every mechanically stable tensor, Quantas calculates extrema of:

- Young's modulus;
- linear compressibility;
- shear modulus;
- Poisson's ratio.

Young's modulus and linear compressibility depend on the longitudinal direction
:math:`\mathbf n`.  Shear modulus and Poisson's ratio also depend on a
transverse measurement direction :math:`\mathbf m\perp\mathbf n`.

The global-extrema report uses a fixed historical search strategy:

- two-angle properties use a ``25 × 25`` brute-force angular search;
- three-angle properties use a ``10 × 10 × 10`` brute-force search;
- the best grid point is refined with a local Nelder--Mead minimization;
- maxima are found by minimizing the negative property.

This search is independent of the ``--ntheta`` and ``--nphi`` options used for
plotting fields.  It is designed to locate robust extrema without constructing
a dense global surface, but it remains a numerical global-search approximation.
For exceptionally sharp or nearly degenerate extrema, verify the result against
a denser 3D field or against symmetry directions known from the crystal.

The reported anisotropy is

.. math::

   A=\frac{x_{\max}}{x_{\min}}

when the minimum is positive.  If a property can cross zero, such as linear
compressibility or Poisson's ratio, this ratio is not a useful complete measure
and may be infinite.  Inspect the signed extrema and their directions instead.

Principal-plane 2D fields
-------------------------

With ``--2d``, Quantas samples the ``xy``, ``xz``, and ``yz`` Cartesian planes.
The default grid has **361 points** on each closed curve, corresponding to a
one-degree angular interval including both endpoints.

All four property families are calculated on every plane in a single batched
Cartesian-compliance evaluation.  For each longitudinal direction:

- Young's modulus and linear compressibility are direct contractions;
- the minimum and maximum shear moduli are obtained algebraically from the
  eigenvalues of a projected symmetric ``2 × 2`` compliance form;
- the minimum and maximum Poisson ratios are obtained by the analogous exact
  transverse-plane eigensystem.

Thus the transverse extrema on each sampled direction do **not** depend on a
local optimizer or on an auxiliary sampling of the transverse angle.  The only
sampling approximation is the spacing of the longitudinal directions around
the selected plane.

For signed quantities, separate plotting branches are retained:

- positive and negative linear compressibility;
- negative Poisson ratio;
- positive minimum Poisson ratio;
- maximum Poisson ratio.

This prevents a negative radius from being confused with a direction reversal
in a polar plot.

Choosing the 2D resolution
~~~~~~~~~~~~~~~~~~~~~~~~~~

``ntheta=361`` is appropriate for publication-quality smooth curves in most
materials.  Lower values are useful for rapid checks, while very narrow lobes
or near-zero sign changes may require refinement.

A practical test is:

.. list-table:: Suggested 2D convergence check
   :header-rows: 1
   :widths: 22 20 20 20 18

   * - Angular points
     - Plane maximum
     - Plane minimum
     - Zero crossings
     - Time
   * - 181
     - …
     - …
     - …
     - …
   * - 361
     - …
     - …
     - …
     - …
   * - 721
     - …
     - …
     - …
     - …

Refine the grid when extrema or sign-change locations move by more than the
scientific precision required for the study.  Increasing the number of points
only makes the sampled curve denser; it does not improve the source tensor.

Three-dimensional directional fields
-------------------------------------

With ``--3d``, Quantas builds a regular full-sphere grid.  The default is:

.. code-block:: text

   ntheta = 61
   nphi   = 121
   directions = 61 × 121 = 7381 grid positions

The polar grid includes both poles.  The azimuthal grid excludes ``2π`` because
that meridian duplicates ``0``.  The default relation used by the CLI is

.. math::

   N_\phi=2N_\theta-1,

which gives similar nominal polar and azimuthal spacing away from the poles.

Selected properties are evaluated with vectorized NumPy contractions.  As in
the 2D workflow, the shear and Poisson transverse extrema are solved
algebraically for every longitudinal direction.  No interpolation or local
optimizer is used to construct the surface fields.

Available surface branches are:

- Young's modulus;
- positive and negative linear compressibility;
- minimum and maximum shear modulus;
- negative, positive-minimum, and maximum Poisson ratio.

A branch that is absent everywhere is omitted and recorded as a non-fatal
surface warning.

Persisted and transient surfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two valid workflows.

**Persist during the calculation**
   Use ``elasticity run --3d``.  The numerical surfaces and their sampling
   metadata are stored in HDF5 and can be plotted later without recomputation.

**Calculate only for a plot**
   Run the basic analysis first, then use ``elasticity plot --3d``.  If the HDF5
   file does not contain the requested surfaces, Quantas reconstructs them
   transiently from the stored stiffness tensor.  Explicit plot-grid options
   always request a fresh transient calculation.

Persist surfaces when they will be reused, exported, or compared numerically.
Use transient plotting when only a few final figures are required and a smaller
HDF5 file is preferred.

Physical and unit-sphere geometry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``physical``
   The plotted radius equals the property magnitude.  Lobes and indentations
   show the directional response directly.

``unit-sphere``
   The geometry remains spherical and the property is represented by color.
   This is useful when sign, relative variation, or several branches are easier
   to compare without large geometric distortion.

The geometry changes only the plot specification.  It does not change the
stored physical values.

Three-dimensional convergence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cost scales approximately as

.. math::

   N_\theta N_\phi N_{\mathrm{selected\ properties}}.

A useful convergence sequence is:

.. list-table:: Suggested 3D convergence check
   :header-rows: 1
   :widths: 18 18 18 20 26

   * - :math:`N_\theta`
     - :math:`N_\phi`
     - Directions
     - Surface extrema
     - Visual/sign topology
   * - 31
     - 61
     - 1891
     - …
     - …
   * - 61
     - 121
     - 7381
     - …
     - …
   * - 121
     - 241
     - 29161
     - …
     - …

The fixed global extrema in the report do not change when this grid is changed.
Use the comparison to assess the sampled surfaces, the visibility of narrow
features, and agreement with the independently optimized extrema.

Batching and memory
-------------------

The public 3D surface API uses a default ``batch_size`` of **65536** directions.
For the default 7381-point surface this normally evaluates the complete grid in
one batch.  The batch size controls temporary NumPy-array memory and call
overhead; it does **not** alter the angular grid, formulas, or numerical
accuracy.

- decrease it if a very dense grid or limited-memory environment causes high
  peak memory;
- increase it only when a dense calculation is split into many small batches
  and sufficient memory is available;
- do not use it as a convergence parameter.

The CLI currently exposes surface resolution but not this advanced batching
control.  Python users can set it through
:class:`quantas.api.elasticity.SurfaceOptions`.

Performance and practical acceleration
--------------------------------------

The inexpensive stages are matrix validation, inversion, symmetry detection,
VRH averaging, and stability.  The fixed global-extrema search is more costly
but independent of requested plot resolution.  Optional sampled fields dominate
runtime and HDF5 size when dense grids are selected.

Recommended acceleration strategy:

1. run the basic workflow first without ``--2d`` or ``--3d``;
2. inspect symmetry, stability, aggregate moduli, and global extrema;
3. add ``--2d`` for a low-cost directional overview;
4. request only the 3D properties needed for the scientific question;
5. start from ``31 × 61`` or the default ``61 × 121`` grid;
6. refine only if extrema, sign branches, or surface topology are unresolved;
7. persist 3D fields only when they will be reused.

Selecting only ``young`` is substantially cheaper and smaller than storing all
Young, compressibility, shear, and Poisson branches.  Plot resolution and image
DPI are renderer concerns and should not be confused with the numerical sphere
resolution.

Defaults and rationale
----------------------

.. list-table:: Elasticity defaults
   :header-rows: 1
   :widths: 27 24 49

   * - Control
     - Default
     - Rationale
   * - 2D calculation
     - disabled
     - Basic scientific results do not require large plane arrays.
   * - 2D angular points
     - 361
     - One-degree closed curves are smooth for most figures.
   * - 3D calculation
     - disabled
     - Complete surfaces increase runtime and HDF5 size.
   * - 3D grid
     - 61 × 121
     - Balanced first production grid for full-sphere visualization.
   * - 3D properties
     - all four families
     - Complete API default; CLI users can select only required properties.
   * - Surface batch size
     - 65536
     - Vectorizes ordinary grids in one batch while remaining configurable.
   * - Surface geometry
     - physical
     - Directly represents property magnitude as radius.
   * - Mechanical-stability tolerance
     - 0 GPa
     - Tests strict positive definiteness of the supplied representation.

.. note::

   The value ``512`` is **not** an Elasticity convergence or batch default.  In
   the current source, the relevant 3D Elasticity batch default is 65536.  The
   value 512 belongs to the SEISMIC batching workflow described in
   :doc:`seismic`.

Warnings and failure modes
--------------------------

The calculation stops when:

- the stiffness matrix has the wrong shape or non-finite entries;
- the matrix is asymmetric beyond tolerance;
- the matrix is singular;
- a rotation matrix is not a proper orthogonal transformation;
- a requested sampling grid or property identifier is invalid;
- a required compliance denominator is non-positive or non-finite.

A mechanically unstable but invertible tensor produces a result and warning,
but directional analysis is intentionally skipped.

Surface warnings can indicate that a signed branch is absent—for example, no
negative compressibility or no auxetic branch.  This is a physical result, not
an execution failure.

Results, export, and interoperability
-------------------------------------

The native result can contain:

- analysis-frame stiffness and compliance;
- detected symmetry;
- Voigt, Reuss, and Hill estimates;
- stability eigenvalues;
- global directional extrema and axes;
- optional principal-plane fields;
- optional 3D surfaces and sampling diagnostics;
- tensor-frame provenance.

The tabular exporter operates on stored 2D fields.  A stiffness tensor and
optional density can also be reused by SEISMIC.  Thermoelastic Point analyses
can write the same text format, allowing a tensor reconstructed at one
:math:`(P,T)` state to enter this workflow without a private conversion step.

Decision guide
--------------

.. list-table:: Practical Elasticity decisions
   :header-rows: 1
   :widths: 36 28 36

   * - Situation
     - Recommended choice
     - Reason
   * - First check of a new tensor
     - basic run only
     - Validate matrix, symmetry, stability, VRH values, and extrema quickly.
   * - Publication polar plots
     - ``--2d`` with 361 points
     - Smooth principal-plane curves at modest cost.
   * - Search for auxeticity or negative compressibility
     - 2D first, then selected 3D property
     - Signed branches are easier to diagnose progressively.
   * - Reusable 3D dataset
     - persist ``--3d`` fields
     - Avoid repeated sampling and preserve numerical provenance.
   * - One final 3D figure only
     - transient plot from basic HDF5
     - Keeps the native result compact.
   * - Sharp surface features
     - compare 61 × 121 and 121 × 241
     - Demonstrates angular convergence.
   * - Very dense grid causes high memory use
     - reduce API ``batch_size``
     - Changes memory batching, not scientific resolution.
   * - Tensor at finite hydrostatic pressure
     - verify Wallace convention first
     - Positive-definiteness and waves require stress-corrected coefficients.
   * - Different laboratory/crystal frames
     - apply an explicit rotation
     - Maintains provenance and transforms directions consistently.

Related documentation
---------------------

- Scientific theory: :doc:`../theory/elasticity`
- Input format: :doc:`../formats/elasticity_input`
- Detailed CLI/API example: :doc:`../tutorials/elasticity`
- Public Python API: :doc:`../api/elasticity`
- Command reference: :doc:`../cli/elasticity`
- Acoustic propagation: :doc:`seismic`
