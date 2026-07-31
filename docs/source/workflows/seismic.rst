Seismic-wave analysis: implementation and workflow
==================================================

Purpose and scope
-----------------

The SEISMIC workflow evaluates acoustic-wave propagation through an anisotropic
elastic medium.  Its scientific input is a mechanically stable stiffness tensor
and a positive mass density.  For each sampled wave-normal direction, Quantas
can calculate phase velocities and polarization axes, analytical group
velocities and ray directions, and curvature-derived acoustic enhancement.

The workflow addresses questions such as:

- how do :math:`V_P`, :math:`V_{S1}`, and :math:`V_{S2}` vary with direction?;
- where is shear-wave splitting largest or where are the two shear modes
  degenerate?;
- how far does the energy-flow direction deviate from the wave normal?;
- how do polarization axes evolve across the sphere?;
- where does the ray mapping show strong focusing or possible caustic behavior?

The physical derivation is given in :doc:`../theory/seismic`.  This page
explains how Quantas samples and organizes that theory, how the numerical
controls should be interpreted, and which approximations remain in a finite
spherical grid.

Computational pipeline
----------------------

.. code-block:: text

   stiffness C + density rho
       │
       ├─ validate matrix, density, and numerical tolerances
       ├─ optionally transform tensor components to another frame
       ├─ require positive-definite stiffness
       ├─ calculate VRH isotropic reference velocities
       ├─ construct a regular spherical wave-normal grid
       ├─ solve the Christoffel eigensystem in vectorized batches
       ├─ optionally calculate analytical group quantities
       ├─ optionally calculate analytical curvature and enhancement
       ├─ optionally track polarization continuity across the grid
       └─ build ResultData → HDF5 / report / CSV / 2D and 3D plots

Phase results are always calculated.  ``group`` includes the phase stage, and
``enhancement`` includes both phase and group stages.

Input and physical preconditions
--------------------------------

A shared Quantas elastic input can be generated directly from supported external
outputs with either the public API or the CLI::

   quantas seismic inpgen OUTPUT --interface crystal --output material.dat
   quantas seismic inpgen OUTCAR --interface vasp --output material.dat

Both entry points use ``quantas.api.seismic.create_input`` and require density
metadata in addition to the stiffness tensor.

The stiffness matrix must be finite, symmetric, and expressed in GPa.  SEISMIC
uses the historical element-wise symmetry criterion

.. math::

   \max_{IJ}|C_{IJ}-C_{JI}|\leq10^{-8}\ \mathrm{GPa},

which is stricter than the aggregate Frobenius criterion used by the standalone
Elasticity workflow.  The density must be finite, positive, and expressed in
``kg m^-3``.

Positive definiteness is mandatory.  Unlike Elasticity, which can still report
basic diagnostics for an unstable tensor, SEISMIC stops because negative
Christoffel eigenvalues would make acoustic phase speeds physically undefined
for at least some directions.

At finite hydrostatic pressure, the tensor must contain the appropriate
stress--strain coefficients for wave propagation, normally the Wallace form.
Passing an uncorrected energy Hessian can produce apparently precise but
physically inconsistent velocities.

Tensor rotation
---------------

Rotation is applied before constructing the elastic medium.  It changes tensor
components, wave-normal coordinates, polarization coordinates, and reported
extremal directions, while leaving scalar physical predictions invariant under
a consistent frame transformation.

Use the same frame convention for:

- stiffness components;
- density-associated structure metadata;
- crystallographic axis labels;
- comparison with experimental or geological directions.

A plot labelled ``[100]``, ``[010]``, and ``[001]`` refers to the analysis
frame represented by the stored tensor, not automatically to the original
source-code cell.

Isotropic reference velocities
------------------------------

Before directional sampling, Quantas computes isotropic velocities from the
Hill bulk and shear moduli and the supplied density.  These values are useful
references for an equivalent isotropic aggregate.

They are **not**:

- averages of all sampled single-crystal phase velocities;
- replacements for directional extrema;
- values inferred from the most symmetric propagation directions.

A large difference between directional fields and the Hill reference is a
measure of the importance of single-crystal anisotropy.

Spherical sampling domain
-------------------------

Quantas uses a regular tensor-product grid in polar angle :math:`\theta` and
azimuth :math:`\phi`.  The azimuthal endpoint at :math:`2\pi` is excluded so
that the first meridian is not duplicated.

The supported domains are:

``upper`` — default
   :math:`0\leq\theta\leq\pi/2`.

``lower``
   :math:`\pi/2\leq\theta\leq\pi`.

``full``
   :math:`0\leq\theta\leq\pi`.

For an ordinary elastic medium, phase-speed magnitudes are antipodally
symmetric: :math:`\mathbf n` and :math:`-\mathbf n` describe the same wave-normal
axis.  The upper hemisphere is therefore sufficient for most scalar maps and
can be completed antipodally for 3D plotting.  It also avoids storing a second
physical copy of the same directional information.

The meaning of ``ntheta`` depends on the selected domain.  With the default
``upper`` hemisphere, ``ntheta=91`` gives approximately one-degree polar
spacing.  The same ``ntheta`` on the ``full`` sphere spans twice the polar range
and therefore gives approximately two-degree spacing.  To retain similar polar
resolution when switching from upper to full, approximately double
``ntheta``.

Default grid and sampled extrema
--------------------------------

The default calculation uses:

.. code-block:: text

   hemisphere = upper
   ntheta     = 91
   nphi       = 181
   positions  = 91 × 181 = 16471

The polar spacing is about one degree and the azimuthal spacing about two
degrees.  Extrema, acoustic-axis candidates, splitting maxima, and caustic
candidates in the report are selected from this finite grid.  They are not
continuous global optimizations.

Consequently:

- smooth phase-velocity extrema often converge rapidly;
- a narrow degeneracy or caustic can lie between grid points;
- enhancement and area factors usually require finer convergence checks than
  phase speed;
- a count of zero candidates means none were found on the selected grid, not
  that the continuous surface contains none.

Sampling levels
---------------

``phase``
~~~~~~~~~

For every direction, Quantas constructs the symmetric Christoffel matrix and
uses ``numpy.linalg.eigh`` to obtain three ordered eigenpairs.  The eigenvalues
are stored in ascending order as:

.. code-block:: text

   V_S2  slow quasi-shear
   V_S1  fast quasi-shear
   V_P   quasi-longitudinal

This level stores:

- eigenvalues;
- phase speeds;
- local polarization axes;
- absolute and relative eigenvalue gaps;
- validity and clamping masks;
- degeneracy masks.

Choose ``phase`` when only directional phase speeds, splitting, and local
polarizations are needed.  It is the fastest and smallest output level.

``group``
~~~~~~~~~

The group level additionally evaluates analytical first derivatives of the
Christoffel eigenvalues with respect to the wave-normal coordinates.  From
these derivatives Quantas obtains:

- group-velocity vectors and magnitudes;
- ray directions;
- power-flow angles between wave normal and ray direction.

No finite angular differences are used.  Nevertheless, group quantities are
marked unresolved inside numerically degenerate eigenspaces because an
individual eigenvector derivative is not unique there.

Choose ``group`` when ray geometry and energy-flow deviation are required but
curvature focusing is not.

``enhancement`` — default
~~~~~~~~~~~~~~~~~~~~~~~~~

The enhancement level also calculates analytical eigenvalue Hessians.  It uses
those Hessians to differentiate the ray-direction mapping, determine an area
factor, and evaluate acoustic enhancement:

.. math::

   A=\frac{1}{J},

where :math:`J` is the local area factor of the wave-normal-to-ray mapping.
The logarithmic enhancement and caustic-candidate mask are stored as well.

This is the most expensive level because it allocates and evaluates several
``3 × 3`` tensors for every direction and mode.  It should be selected when
phonon focusing, enhancement maps, or caustic searches are part of the
scientific objective.  It need not be the production default for a study that
only reports velocities.

.. list-table:: Choosing a sampling level
   :header-rows: 1
   :widths: 28 24 24 24

   * - Required result
     - ``phase``
     - ``group``
     - ``enhancement``
   * - Phase speeds and splitting
     - yes
     - yes
     - yes
   * - Local polarizations
     - yes
     - yes
     - yes
   * - Group velocity and ray direction
     - no
     - yes
     - yes
   * - Power-flow angle
     - no
     - yes
     - yes
   * - Area factor and enhancement
     - no
     - no
     - yes
   * - Lowest runtime and file size
     - best
     - intermediate
     - highest cost

The CLI default remains ``enhancement`` because it preserves the complete
historical SEISMIC result in one run.  For large convergence studies, choosing
the lowest level that answers the scientific question is the most effective
acceleration.

Christoffel eigenvalue validity
-------------------------------

The reduced Christoffel eigenvalues have units of ``km² s⁻²``.  Numerical
rounding can occasionally produce a very small negative value near a physically
zero eigenvalue.  Quantas defines a direction-dependent threshold

.. math::

   \epsilon_\lambda
   =\epsilon_{\mathrm{abs}}
   +\epsilon_{\mathrm{rel}}\max_i|\lambda_i|.

The defaults are:

.. code-block:: text

   eigenvalue_rtol = 1e-10
   eigenvalue_atol = 1e-12 km² s⁻²

An eigenvalue in :math:`[-\epsilon_\lambda,0)` is accepted and clamped to zero
before taking its square root.  A more negative value is invalid and produces
``NaN`` for the corresponding phase speed.

These tolerances distinguish floating-point noise from a genuinely invalid
Christoffel solution.  They must not be enlarged merely to hide an unstable or
incorrect stiffness tensor.  The clamped and invalid masks are preserved for
diagnostics and export.

Degeneracy detection
--------------------

Adjacent eigenvalue gaps are compared with

.. math::

   \epsilon_{\mathrm{deg}}
   =\delta_{\mathrm{abs}}
   +\delta_{\mathrm{rel}}\max_i|\lambda_i|,

using defaults:

.. code-block:: text

   degeneracy_rtol = 1e-8
   degeneracy_atol = 1e-10 km² s⁻²

A mode whose smallest adjacent gap falls below this threshold belongs to a
numerically degenerate eigenspace.  Quantas retains the phase speeds, but:

- individual polarization vectors within that eigenspace are not unique;
- group and enhancement quantities can be valid but not uniquely resolved by
  branch;
- tracking may rotate a basis inside the degenerate shear subspace;
- candidate counts depend on the selected tolerance and grid.

The tolerance should reflect numerical precision and the scientific resolution
of near-degeneracies.  It is an advanced parameter, not a generic smoothing
control.

Polarization tracking
---------------------

Eigenvectors are axial: :math:`\mathbf p` and :math:`-\mathbf p` describe the
same polarization.  Moreover, the two locally ordered shear eigenvectors can
exchange identities when their speeds approach one another.

When tracking is enabled, Quantas traverses the regular grid in a deterministic
serpentine path and builds continuity branches ``shear_a``, ``shear_b``, and
``p``.  It can:

- flip an eigenvector sign to maximize continuity;
- exchange the two local shear assignments;
- identify ambiguous permutations;
- rotate the basis inside a degenerate shear subspace;
- start a new continuity segment when no unique predecessor exists.

Tracking does not change eigenvalues, locally ordered phase speeds, or the
physical eigenspace.  It provides a continuous representation for plotting and
for following polarization branches.

Disable tracking when:

- only scalar velocities are needed;
- minimizing output size and post-processing time is important;
- local eigenvectors are sufficient and no continuous polarization plot will be
  produced.

Keep it enabled when shear polarization maps, branch continuity, or mode
exchange are scientifically relevant.

Analytical group quantities
---------------------------

Group velocities are calculated from analytical Christoffel eigenvalue
gradients, not from finite differences between neighboring spherical samples.
Therefore increasing ``ntheta`` and ``nphi`` does not improve the derivative at
an already sampled direction; it improves only the angular coverage and the
chance of resolving spatial extrema or narrow features between directions.

The group velocity vector is normalized to obtain the ray direction.  The
power-flow angle is the angle between this ray direction and the original wave
normal.  In an isotropic medium it vanishes; in an anisotropic crystal it can be
substantial.

Near a degeneracy, individual group branches are marked unresolved because the
selected eigenvector basis is not unique.  This is a physical/numerical
property of the eigenspace, not a failure of the analytical derivative.

Enhancement, pseudoinverse, and caustic candidates
--------------------------------------------------

The analytical eigenvalue Hessian contains a pseudoinverse of a shifted
Christoffel matrix.  Its relative singular-value cutoff is controlled by
``pseudoinverse_rcond``, default ``1e-10``.

A smaller cutoff retains more nearly singular information but can amplify
floating-point noise near degeneracies.  A larger cutoff regularizes more
aggressively but may suppress legitimate curvature.  Change it only during a
focused convergence study and always inspect resolved masks and degeneracy
diagnostics.

The caustic threshold is scaled locally:

.. math::

   \epsilon_J
   =J_{\mathrm{abs}}+J_{\mathrm{rel}}\lVert\operatorname{cof}(D\mathbf r)\rVert,

with defaults:

.. code-block:: text

   caustic_rtol = 1e-10
   caustic_atol = 1e-12

A sampled point is labelled a **caustic candidate** when its calculated area
factor is no larger than this threshold.  The term candidate is intentional:

- the result is grid dependent;
- nearby unresolved degeneracies can invalidate curvature branches;
- a true continuous zero may fall between samples;
- finite but very large enhancement is not identical to a zero area factor.

A robust caustic claim requires grid refinement, tolerance checks, and
inspection of neighboring values—not only a non-zero candidate count.

The meaning of ``batch_size=512``
---------------------------------

SEISMIC evaluates independent wave-normal directions in vectorized NumPy
batches.  The default is:

.. code-block:: text

   batch_size = 512 directions

For the default 16471-position grid, this corresponds to 33 batches.  The value
is a compromise between:

- reducing Python-loop and small-kernel overhead;
- bounding temporary arrays used by eigenvectors, gradients, Hessians, and
  pseudoinverses;
- providing frequent enough progress updates for interactive CLI use.

**Batch size is not an accuracy or convergence parameter.**  The same direction
is solved by the same vectorized formulas regardless of which batch contains
it.  Changing the batch size should preserve numerical arrays to floating-point
reproducibility.

Change it when:

- **lower values** are needed to reduce peak temporary memory on constrained
  machines or extremely dense enhancement grids;
- **higher values** improve throughput on a machine with sufficient memory and
  the default calculation is split into many batches.

Do not increase it to obtain more accurate velocities or caustics.  Increase
``ntheta`` and ``nphi`` to refine the angular grid.

A useful performance-only benchmark is:

.. list-table:: Batch-size benchmark template
   :header-rows: 1
   :widths: 22 24 26 28

   * - Batch size
     - Runtime
     - Peak memory
     - Maximum array difference
   * - 128
     - …
     - …
     - reference
   * - 512
     - …
     - …
     - …
   * - 2048
     - …
     - …
     - …

Choose the fastest value that fits comfortably in memory.  The maximum array
difference should be zero or at ordinary floating-point round-off.

Angular convergence
-------------------

The total number of grid positions is

.. math::

   N=N_\theta N_\phi.

Runtime is approximately linear in :math:`N`, with a much larger per-point
constant for ``enhancement`` than for ``phase``.  HDF5 size also grows linearly
with :math:`N`, while individual modes and tensor fields add fixed multiples of
that size.

A staged convergence study is preferable to starting from an extremely dense
grid:

.. list-table:: Suggested upper-hemisphere convergence sequence
   :header-rows: 1
   :widths: 18 18 18 23 23

   * - :math:`N_\theta`
     - :math:`N_\phi`
     - Positions
     - Phase/splitting extrema
     - Enhancement/candidates
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
   * - 91
     - 181
     - 16471
     - …
     - …
   * - 121
     - 241
     - 29161
     - …
     - …

Recommended procedure:

1. converge phase extrema and maximum shear splitting using ``level=phase``;
2. converge power-flow extrema using ``level=group``;
3. only then run ``level=enhancement`` on the resolutions required for
   curvature features;
4. compare candidate locations and neighboring area factors, not only counts;
5. use the same hemisphere definition when comparing grids.

If only a few specific crystallographic directions are required, the public
low-level physical objects can solve those directions directly, but the
persisted workflow is intentionally organized around regular fields for maps,
reports, export, and GUI reuse.

Performance and memory strategy
-------------------------------

The principal costs are:

``phase``
   Batched ``3 × 3`` eigensystems and polarization storage.

``group``
   Phase cost plus analytical gradient tensors and group vectors.

``enhancement``
   Group cost plus Hessians, pseudoinverses, ray gradients, cofactors, and area
   factors.

``polarization tracking``
   A separate deterministic traversal after phase sampling.  Its cost is
   usually lower than enhancement but grows with every grid position.

Practical acceleration sequence:

1. validate a new tensor with Elasticity first;
2. start with ``level=phase`` and a moderate grid;
3. disable tracking if no polarization output is needed;
4. refine angular resolution before enabling enhancement;
5. select the upper hemisphere unless a lower/full domain is specifically
   required;
6. use ``batch_size`` only to balance memory and throughput;
7. render only the maps and surfaces required for the study.

The number of contour levels, image DPI, colormap, polarization stride, and 3D
plot geometry affect rendering only.  They do not refine the stored acoustic
field.

Defaults and rationale
----------------------

.. list-table:: SEISMIC defaults
   :header-rows: 1
   :widths: 28 23 49

   * - Control
     - Default
     - Rationale
   * - Hemisphere
     - upper
     - Uses antipodal symmetry without duplicating the physical axis field.
   * - Sampling level
     - enhancement
     - Produces the complete persisted acoustic dataset in one run.
   * - Polar grid
     - 91
     - Approximately one-degree spacing over the upper hemisphere.
   * - Azimuthal grid
     - 181
     - Approximately two-degree spacing without a duplicated seam.
   * - Batch size
     - 512
     - Balances vectorization, temporary memory, and progress cadence.
   * - Polarization tracking
     - enabled
     - Produces continuous axes suitable for maps and branch analysis.
   * - Eigenvalue tolerance
     - ``1e-10`` relative, ``1e-12`` absolute
     - Clamps only tiny negative numerical eigenvalues.
   * - Degeneracy tolerance
     - ``1e-8`` relative, ``1e-10`` absolute
     - Marks near-equal eigenspaces without requiring exact equality.
   * - Pseudoinverse cutoff
     - ``1e-10``
     - Regularizes analytical Hessians near singular shifted systems.
   * - Caustic tolerance
     - ``1e-10`` relative, ``1e-12`` absolute
     - Identifies sampled area factors numerically consistent with zero.

Warnings and diagnostic masks
-----------------------------

The workflow stops when:

- stiffness or density is invalid;
- the stiffness is not positive definite;
- spherical-grid dimensions are invalid;
- batch size or advanced tolerances are invalid;
- an explicit tensor rotation is invalid.

Non-fatal diagnostics are preserved in arrays and metadata:

- invalid phase-mode values;
- small negative eigenvalues clamped to zero;
- mode and pair degeneracies;
- shear acoustic-axis candidates;
- unresolved group/enhancement branches;
- non-finite enhancement values;
- caustic candidates;
- polarization sign flips, swaps, ambiguous assignments, and subspace rotations.

A warning count should be interpreted together with its mask and location.  For
example, unresolved shear polarizations on a symmetry axis are expected and do
not invalidate the phase speeds there.

During sampling, the calculator emits one operational ``PROGRESS`` event after
each completed batch.  Its ``current`` and ``total`` fields are monotonic and
the final progress value is one.  Progress events are transport state for a live
frontend and are intentionally not duplicated in the persisted event history.
Settings, input, isotropic references, field completion, warnings, final
completion, and errors use the same frontend-neutral event model; non-progress
events are retained in the result envelope and native HDF5 file.

Native HDF5 writing and reopening recheck finite positive density, stiffness
symmetry, positive definiteness, and consistency of the stored stability
eigenvalues.  A non-propagating medium therefore cannot be published or reopened
as an apparently valid native SEISMIC result.

Reports, export, and plots
--------------------------

The HDF5 result preserves the grid and every field available at the selected
level.  The report summarizes:

- tensor and density;
- VRH isotropic reference velocities;
- sampled extrema and multiplicities;
- shear splitting and anisotropy;
- degeneracy and acoustic-axis diagnostics;
- group and enhancement diagnostics when available.

CSV export writes one row per sampled grid position and acoustic mode, including
wave-normal coordinates, velocities, polarization, gaps, group quantities,
enhancement, and diagnostic masks where available.

Plotting is a post-processing stage.  It can produce:

- equal-area or stereographic 2D maps;
- a six-panel summary;
- phase, slowness, or group 3D surfaces;
- unit-sphere or physical geometry;
- optional tracked-polarization overlays.

Plot ``levels`` control contour rendering only.  ``polarization_stride``
reduces visual clutter by plotting a subset of already calculated axes; it does
not change tracking or field resolution.  ``complete-surface`` uses antipodal
symmetry to display a full surface from a hemispherical result.

Interoperability
----------------

SEISMIC consumes the same stiffness-plus-density text format produced by a
Thermoelasticity Point analysis.  This makes it possible to reconstruct an
adiabatic tensor at one pressure and temperature and immediately calculate its
directional wave field.

For comparison with the standalone Elasticity result, remember:

- VRH moduli and stability should agree for the same tensor;
- SEISMIC additionally requires density;
- isotropic velocities come from the Hill aggregate;
- directional velocities come from the full Christoffel tensor;
- tensor rotations must be applied consistently in both workflows.

Decision guide
--------------

.. list-table:: Practical SEISMIC decisions
   :header-rows: 1
   :widths: 37 27 36

   * - Scientific objective
     - Recommended settings
     - Reason
   * - Phase velocities and splitting only
     - ``level=phase``
     - Lowest cost and smallest HDF5 result.
   * - Ray paths or power-flow angles
     - ``level=group``
     - Adds analytical first derivatives without Hessian cost.
   * - Focusing or caustic search
     - ``level=enhancement``
     - Required for area factor and enhancement.
   * - First check of a new material
     - upper, 31 × 61 or 61 × 121, phase
     - Rapidly validates tensor, density, and broad anisotropy.
   * - Publication velocity maps
     - converge at least 61 × 121 against 91 × 181
     - Demonstrates sampled extrema and map stability.
   * - Caustic candidates
     - refine beyond the phase-converged grid
     - Curvature features converge more slowly than velocities.
   * - Polarization plots
     - tracking enabled
     - Preserves axial continuity and records branch swaps.
   * - Scalar fields only
     - tracking disabled
     - Saves post-processing and storage without changing speeds.
   * - Limited memory
     - lower ``batch_size``
     - Bounds temporary arrays without changing resolution.
   * - Faster machine with ample memory
     - benchmark larger batches
     - May reduce overhead; results should remain unchanged.
   * - Full-sphere calculation
     - double ``ntheta`` relative to upper for similar polar spacing
     - Same point count over a wider domain otherwise lowers resolution.
   * - Unexpected invalid eigenvalues
     - inspect tensor/stability; do not loosen tolerance first
     - Large negative values are physical/input problems, not rounding noise.

Related documentation
---------------------

- Scientific theory: :doc:`../theory/seismic`
- Input format: :doc:`../formats/elasticity_input`
- Detailed CLI/API example: :doc:`../tutorials/seismic`
- Public Python API: :doc:`../api/seismic`
- Command reference: :doc:`../cli/seismic`
- Tensor analysis: :doc:`elasticity`
- Thermoelastic Point interoperability: :doc:`thermoelasticity`
