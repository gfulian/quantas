Elasticity and SEISMIC native HDF5 payloads
===========================================

Elasticity and SEISMIC use the common Quantas envelope described in
:doc:`hdf5`.  Both store the normalized stiffness tensor under ``/results``;
SEISMIC additionally stores density and sampled wave fields.

Elasticity payload
------------------

.. code-block:: text

   /results
       attributes: jobname, crystal_system
       /stiffness
       /compliance
       /averages/{voigt,reuss,hill}
       /stability/eigenvalues
       /variations/<property>/...
       /properties_2d/{xy,xz,yz}/...
       /properties_3d/...                     # optional
       /metadata

Core tensors and isotropic estimates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 32 20 48

   * - Path
     - Shape
     - Meaning
   * - ``/results/stiffness``
     - ``(6, 6)``
     - Stiffness matrix in GPa and engineering Voigt convention.
   * - ``/results/compliance``
     - ``(6, 6)``
     - Inverse stiffness in reciprocal GPa.
   * - ``/results/averages/voigt``
     - ``(4,)``
     - ``K, E, G, nu`` in that order.
   * - ``/results/averages/reuss``
     - ``(4,)``
     - Reuss ``K, E, G, nu``.
   * - ``/results/averages/hill``
     - ``(4,)``
     - Hill ``K, E, G, nu``.
   * - ``/results/stability/eigenvalues``
     - ``(6,)``
     - Ordered stiffness eigenvalues in GPa.

The stability group also stores ``is_stable`` and the applied ``tolerance`` as
attributes.

Directional extrema
~~~~~~~~~~~~~~~~~~~

Each child of ``/results/variations`` corresponds to a directional property.
The group stores scalar attributes:

- ``minimum``;
- ``maximum``;
- ``anisotropy``;

and Cartesian direction datasets:

- ``minimum_axis``;
- ``maximum_axis``;
- optional transverse ``minimum_measurement_axis`` and
  ``maximum_measurement_axis``.

Longitudinal axes have shape ``(3,)``.  Measurement axes are present for
properties such as shear modulus and Poisson ratio that depend on a second
orthogonal direction.

Principal-plane fields
~~~~~~~~~~~~~~~~~~~~~~

When 2D calculation was requested, ``/results/properties_2d`` contains ``xy``,
``xz``, and ``yz`` groups.  Each plane stores:

``theta`` and ``phi``
   Angular coordinates with shape ``(nangle,)``.

``young_modulus``
   One value per angle, shape ``(nangle,)``.

``linear_compressibility``
   Two branches with shape ``(nangle, 2)``: positive and magnitude of negative
   compressibility.

``shear_modulus``
   Minimum and maximum transverse branches, shape ``(nangle, 2)``.

``poisson_ratio``
   Negative minimum, positive minimum, and maximum branches, shape
   ``(nangle, 3)``.

The branch order is part of the current payload contract and should be retained
in any derived table.

Persisted 3D surfaces
~~~~~~~~~~~~~~~~~~~~~

Three-dimensional fields are optional.  Plot-only calculations may reconstruct
them transiently without storing them.

.. code-block:: text

   /results/properties_3d
       attribute: schema_version = "1.0"
       /theta                         # (ntheta, nphi), rad
       /phi                           # (ntheta, nphi), rad
       /warnings
       /metadata
       /surfaces/<key>
           attributes: property_name, branch, unit
           /radius                   # (ntheta, nphi)
           /values                   # (ntheta, nphi)
           /metadata

All stored surfaces share the same angular grid.  Cartesian ``x``, ``y``, and
``z`` coordinates are reconstructed by the reader from ``radius``, ``theta``,
and ``phi``; they are not redundantly stored.

SEISMIC payload
---------------

.. code-block:: text

   /results
       attributes: jobname, density, density_unit, level, batch_size
       /stiffness
       /stability
       /averages
       /isotropic_reference
       /grid
       /fields/phase
       /fields/group                         # level >= group
       /fields/enhancement                   # level = enhancement
       /fields/tracking                      # when tracking enabled

   /diagnostics/seismic

The mode order is stored as an attribute and is currently:

.. code-block:: text

   V_S2, V_S1, V_P

Do not assume a different P/S ordering when slicing the final axis.

Base tensor and grid
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 36 22 42

   * - Path
     - Shape
     - Meaning
   * - ``/results/stiffness``
     - ``(6, 6)``
     - Wallace stiffness in GPa.
   * - ``/results/stability/eigenvalues``
     - ``(6,)``
     - Positive-definiteness diagnostic.
   * - ``/results/averages/{voigt,reuss,hill}``
     - ``(4,)``
     - ``K, E, G, nu``.
   * - ``/results/grid/theta``
     - ``(ntheta, nphi)``
     - Polar angle in radians.
   * - ``/results/grid/phi``
     - ``(ntheta, nphi)``
     - Azimuth in radians, without duplicated ``2*pi`` endpoint.
   * - ``/results/grid/directions``
     - ``(ntheta, nphi, 3)``
     - Unit wave-normal directions.

The grid group stores ``hemisphere``, ``ntheta``, ``nphi``, and
``azimuth_endpoint_included`` as attributes.

Phase fields
~~~~~~~~~~~~

``/results/fields/phase`` contains:

.. list-table::
   :header-rows: 1
   :widths: 38 24 38

   * - Dataset
     - Shape
     - Meaning
   * - ``eigenvalues``
     - ``(ntheta, nphi, 3)``
     - Christoffel eigenvalues in km\ :sup:`2` s\ :sup:`-2`.
   * - ``phase_speeds``
     - ``(ntheta, nphi, 3)``
     - Phase speeds in km s\ :sup:`-1`.
   * - ``polarizations``
     - ``(ntheta, nphi, 3, 3)``
     - Row-oriented polarization axes for every mode.
   * - ``mode_eigenvalue_gaps``
     - ``(ntheta, nphi, 3)``
     - Smallest adjacent gap associated with each mode.
   * - ``pair_eigenvalue_gaps``
     - ``(ntheta, nphi, 2)``
     - Adjacent pair gaps for ``S2-S1`` and ``S1-P``.
   * - relative gap arrays
     - matching mode/pair shape
     - Gaps normalized by the largest directional eigenvalue.
   * - ``valid_mask``
     - ``(ntheta, nphi, 3)``
     - Numerically valid modes.
   * - ``clamped_mask``
     - ``(ntheta, nphi, 3)``
     - Small negative eigenvalues clamped within tolerance.
   * - ``degeneracy_mask``
     - ``(ntheta, nphi, 3)``
     - Per-mode degeneracy state.
   * - ``pair_degeneracy_mask``
     - ``(ntheta, nphi, 2)``
     - Adjacent-pair acoustic-axis candidates.

Threshold arrays are stored explicitly so that external analyses can reproduce
the numerical classification applied at every point.

Group fields
~~~~~~~~~~~~

When ``level`` is ``group`` or ``enhancement``, ``/fields/group`` contains:

- ``eigenvalue_gradients`` with shape ``(..., 3, 3)``;
- ``group_velocities`` with shape ``(..., 3, 3)``;
- ``group_speeds`` with shape ``(..., 3)``;
- ``ray_directions`` with shape ``(..., 3, 3)``;
- ``power_flow_angles`` with shape ``(..., 3)``;
- per-mode ``valid_mask`` and ``resolved_mask``.

Here ``...`` means ``(ntheta, nphi)``.  The first trailing index is mode and the
last index is Cartesian component where present.

Enhancement fields
~~~~~~~~~~~~~~~~~~

At ``enhancement`` level, ``/fields/enhancement`` contains:

- eigenvalue Hessians ``(..., 3, 3, 3)``;
- ray-direction gradients ``(..., 3, 3, 3)``;
- area factors, caustic thresholds, and ``log10_enhancement`` with shape
  ``(..., 3)``;
- per-mode valid, resolved, finite, and caustic-candidate masks.

The stored enhancement representation is logarithmic and carries the attribute
``representation = log10(A)``.  Do not exponentiate invalid or non-finite
entries without applying the corresponding masks.

Polarization tracking
~~~~~~~~~~~~~~~~~~~~~

When tracking is enabled, ``/fields/tracking`` stores branch order
``shear_a, shear_b, p`` and:

- tracked ``polarizations`` with shape ``(..., 3, 3)``;
- ``branch_mode_indices`` with shape ``(..., 3)``;
- sign-flip, resolved, segment-start, and shear-swap masks;
- continuity scores;
- ambiguous-permutation and shear-subspace-rotation masks.

A local mode index of ``-1`` means that no unique local eigenvector can be
assigned inside a degenerate subspace.  It is not a missing phase speed.

Diagnostics metadata
~~~~~~~~~~~~~~~~~~~~

SEISMIC-specific metadata are stored under ``/diagnostics/seismic`` as a
recursive mapping.  This includes sampling and numerical diagnostics not
naturally represented as dense scientific arrays.

Typed reading examples
----------------------

Elasticity:

.. code-block:: python

   from quantas.api import elasticity

   result = elasticity.get_result(
       elasticity.read_result("elasticity.hdf5")
   )
   print(result.stiffness.shape)
   print(result.properties_2d.keys())
   print(result.has_3d_data())

SEISMIC:

.. code-block:: python

   from quantas.api import seismic

   result = seismic.get_result(seismic.read_result("seismic.hdf5"))
   print(result.grid.shape)
   print(result.field.phase.phase_speeds.shape)
   print(result.field.group is not None)
   print(result.field.enhancement is not None)

Reading masks before reductions
-------------------------------

For SEISMIC, a valid global minimum or maximum should normally be computed only
after applying:

- the phase ``valid_mask``;
- the group or enhancement ``resolved_mask`` when branch uniqueness matters;
- the enhancement ``finite_mask``;
- any scientific restriction on hemisphere or mode.

For Elasticity, check mechanical stability before interpreting directional
fields.  Persisted arrays from an unstable input may be incomplete by design.

Related pages
-------------

- :doc:`elasticity_input`
- :doc:`hdf5`
- :doc:`hdf5_inspection`
- :doc:`tabular_outputs`
