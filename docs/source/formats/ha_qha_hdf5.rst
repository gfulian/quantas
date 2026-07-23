HA and QHA native HDF5 payloads
===============================

HA and QHA use the common Quantas result envelope described in :doc:`hdf5`.
Their module-specific arrays are stored under ``/results``.  The current common
envelope schema is ``2.0``.

The numerical basis stored by the thermodynamic workflows is the normalization
cell declared by the input.  Energies are not silently converted to a molar or
per-atom basis in native HDF5.  Formula-unit and atom normalization metadata are
stored with the result so that renderers and external scripts can perform
explicit conversions.

HA payload
----------

The HA result tree is compact:

.. code-block:: text

   /results
       attributes: jobname
       /metadata
       /volume
       /temperature
       /static_energy
       /zero_point_energy
       /thermal_energy
       /internal_energy
       /entropy
       /vibrational_free_energy
       /free_energy
       /isochoric_heat_capacity

Let ``nV`` be the number of sampled volumes and ``nT`` the number of requested
temperatures.

.. list-table:: HA datasets
   :header-rows: 1
   :widths: 30 20 50

   * - Dataset
     - Shape
     - Meaning
   * - ``volume``
     - ``(nV,)``
     - Sampled normalization-cell volumes.
   * - ``temperature``
     - ``(nT,)``
     - Temperature grid in K.
   * - ``static_energy``
     - ``(nV,)``
     - Static electronic energies.
   * - ``zero_point_energy``
     - ``(1, nV)``
     - Temperature-independent zero-point contribution.  The leading length-one
       axis is retained deliberately rather than broadcasting stored values.
   * - ``thermal_energy``
     - ``(nT, nV)``
     - Thermal vibrational contribution.
   * - ``internal_energy``
     - ``(nT, nV)``
     - Total internal energy.
   * - ``entropy``
     - ``(nT, nV)``
     - Vibrational entropy in native energy per cell per K.
   * - ``vibrational_free_energy``
     - ``(nT, nV)``
     - Vibrational Helmholtz free energy.
   * - ``free_energy``
     - ``(nT, nV)``
     - Total Helmholtz free energy.
   * - ``isochoric_heat_capacity``
     - ``(nT, nV)``
     - Constant-volume heat capacity in native energy per cell per K.

The ``/results/metadata/units`` mapping records the actual energy, volume,
frequency, entropy, and heat-capacity units used by the run.  HA payload arrays
from the current writer do not all carry local ``unit`` attributes, so external
scripts should read this metadata mapping.

QHA payload
-----------

The QHA payload adds pressure, equilibrium fields, structural data, uncertainty
arrays, and fit diagnostics:

.. code-block:: text

   /results
       attributes: jobname, completed
       /metadata
       /volume
       /static_energy
       /temperature
       /pressure
       /equilibrium_volume
       ... thermodynamic grid arrays
       /valid_mask
       /uncertainties/<property>

   /diagnostics
       /fit_records/0, /1, ...
       /failed_points/0, /1, ...

Let ``nP`` be the number of pressures.  Unless noted otherwise, equilibrium
thermodynamic fields use axis order:

.. code-block:: text

   (temperature, pressure)

Core coordinates and state fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 31 20 49

   * - Dataset
     - Shape
     - Meaning
   * - ``volume``
     - ``(nV,)``
     - Sampled static volumes used to construct the free-energy surface.
   * - ``static_energy``
     - ``(nV,)``
     - Static energies at the sampled volumes.
   * - ``temperature``
     - ``(nT,)``
     - Temperature axis.
   * - ``pressure``
     - ``(nP,)``
     - Pressure axis in GPa.
   * - ``equilibrium_volume``
     - ``(nT, nP)``
     - Equilibrium volume at each state.
   * - ``valid_mask``
     - ``(nT, nP)``
     - Whether the state was completed successfully.

Thermodynamic grid arrays
~~~~~~~~~~~~~~~~~~~~~~~~~

When calculated, the following fields normally have shape ``(nT, nP)``:

- ``zero_point_energy``;
- ``thermal_energy``;
- ``internal_energy``;
- ``entropy``;
- ``vibrational_free_energy``;
- ``free_energy``;
- ``isochoric_heat_capacity``;
- ``isobaric_heat_capacity``;
- ``heat_capacity_difference``;
- ``isothermal_bulk_modulus``;
- ``adiabatic_bulk_modulus``;
- ``bulk_modulus_derivative``;
- ``thermal_expansion``;
- ``thermal_expansion_mixed``;
- ``thermal_expansion_mode``;
- ``thermal_expansion_numerical``;
- ``enthalpy``;
- ``gibbs_free_energy``;
- ``gruneisen``;
- ``mode_weighted_gruneisen``.

``thermal_expansion_source`` has the same grid shape but stores the
machine-readable source code used for the selected alpha value at each state.
Interpret the code through QHA metadata rather than assuming an undocumented
integer mapping.

Mode-resolved Gruneisen data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``mode_gruneisen`` is not a pressure--temperature grid.  It stores
mode-resolved values on the sampled volume path.  Its exact leading dimensions
depend on the available q-point and band representation; the current result
model preserves the native sampled mode grid rather than forcing it into the
``(nT, nP)`` convention.

Always inspect its shape and the stored input metadata before reducing or
plotting it.

Structural fields
~~~~~~~~~~~~~~~~~

When the input contains a structural volume path, QHA may store:

.. list-table::
   :header-rows: 1
   :widths: 35 25 40

   * - Dataset
     - Shape
     - Meaning
   * - ``equilibrium_lattice``
     - ``(nT, nP, 3, 3)``
     - Direct-lattice vectors by rows.
   * - ``lattice_parameters``
     - ``(nT, nP, 6)``
     - ``a, b, c, alpha, beta, gamma``.
   * - ``lattice_parameter_derivatives``
     - ``(nT, nP, 6)``
     - Temperature derivatives of lengths and angles.
   * - ``axial_thermal_expansion``
     - ``(nT, nP, 3)``
     - ``alpha_a, alpha_b, alpha_c``.
   * - ``thermal_expansion_tensor``
     - ``(nT, nP, 3, 3)``
     - Symmetric Cartesian expansion tensor.
   * - ``structural_extrapolation_mask``
     - ``(nT, nP)``
     - Equilibrium states outside the sampled structural path.

An absent structural dataset means that the source YAML lacked the required
path or that the quantity was not requested.  It must not be replaced by an
isotropic guess.

Uncertainties
~~~~~~~~~~~~~

``/results/uncertainties`` contains only properties for which uncertainty
propagation was performed.  Child names are stable property-oriented keys such
as ``sigma_VT`` or ``sigma_KT``.  Their shapes follow the corresponding result
arrays.

An empty uncertainty group is valid and means that no uncertainty arrays were
available.  It does not imply zero uncertainty.

Fit records and failed points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``/diagnostics/fit_records`` stores structured records for local polynomial or
EOS fits.  Each numbered child can contain:

- quantity and method;
- pressure and/or temperature coordinates;
- success flag and message;
- fitted parameters, covariance, residual diagnostics, and quality;
- additional metadata.

``/diagnostics/failed_points`` stores pressure--temperature states that could
not be completed.  The grid arrays preserve their rectangular shape and the
``valid_mask`` identifies usable states.

Completion state
~~~~~~~~~~~~~~~~

The ``completed`` attribute tells whether the workflow reached the end of the
requested grid.  A completed result may still contain invalid individual
states; an incomplete result may contain scientifically usable earlier states.
Read ``completed``, ``valid_mask``, warnings, and failed-point records together.

Typed reading examples
----------------------

HA:

.. code-block:: python

   from quantas.api import ha

   envelope = ha.read_result("ha.hdf5")
   data = ha.get_result(envelope)

   print(data.temperature.shape)
   print(data.zero_point_energy.shape)
   print(data.metadata["units"])

QHA:

.. code-block:: python

   from quantas.api import qha

   envelope = qha.read_result("qha.hdf5")
   data = qha.get_result(envelope)

   print(data.equilibrium_volume.shape)
   print(data.valid_mask.shape)
   print(data.uncertainties.keys())
   print(len(data.fit_records), len(data.failed_points))

Historical unit migration
-------------------------

Selected schema-1.0 pre-release HA/QHA files stored entropy and heat capacities
with an implicit milli-energy numerical scale while reporting the unprefixed
energy unit.  Current module readers apply an explicit migration to true native
energy per cell per K and append a warning.  A direct h5py reader does not apply
this correction automatically.

For archival analysis of historical files, use the current Quantas reader and
preserve the migration warning in derived provenance.

Related pages
-------------

- :doc:`phonon_yaml`
- :doc:`hdf5`
- :doc:`hdf5_inspection`
- :doc:`tabular_outputs`
