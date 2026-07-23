Thermoelastic HDF5 payload
==========================

Thermoelasticity uses the common Quantas envelope and stores its scientific
payload under ``/results``.  The module-local payload schema is currently
``1.0`` and is recorded as the ``thermoelastic_schema_version`` attribute.

A thermoelastic file can represent:

- a calibration result containing the cold reference EOS and component fits;
- a reconstructed pressure--temperature grid;
- one or more named geological profiles;
- any combination of the above when later analyses have been appended to the
  result object before writing a new file.

Root payload
------------

.. code-block:: text

   /results
       attributes: thermoelastic_schema_version, jobname, completed
       /independent_labels
       /temperature
       /pressure
       /equilibrium_volume
       /density
       /elastic_extrapolation
       /qha_extrapolation
       ... optional grid arrays
       /reference_eos
       /component_fits/<label>
       /profiles/<name>
       /mechanical_stability                 # optional
       /metadata

Grid axis order
---------------

Let ``nT`` be the number of temperatures, ``nP`` the number of pressures, and
``nC`` the number of independent elastic components.  Grid arrays use axis
order:

.. code-block:: text

   (temperature, pressure, ...)

The coordinate arrays are one-dimensional.  A property evaluated at all states
uses shape ``(nT, nP)`` followed by component or tensor axes.

Core grid datasets
------------------

.. list-table::
   :header-rows: 1
   :widths: 39 25 36

   * - Dataset
     - Shape
     - Unit and meaning
   * - ``temperature``
     - ``(nT,)``
     - K; archived temperature axis.
   * - ``pressure``
     - ``(nP,)``
     - GPa; archived pressure axis.
   * - ``equilibrium_volume``
     - ``(nT, nP)``
     - Angstrom cubed; QHA equilibrium volume.
   * - ``sigma_equilibrium_volume``
     - ``(nT, nP)`` optional
     - One-sigma volume uncertainty.
   * - ``density``
     - ``(nT, nP)``
     - kg m\ :sup:`-3`.
   * - ``elastic_extrapolation``
     - ``(nT, nP)`` boolean
     - Equilibrium volume outside the calibrated elastic-volume interval.
   * - ``qha_extrapolation``
     - ``(nT, nP)`` boolean
     - Requested state outside the archived QHA coordinate grid.

The two extrapolation masks describe different support failures and must not be
combined without preserving their meaning.

Independent and full stiffness fields
-------------------------------------

``independent_labels`` is a UTF-8 array of length ``nC`` defining the last axis
of independent-component arrays.

.. list-table::
   :header-rows: 1
   :widths: 42 28 30

   * - Dataset
     - Shape
     - Meaning
   * - ``independent_stiffness``
     - ``(nT, nP, nC)``
     - Independent isothermal components in GPa.
   * - ``sigma_independent_stiffness``
     - ``(nT, nP, nC)``
     - One-sigma estimates in GPa.
   * - ``independent_stiffness_covariance``
     - ``(nT, nP, nC, nC)``
     - Full covariance among independent components in GPa squared.
   * - ``stiffness_isothermal``
     - ``(nT, nP, 6, 6)``
     - Full symmetry-constrained Wallace tensor in GPa.
   * - ``sigma_stiffness_isothermal``
     - ``(nT, nP, 6, 6)``
     - Element-wise one-sigma estimates in GPa.

A missing optional array means that the file represents calibration only or
that the requested analysis did not calculate that quantity.  Absence does not
mean zero.

Reference EOS
-------------

``/results/reference_eos`` stores the shared cold EOS used by all component
fits.

Attributes include:

- ``eos``;
- ``reference_volume_A3``;
- ``bulk_modulus_GPa``;
- ``bulk_modulus_derivative``;
- ``bulk_modulus_second_derivative_GPa-1``.

Optional ``covariance_V0_K0_Kprime`` stores the covariance of the parameters
actually propagated by the component model.  The ``fit`` child contains the
complete generic fit result and ``metadata`` stores scientific provenance and
normalization details.

Component fits
--------------

Every independent component appears below:

.. code-block:: text

   /results/component_fits/C11
       attributes: active, zero_by_tolerance, wallace_delta
       /entries
       /volumes
       /pressures
       /observed
       /fitted
       /residuals
       /relative_residuals
       /symmetry_spread
       /fit                              # absent for exact-zero components
       /quality                          # when classified
       /metadata

``entries`` records the tensor positions and signs contributing to the
symmetry-reduced component.  ``observed`` and ``fitted`` are in GPa.
``volumes`` and ``pressures`` preserve the sampled static states.

``zero_by_tolerance = True`` means that the component was represented as exactly
zero rather than passed to a numerical fit.  It is not a failed fit.

Fit support and quality
-----------------------

``quality`` may contain:

- classification and machine-readable issue identifiers;
- finite-strain span and reference-volume bracketing;
- design rank and normalized condition number;
- point leverage;
- symmetry spread;
- leave-one-out sensitivity;
- order sensitivity;
- thresholds used for classification.

These fields describe support.  They do not alter the fitted parameters.
External analyses should retain both the numerical fit and its quality
classification.

Adiabatic fields
----------------

When the QHA source supplies all required thermodynamic fields and adiabatic
conversion is enabled, the root payload may contain:

.. list-table::
   :header-rows: 1
   :widths: 43 27 30

   * - Dataset
     - Shape
     - Meaning
   * - ``isochoric_heat_capacity_cell``
     - ``(nT, nP)``
     - J cell\ :sup:`-1` K\ :sup:`-1`.
   * - ``sigma_isochoric_heat_capacity_cell``
     - ``(nT, nP)``
     - One-sigma heat-capacity estimate.
   * - ``thermal_expansion_tensor``
     - ``(nT, nP, 3, 3)``
     - Cartesian expansion tensor in K\ :sup:`-1`.
   * - ``sigma_thermal_expansion_tensor``
     - ``(nT, nP, 3, 3)``
     - One-sigma tensor estimate.
   * - ``stiffness_adiabatic``
     - ``(nT, nP, 6, 6)``
     - Adiabatic stiffness in GPa.
   * - ``sigma_stiffness_adiabatic``
     - ``(nT, nP, 6, 6)``
     - One-sigma estimate.
   * - ``adiabatic_correction``
     - ``(nT, nP, 6, 6)``
     - Stored ``C_S - C_T`` correction in GPa.
   * - ``adiabatic_thermal_stress``
     - ``(nT, nP, 6)``
     - Voigt vector ``lambda = C_T alpha`` in GPa K\ :sup:`-1`.
   * - ``adiabatic_valid_mask``
     - ``(nT, nP)`` boolean
     - Validity of the conversion.

At non-zero temperature, invalid adiabatic states are represented by ``False``
in the validity mask and NaN-valued arrays.  They are not silently replaced by
isothermal values.

Mechanical stability
--------------------

A reconstructed grid may contain ``/results/mechanical_stability`` with:

- ``eigenvalues`` of shape ``(nT, nP, 6)`` in GPa;
- ``minimum_eigenvalue`` of shape ``(nT, nP)``;
- ``stable_mask``;
- ``indeterminate_mask``;
- ``criterion`` and ``tolerance_GPa`` attributes;
- a ``metadata`` mapping describing tensor convention and references.

The stability field applies to the reconstructed isothermal Wallace tensor.
Read it together with extrapolation masks.

Depth-profile payloads
----------------------

Named profiles are stored under ``/results/profiles/<name>``.  A slash is not
allowed in the profile name because it would change the HDF5 hierarchy.

For ``nD`` depth points, required arrays are:

.. list-table::
   :header-rows: 1
   :widths: 43 26 31

   * - Dataset
     - Shape
     - Meaning
   * - ``depth``
     - ``(nD,)``
     - Depth in km.
   * - ``pressure``
     - ``(nD,)``
     - Pressure in GPa.
   * - ``temperature``
     - ``(nD,)``
     - Temperature in K.
   * - ``volume`` and ``density``
     - ``(nD,)``
     - Equilibrium volume and density.
   * - ``independent_stiffness``
     - ``(nD, nC)``
     - Independent components.
   * - ``independent_stiffness_covariance``
     - ``(nD, nC, nC)``
     - Component covariance.
   * - ``stiffness_isothermal``
     - ``(nD, 6, 6)``
     - Full tensor along the path.
   * - extrapolation masks
     - ``(nD,)``
     - QHA-coordinate and elastic-volume support.

Adiabatic and mechanical-stability groups may appear inside each profile using
the corresponding one-dimensional leading shape.

Profile ``metadata`` stores model names, citations, joins, offsets, interpolation
choices, and other provenance from the profile specification.

Explicit absence markers
------------------------

Optional thermoelastic arrays use explicit ``<name>_is_none`` attributes when
absent.  This distinguishes a quantity that was not available from an empty
array.  The public reader reconstructs these markers as ``None``.

Typed reading
-------------

.. code-block:: python

   import numpy as np
   from quantas.api import thermoelasticity

   envelope = thermoelasticity.read_result("thermoelastic.hdf5")
   data = thermoelasticity.get_result(envelope)

   print(data.reference_eos.eos)
   print(data.independent_labels)
   print(data.component_fits.keys())

   if data.stiffness_isothermal is not None:
       c11 = np.asarray(data.stiffness_isothermal)[..., 0, 0]
       supported = ~np.asarray(data.extrapolation_mask)
       print(np.nanmin(np.where(supported, c11, np.nan)))

   for name, profile in data.profiles.items():
       print(name, profile.depth.shape)

Related pages
-------------

- :doc:`thermoelastic_input`
- :doc:`earth_profile_spec`
- :doc:`hdf5`
- :doc:`hdf5_inspection`
- :doc:`tabular_outputs`
