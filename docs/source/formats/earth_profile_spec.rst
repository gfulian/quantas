Earth pressure-temperature profile specification
================================================

Thermoelastic properties can be evaluated along a geological depth path rather
than over a rectangular pressure-temperature grid.  Quantas represents such a
path as three aligned one-dimensional arrays:

.. math::

   z_i,\qquad P(z_i),\qquad T(z_i).

Pressure and temperature may be supplied directly in a table or generated from
independent physical models in a YAML specification.  The specification format
is useful when provenance, layers, boundaries, and joining policies must remain
explicit and reproducible.

Generate a fully commented starting point with::

   quantas thermoelasticity profile-template my_profile.yaml

The generated numerical values are examples, not a site-specific geotherm.
Review every density, conductivity, heat-production value, boundary condition,
and citation before quantitative use.

Units and general validation
----------------------------

Earth-profile inputs use fixed physical units:

=============================  ==================
Quantity                       Unit
=============================  ==================
Depth                          km
Pressure                       GPa
Temperature                    K
Density                        kg m\ :sup:`-3`
Gravity                        m s\ :sup:`-2`
Thermal conductivity           W m\ :sup:`-1` K\ :sup:`-1`
Volumetric heat production     microW m\ :sup:`-3`
Surface heat flow              mW m\ :sup:`-2`
Thermal diffusivity            m\ :sup:`2` s\ :sup:`-1`
Plate age                      Ma
=============================  ==================

Depth values must be finite, non-negative, unique, and strictly increasing after
sorting.  Pressure and absolute temperature must be finite and non-negative.
Every evaluated pressure and temperature array must have the same length as the
depth array.

A profile is a prescribed path, not a mineral-stability calculation.  Quantas
does not infer phase transitions, chemical reactions, or whether the selected
material remains stable along the full path.

Complete depth-pressure-temperature table
-----------------------------------------

The ``--profile`` option accepts comma-separated or whitespace-separated data.
Blank lines and lines beginning with ``#`` are ignored.  A header is required.
The canonical fields are::

   depth_km  P_GPa  T_K
   0         0.000  288.15
   35        1.05   800.0
   100       3.10   1500.0

Accepted aliases are:

================  ================================================
Field             Accepted names
================  ================================================
Depth             ``depth_km``, ``depth``, ``z_km``
Pressure          ``P_GPa``, ``pressure_GPa``, ``pressure``, ``P``
Temperature       ``T_K``, ``temperature_K``, ``temperature``, ``T``
================  ================================================

Rows are sorted by depth using a stable ordering.  Duplicate depths remain an
error because pressure and temperature would otherwise be ambiguous at the same
coordinate.  Additional columns are ignored by the simple profile reader; keep
ancillary geological labels in a separate provenance file or in the YAML model
metadata.

Composed YAML schema
--------------------

``--profile-spec`` accepts schema version 1.  The top-level object is a mapping
with four sections:

.. code-block:: yaml

   schema_version: 1
   name: prem-custom-temperature

   depth:
     min_km: 0.0
     max_km: 100.0
     step_km: 1.0
     include_critical_depths: true

   pressure:
     model: prem

   temperature:
     source: table
     file: custom_temperature.dat
     interpolation: pchip
     citation: >-
       Complete bibliographic citation or description of the user model.

``name`` identifies the evaluated profile in HDF5 and tabular exports.  It must
not be empty.  Relative table paths are resolved from the directory containing
the YAML specification, making a profile directory relocatable as a unit.

Depth grid
~~~~~~~~~~

The ``depth`` mapping defines a regular base grid:

``min_km``
   Lower depth bound.  If omitted, the common lower bound of the selected
   pressure and temperature models is used.

``max_km``
   Upper depth bound.  If omitted, the common upper bound is used.

``step_km``
   Positive regular spacing.  A smaller value increases the number of evaluated
   states but does not add physical information to a coarse tabulated source.

``include_critical_depths``
   When true, model boundaries and original tabulated knots are inserted in
   addition to the regular grid.  This prevents a regular step from skipping a
   layer boundary, piecewise join, or source observation.

Pressure and temperature are built independently and then evaluated on this
shared grid.  Their model domains must overlap over the requested depth range.

Tabulated interpolation
~~~~~~~~~~~~~~~~~~~~~~~

Pressure and temperature tables support:

``linear``
   Piecewise-linear interpolation.  It is local and does not overshoot the
   interval endpoints, but has discontinuous first derivatives at knots.

``pchip``
   Shape-preserving piecewise cubic Hermite interpolation.  This is the default
   for tabulated profile components.  It preserves monotonic trends better than
   an unconstrained cubic spline and provides smoother first derivatives.

Interpolation is not extrapolation.  The requested profile depth range must be
inside the tabulated model domain.  Adding a denser destination grid changes
only sampling of the interpolant.

Pressure models
---------------

PREM pressure
~~~~~~~~~~~~~

.. code-block:: yaml

   pressure:
     model: prem
     max_depth_km: 2891.0
     integration_step_km: 0.25

The model integrates hydrostatic pressure using PREM density.  The
``integration_step_km`` controls the internal numerical integration, not the
output depth spacing.  Reduce it only after a convergence comparison of
pressure at representative depths.  ``max_depth_km`` cannot exceed the model's
supported mantle interval.

Layered lithostatic pressure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   pressure:
     model: layered-lithostatic
     name: local-crustal-column
     gravity_m_s2: 9.80665
     citation: Complete source for the selected densities.
     layers:
       - name: upper crust
         thickness_km: 15.0
         density_kg_m3: 2700.0
       - name: lower crust
         thickness_km: 20.0
         density_kg_m3: 2950.0

Layers are contiguous and begin at zero depth.  Pressure is integrated using the
specified constant density in each layer and constant gravity for the complete
column.  Layer thicknesses must be positive and densities must be positive.
This model is appropriate for an explicitly parameterized local column, not as
a replacement for a self-consistent spherical Earth model at mantle depth.

Tabulated pressure
~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   pressure:
     source: table
     name: experimental-pressure-calibration
     file: pressure.dat
     interpolation: pchip
     citation: Complete source for the pressure table.

The table must contain depth and one accepted pressure column.  The original
file path, interpolation kind, model name, and optional citation are retained in
profile metadata.

Temperature models
------------------

Continental conductive geotherm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   temperature:
     model: continental-conductive
     name: local-continental-model
     surface_temperature_K: 288.15
     surface_heat_flow_mW_m2: 55.0
     citation: Complete source for the selected thermal parameters.
     layers:
       - name: upper crust
         thickness_km: 15.0
         conductivity_W_mK: 2.5
         heat_production_uW_m3: 0.8
       - name: lower crust
         thickness_km: 20.0
         conductivity_W_mK: 2.5
         heat_production_uW_m3: 0.4
       - name: lithospheric mantle
         thickness_km: 85.0
         conductivity_W_mK: 3.3
         heat_production_uW_m3: 0.02

The model solves one-dimensional steady conduction through contiguous layers.
Conductivity and thickness must be positive; heat production must be
non-negative.  The surface heat-flow boundary condition is applied at zero
depth.  The total model domain is the sum of the layer thicknesses.

Oceanic half-space cooling
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   temperature:
     model: oceanic-half-space
     age_Ma: 50.0
     surface_temperature_K: 273.15
     mantle_temperature_K: 1623.15
     diffusivity_m2_s: 1.0e-6
     max_depth_km: 300.0

This model treats the lithosphere as a cooling semi-infinite half-space.  Age,
diffusivity, and maximum depth must be positive.  It does not include a finite
plate thickness or basal boundary condition.

Oceanic plate cooling
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   temperature:
     model: oceanic-plate
     age_Ma: 50.0
     plate_thickness_km: 125.0
     surface_temperature_K: 273.15
     mantle_temperature_K: 1623.15
     diffusivity_m2_s: 1.0e-6
     series_terms: 200

The finite-plate solution adds a basal-temperature condition.  ``series_terms``
controls truncation of the analytical series.  It changes numerical convergence,
not output depth resolution.  Compare representative temperatures before
increasing it substantially.  ``max_depth_km`` may be supplied explicitly; if
omitted, the model uses its physical plate domain.

Katsura mantle adiabat
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   temperature:
     model: katsura-2022
     transition_epsilon_km: 0.001

The accepted model domain is 50--2800 km.  ``transition_epsilon_km`` controls
sampling immediately around internal transition depths so that the chosen side
of a discontinuity is unambiguous.  It is not a smoothing width.

Parametric basal thermal boundary layer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   temperature:
     model: linear-boundary-layer
     name: user-cmb-boundary-layer
     depth_top_km: 2800.0
     depth_bottom_km: 2891.0
     temperature_top_K: 2587.0
     temperature_bottom_K: 4000.0
     exponent: 1.0
     citation: Complete source for the chosen CMB temperature and shape.

This is a user-parameterized transition between two boundary temperatures.
Quantas provides no hidden canonical core-mantle-boundary temperature.  The
bottom depth must exceed the top depth; temperatures must be non-negative; the
exponent must be positive.

Tabulated temperature
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   temperature:
     source: table
     file: temperature.dat
     interpolation: linear
     citation: Complete source for the temperature table.

The table must contain depth and one accepted temperature column.  Use this
format for an externally calculated geotherm whose numerical values must be
preserved independently of Quantas model parameters.

Piecewise temperature profiles
------------------------------

A piecewise model combines any supported temperature models.  Segments must be
ordered, contiguous, non-overlapping, and have positive depth extent.  The first
segment is evaluated directly.  Every later segment declares how it joins the
previous one.

.. code-block:: yaml

   temperature:
     model: piecewise
     name: continental-to-mantle
     segments:
       - depth_min_km: 0.0
         depth_max_km: 80.0
         source: table
         file: lithosphere_temperature.dat
         interpolation: pchip
         citation: Complete lithosphere-profile source.

       - depth_min_km: 80.0
         depth_max_km: 2800.0
         model: katsura-2022
         join:
           mode: continuous-offset

The allowed join modes are:

``direct``
   Preserve both model values without an offset.  A real temperature jump may
   therefore remain at the boundary.

``continuous-offset``
   Add a constant offset to the new segment so that its first value equals the
   preceding segment at the boundary.  The applied offset is retained in
   metadata.

``blend``
   Blend the two segments across a finite interval.  ``width_km`` is required
   and must fit inside the neighboring segment domains.  The blend interval is
   recorded in metadata.

No join transformation is silent.  Do not use continuity or blending merely to
hide a physically meaningful thermal or phase boundary.

Provenance and citations
------------------------

User-defined layered and tabulated models should include a complete ``citation``
or an equally precise description of their origin.  Built-in scientific models
retain their own canonical references.  The evaluated profile metadata records:

* profile schema and name;
* source specification path;
* selected pressure and temperature model kinds;
* table paths and interpolation methods;
* layer parameters and boundary conditions;
* critical depths and join transformations;
* user-supplied citations.

Preserve the YAML file and every referenced table together with the
thermoelastic HDF5 result.  A profile exported only as three numerical columns
cannot reproduce model choices such as series truncation, internal integration
step, or piecewise joining.

CLI examples
------------

List the built-in scientific presets::

   quantas thermoelasticity analysis profile fit.hdf5 --list-presets

Evaluate a scientific preset::

   quantas thermoelasticity analysis profile fit.hdf5 \
      --preset mantle-katsura-2022 \
      --extrapolation warn \
      --output mantle_path.hdf5

Evaluate a complete table::

   quantas thermoelasticity analysis profile fit.hdf5 \
      --profile profile.dat \
      --extrapolation warn \
      --output tabulated_path.hdf5

Evaluate a composed user profile::

   quantas thermoelasticity analysis profile fit.hdf5 \
      --profile-spec profile.yaml \
      --extrapolation warn \
      --output custom_path.hdf5

Before using the result, inspect both QHA and elastic extrapolation masks.  A
geologically plausible path can still leave the numerical support of the
calibrated material model.
