Thermoelasticity post-fit analysis and export API
=================================================

Post-fit analysis evaluates the calibrated model at requested states.  A grid
analysis and a profile-only analysis are separate operations because a
generic P--T rectangle should not be constructed when only a geological path
is required.

Grid construction and analysis
------------------------------

.. autofunction:: quantas.api.thermoelasticity.regular_grid

.. autofunction:: quantas.api.thermoelasticity.analyze_grid

.. autofunction:: quantas.api.thermoelasticity.analyze_profiles

Profiles
--------

.. autofunction:: quantas.api.thermoelasticity.profile_presets

.. autofunction:: quantas.api.thermoelasticity.build_profile_preset

.. autofunction:: quantas.api.thermoelasticity.read_depth_profile

.. autofunction:: quantas.api.thermoelasticity.read_profile_spec

.. autofunction:: quantas.api.thermoelasticity.write_profile_template

Tensor selection and report tables
----------------------------------

.. autofunction:: quantas.api.thermoelasticity.select_stiffness_tensor

.. autofunction:: quantas.api.thermoelasticity.point_table

.. autofunction:: quantas.api.thermoelasticity.grid_info_table

Tabular and interoperable exports
---------------------------------

.. autofunction:: quantas.api.thermoelasticity.write_grid_table

.. autofunction:: quantas.api.thermoelasticity.write_profile_table

.. autofunction:: quantas.api.thermoelasticity.write_state_input

.. autofunction:: quantas.api.thermoelasticity.write_tensor_export
