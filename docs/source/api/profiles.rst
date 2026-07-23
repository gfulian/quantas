Terrestrial profile API
=======================

:mod:`quantas.api.profiles` exposes passive pressure--temperature--depth
profiles independently of Thermoelasticity.  Profiles can be built from a
scientific preset, a mapping, or an Earth-profile YAML specification and then
passed to a workflow that consumes depth paths.

.. code-block:: python

   from quantas.api import profiles

   print(profiles.preset_names())
   profile = profiles.build_preset("continental-reference")
   print(profile.depth, profile.pressure, profile.temperature)

Profile contracts
-----------------

.. autoclass:: quantas.api.profiles.DepthProfile
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.profiles.Model
   :members:
   :show-inheritance:

.. autoclass:: quantas.api.profiles.Preset
   :members:
   :show-inheritance:

Preset discovery and construction
---------------------------------

.. autofunction:: quantas.api.profiles.preset_names

.. autofunction:: quantas.api.profiles.presets

.. autofunction:: quantas.api.profiles.build_preset

Custom specifications
---------------------

.. autofunction:: quantas.api.profiles.from_mapping

.. autofunction:: quantas.api.profiles.read_spec

See also
--------

- :doc:`../theory/earth_profiles`
- :doc:`../formats/earth_profile_spec`
- :doc:`../workflows/thermoelasticity`
- :doc:`../tutorials/thermoelasticity`
