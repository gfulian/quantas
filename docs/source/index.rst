Quantas documentation
=====================

.. image:: _static/branding/Quantas_banner-1600.png
   :alt: Quantas banner
   :class: hero-banner

**Quantas** stands for **Quantitative Analysis of Solids**. It is a Python
library and scientific application for the analysis of solid-state
thermodynamics, elasticity, seismic-wave propagation, equations of state, and
thermoelastic properties.

Quantas 2.0 exposes the same scientific workflows through a reusable Python API
and a Click/Rich command-line interface. Frontend-neutral results, reports, and
plot descriptions are designed to support Quantas GUI and other graphical
frontends without changing the numerical core.

Start here
----------

- New users: :doc:`introduction/overview` and
  :doc:`getting_started/installation`.
- Researchers choosing a method: :doc:`introduction/capabilities` and
  :doc:`theory/index`.
- Command-line users: :doc:`getting_started/command_line`.
- Python users: :doc:`getting_started/python_api`.
- Contributors: :doc:`developer/index`.

.. toctree::
   :maxdepth: 2
   :caption: INTRODUCTION
   :hidden:

   introduction/overview
   introduction/capabilities
   introduction/changes_from_0_9
   introduction/citing_quantas

.. toctree::
   :maxdepth: 2
   :caption: GETTING STARTED
   :hidden:

   getting_started/installation
   getting_started/first_calculation
   getting_started/python_api
   getting_started/command_line
   getting_started/results

.. toctree::
   :maxdepth: 2
   :caption: SCIENTIFIC BACKGROUND
   :hidden:

   theory/ha
   theory/qha
   theory/elasticity
   theory/seismic
   theory/eos
   theory/thermoelasticity
   theory/earth_profiles

.. toctree::
   :maxdepth: 3
   :caption: IMPLEMENTATION AND WORKFLOWS
   :hidden:

   workflows/concepts/index
   workflows/phonon_input_generation
   workflows/ha
   workflows/qha
   workflows/elasticity
   workflows/seismic
   workflows/eos
   workflows/thermoelasticity
   workflows/interoperability

.. toctree::
   :maxdepth: 3
   :caption: TUTORIALS
   :hidden:

   tutorials/ha
   tutorials/qha
   tutorials/elasticity
   tutorials/seismic
   tutorials/eos
   tutorials/thermoelasticity

.. toctree::
   :maxdepth: 2
   :caption: INPUT AND OUTPUT FORMATS
   :hidden:

   formats/index
   formats/elasticity_input
   formats/phonon_yaml
   formats/eos_input
   formats/eos_spec
   formats/thermoelastic_input
   formats/earth_profile_spec
   formats/hdf5
   formats/ha_qha_hdf5
   formats/elasticity_seismic_hdf5
   formats/eos_hdf5
   formats/thermoelastic_hdf5
   formats/hdf5_inspection
   formats/tabular_outputs

.. toctree::
   :maxdepth: 2
   :caption: COMMAND REFERENCE
   :hidden:

   cli/index
   cli/conventions
   cli/ha
   cli/qha
   cli/elasticity
   cli/seismic
   cli/eos
   cli/thermoelasticity

.. toctree::
   :maxdepth: 2
   :caption: API REFERENCE
   :hidden:

   api/index
   api/common
   api/plotting
   api/ha
   api/qha
   api/elasticity
   api/seismic
   api/eos
   api/thermoelasticity
   api/rendering
   api/profiles
   api/interoperability
   api/registry

.. toctree::
   :maxdepth: 2
   :caption: SCIENTIFIC VALIDATION
   :hidden:

   validation/strategy
   validation/matrix
   validation/ha_qha
   validation/elasticity
   validation/seismic
   validation/eos
   validation/thermoelasticity
   validation/precision

.. toctree::
   :maxdepth: 2
   :caption: DEVELOPMENT GUIDE
   :hidden:

   developer/getting_started
   developer/architecture
   developer/module_anatomy
   developer/numerical_methods
   developer/public_api
   developer/change_recipes
   developer/extending
   developer/interfaces
   developer/events
   developer/persistence
   developer/rendering_frontends
   developer/citation_registry
   developer/testing
   developer/documentation
   developer/thermoelastic_architecture
   developer/packaging
   developer/review_checklist
