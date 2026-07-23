Thermoelastic architecture and maintenance
===========================================

Worked architectural case study
---------------------------------

This page applies the general rules from :doc:`architecture`,
:doc:`module_anatomy`, :doc:`persistence`, and
:doc:`rendering_frontends` to the most structurally complex current Quantas
workflow. It is not a substitute for those general contracts.

The thermoelastic workflow is split into four scientific stages. Every
frontend must call the same services and therefore obtain identical float64
results.

.. code-block:: text

   normalized CRYSTAL/QHA inputs
              |
              v
   calibration calculator  ---> calibration HDF5
              |
              v
   analysis engine: point | grid | profile
              |
              v
   report | table | tensor export | plot specification

Layer responsibilities
----------------------

``quantas.core``
   Finite strain, elasticity, adiabatic conversion, numerical interpolation,
   Earth-profile physics, and stability criteria.  No workflow or frontend
   dependency is permitted.

``quantas.modules.thermoelasticity``
   Calibration and analysis orchestration, passive contracts, scientific
   policies, HDF5 payload adapters, reports, exports, and frontend-neutral plot
   specifications.

``quantas.cli``
   Click options, prompts, dispatch, Rich presentation, and output naming only.

Calibration contract
--------------------

``run`` always calibrates the static reference EOS and the independent
:math:`C_{IJ}(V)` models.  It never reconstructs a P--T grid or a depth profile.
The EOS parameters :math:`V_0`, :math:`K_0`, and :math:`K'_0` come from the QHA
``static_energy`` versus volume field.  They represent the static 0 K,
zero-pressure reference and exclude thermal and zero-point contributions.

CRYSTAL calculations must use ``PRESSURE``.  The sampled matrices are already
Wallace stress--strain coefficients.  ``wallace_delta`` in the third-order
finite-strain polynomial is an analytical expansion coefficient; it is not a
second pressure correction applied to the observations.

Analysis contract
-----------------

``ThermoelasticAnalysisEngine`` owns point, rectangular-grid, geological-profile,
and constant-pressure/temperature-section evaluation.  All QHA fields use the
shared ``RectilinearFieldInterpolator``.  It supports nonuniform and singleton
axes, trailing tensor dimensions, float64 arithmetic, and explicit
extrapolation masks.

Passive data shapes
-------------------

.. list-table::
   :header-rows: 1

   * - Quantity
     - Shape
     - Unit
   * - QHA temperature, pressure
     - ``(nT,)``, ``(nP,)``
     - K, GPa
   * - equilibrium volume, density
     - ``(nT, nP)``
     - Å³ per normalized cell, kg m⁻³
   * - independent stiffness
     - ``(nT, nP, nC)``
     - GPa
   * - full stiffness/covariance-derived sigma
     - ``(nT, nP, 6, 6)``
     - GPa
   * - profile fields
     - leading dimension ``(nDepth,)``
     - as above
   * - thermal expansion tensor
     - ``(..., 3, 3)``
     - K⁻¹
   * - isochoric heat capacity
     - ``(...)``
     - J cell⁻¹ K⁻¹

Extension rules
---------------

* Add new QHA scalar/tensor fields to the passive result contract, HDF5 helper,
  shared interpolation service, and round-trip tests together.
* Add a new analysis operation to ``ThermoelasticAnalysisEngine`` before
  exposing it in CLI or GUI adapters.
* Keep component semantics in ``thermoelasticity.components``; plotting styles
  must not be imported by I/O.
* Keep complete numerical values in HDF5 and tables. Display precision belongs
  to report/render profiles only.
* Never clip unstable tensors, fabricate missing covariances, or silently
  replace invalid adiabatic states.

Schemas and compatibility
-------------------------

The current pre-release thermoelastic YAML input schema is 1.0 and requires
complete common-frame provenance.  Earlier unnormalized schemas are rejected.
The HDF5 payload has an independent ``thermoelastic_schema_version`` attribute,
currently 1.0.  Readers may retain older HDF5 support when interpretation is
unambiguous, but new writers emit only the current schema.

Required regression matrix
--------------------------

Every structural change must preserve:

* CLI/API calibration equivalence;
* point/grid/profile equivalence at identical states;
* frozen MgO and dolomite reference values;
* exact 0 K equality :math:`C^S=C^T`;
* frame-rotation invariants and stability masks;
* HDF5 read/write/read equality;
* wide-table and full-tensor export semantics;
* all staged architecture, CLI, plotting, packaging, and documentation tests.
