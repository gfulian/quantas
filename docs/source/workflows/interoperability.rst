Interoperability between workflows
==================================

Purpose and scope
-----------------

Interoperability in Quantas is a scientific transformation between typed
results, not a shortcut for copying similarly shaped arrays. A downstream
workflow must receive the units, tensor condition, normalization, masks, and
provenance that define the upstream result.

The currently supported chain is:

.. code-block:: text

   QHA input or result
          |
          v
   validated QHA--elasticity context
          |
          v
   thermoelastic calibration Cij(V)
          |
          v
   one reconstructed state Cij(P,T), rho(P,T)
          |                         |
          v                         v
      Elasticity                 SEISMIC

The physical approximations belong to the individual workflows. This page is
concerned with the **contracts between them**: what is transferred, what is
validated, and what remains the user's scientific responsibility.

Supported transformations
-------------------------

.. list-table:: Public interoperability boundaries
   :header-rows: 1
   :widths: 23 29 28 20

   * - Producer
     - Contract checked at the boundary
     - Consumer
     - Public interface
   * - QHA input, result, or native HDF5
     - Volume normalization, P--T grid, equilibrium-volume field, structural
       identity, and elastic-volume coverage
     - Thermoelastic calibration context
     - :func:`quantas.api.interop.qha_to_thermoelastic_context`
   * - Thermoelastic fit result or native HDF5
     - Requested P, T, tensor condition, extrapolation policy, stiffness, and
       QHA-consistent density
     - SEISMIC input
     - :func:`quantas.api.interop.thermoelastic_to_seismic`
   * - Thermoelastic point analysis
     - One selected stiffness tensor and density in the shared text contract
     - Elasticity or SEISMIC CLI
     - ``quantas thermoelasticity analysis point``

These are deliberate, narrow transformations. Quantas does not infer a general
workflow graph from file names or array shapes.

QHA to Thermoelasticity
-----------------------

QHA provides the thermodynamic path

.. math::

   (P,T) \longmapsto V(P,T),

while the elastic-volume series provides the static response

.. math::

   V \longmapsto C_{IJ}^{\mathrm{cold}}(V).

The interoperability context checks that these two datasets can be combined
before the QSA fit starts. In particular it verifies the equilibrium-volume
field, converts its volume unit, compares it with the elastic calibration
interval, records extrapolated states, and checks primitive atomic identity
when that information is available.

The equilibrium volume is mandatory. Heat capacity and the Cartesian thermal-
expansion tensor become mandatory only when an adiabatic conversion is
requested. The detailed QSA requirements are described in
:doc:`thermoelasticity`; the QHA result contract is described in
:doc:`../formats/ha_qha_hdf5`.

For production calculations, a native QHA HDF5 file is normally the clearest
boundary because it freezes the QHA options and can be inspected independently.
The Python API may instead pass a completed ``ResultData`` object directly; the
numerical transformation is the same.

Thermoelasticity to a material state
------------------------------------

The reusable thermoelastic archive stores a calibrated model, not a precomputed
P--T tensor grid. A downstream calculation requests one concrete state:

.. math::

   (P,T,\mathrm{condition})
   \longmapsto [C_{IJ}(P,T),\rho(P,T)].

The condition is ``isothermal`` or ``adiabatic`` and is part of the scientific
state. SEISMIC normally uses the adiabatic tensor when it is available because
elastic-wave propagation is adiabatic on the relevant timescale. Elasticity may
analyze either condition as long as it is identified correctly.

Two extrapolation questions remain separate throughout the transformation:

- whether P and T lie outside the archived QHA coordinate grid;
- whether the reconstructed equilibrium volume lies outside the static elastic
  calibration interval.

The ``fail``, ``warn``, and ``allow`` policies control how a requested state is
handled; the masks themselves remain part of the result.

Portable state files and in-memory states
-----------------------------------------

``thermoelasticity analysis point`` writes a compact stiffness-plus-density
text file that both Elasticity and SEISMIC understand. It is useful when the
selected state itself should be inspected, archived, or moved between
machines. The native thermoelastic HDF5 should still be retained because the
text file contains only one state and no calibration history.

The text representation is formatted, while the in-memory API keeps the
original ``float64`` values. Downstream CLI and API calculations should
therefore agree to the precision of the shared text contract, but file-mediated
results need not be bitwise identical to an in-memory calculation.

Frontend boundaries
-------------------

The CLI and API expose the same scientific operations with different artifact
boundaries:

.. list-table:: Equivalent workflow stages
   :header-rows: 1
   :widths: 24 36 40

   * - Stage
     - CLI
     - Python API
   * - QHA
     - ``quantas qha run`` writes native HDF5
     - :func:`quantas.api.qha.run` returns ``ResultData``
   * - Coupling
     - ``thermoelasticity run`` reads QHA HDF5 or YAML
     - :func:`quantas.api.interop.qha_to_thermoelastic_context`
   * - Calibration
     - ``quantas thermoelasticity run``
     - :func:`quantas.api.thermoelasticity.run_context`
   * - State reconstruction
     - ``thermoelasticity analysis point``
     - :func:`quantas.api.interop.thermoelastic_to_seismic`
   * - Downstream analysis
     - shared state file passed to Elasticity or SEISMIC
     - typed ``elasticity.Input`` or ``seismic.Input``

Executable end-to-end examples are distributed separately so this workflow page
does not duplicate a tutorial:

- :download:`CLI workflow <../_downloads/interoperability/workflow_cli.sh>`
- :download:`Python API workflow <../_downloads/interoperability/workflow_api.py>`

Run them from the project root when a full reproducibility check is required.

Provenance and persistence
--------------------------

The native archives remain the authoritative reusable checkpoints:

``QHA HDF5``
   Stores the resolved thermodynamic surface, options, diagnostics, warnings,
   and numerical precision metadata.

``Thermoelastic fit HDF5``
   Stores the reference EOS, component fits, covariance, source fields, tensor
   provenance, and calibration diagnostics.

``Elasticity / SEISMIC HDF5``
   Stores the downstream analysis and identifies the state supplied to that
   workflow.

The in-memory QHA--thermoelastic context is a validation object rather than a
new persistent file format.

Failure modes worth distinguishing
----------------------------------

``structural mismatch``
   Different primitive composition or atom ordering is rejected rather than
   reconciled automatically.

``missing QHA fields``
   Missing equilibrium volume prevents QSA. Missing heat capacity or thermal
   expansion may still permit an isothermal calibration.

``unsupported state``
   Coordinate and elastic-volume extrapolation are reported independently.

``invalid density``
   SEISMIC requires finite positive density and rejects a state for which it
   cannot be reconstructed.

``unstable stiffness``
   The state may still be meaningful to the Elasticity stability analysis, but
   SEISMIC cannot propagate waves through a non-positive-definite stiffness
   matrix.

``wrong tensor condition``
   File compatibility never changes an isothermal tensor into an adiabatic one.
   The intended experiment determines which condition is appropriate.

Current boundaries
------------------

The interoperability layer does not perform phase selection, reconcile
different compositions or primitive cells, create multiphase aggregates, infer
the appropriate thermodynamic tensor condition, or act as a workflow scheduler.
New transformations should be added only when their scientific contract can be
stated and tested explicitly.

Related documentation
---------------------

- QHA implementation: :doc:`qha`
- Thermoelastic implementation: :doc:`thermoelasticity`
- Elasticity implementation: :doc:`elasticity`
- SEISMIC implementation: :doc:`seismic`
- Public interoperability API: :doc:`../api/interoperability`
- Native HDF5 contracts: :doc:`../formats/hdf5`
