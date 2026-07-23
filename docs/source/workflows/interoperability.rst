Interoperability between workflows
==================================

Quantas treats interoperability as an explicit scientific transformation, not
as accidental compatibility between file formats.  A downstream module must
receive the physical quantities, units, tensor convention, thermodynamic
condition, masks, and provenance that it requires.  Moving arrays manually may
appear simple, but it can silently lose exactly the information that makes a
cross-module result reproducible.

The currently supported public path is:

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

Two public transformations in :mod:`quantas.api.interop` implement the
non-trivial parts of this path:

- :func:`quantas.api.interop.qha_to_thermoelastic_context` validates and couples
  QHA data with a static elastic-volume series;
- :func:`quantas.api.interop.thermoelastic_to_seismic` reconstructs one
  thermoelastic state and returns a validated SEISMIC input.

The CLI represents the same boundaries with native HDF5 files and a shared
Elasticity/SEISMIC text input.  The Python API may keep the same typed contracts
in memory and only persist the artifacts that the user wants to retain.

.. important::

   Interoperability does not make two scientifically incompatible calculations
   compatible.  It verifies the contracts that Quantas can inspect, reports
   extrapolation and missing fields, and refuses invalid transformations.  The
   user remains responsible for phase stability, electronic-structure quality,
   phonon convergence, and the scientific suitability of the selected
   approximation.

Supported transformations
-------------------------

The current public interoperability surface is intentionally small.

.. list-table:: Supported cross-module paths
   :header-rows: 1
   :widths: 25 24 27 24

   * - Source
     - Transformation
     - Destination
     - Public interface
   * - QHA input, result object, or native HDF5
     - Validate volume coverage, atomic normalization, units, and required QHA
       fields
     - :class:`quantas.api.thermoelasticity.Context`
     - :func:`quantas.api.interop.qha_to_thermoelastic_context`
   * - Thermoelastic calibration result or native HDF5
     - Reconstruct one P--T state, select :math:`C^T` or :math:`C^S`, and attach
       QHA-consistent density
     - :class:`quantas.api.seismic.Input`
     - :func:`quantas.api.interop.thermoelastic_to_seismic`
   * - Thermoelastic point analysis
     - Write the selected stiffness tensor and density using the shared text
       contract
     - Elasticity and SEISMIC CLI input
     - ``quantas thermoelasticity analysis point``

There is no generic workflow graph engine and no automatic chaining based on
file names.  Unsupported or future transformations are not inferred from
similar-looking arrays.

Why the CLI and API look different
----------------------------------

The CLI communicates between independent processes.  Persistent files are
therefore the natural boundary:

.. code-block:: text

   qha run
   -> QHA HDF5
   -> thermoelasticity run
   -> thermoelastic fit HDF5
   -> thermoelasticity analysis point
   -> shared text input
   -> elasticity run / seismic run

This design provides several advantages:

- every expensive stage can be resumed independently;
- intermediate results can be inspected before they are consumed;
- the exact source archive remains available for provenance;
- downstream calculations do not repeat QHA or thermoelastic calibration;
- shell scripts and batch systems can schedule each stage separately.

The Python API can pass :class:`~quantas.api.common.ResultData` and passive
input contracts directly:

.. code-block:: text

   QHA ResultData
   -> Thermoelastic Context
   -> thermoelastic ResultData
   -> Seismic Input
   -> seismic ResultData

This avoids temporary files while preserving the same numerical path.  Writing
HDF5 remains recommended at scientific checkpoints, even in an in-memory
workflow.

General interoperability rules
------------------------------

Use typed transformations
~~~~~~~~~~~~~~~~~~~~~~~~~

Prefer :mod:`quantas.api.interop` or a documented module export function to
manual array copying.  The transformation should own:

- unit conversion;
- shape validation;
- tensor-condition selection;
- density propagation;
- coverage and extrapolation checks;
- metadata and source provenance;
- clear failure messages.

Freeze expensive upstream results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For exploratory notebooks, passing a QHA result object directly is convenient.
For publication calculations, a native QHA HDF5 is normally the better source:

- it fixes the exact QHA options and versioned result;
- it avoids an accidental recalculation with changed defaults;
- it can be inspected independently;
- it makes the downstream calibration reproducible without the original
  Python process.

The same principle applies to the reusable thermoelastic fit archive.

Do not infer units from array magnitude
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The QHA result envelope records the volume unit.  The interoperability layer
converts equilibrium volumes to the angstrom-based convention used by the
thermoelastic elastic-volume series.  Users should not reproduce this logic by
assuming that a numerical range "looks like" cubic angstroms.

Preserve tensor condition
~~~~~~~~~~~~~~~~~~~~~~~~~

Thermoelastic results may contain both isothermal and adiabatic stiffness
coefficients.  The selected condition is part of the downstream input, not a
plotting preference.

- Elasticity may analyze either condition if it is named and interpreted
  correctly.
- SEISMIC normally uses the adiabatic tensor because elastic-wave propagation
  is adiabatic on the relevant timescale.
- Requesting ``adiabatic`` is valid only when the upstream result contains the
  data needed for the conversion.

Treat extrapolation as a source property
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The transformation distinguishes two limitations:

- QHA coordinate extrapolation outside the stored P--T grid;
- elastic-volume extrapolation outside the calibrated static :math:`C_{ij}(V)`
  range.

A downstream result must not erase these distinctions.  The point
transformation accepts ``fail``, ``warn``, or ``allow`` so the policy is an
explicit user decision.

QHA to Thermoelasticity
-----------------------

Scientific purpose
~~~~~~~~~~~~~~~~~~

The QHA result supplies the thermodynamic path:

.. math::

   (P,T) \longmapsto V(P,T),

while the static elastic series supplies the cold response:

.. math::

   V \longmapsto C_{IJ}^{\mathrm{cold}}(V).

The coupling context verifies that these datasets can be used together before
any finite-strain calibration is attempted.  The quasi-static approximation is
explained in :doc:`thermoelasticity` and its thermodynamic basis in
:doc:`../theory/thermoelasticity`.

Accepted QHA sources
~~~~~~~~~~~~~~~~~~~~

:func:`quantas.api.interop.load_qha_result` accepts:

- an existing QHA :class:`~quantas.api.common.ResultData`;
- a native QHA HDF5 path;
- a QHA or phonon input object;
- a YAML path that must first be calculated.

When a calculation is required, :class:`quantas.api.qha.Options` must describe
the desired P--T domain and numerical strategy.  Options are ignored when the
source is already a completed QHA result.

Context checks
~~~~~~~~~~~~~~

The public transformation performs the following checks and normalizations:

1. normalize the thermoelastic elastic-volume input;
2. load or calculate a valid QHA result envelope;
3. require ``equilibrium_volume`` with shape ``(nT, nP)``;
4. convert the QHA volume unit to angstroms;
5. compare every equilibrium volume with the sampled elastic-volume interval;
6. list QHA fields that are absent;
7. compare primitive atomic numbers and ordering when the QHA input stores
   them;
8. record source mode, volume bounds, and the number of extrapolated states.

The returned context exposes:

``input_data``
   Normalized thermoelastic input.

``qha_result_data``
   Complete QHA result envelope.

``qha``
   Typed QHA payload.

``extrapolation_mask``
   Boolean ``(nT, nP)`` mask for QHA equilibrium volumes outside the elastic
   calibration interval.

``missing_qha_fields``
   Names of optional or required downstream fields absent from the QHA result.

``metadata``
   Elastic and QHA volume bounds, point counts, and source mode.

Missing fields do not all have the same consequence.  Equilibrium volume is
mandatory for the QSA.  Heat capacity and the Cartesian thermal-expansion
tensor are additionally required when the adiabatic correction is requested.
Structural fields may be absent without preventing an isothermal tensor
calibration.

CLI: explicit two-stage calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most reproducible CLI path first freezes QHA:

.. code-block:: console

   quantas qha run examples/qha/crystal-phonons/dol_pbe0.yaml \
       --scheme freq \
       --minimization poly \
       --thermal-expansion mixed_derivative \
       --temperature 300 1800 100 \
       --pressure 0 10 2 \
       --energy-degree 3 \
       --frequency-degree 3 \
       --output dolomite_qha.hdf5 \
       --report dolomite_qha.log \
       --no-progress \
       --force

The native result is then consumed by the thermoelastic calibration:

.. code-block:: console

   quantas thermoelasticity run \
       examples/thermoelasticity/dol_pbe0_thermoelastic.yaml \
       dolomite_qha.hdf5 \
       --reference-eos BM3 \
       --finite-strain-order 3 \
       --adiabatic auto \
       --validation standard \
       --output dolomite_fit.hdf5 \
       --report dolomite_fit.log \
       --no-progress \
       --force

The second command reads the QHA archive.  It does not rerun QHA.

CLI: direct YAML source
~~~~~~~~~~~~~~~~~~~~~~~

``thermoelasticity run`` may also receive a QHA YAML input directly:

.. code-block:: console

   quantas thermoelasticity run \
       examples/thermoelasticity/dol_pbe0_thermoelastic.yaml \
       examples/qha/crystal-phonons/dol_pbe0.yaml \
       --qha-temperature 300 1800 100 \
       --qha-pressure 0 10 2 \
       --qha-scheme freq \
       --qha-minimization poly \
       --qha-degree 3 \
       --reference-eos BM3 \
       --finite-strain-order 3 \
       --output dolomite_fit.hdf5 \
       --no-progress \
       --force

This is convenient for a compact calculation, but the intermediate QHA result
is less visible as an independent scientific checkpoint.  Use the explicit
two-stage path when the QHA result must be reviewed, reused, or cited
separately.

API: inspect the coupling before fitting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The API makes the prevalidated context directly available:

.. code-block:: python

   from quantas.api import interop, qha, thermoelasticity

   qha_result = qha.run(
       "examples/qha/crystal-phonons/dol_pbe0.yaml",
       options=qha.Options(
           temperature_min=300.0,
           temperature_max=1800.0,
           temperature_step=100.0,
           pressure_min=0.0,
           pressure_max=10.0,
           pressure_step=2.0,
           scheme="freq",
           minimization="poly",
           thermal_expansion_method="mixed_derivative",
       ),
   )

   context = interop.qha_to_thermoelastic_context(
       "examples/thermoelasticity/dol_pbe0_thermoelastic.yaml",
       qha_result,
   )

   print(context.metadata)
   print(context.missing_qha_fields)
   print(context.extrapolation_mask.sum())

   fit_result = thermoelasticity.run_context(
       context,
       options=thermoelasticity.Options(
           reference_eos="BM3",
           finite_strain_order=3,
           adiabatic_mode="auto",
       ),
   )

The convenience function :func:`quantas.api.thermoelasticity.run` performs the
context preparation and calibration in one call.  Use the explicit context
when coverage or completeness must be inspected before fitting.

Thermoelasticity to one material state
--------------------------------------

A thermoelastic fit archive stores a reusable model.  A downstream Elasticity
or SEISMIC calculation requires one concrete state:

.. math::

   (P,T,\mathrm{condition})
   \longmapsto
   \left[C_{IJ}(P,T),\rho(P,T)\right].

The state transformation therefore needs:

- pressure in GPa;
- temperature in K;
- ``isothermal`` or ``adiabatic`` stiffness condition;
- an extrapolation policy.

It reconstructs a one-point grid through the same post-fit engine used by the
thermoelastic CLI, selects the requested tensor, validates density, and creates
a descriptive job name containing P, T, and tensor condition.

CLI: write a shared state input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: console

   quantas thermoelasticity analysis point \
       dolomite_fit.hdf5 5 800 \
       --tensor-condition adiabatic \
       --extrapolation fail \
       --output dolomite_5GPa_800K.dat \
       --force

The output contains:

- the complete symmetric stiffness matrix in GPa;
- QHA-consistent density in kg m\ :sup:`-3`;
- P, T, and tensor condition in the title;
- the common text structure understood by both downstream modules.

The file is intentionally a portable boundary.  It can be inspected, archived,
or passed to another machine without the thermoelastic HDF5 reader.

The portable text representation writes stiffness values with six digits after
the decimal point and density with eight.  An in-memory API transformation
retains the original ``float64`` values.  Downstream CLI and API results are
therefore equivalent to the precision of the shared text contract, but they
need not be bitwise identical after the CLI has re-read that file.  Retain the
native thermoelastic HDF5 when full stored precision is required.

State to Elasticity
-------------------

CLI
~~~

The shared state file can be analyzed directly:

.. code-block:: console

   quantas elasticity run dolomite_5GPa_800K.dat \
       --output dolomite_state_elasticity.hdf5 \
       --report dolomite_state_elasticity.log \
       --no-progress \
       --force

Elasticity uses the stiffness tensor.  Density is retained in the input file
for compatibility with SEISMIC but is not needed for static directional
elastic properties.

API
~~~

A dedicated ``thermoelastic_to_elasticity`` helper is not required because the
Elasticity public input contains only a title and a stiffness tensor.  The same
state returned for SEISMIC can be adapted without unit conversion or private
imports:

.. code-block:: python

   import numpy as np

   from quantas.api import elasticity, interop

   state = interop.thermoelastic_to_seismic(
       "dolomite_fit.hdf5",
       pressure=5.0,
       temperature=800.0,
       tensor_condition="adiabatic",
       extrapolation_policy="fail",
   )

   elastic_input = elasticity.Input(
       jobname=state.jobname,
       stiffness=np.asarray(state.stiffness, dtype=np.float64),
       source=state.source,
   )
   elastic_result = elasticity.run(elastic_input)

Alternatively, write the shared state input with
:func:`quantas.api.thermoelasticity.write_state_input` and pass its path to
:func:`quantas.api.elasticity.run`.  The file-mediated route is useful when the
state itself is an artifact that must be retained.

State to SEISMIC
----------------

CLI
~~~

.. code-block:: console

   quantas seismic run dolomite_5GPa_800K.dat \
       --level phase \
       --ntheta 31 \
       --nphi 61 \
       --output dolomite_state_seismic.hdf5 \
       --report dolomite_state_seismic.log \
       --no-progress \
       --force

The reduced phase-only grid above is suitable for an interoperability check.
A final acoustic study should choose sampling and property level according to
:doc:`seismic` and :doc:`../tutorials/seismic`.

API
~~~

The direct API path returns a validated :class:`quantas.api.seismic.Input`:

.. code-block:: python

   from quantas.api import interop, seismic

   seismic_input = interop.thermoelastic_to_seismic(
       "dolomite_fit.hdf5",
       pressure=5.0,
       temperature=800.0,
       tensor_condition="adiabatic",
       extrapolation_policy="fail",
   )

   seismic_result = seismic.run(
       seismic_input,
       options=seismic.Options(
           ntheta=31,
           nphi=61,
           hemisphere="upper",
           level="phase",
       ),
   )

The transformation centralizes reconstruction, condition selection, density,
and provenance.  Do not extract ``stiffness_adiabatic[0, 0]`` and density by
hand unless reproducing the same validation logic is an explicit research
requirement.

Complete CLI example
--------------------

The complete shell workflow is distributed as:

:download:`workflow_cli.sh <../_downloads/interoperability/workflow_cli.sh>`

Run it from the project root:

.. code-block:: console

   bash examples/interoperability/workflow_cli.sh interoperability_cli

It writes:

.. code-block:: text

   interoperability_cli/
   |-- dolomite_qha.hdf5
   |-- dolomite_qha.log
   |-- dolomite_fit.hdf5
   |-- dolomite_fit.log
   |-- dolomite_5GPa_800K.dat
   |-- dolomite_state_elasticity.hdf5
   |-- dolomite_state_elasticity.log
   |-- dolomite_state_seismic.hdf5
   `-- dolomite_state_seismic.log

Complete Python API example
---------------------------

The equivalent public-API workflow is available as:

:download:`workflow_api.py <../_downloads/interoperability/workflow_api.py>`

Run it from the project root:

.. code-block:: console

   python examples/interoperability/workflow_api.py \
       --output-dir interoperability_api

The script demonstrates both styles:

- in-memory typed transformations for the numerical chain;
- native HDF5 and shared text artifacts at reusable checkpoints.

For the distributed dolomite example at 5 GPa and 800 K, the verified output is:

.. code-block:: text

   Context extrapolated states: 0
   State density / kg m^-3: 2908.736749
   State C11_S / GPa: 230.296285
   Hill bulk modulus / GPa: 109.557211
   Sampled V_P range / km s^-1: 6.600439 9.171102

The V\ :sub:`P` range belongs to the deliberately coarse ``31 x 61`` phase
sampling used by the example.  It is a reproducibility checkpoint, not a
converged directional study.

CLI/API equivalence
-------------------

The two frontends use the same public numerical operations but expose different
artifact boundaries.

.. list-table:: Equivalent stages
   :header-rows: 1
   :widths: 23 36 41

   * - Stage
     - CLI
     - Python API
   * - QHA
     - ``quantas qha run`` writes HDF5
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
   * - Elasticity input
     - shared ``.dat`` file
     - :class:`quantas.api.elasticity.Input`
   * - SEISMIC input
     - shared ``.dat`` file
     - :class:`quantas.api.seismic.Input`

Numerical equivalence should be checked at the boundary quantities before
comparing derived reports or sampled extrema:

- pressure and temperature;
- selected tensor condition;
- density;
- all 21 independent entries of the symmetric stiffness matrix;
- downstream options such as sphere sampling.

The QHA and thermoelastic HDF5 paths preserve native arrays and can be compared
exactly when options are identical.  The portable state text file is a
formatted export; compare downstream results with a tolerance consistent with
its six-decimal stiffness representation rather than demanding bitwise
equality.  For the distributed example, the resulting phase speeds differ by
less than :math:`2	imes10^{-8}` km s\ :sup:`-1` between the in-memory API path
and the file-mediated CLI path.

Provenance and persistence
--------------------------

QHA result
~~~~~~~~~~

The QHA HDF5 stores normalized input, resolved options, equilibrium properties,
fit diagnostics, warnings, and numerical precision metadata.

Thermoelastic context
~~~~~~~~~~~~~~~~~~~~~

The context is an in-memory validation object.  Its diagnostics become part of
the calibration result and report; it is not a separate native file format.

Thermoelastic fit archive
~~~~~~~~~~~~~~~~~~~~~~~~~

The fit HDF5 is the reusable source for all later Point, Grid, and Profile
analyses.  It stores the reference EOS, component fits, covariance,
classification, QHA source fields, and calibration metadata.

Shared state input
~~~~~~~~~~~~~~~~~~

The text state file is a compact interoperability artifact.  It does not
replace the full thermoelastic archive because it contains only one selected
state and no fit history, uncertainty field, or source grid.

Downstream HDF5 files
~~~~~~~~~~~~~~~~~~~~~

Elasticity and SEISMIC each write their own native result envelopes.  Their
metadata identify the downstream workflow, while the input title and source
retain the state origin where available.

Failure modes and diagnostics
-----------------------------

Different primitive atomic normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the QHA input stores atomic numbers, their sequence must match the
thermoelastic reference structure.  A mismatch may represent a different
primitive cell, atom ordering, composition, or phase and is rejected.

Missing equilibrium volume
~~~~~~~~~~~~~~~~~~~~~~~~~~

Thermoelastic QSA cannot proceed without ``equilibrium_volume``.  A harmonic
result or incomplete QHA payload is not a substitute.

Incomplete adiabatic inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~

An isothermal calibration may be valid even when QHA heat capacity or Cartesian
thermal expansion is absent.  Requesting an adiabatic state then fails or is
reported according to the selected thermoelastic adiabatic mode.  Quantas does
not silently label :math:`C^T` as :math:`C^S`.

QHA volume outside elastic support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The coupling context counts these states before fitting.  A point
transformation with ``extrapolation_policy="fail"`` refuses them.  ``warn`` or
``allow`` should be used only with a documented scientific justification.

Requested P or T outside the stored QHA grid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is coordinate extrapolation and is distinct from elastic-volume
extrapolation.  A state can trigger either mask, both, or neither.

Invalid density
~~~~~~~~~~~~~~~

SEISMIC requires finite positive density.  The direct transformation rejects a
state whose QHA-consistent density cannot be reconstructed.

Unstable stiffness
~~~~~~~~~~~~~~~~~~

The transformation can construct the state, but SEISMIC rejects a stiffness
matrix that is not positive definite.  Elasticity may report instability and
skip directional properties according to its own workflow contract.

Wrong tensor condition
~~~~~~~~~~~~~~~~~~~~~~

A state exported as isothermal remains isothermal.  Downstream file compatibility
does not imply that the tensor condition is appropriate for the intended
experiment or wave-propagation calculation.

Repeated upstream calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Passing a YAML source to ``load_qha_result`` or ``thermoelasticity run`` may
calculate QHA.  Passing a QHA ``ResultData`` or HDF5 reads the completed result.
Use a frozen HDF5 when recomputation would be undesirable.

Performance and reuse
---------------------

Interoperability itself is inexpensive.  The dominant costs are the upstream
and downstream scientific calculations:

- QHA over the selected P--T grid;
- thermoelastic component calibration;
- Elasticity 2D/3D fields when requested;
- SEISMIC spherical sampling, especially group and enhancement levels.

A productive workflow is therefore:

1. calculate and inspect QHA once;
2. calibrate and inspect the thermoelastic fit once;
3. reconstruct many inexpensive Point states from the fit archive;
4. run coarse downstream checks;
5. increase Elasticity or SEISMIC sampling only for selected states.

Do not rerun QHA for every seismic state.  Do not recalibrate
:math:`C_{IJ}(V)` for every point on a profile.  The reusable fit archive exists
precisely to separate those costs.

Decision guide
--------------

.. list-table:: Recommended interoperability choices
   :header-rows: 1
   :widths: 32 32 36

   * - Situation
     - Recommended path
     - Reason
   * - One exploratory notebook
     - Pass ``ResultData`` in memory
     - Minimal temporary I/O while retaining typed contracts
   * - Publication or shared calculation
     - Freeze QHA and thermoelastic HDF5 checkpoints
     - Reproducibility and independent inspection
   * - Need to review coverage before fitting
     - Build an explicit interoperability context
     - Exposes missing fields and extrapolation mask
   * - One reusable P--T state
     - Write the shared text input
     - Portable boundary for both downstream CLIs
   * - Direct Python SEISMIC calculation
     - Use ``thermoelastic_to_seismic``
     - Centralized density, tensor condition, and validation
   * - Direct Python Elasticity calculation
     - Construct ``elasticity.Input`` from the validated state stiffness
     - Elasticity requires no density-specific transformation
   * - Many states along a profile
     - Analyze the profile once and select states deliberately
     - Avoid repeated calibration and hidden gridding
   * - State outside support
     - Keep ``fail`` until the extrapolation is scientifically reviewed
     - Prevents silent use of unsupported tensors

Boundaries of the current implementation
----------------------------------------

The current interoperability layer does not:

- perform phase-equilibrium selection;
- transform between different chemical compositions;
- reconcile different primitive cells or atom orderings automatically;
- derive QHA q-point weights;
- add explicit phonon-strain elastic terms absent from the QSA;
- infer whether isothermal or adiabatic coefficients are appropriate;
- build multiphase aggregate properties;
- chain EOS fitting into QHA or thermoelasticity implicitly;
- create a general-purpose workflow scheduler.

New transformations should be added only when the scientific contract is
well-defined and testable.  The development procedure is described in
:doc:`../developer/extending`, :doc:`../developer/public_api`, and
:doc:`../developer/change_recipes`.

Related documentation
---------------------

- QHA implementation: :doc:`qha`
- Thermoelastic implementation: :doc:`thermoelasticity`
- Elasticity implementation: :doc:`elasticity`
- SEISMIC implementation: :doc:`seismic`
- Complete thermoelastic tutorial: :doc:`../tutorials/thermoelasticity`
- Public interoperability functions: :doc:`../api/interoperability`
- QHA HDF5 schema: :doc:`../formats/ha_qha_hdf5`
- Thermoelastic HDF5 schema: :doc:`../formats/thermoelastic_hdf5`
- Shared Elasticity/SEISMIC input: :doc:`../formats/elasticity_input`
