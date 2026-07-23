Your first calculation
======================

This page runs one small calculation twice: first from the command line and
then through the public Python API.  Both paths use the same calcite stiffness
tensor and must reproduce the same scientific result.

The example is intentionally short.  Its purpose is to show the complete
Quantas pattern:

.. code-block:: text

   input file
   → validated scientific calculation
   → typed result
   → deterministic report
   → native HDF5 file

The detailed interpretation of the elastic properties is developed in the
:doc:`../tutorials/elasticity` tutorial.

Prepare a working directory
---------------------------

Create an empty directory and place the example input inside it:

.. code-block:: console

   mkdir quantas-first-calculation
   cd quantas-first-calculation

:download:`Download calcite.dat <../_downloads/calcite.dat>`

The file contains a title, the upper-triangular part of a stiffness matrix in
Voigt notation, and an optional density.  Elastic constants are expressed in
GPa.  The Elasticity reader validates the matrix and expands it to a symmetric
:math:`6\times6` tensor before the calculation starts.

Run the command-line calculation
--------------------------------

Execute:

.. code-block:: console

   quantas elasticity run calcite.dat \
       --output calcite_elasticity.hdf5 \
       --report calcite_elasticity.log \
       --no-progress \
       --force

The command performs four main operations:

#. it reads and validates the input;
#. it calculates the compliance tensor, Voigt--Reuss--Hill averages,
   mechanical stability, and global directional extrema;
#. it writes a deterministic plain-text report;
#. it stores the complete result in native HDF5 format.

``--no-progress`` only disables the transient terminal progress display.  It
does not change the calculation.  ``--force`` permits the two requested files
to be replaced when the example is repeated.

Check the terminal or report output
-----------------------------------

The report should contain the following checkpoints to the displayed
precision:

.. code-block:: text

   Crystal system: trigonal_low
   Mechanically stable: True
   Hill bulk modulus / GPa: 83.099307
   Hill shear modulus / GPa: 35.241563
   Young minimum / GPa: 57.603207
   Young maximum / GPa: 154.060777

The Hill values describe an equivalent isotropic polycrystalline aggregate.
The Young-modulus extrema describe the directional response of the single
crystal.  They answer different physical questions and should not be compared
as if they were the same type of average.

Two files should now be present:

``calcite_elasticity.log``
   Human-readable plain text.  It is useful for inspection, archiving, and
   attaching a stable record to a calculation.

``calcite_elasticity.hdf5``
   Native structured result.  It retains numerical arrays at full precision,
   metadata, normalized inputs and options, warnings, and the scientific
   payload needed by later plotting or export commands.

Repeat the calculation through Python
-------------------------------------

The equivalent public-API program is available here:

:download:`Download first_elasticity_api.py <../_downloads/getting_started/first_elasticity_api.py>`

Place the script beside ``calcite.dat`` and run:

.. code-block:: console

   python first_elasticity_api.py

The essential code is:

.. literalinclude:: ../_downloads/getting_started/first_elasticity_api.py
   :language: python
   :linenos:

The script uses only supported public namespaces:

.. code-block:: python

   from quantas.api import elasticity, rendering

``elasticity.run`` returns a complete
:class:`quantas.api.common.ResultData` envelope.  ``elasticity.get_result``
extracts and validates the typed Elasticity payload.  The report builder
returns frontend-neutral tables, and :func:`quantas.api.rendering.render_tables`
converts those tables to deterministic plain text.  The scientific calculation
never depends on terminal or plotting objects.

Compare the two frontends
-------------------------

The API script creates:

.. code-block:: text

   calcite_elasticity_api.hdf5
   calcite_elasticity_api.log

The displayed checkpoints must agree with the command-line calculation.  A
strict numerical comparison can be made directly from the two HDF5 files:

.. code-block:: python

   import numpy as np
   from quantas.api import elasticity

   cli = elasticity.get_result(
       elasticity.read_result("calcite_elasticity.hdf5")
   )
   api = elasticity.get_result(
       elasticity.read_result("calcite_elasticity_api.hdf5")
   )

   np.testing.assert_allclose(cli.stiffness, api.stiffness)
   np.testing.assert_allclose(cli.compliance, api.compliance)
   assert cli.averages == api.averages
   assert cli.stability.is_stable == api.stability.is_stable

This equality is a central Quantas design rule: changing the frontend must not
change the scientific result.

Where to continue
-----------------

- Read :doc:`python_api` to understand result envelopes, typed payloads,
  observers, and module discovery.
- Read :doc:`command_line` for command organization, common options, progress,
  reports, and exit behavior.
- Read :doc:`results` for the roles of HDF5, reports, exports, and plots.
- Continue with :doc:`../tutorials/index` for complete scientific examples.
