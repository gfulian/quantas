External-code interfaces and parsers
====================================

Code-specific parsers live under :mod:`quantas.interfaces`. Their purpose is to
translate external output syntax into explicit physical data while preserving
provenance. They do not run Quantas calculations.

Separation of responsibilities
------------------------------

.. code-block:: text

   CRYSTAL / VASP / Phonopy output
                |
                v
   quantas.interfaces.<code> parser
                |
                v
   code-neutral parsed quantities
                |
                v
   module input generator / normalizer
                |
                v
   Quantas Input contract

This separation allows a second code to provide the same normalized module
input without duplicating the scientific workflow.

Reader lifecycle
----------------

Concrete parsers usually derive from :class:`quantas.models.BasicReader` and
maintain:

``completed``
   ``True`` only after all required data have been parsed and validated.

``error``
   A concise contextual message when the reader could not complete.

``load(path)``
   Reset state, identify the file, check completion, parse fields, validate
   semantics, and set ``completed``.

A parser may expose typed properties for stiffness, density, volume, energy,
pressure, structure, q-point data, or other source quantities. These properties
must state units and conventions.

File recognition and completion
-------------------------------

Do not parse every text file optimistically. Use stable markers to distinguish:

* correct program and calculation type;
* completed calculation;
* required output section;
* code version or format variant when relevant.

A file can be recognized but incomplete. Report these conditions separately so
the user knows whether the wrong file was selected or the external calculation
failed.

Units and conventions
---------------------

Convert units at the parser or normalization boundary and document the target.
Examples include:

* VASP elastic constants from kbar to GPa;
* crystal density to kg m\ :sup:`-3`;
* static energy to hartree or another declared native input unit;
* q-point weights exactly as supplied upstream;
* Voigt order conversion between external and Quantas conventions.

Never assume that two codes use the same shear ordering, stress sign,
crystallographic cell, or pre-stress correction.

Structures and symmetry
-----------------------

When a parser supplies a structure, preserve:

* lattice vectors;
* fractional or Cartesian coordinate convention;
* atomic numbers and atom order;
* primitive/conventional/supercell basis;
* transformation or repetition matrix;
* source symmetry information;
* Quantas/spglib analysis settings and results when calculated.

Do not reorder atoms without recording and testing the mapping. Thermoelastic
co-rotation and multi-volume structural paths depend on consistent atom
identity.

Phonon eigenvectors and mode continuity
---------------------------------------

Phonon parsers follow the same separation of responsibilities, with one extra
boundary:

.. code-block:: text

   external phonon output
           |
           v
   code-specific parser
           |
           +-- frequencies
           +-- normalized eigenvectors
           +-- atom ordering and q metadata
           |
           v
   PhononModeData
           |
           v
   core.numerics phonon tracker
           |
           v
   normalized PhononInputData / QHA input

The parser is responsible for reconstructing the eigenvector representation of
the external code and documenting its normalization.  It must not decide that a
QHA branch correspondence is acceptable merely because two raw mode indices
match.

The backend-neutral tracker receives ``float64`` frequencies and ``complex128``
unit-norm eigenvectors.  It knows nothing about CRYSTAL markers, Phonopy YAML,
or future VASP/QE syntax.  Conversely, the CRYSTAL parser does not know QHA
failure policy or CLI rendering.

For CRYSTAL, general-q vectors are reconstructed from in-phase and anti-phase
components and converted to unit-norm mass-weighted directions before they
leave the interface layer.  Degenerate-subspace matching, Hungarian assignment,
ambiguity classification, and leave-one-out validation belong to the numerical
tracking layer.

.. important::

   Do not move mode-continuity policy into a code-specific parser.  A future
   interface must be able to supply the same neutral ``PhononModeData`` and
   obtain the same tracking result from the same normalized arrays.

When an external workflow has already established continuity, preserve that
fact explicitly as source provenance rather than relabelling it as a Quantas
tracking result.  The ``crystal-qha`` path is the current example.

Provenance
----------

Store the information required to understand the parsed quantity later:

* source path and raw text when appropriate;
* external code and calculation type;
* relevant keyword values;
* pressure correction or relaxation state;
* cell normalization;
* parser version or schema label;
* warnings about missing optional fields.

For CRYSTAL elasticity, for example, the ``PRESSURE`` keyword and reported
elastic pressure are scientifically relevant because they establish whether
the output contains the stress-corrected coefficients required under
hydrostatic pre-stress.

CRYSTAL elastic volume series
-----------------------------

:func:`quantas.interfaces.crystal.read_crystal_elastic_series` composes a set
of completed ELASTCON or ELAPIEZO outputs into the backend-neutral
:class:`quantas.models.elastic_states.ElasticStateSeries` contract.  This is
the interface boundary used by Kieffer and available to other multi-volume
elastic workflows; it does not create a Kieffer-specific elastic format.

The importer requires finite volume, density, static energy, and stiffness at
every state.  It sorts the resulting states by increasing volume and selects
the minimum-static-energy state as the reference.  Tensor axes remain in the
CRYSTAL Cartesian frame.

Pressure selection is explicit:

``auto``
   Preserve tensors already corrected by CRYSTAL when the ``PRESSURE`` keyword
   is present.  For raw tensors, use the pressure printed for the unstrained
   stress tensor.

``output_stress``
   Require raw tensors and explicitly use their reported unstrained-stress
   pressure.

``manual``
   Require one finite pressure in GPa per input file.  Values follow input-file
   order before volume sorting.  Positive pressure denotes compression.

``deferred``
   Retain a raw tensor without attaching pressure. This adapter-level route
   requires ``apply_prestress_correction=False`` and exists so a composing
   workflow can attach independently derived pressure provenance before a
   separate correction. It is not exposed as a user-facing ``add-kieffer``
   pressure source.

By default, raw energy--strain tensors are converted once to Wallace
hydrostatic coefficients.  The pressure value, source, method, source tensor
kind, and software applying the correction are retained in each state.
Passing a non-auto pressure policy for a tensor already corrected by CRYSTAL is
an error, preventing an accidental second correction.

.. code-block:: python

   from quantas.interfaces.crystal import read_crystal_elastic_series

   series = read_crystal_elastic_series(
       ["state_01.out", "state_02.out", "state_03.out"],
       pressure_policy="output_stress",
   )

If the structural block reconstructed from an output does not have the same
volume as the final elastic scalar, the final elastic volume remains
authoritative.  The inconsistent lattice is not attached to the state; its
volume and the failed consistency check are recorded in metadata.  This avoids
silently coupling a tensor to a stale geometry block while retaining the
diagnostic needed to inspect the source output.

Kieffer input enrichment
------------------------

The public HA and QHA APIs expose ``add_kieffer_input``.  Their shared
implementation reads the phonon input and the CRYSTAL elastic volume series,
builds the anisotropic acoustic averages, validates the appropriate HA or QHA
applicability contract, and writes a new YAML file.  The corresponding command
is registered under both workflows:

.. code-block:: console

   quantas ha add-kieffer ha.yaml state.out -o ha-kieffer.yaml
   quantas qha add-kieffer qha.yaml --elastic-list elastic-files.txt \
       --interface crystal -o qha-kieffer.yaml

Paths inside ``elastic-files.txt`` are resolved relative to the list file. Blank
lines and lines beginning with ``#`` are ignored. This makes the list portable
when the complete calculation directory is moved.

The default ``--pressure-source auto`` preserves tensors corrected by CRYSTAL's
``PRESSURE`` keyword and otherwise uses pressure from the unstrained stress.
Manual pressure values can be supplied in input-file order:

.. code-block:: console

   quantas qha add-kieffer qha.yaml --elastic-list elastic-files.txt \
       --pressure-source manual \
       --pressure 11.53 --pressure 8.718 --pressure 6.069 \
       -o qha-kieffer.yaml

For multi-volume QHA inputs, pressure may instead be evaluated from the static
energy-volume arrays already present in the phonon input:

.. code-block:: console

   quantas qha add-kieffer qha.yaml --elastic-list elastic-files.txt \
       --interface crystal --pressure-source energy-eos --eos BM3 \
       -o qha-kieffer.yaml

   quantas qha add-kieffer qha.yaml --elastic-list elastic-files.txt \
       --interface crystal --pressure-source energy-polynomial --degree 3 \
       -o qha-kieffer.yaml

The reusable fit operations live in :mod:`quantas.core.physics.eos`; pressure
assignment and hydrostatic correction remain separate operations in
:mod:`quantas.core.physics.elasticity`. This boundary lets tests verify that
``P(V)`` is attached to an unmodified raw tensor before the tensor is corrected
exactly once.

The destination defaults to ``<input-stem>-kieffer.yaml`` and must differ from
the source path. An existing Kieffer block is never replaced silently.  The
generated top-level ``kieffer`` mapping identifies the sine-wave method and its
additive composition, declares canonical units, and stores one state per
volume with:

* cutoff frequencies in Hz;
* effective slow-shear, fast-shear, and longitudinal velocities in km/s;
* pressure and tensor convention;
* elastic-state association;
* spherical-quadrature diagnostics;
* original phonon and elastic source paths.

The public ``read_kieffer_input`` operation restores this block as a validated
:class:`quantas.models.kieffer.KiefferVolumeSeries`.  It can therefore be
passed explicitly to the HA/QHA calculation APIs without reconstructing the
elastic calculation.

Error handling
--------------

Catch only exceptions that can be converted into a useful parser error. Avoid a
broad ``except Exception`` that turns programming errors into misleading input
messages.

A useful error identifies:

* the file;
* the section or quantity;
* the expected marker or shape;
* the observed problem.

Do not return a zero-filled scientific array after a required section failed to
parse.

Fixture strategy
----------------

Parser tests should include:

* a small complete real output;
* incomplete output;
* wrong calculation type;
* malformed numerical section;
* optional section absent;
* supported code-version variants;
* unit and convention conversion;
* structure and atom ordering;
* integration through the normalized module input generator.

Keep fixtures scientifically recognizable but small enough for the repository.
When a full external output is too large, retain the required sections without
inventing numerical values.

Adding support for a new code
-----------------------------

#. Start with one mature workflow and one clearly defined output quantity.
#. Implement the parser under ``quantas.interfaces.<code>``.
#. Add a code-neutral adapter to the existing module input contract.
#. Add an interface selector only after the parser is tested.
#. Document the supported code version, calculation settings, and limitations.
#. Compare normalized outputs from two codes when equivalent datasets are
   available.
#. Generalize shared parsing concepts only after two real implementations show
   the common abstraction.

Do not create a large universal parser hierarchy before the external formats
have demonstrated a stable common structure.
