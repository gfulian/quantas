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
