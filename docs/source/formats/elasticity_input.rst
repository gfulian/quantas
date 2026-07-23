Elasticity and SEISMIC text input
=================================

Elasticity and SEISMIC use the same compact text representation for a
second-order stiffness tensor.  SEISMIC adds one required density value after
the tensor.  The format is intentionally small so that a tensor reconstructed
by Thermoelasticity, generated from an external code, or written by hand can be
used without a code-specific wrapper.

Canonical example
-----------------

.. literalinclude:: ../_downloads/hydroxylapatite.dat
   :language: text

The logical order is:

#. an optional one-line job description;
#. six rows defining the stiffness matrix;
#. for SEISMIC only, one positive density value.

The input does not declare units.  They are fixed by the format:

================= =====================
Quantity          Required unit
================= =====================
Stiffness matrix  GPa
Density           kg m\ :sup:`-3`
================= =====================

Job description
---------------

The first non-empty line is interpreted as a job description when its first
token is not numeric.  When the first token is numeric, Quantas assumes that
the matrix starts immediately and uses ``Unknown`` as the job name.

A description is free text and is stored as provenance.  It does not affect
symmetry detection, rotations, stability, or any numerical result.

Stiffness matrix
----------------

Quantas accepts three equivalent layouts.

Full symmetric matrix
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   C11 C12 C13 C14 C15 C16
   C21 C22 C23 C24 C25 C26
   C31 C32 C33 C34 C35 C36
   C41 C42 C43 C44 C45 C46
   C51 C52 C53 C54 C55 C56
   C61 C62 C63 C64 C65 C66

Upper triangular matrix
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   C11 C12 C13 C14 C15 C16
       C22 C23 C24 C25 C26
           C33 C34 C35 C36
               C44 C45 C46
                   C55 C56
                       C66

In an actual file the indentation is optional; the reader recognizes the row
lengths ``6, 5, 4, 3, 2, 1``.

Lower triangular matrix
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

                       C11
                   C21 C22
               C31 C32 C33
           C41 C42 C43 C44
       C51 C52 C53 C54 C55
   C61 C62 C63 C64 C65 C66

The recognized row lengths are ``1, 2, 3, 4, 5, 6``.  Triangular inputs are
expanded by exact reflection across the diagonal.

Voigt convention
----------------

The matrix uses the engineering Voigt order

.. code-block:: text

   1 -> 11
   2 -> 22
   3 -> 33
   4 -> 23
   5 -> 13
   6 -> 12

and the corresponding engineering shear-strain convention.  A matrix written
with a tensorial-shear convention must be converted before use.  Quantas does
not infer or repair a convention from the numerical values.

Symmetry requirements
---------------------

All entries must be finite floating-point values.  A full matrix must be
symmetric within the reader used by the selected workflow:

- Elasticity preserves the historical Frobenius-norm criterion
  ``||C - C.T|| <= 1e-3 GPa``;
- SEISMIC uses the stricter element-wise criterion
  ``max(abs(Cij - Cji)) <= 1e-8 GPa``.

The format therefore permits a file accepted by Elasticity to be rejected by
SEISMIC when a nominally symmetric full matrix contains small mismatches.
Writing a triangular matrix avoids this ambiguity because the missing half is
constructed exactly.

Density for SEISMIC
-------------------

SEISMIC requires one additional line after the six matrix rows:

.. code-block:: text

   3174.0

Only the first token is read.  The value must be finite and strictly positive.
Elasticity ignores no trailing density line by design; use the same shared file
when possible, but treat the density as part of the SEISMIC contract.

Comments and blank lines
------------------------

The SEISMIC reader ignores blank lines and full-line comments beginning with
``#``.  To keep one file portable between both workflows, avoid inserting
comments inside the six matrix rows and prefer a single descriptive first
line.

Generating inputs from external codes
--------------------------------------

Supported interfaces can generate the normalized text form:

.. code-block:: console

   quantas elasticity inpgen calculation.out --interface crystal
   quantas elasticity inpgen OUTCAR --interface vasp

The parser preserves the Cartesian frame reported by the source code.  User
rotations belong to the scientific workflow and are recorded in the result;
they should not be hidden by manually rotating only part of an input.

Thermoelastic Point analysis can also export a state directly in this format,
including density, so that the same file can be passed to both downstream
modules.

Validation checklist
--------------------

Before using a hand-written file, verify:

- exactly six matrix rows are present;
- the rows form one accepted full or triangular layout;
- values are in GPa;
- the Voigt and shear conventions are correct;
- a full matrix is symmetric to the stricter SEISMIC tolerance when the file
  will be shared;
- density is in kg m\ :sup:`-3` and corresponds to the same pressure,
  temperature, composition, and normalization cell as the stiffness tensor;
- coefficients under hydrostatic pre-stress are Wallace stress--strain
  coefficients rather than an uncorrected energy Hessian.

Related pages
-------------

- :doc:`../workflows/elasticity`
- :doc:`../workflows/seismic`
- :doc:`elasticity_seismic_hdf5`
