Thermoelastic input
===================

Thermoelasticity uses a schema-versioned YAML file describing a series of
hydrostatically pre-stressed elastic tensors at different volumes.  The file is
normally generated from CRYSTAL outputs and then coupled to a separate native
QHA result during calibration.

The current accepted schema is:

.. code-block:: yaml

   schema:
     name: quantas-thermoelastic-input
     version: '1.0'

Earlier unnormalized pre-release layouts are rejected.  Regenerate them with
``quantas thermoelasticity inpgen`` so that tensor-frame provenance is explicit.

Complete example
----------------

:download:`Download the dolomite example
<../_downloads/tutorials/thermoelasticity/dol_pbe0_thermoelastic.yaml>`

Top-level structure
-------------------

A current file contains:

.. code-block:: yaml

   schema: {...}
   job: Dolomite PBE0 QSA thermoelasticity
   method: quasistatic
   interface: crystal
   conventions: {...}
   units: {...}
   reference: {...}
   elastic_data:
     - {...}
     - {...}

``method`` must currently be ``quasistatic``.  ``interface`` and ``job`` are
provenance.  The scientific values are carried by ``reference`` and
``elastic_data``.

Conventions and units
---------------------

The generated file records the conventions used to interpret every tensor:

.. code-block:: yaml

   conventions:
     strain: eulerian-finite-strain
     stiffness: wallace-hydrostatic-stress-strain
     prestress: applied-by-crystal-pressure-keyword
     tensor_orientation: crystal
     voigt_order: [11, 22, 33, 23, 13, 12]

The normative units are:

================= =================
Quantity          Unit
================= =================
Pressure          GPa
Stiffness         GPa
Volume            Angstrom cubed
Density           kg m\ :sup:`-3`
Static energy     Hartree
Lattice vectors   Angstrom
================= =================

The reader validates numerical values according to these conventions; unit
labels are not used as permission to mix unit systems inside one file.

Reference block
---------------

``reference.index`` identifies one element of the sorted elastic series used to
define the common Cartesian frame.  The generator normally chooses the state
closest to zero pressure.  This frame reference is not necessarily the
thermodynamic EOS volume ``V0``.

``reference.structure`` stores:

- ``natom``;
- a ``3 x 3`` direct lattice with vectors by rows;
- ordered atomic numbers;
- ordered fractional coordinates.

``reference.symmetry`` stores the crystallographic interpretation used during
normalization:

- space-group number and international symbol;
- Hall number, Hall symbol, and setting choice;
- point group;
- detected elastic system;
- ``symprec`` and angular tolerance.

Atom order is part of the contract.  The volume series is not treated as a set
of unrelated cells that may be reordered independently.

Elastic data records
--------------------

Each entry in ``elastic_data`` represents one static state:

.. code-block:: yaml

   - source: calculation.out
     pressure: 5.0
     stress_pressure: 4.998
     volume: 104.2
     density: 2930.0
     energy: -1405.123456789
     lattice:
       - [a11, a12, a13]
       - [a21, a22, a23]
       - [a31, a32, a33]
     stiffness:
       - [C11, C12, C13, C14, C15, C16]
       - [...]
     frame: {...}

``source``
   Original electronic-structure output used to construct the record.

``pressure``
   Hydrostatic pressure requested through CRYSTAL's ``PRESSURE`` keyword and
   used for the pre-stress correction.

``stress_pressure``
   Pressure reconstructed from the final unstrained stress tensor.  A mismatch
   with ``pressure`` is a validation diagnostic.

``volume``
   Primitive normalization-cell volume.  It must be finite, positive, unique,
   and consistent with the determinant of ``lattice``.

``density``
   Density of the same primitive cell.  Across the series, ``density * volume``
   should reconstruct one constant cell mass within numerical tolerance.

``energy``
   Static DFT energy of the same primitive cell.  The reference EOS used by the
   thermoelastic calibration is taken from the QHA archive, but the SOEC energy
   remains valuable provenance and a consistency check.

``stiffness``
   Full symmetric ``6 x 6`` Wallace stress--strain tensor in GPa, using the
   engineering Voigt convention.

All records must have been calculated with the hydrostatic pre-stress
correction.  Quantas does not apply that correction a second time and does not
accept a series containing uncorrected points.

Frame metadata
--------------

Schema 1.0 requires a ``frame`` mapping for every state.  At minimum it contains:

``status``
   Normalization status, normally ``normalized``.

``method``
   Current generator value
   ``right_polar_decomposition_corotation``.

``rotation_to_reference``
   ``3 x 3`` orthogonal matrix used to co-rotate lattice and stiffness into the
   reference Cartesian frame.

``principal_logarithmic_strain``
   Three principal logarithmic strains after rigid rotation has been removed.

``source_lattice``
   Original lattice before co-rotation.

Generated records may additionally store removed rotation angle, maximum
ordered-atom displacement, and other diagnostics.  These values are provenance
and validation evidence; they are not extra fit parameters.

Series-level requirements
-------------------------

After parsing, records are ordered by increasing volume.  The series must
satisfy:

- at least one record is present;
- volumes are finite, positive, unique, and strictly increasing after sorting;
- every stiffness matrix is finite and symmetric;
- every density is finite and positive;
- lattice determinants agree with reported volumes;
- all points share species, ordered atoms, symmetry interpretation, and one
  co-rotated Cartesian frame;
- ``reference.index`` is valid;
- every point records applied hydrostatic pre-stress.

A large pressure range does not compensate for poor volume coverage around the
reference EOS minimum.  Fit support is diagnosed later and stored in the native
HDF5 result.

Relationship to the QHA file
----------------------------

The YAML contains the static elastic series only.  The QHA result supplied to
``thermoelasticity run`` provides:

- static volumes and energies used for the cold reference EOS;
- temperature and pressure grids;
- equilibrium volume ``V(P,T)``;
- optional volume uncertainty;
- heat capacity and Cartesian thermal-expansion tensor required for adiabatic
  conversion.

The normalization cells of the elastic YAML and QHA archive must be physically
compatible.  A numerically plausible tensor coupled to a differently normalized
volume is not a valid calculation.

Generating the file
-------------------

For a list of CRYSTAL outputs:

.. code-block:: console

   quantas thermoelasticity inpgen elastic_outputs.txt \
      --list --jobname "My thermoelastic series" \
      --output thermoelastic.yaml

The generator performs code-specific parsing, pressure checks, atom
correspondence, symmetry analysis, and frame co-rotation.  Hand editing should
be limited to provenance fields unless the user is prepared to reproduce all
validation steps.

Inspection from Python
----------------------

.. code-block:: python

   from quantas.api import thermoelasticity

   data = thermoelasticity.read_input("thermoelastic.yaml")
   series = data.elastic_series

   print(series.npoints)
   print(series.volume_bounds)
   print(series.pressures)
   print(series.stiffness.shape)   # (npoint, 6, 6)
   print(series.reference_index)

Validation checklist
--------------------

- The source calculations all use CRYSTAL ``PRESSURE``.
- Volumes bracket or closely surround the intended reference state.
- The reference atom ordering is preserved at every state.
- Frame metadata is present for every record.
- The same primitive-cell normalization is used by volume, energy, density,
  lattice, stiffness, and QHA.
- No rotation or pressure correction has been applied twice.
- The file was generated with the current schema rather than copied from an
  earlier pre-release snapshot.

Related pages
-------------

- :doc:`../workflows/thermoelasticity`
- :doc:`thermoelastic_hdf5`
- :doc:`earth_profile_spec`
