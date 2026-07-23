HA and QHA phonon YAML
======================

HA and QHA share one normalized YAML contract for static energies,
volume-dependent phonon frequencies, q-point integration data, and optional
structural information.  The same file can be used by both modules; QHA adds a
mode-continuity declaration because mode-resolved interpolation has stricter
requirements than harmonic summation.

A complete current example is available here:

:download:`Download the MgO example <../_downloads/mgo_b3lyp.yaml>`

The public readers return :class:`quantas.api.ha.Input` or
:class:`quantas.api.qha.Input`.  Arrays are validated and converted to
``float64`` before calculation.

Minimal contract
----------------

A compact multi-volume file has this logical form:

.. code-block:: yaml

   job: Example material
   natom: 2
   formula_units: 1
   qpoints: 2
   supercell:
     - [2, 0, 0]
     - [0, 2, 0]
     - [0, 0, 2]

   volume: [18.0, 19.0, 20.0]
   energy: [-100.01, -100.03, -100.02]
   mode_continuity: assumed

   phonon:
     - q-position: [0.0, 0.0, 0.0]
       weight: 1
       band:
         - frequency: [0.0, 0.0, 0.0]
         - frequency: [0.0, 0.0, 0.0]
         - frequency: [0.0, 0.0, 0.0]
         - frequency: [200.0, 190.0, 180.0]
         - frequency: [300.0, 285.0, 270.0]
         - frequency: [500.0, 480.0, 460.0]
     - q-position: [0.5, 0.0, 0.0]
       weight: 3
       band:
         # six frequency records, each with three volume values

Required top-level fields
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 22 54

   * - Field
     - Normalized shape
     - Meaning
   * - ``natom``
     - scalar integer
     - Number of atoms in the thermodynamic normalization cell.  The number of
       branches is exactly ``3 * natom``.
   * - ``formula_units``
     - scalar integer
     - Number of chemical formula units in that cell.  Historical aliases
       ``z`` and ``Z`` are accepted; omission defaults to one.
   * - ``qpoints``
     - scalar integer
     - Number of entries in ``phonon``.
   * - ``volume``
     - ``(nvol,)``
     - Volumes of the thermodynamic normalization cell.
   * - ``energy``
     - ``(nvol,)``
     - Static electronic energies for the same volume sequence.
   * - ``phonon``
     - sequence of ``qpoints`` mappings
     - Q-point weights, optional coordinates, and exactly ``3 * natom`` band
       records per q-point.

``job`` and ``supercell`` are strongly recommended but are not required by the
basic reader.  ``job`` is provenance.  ``supercell`` records the expansion used
by the upstream phonon calculation and is also useful for structural
normalization.

Units
-----

The YAML values do not carry per-field unit tags.  Their units are supplied by
the HA/QHA calculation options.  The public defaults are:

================= =========================
Quantity          Default interpretation
================= =========================
Energy            Hartree per cell
Length/volume     Angstrom / Angstrom cubed
Frequency         cm\ :sup:`-1`
================= =========================

Changing an input-unit option changes how the same numerical values are
interpreted; it is not a display-only conversion.  All volumes, energies, and
frequencies in one file must use one consistent basis and unit system.

Volume and energy arrays
------------------------

``volume`` and ``energy`` must describe the same ordered states.  A scalar is
accepted for single-volume HA and is normalized to a length-one array.
Multi-volume HA retains one independent harmonic result for every sampled
volume.  QHA requires enough distinct volumes to support the selected
interpolation and minimization models.

The cell basis must be consistent across all fields:

- ``natom`` and ``formula_units`` describe the same cell as ``volume``;
- ``energy`` is the static energy of that cell;
- the number of modes is ``3 * natom``;
- optional structural lattices must have determinants equal to the same
  volumes, within their validation tolerance.

Q-point records
---------------

Each item under ``phonon`` has:

``weight``
   Required integration weight.  It must be finite and the sum of all weights
   must be positive.

``q-position``
   Optional fractional q-point coordinates with three components.  Coordinates
   are stored for provenance and mode analysis; they are not used to infer
   missing weights.

``band``
   A sequence of exactly ``3 * natom`` mappings.  Each mapping contains
   ``frequency`` as either a scalar for one volume or a sequence of length
   ``nvol``.

After reading, frequencies have shape:

.. code-block:: text

   (nqpoint, nmode, nvolume)

and weights have shape ``(nqpoint,)``.

Q-point weights
---------------

Quantas normalizes the supplied weights by their sum before evaluating the
thermodynamic sums.  It does not currently derive irreducible-Brillouin-zone
multiplicities from a separate crystal-symmetry analysis.

Consequently:

- a full uniform mesh may use equal weights;
- an irreducible mesh requires the symmetry multiplicities supplied by the
  upstream calculation;
- multiplicities must not be guessed from q-point coordinates alone;
- a common multiplicative factor is harmless because the weights are
  normalized, but relative weights are scientifically significant;
- Phonopy mappings may carry their native ``weight`` values into the Quantas
  representation.

A weight error affects every extensive vibrational quantity while leaving the
file syntactically valid.  Weight provenance should therefore be retained with
the calculation.

Frequency records and non-positive modes
----------------------------------------

A frequency series occupies one stable band index across all volumes.  Zero or
negative values are preserved by the input contract; their scientific handling
belongs to the HA/QHA workflow and is reported there.  Do not delete acoustic
or imaginary modes merely to make the YAML rectangular.  The file should
represent the upstream calculation faithfully, and the calculation report
should document any excluded contributions.

Mode continuity for QHA
-----------------------

``mode_continuity`` accepts four values:

``verified``
   Branch correspondence was established by a documented tracking or mapping
   procedure.

``assumed``
   The producer assumes that mode indices remain physically corresponding, but
   no explicit verification was stored.

``unknown``
   Continuity has not been assessed.

``unreliable``
   Known crossings, reordering, or reconstruction failures make mode-resolved
   interpolation unsafe.

The ``freq`` QHA scheme and mode Gruneisen properties require meaningful
mode-by-mode correspondence.  The ``td`` scheme interpolates already integrated
thermodynamic quantities and does not require branch continuity.

Optional structural path
------------------------

Current inputs may include a ``structure`` mapping.  This block allows QHA to
reconstruct equilibrium lattice parameters, axial expansion, and a Cartesian
thermal-expansion tensor along the sampled volume path.

The principal children are:

``reference``
   Reference primitive lattice, fractional positions, atomic numbers, label,
   and provenance metadata.

``volume_series``
   ``volume``, ``lattice``, and ``fractional_positions`` arrays for every
   sampled state.  Lattices have shape ``(nvol, 3, 3)`` and fractional
   positions have shape ``(nvol, natom, 3)``.

``normalization``
   Basis, source basis, integer expansion matrix, number of repetitions, source
   atom count, and normalized atom count.

``symmetry``
   Optional space-group and Hall metadata, symmetry tolerance, equivalent-atom
   mapping, and transformations.

``reconstruction``
   Optional diagnostics describing how a compact primitive path was recovered
   from an upstream supercell representation.

``transformations``
   Optional matrices such as ``primitive_to_crystallographic``.

``orientation`` and ``reference_index``
   Cartesian convention and the sampled state used as structural reference.

A historical YAML without ``structure`` remains valid for HA and scalar QHA
thermodynamics; structural and axial outputs are then unavailable rather than
invented.

Input generation
----------------

Use module-specific input generators instead of hand-copying long frequency
arrays.  They preserve upstream metadata and apply the supported normalization
path.  Typical commands are documented in the HA and QHA command references
and tutorials.

Validation checklist
--------------------

Before a production run, verify:

- ``len(volume) == len(energy)``;
- every band frequency record has the same number of volume values;
- every q-point has exactly ``3 * natom`` bands;
- ``qpoints == len(phonon)``;
- all weights are finite and their sum is positive;
- the cell normalization is consistent across energy, volume, modes, and
  structure;
- the selected input units match the numbers in the file;
- ``mode_continuity`` is honest for the intended QHA scheme;
- structural lattices and volumes describe the same cell.

Related pages
-------------

- :doc:`../workflows/ha`
- :doc:`../workflows/qha`
- :doc:`ha_qha_hdf5`
