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

For CRYSTAL, complex general-q vectors are reconstructed from in-phase and
anti-phase components and converted to unit-norm mass-weighted directions before they
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

VASP run documents
------------------

The generic VASP interface treats one calculation directory as one run source.
:func:`quantas.interfaces.vasp.resolve_vasp_run_source` accepts the directory
itself, ``vasprun.xml``, or ``OUTCAR`` and resolves the sibling files.  During
the b13 interface-maintenance tranche, ``vasprun.xml`` is the required primary
structured document and ``OUTCAR`` is optional complementary evidence.  An
``OUTCAR`` without the sibling XML file is therefore not silently promoted to a
complete run source.

:class:`quantas.interfaces.vasp.document.VaspRunDocument` owns XML/text syntax.
It exposes generator metadata, explicitly recorded INCAR values, effective
scalar parameters, atom ordering, and ionic-state containers without deciding
which state or energy should feed a Quantas scientific workflow.
:class:`quantas.interfaces.vasp.output.VaspOutputParser` converts these records
to canonical :class:`quantas.models.structures.CrystalStructure` objects and
VASP-specific ionic-state records containing:

* lattice vectors in angstrom and fractional coordinates as ``float64``;
* atomic numbers in the exact VASP atom order;
* ``e_fr_energy``, ``e_wo_entrp``, and ``e_0_energy`` as separate values in eV;
* forces in eV/angstrom;
* stress tensors in kbar, without a premature sign or pressure conversion;
* electronic/ionic convergence facts that can be established explicitly;
* source and resolution provenance.

The parser accepts the older layout in which each ionic state is enclosed by a
``<calculation>`` element and the current documented flat ionic-state layout.
Support for a layout means that Quantas understands its structure; scientific
validation against a specific VASP version still requires a real reference
output for that version.

VASP energy semantics require special care.  In VASP 5.4.4 the outer
``calculation/energy`` record is affected by a documented output bug: the
``e_wo_entrp`` and ``e_0_energy`` tags can contain the extrapolated energy and
electronic entropy term, respectively, instead of their nominal quantities.
The b13 parser does not branch on a hard-coded version string.  For
``<calculation>``-style output it instead takes the relative values of ``F``,
``E``, and ``E0`` from the final electronic ``scstep`` and transfers only the
shift in the outer ``e_fr_energy``.  This preserves an additive correction that
is present only in the ionic-state total while avoiding the mislabeled outer
tags.  The raw outer values, the applied shift, and whether the known VASP-5
pattern was observed remain in metadata.  When an ``OUTCAR`` can be paired
unambiguously, the resolved XML energies are checked against its final
``TOTEN``, ``energy without entropy``, and ``energy(sigma->0)`` values.

This resolution policy follows the `VASP developers' description of the VASP
5.4.4 XML issue <https://vasp.at/forum/viewtopic.php?t=17839>`_ and its
correction in VASP 6.  It is an interface-level source correction; workflow
adapters select scientific quantities only after the three VASP energy values
have been resolved.

VASP elasticity adaptation
~~~~~~~~~~~~~~~~~~~~~~~~~~

The VASP elasticity reader accepts either a calculation directory, ``OUTCAR``,
or ``vasprun.xml`` with a sibling ``OUTCAR``.  Elastic moduli are read from the
human-readable OUTCAR because that is where VASP 5.4.4 reports the complete
finite-difference elasticity decomposition.  Quantas preserves separately:

* ``SYMMETRIZED ELASTIC MODULI`` (clamped-ion);
* ``ELASTIC MODULI CONTR FROM IONIC RELAXATION`` when present;
* ``TOTAL ELASTIC MODULI`` (relaxed-ion), selected by default when available.

VASP labels the six components ``XX YY ZZ XY YZ ZX``.  The parser maps both
matrix axes explicitly by label to Quantas' ``11 22 33 23 13 12`` convention;
it does not use positional shear swaps.  The first stress block is retained as
the unstrained reference stress, with VASP's positive-compression sign
convention, and the first reported cell volume is used for density because
later ``IBRION=6`` records contain trial lattice distortions.

The interface also exposes :func:`quantas.interfaces.vasp.read_vasp_elastic_series`
to collect several VASP elastic calculations into the shared
:class:`~quantas.models.elastic_states.ElasticStateSeries` contract.  This raw
adapter remains factual: it sorts states by primitive-cell volume, preserves
the selected VASP stiffness, density, and reference stress/pressure provenance,
and applies **no** finite-prestress correction.  Raw VASP states are classified
as ``raw_stress_strain`` and therefore continue to fail the shared incremental
stiffness gate.

For a genuinely hydrostatic VASP reference stress, Quantas now provides an
explicit second step via
:func:`quantas.interfaces.vasp.convert_vasp_hydrostatic_elastic_series`.  The
conversion follows Appendix A of Singh et al., *MechElastic* (Computer Physics
Communications 267, 108068, 2021): pressure is subtracted from all
six Voigt diagonal terms and added to ``C12``, ``C13`` and ``C23`` (and their
symmetric partners).  Pressure is positive in compression.  The complete
unstrained stress tensor is checked against ``P I`` before this scalar
hydrostatic adjustment is allowed; appreciable deviatoric stress is rejected.
The corrected state is then labelled ``wallace_hydrostatic`` and records the
raw ``raw_stress_strain`` source kind and correction DOI in provenance.

This is a VASP-specific ingestion rule.  The CRYSTAL Erba/Barron--Klein
energy--strain transformation and Quantas' generic Eulerian finite-strain
operator are not reused.  Pressure selection remains a separate operation from
tensor conversion.  :func:`quantas.interfaces.vasp.assign_vasp_manual_pressures`
can replace the raw output-stress pressure with explicit hydrostatic values, and
:func:`quantas.interfaces.vasp.resolve_vasp_energy_derived_pressures` reuses the
backend-neutral E(V) pressure fitter and volume matcher.  Both operations leave
the VASP stiffness matrix raw; only the explicit hydrostatic conversion marks
the tensor incremental.  The original VASP reference stress is still retained
and must itself be hydrostatic, while the selected correction pressure and its
difference from the VASP output pressure are recorded in provenance.  No
Kieffer/HA/QHA coupling is introduced at this checkpoint.

VASP Gamma phonon adaptation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`quantas.interfaces.vasp.phonons.VaspPhononReader` currently exposes
only primitive-cell Gamma phonons from a completed VASP run.  Frequencies and
the real/imaginary mode label are read from the human-readable ``OUTCAR``;
higher-precision normalized eigenvectors are read from
``vasprun.xml/dynmat/eigenvectors``.  The interface stores those vectors as the
backend-neutral unit-norm mass-weighted representation used by
:class:`~quantas.models.phonons.PhononModeData`.

The three Gamma translations are not identified by a fixed frequency cutoff.
Quantas projects every VASP eigenvector onto the three-dimensional
mass-weighted rigid-translation subspace.  Exactly three well-resolved
translations are required.  Their raw VASP frequencies and projection scores
remain in provenance, while the thermodynamic frequencies are set to exactly
zero so numerical acoustic-sum-rule drift is not counted as a physical
harmonic oscillator.  Other imaginary modes retain negative frequencies.

.. important::

   Direct phonon dispersion from VASP output is **not implemented yet**.  The
   current VASP reader requires the calculation cell itself to be primitive and
   exposes a single q-point, ``Gamma = (0, 0, 0)``, with unit weight and an
   identity phonon-supercell matrix.  A reducible source cell is rejected
   rather than treating folded supercell Gamma modes as a primitive-cell
   dispersion.  VASP 6 ``LPHON_DISPERSION``/``QPOINTS`` output is likewise
   outside the current parser contract.

The initial undistorted VASP state supplies the static
``energy(sigma->0)`` value for HA/QHA input generation.  Its source unit remains
``eV``; phonon frequencies are exposed in ``cm^-1`` and structural quantities
in angstrom.  The first validated characterization target is VASP 5.4.4
``IBRION=6`` output for primitive MgO.

VASP Energy EOS adaptation
~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`quantas.interfaces.vasp.energy_volume.VaspEnergyVolumeReader` adapts one
completed VASP run to one backend-neutral
:class:`~quantas.models.computation.StructureEnergyPoint`.  The current b13
scope is a zero-electronic-temperature/static E--V dataset, so the adapter
selects VASP ``e_0_energy`` / ``energy(sigma->0)``.  This is not a statement
that ``E0`` is the appropriate quantity for every VASP workflow: a calculation
that intentionally represents finite electronic temperature has different
thermodynamic semantics.  For accurate bulk total-energy calculations VASP
recommends the tetrahedron method with Blöchl corrections (``ISMEAR=-5``);
when Gaussian or Methfessel--Paxton smearing is used, ``energy(sigma->0)`` is
an extrapolation and convergence with respect to ``SIGMA`` remains the user's
scientific responsibility.  See the `VASP smearing guidance
<https://www.vasp.at/wiki/index.php/Smearing_technique>`_.

Each Energy EOS source must resolve to exactly one ionic state.  An optimization
history is therefore never flattened into the E--V series.  The source cell is
passed through the shared spglib primitive-cell normalization.  If the VASP
cell is already primitive, its lattice basis and orientation are preserved.  If
it contains multiple primitive repetitions, structure and volume are reduced
and the source-cell energy is divided by the same integer multiplicity.  This
produces the same primitive normalization expected by the backend-neutral EOS
collector; ``--crystal-reference crystallographic`` may subsequently scale the
complete series to one fixed crystallographic cell.

Before independent runs are merged, Quantas requires one explicit VASP energy
compatibility signature.  It records the selected quantity, ``ISMEAR`` and
``SIGMA`` where active, relevant exchange-correlation/spin/charge/hybrid/DFT+U
settings, plane-wave precision/cutoff, Brillouin-zone sampling, and the
pseudopotential labels stored in ``vasprun.xml``.  A mismatch is rejected rather
than silently mixing energies computed on different electronic surfaces.  The
check is intentionally conservative; a future workflow may relax individual
fields only with an explicit scientific policy and characterization tests.

The first characterization fixture is MgO/periclase calculated with VASP
``5.4.4.18Apr17-6-g9f103f2a35``.  Full user calculations were used to verify
three consecutive cell optimizations and seven fixed-cell EOS states; compact
fixtures retain the real generator, INCAR, atom, structure, electronic-energy,
force, stress, and convergence records needed by the repository tests.

CRYSTAL static-energy semantics
-------------------------------

CRYSTAL distinguishes the converged electronic SCF energy from the physical
total energy used when a-posteriori corrections are active.  Quantas preserves
both quantities at the interface boundary:

``SCF energy``
   The electronic energy printed on ``SCF ENDED - CONVERGENCE ON ENERGY`` and
   repeated by ``TOTAL ENERGY(DFT)(AU)``.

``total energy``
   The corrected energy printed by CRYSTAL when available.  Recognized forms
   include ``TOTAL ENERGY + DISP (AU)``, ``TOTAL ENERGY + GCP (AU)``, and
   ``TOTAL ENERGY + DISP + GCP (AU)``.  If CRYSTAL prints no corrected total,
   the total energy is identical to the SCF energy.

The generic :class:`quantas.interfaces.crystal.output.CrystalOutputParser`
resolves these values state by state and does not attach a correction printed
for one SCF calculation to a later state.  The corrected total printed by the
backend is authoritative; the interface does not reconstruct it by summing empirical
components.  Correction labels and the difference between total and SCF
energy are retained as provenance.

Scientific workflows consume the resolved **total energy**.  In particular,
CRYSTAL phonon input generation continues to use the ``CENTRAL POINT`` energy,
which is the total energy attached by CRYSTAL to the undisplaced reference
configuration.  When it can be matched to the preceding SCF state, the
uncorrected SCF energy is additionally retained in input provenance.  CRYSTAL
elastic readers likewise expose ``scf_energy`` and ``total_energy`` while the
historical ``energy`` property is an alias for ``total_energy``.

CRYSTAL Energy EOS state series
-------------------------------

:class:`quantas.interfaces.crystal.energy_volume.CrystalEnergyVolumeReader`
normalizes three CRYSTAL output shapes to one structure--energy contract: a
static SCF result contributes one state, a completed ordinary geometry
optimization contributes its final state, and a completed native ``EOS`` run
contributes all states in the final sorted E(V) table.  The reader reuses
:class:`quantas.interfaces.crystal.output.CrystalOutputParser` for authoritative
total-energy resolution and :class:`quantas.interfaces.crystal.geometry.CrystalGeometryParser`
for structures; it does not duplicate correction-specific regular expressions.

For native EOS output, the sorted E(V) table is checked against independently
parsed final optimized geometries and state-resolved total energies.  The
interface returns :class:`quantas.models.computation.StructureEnergySeries`;
source-list flattening and cross-file compatibility belong to the EOS workflow
layer.  CRYSTALpytools may be used externally as an audit reference but is not a
Quantas runtime dependency.

CRYSTAL elastic volume series
-----------------------------

:func:`quantas.interfaces.crystal.read_crystal_elastic_series` composes a set
of completed ELASTCON or ELAPIEZO outputs into the backend-neutral
:class:`quantas.models.elastic_states.ElasticStateSeries` contract.  This is
the interface boundary used by Kieffer and available to other multi-volume
elastic workflows; it does not create a Kieffer-specific elastic format.

The importer requires finite volume, density, total static energy, and stiffness at
every state.  It sorts the resulting states by increasing volume and selects
the minimum-static-energy state as the reference.  Tensor axes remain in the
CRYSTAL Cartesian frame.

For each elastic output, the structural state is the unstrained reference used
to generate the elastic distortions.  Quantas therefore restricts structure,
static energy, density, and output-stress pressure collection to the part of
the CRYSTAL output preceding the first ``STRAIN MATRIX``.  This distinction is
important when ``COORPRT`` causes geometries from later strained or internally
relaxed configurations to be printed.  The selected lattice must also agree
with the primitive-cell volume reported by the elastic module before it can be
attached to an :class:`~quantas.models.elastic_states.ElasticState`.

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

By default, raw CRYSTAL energy--strain tensors are converted once with the
finite-pressure transformation implemented by CRYSTAL itself
[#erba_mahmoud_belmonte_dovesi_2014]_:

.. math::

   B_{ijkl}=C_{ijkl}+\frac{P}{2}
   \left(2\delta_{ij}\delta_{kl}-\delta_{il}\delta_{jk}-\delta_{ik}\delta_{jl}\right).

In CRYSTAL Voigt order this leaves ``C11``, ``C22``, and ``C33`` unchanged,
adds ``+P`` to ``C12``, ``C13``, and ``C23``, and adds ``-P/2`` to the three
shear diagonals.  This interface conversion is deliberately distinct from the
Eulerian finite-strain Wallace term used internally by the QSA model.  The
pressure value, source, method, source tensor kind, and software applying the
correction are retained in each state. Passing a non-auto pressure policy for
a tensor already corrected by CRYSTAL is an error, preventing an accidental
second correction.


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

The reusable fit operations live in :mod:`quantas.core.physics.eos`.  The
backend-neutral
:func:`quantas.core.physics.elasticity.resolve_energy_derived_pressures`
service combines the selected E(V) fit with explicit volume matching and
pressure assignment while leaving the raw stiffness coefficients unchanged.
Both Kieffer enrichment and thermoelastic input generation use this same
pressure-resolution path.  The subsequent hydrostatic tensor correction
remains interface-specific, so CRYSTAL conventions do not leak into the shared
core.  This boundary lets tests verify that ``P(V)`` is attached to an
unmodified raw tensor before the tensor is corrected exactly once.

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


.. include:: ../_generated/references/developer_interfaces.inc
