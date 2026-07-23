Native HDF5 results
===================

HDF5 is the authoritative machine-readable result format for complete Quantas
workflows.  A native file is not a dump of terminal text: it combines a common
result envelope with a module-specific scientific payload.

A reader must inspect stored metadata and units.  It must not infer the module
from the filename or assume that a dataset name has the same shape in every
schema version.

Common envelope
---------------

Single-shot modules use the following root layout:

.. code-block:: text

   /metadata
       attributes: program, module, method, version,
                   schema_version, created_at, created_by
       /numerics
           attributes: working, working_dtype, working_bits,
                       storage, storage_dtype, storage_bits

   /input
       attribute: source                       # optional
       /raw                                     # optional original text
       /data                                    # normalized recursive mapping

   /options                                     # resolved workflow options
   /results                                     # module-specific payload

   /diagnostics
       /warnings                                # optional UTF-8 array
       /report_text                             # optional embedded plain report
       ... module-specific diagnostic groups

   /events
       attribute: count
       /0, /1, ...                              # persistent event records

EOS uses a session-capable specialization documented in :doc:`eos_hdf5`.

Metadata
--------

``/metadata`` identifies the producer and the responsible reader.

``program``
   Normally ``quantas``.

``module``
   Stable identifier such as ``ha``, ``qha``, ``elasticity``, ``seismic``,
   ``eos``, or ``thermoelasticity``.

``method``
   Scientific workflow or approximation used by the module.

``version``
   Quantas version that wrote the file.

``schema_version``
   Version of the common result envelope.  Module payloads may also store an
   additional local schema version.

``created_at`` and ``created_by``
   Creation timestamp and producer identity when available.

External scripts should reject unsupported schema versions explicitly.  A
newer Quantas version may read selected historical layouts through a controlled
migration path, but a generic h5py script has no such guarantee.

Numerical precision
-------------------

Quantas scientific calculations and native floating-point datasets use IEEE-754
``float64``.  Complex scientific values use ``complex128`` when required.  No
runtime precision selector is exposed.

The fixed policy is recorded under ``/metadata/numerics``.  Historical files
without this group are interpreted by Quantas according to the original
double-precision default.

Display precision is independent.  Rounding in a terminal table, report, CSV,
or plot label does not alter the HDF5 value.

Normalized input
----------------

``/input`` may contain:

``source`` attribute
   Original path or source description.

``raw`` dataset
   Original text when the reader retained it.

``data`` group
   Parsed, normalized input represented as nested groups, datasets, and
   attributes.

The normalized input is stored for provenance and reproducibility.  It is not
necessarily byte-identical to the original source and should not be treated as
a replacement for code-specific archival files when those files contain
additional information.

Options
-------

``/options`` stores the resolved scientific and operational options used by the
calculator.  Nested dataclasses, mappings, sequences, enums, paths, and scalar
values are serialized recursively.

Quantas may store scalar strings and numbers as HDF5 attributes rather than
one-element datasets.  A generic reader must therefore inspect both child nodes
and attributes.

Diagnostics and report text
---------------------------

``/diagnostics/warnings`` stores meaningful warnings collected by the workflow.
Module-specific fit records, failed states, or quality metrics may appear as
additional children.

``/diagnostics/report_text`` is optional.  When present, it is a deterministic
plain-text rendering created by a frontend.  It is a derived view, not the
authoritative numerical payload.

Events
------

Persistent workflow events are stored as numbered groups.  Each record carries:

- ``message``;
- ``level``;
- ``timestamp``;
- optional ``progress``;
- structured ``data``.

Operational progress events are observer-only and are normally not persisted.
Warnings, errors, meaningful workflow messages, and structured result events
may be stored.

Recursive mapping representation
--------------------------------

Quantas uses a neutral recursive representation for metadata, options, and
small structured records:

- mappings become groups;
- NumPy arrays become datasets;
- strings and scalar numbers normally become attributes;
- heterogeneous sequences become numbered child nodes;
- ``None`` is represented by an internal kind marker;
- homogeneous sequences may be stored directly as datasets.

For this reason, a script that only walks datasets will miss scalar attributes.
The inspection examples in :doc:`hdf5_inspection` show how to report both.

Units and descriptions
----------------------

Module payload writers attach ``unit`` attributes to many scientific datasets
and may attach a human-readable ``description``.  Some older or recursively
serialized fields keep unit information in a neighboring metadata mapping
instead.

A robust external analysis should:

#. read the file and module schema;
#. read the dataset ``unit`` attribute when present;
#. consult module metadata when the dataset has no local unit;
#. convert units explicitly in the analysis script;
#. preserve the original unit in derived output metadata.

Missing and invalid values
--------------------------

Different states have different meanings:

- an absent optional dataset means the quantity was not calculated or was not
  available from the source workflow;
- a present boolean validity mask describes which stored states are usable;
- ``NaN`` may represent an explicitly invalid numerical state while preserving
  array shape;
- an EOS ``not_processed`` slot is a session state, not a NaN-valued fit.

Do not replace these cases by zeros.  Read the corresponding validity,
extrapolation, completion, or acceptance fields.

Compression and chunking
------------------------

Large module arrays may use gzip compression, chunking, and byte shuffling.
These are storage details and do not change numerical values.  h5py transparently
returns ordinary NumPy arrays.

Preferred public reading path
-----------------------------

For ordinary analysis, use metadata-driven dispatch through the public API:

.. code-block:: python

   from quantas.api import registry

   result_or_archive = registry.open_result("result.hdf5")

For HA, QHA, Elasticity, SEISMIC, and Thermoelasticity the return value is a
:class:`quantas.api.common.ResultData`.  EOS returns an active archive object
that must be closed.

A module-specific reader gives stronger typing:

.. code-block:: python

   from quantas.api import qha

   envelope = qha.read_result("mgo_qha.hdf5")
   result = qha.get_result(envelope)

   print(envelope.metadata.module)
   print(result.temperature.shape)
   print(result.equilibrium_volume.shape)

Direct h5py access
------------------

Direct access is appropriate when:

- Quantas is unavailable in the analysis environment;
- the user needs only a few arrays;
- data will be consumed by another language;
- the schema itself is being audited.

The public module reader remains safer because it validates metadata, units,
shapes, and supported migrations.  See :doc:`hdf5_inspection` for complete
examples of both approaches.

Writing or modifying native files
---------------------------------

Do not edit native Quantas HDF5 files in place with h5py unless developing a
schema migration.  Changing a numerical dataset without updating options,
metadata, masks, diagnostics, and provenance creates an internally inconsistent
archive.

For derived analysis:

- leave the original file unchanged;
- write a new HDF5/NetCDF/CSV/NPZ product;
- record the source file path, checksum, module, schema version, dataset paths,
  units, and transformation applied.

Use the public ``write_result`` function when constructing a new Quantas result
envelope from a valid typed result.

Module payload references
-------------------------

- :doc:`ha_qha_hdf5`
- :doc:`elasticity_seismic_hdf5`
- :doc:`eos_hdf5`
- :doc:`thermoelastic_hdf5`
