Native HDF5 persistence and schema evolution
============================================

HDF5 is the native Quantas result format. Persistence is part of the scientific
contract: a file must preserve the values and metadata needed to reconstruct
and interpret a result independently of the frontend that produced it.

Shared envelope
---------------

Current single-shot results use schema version ``2.0`` with the top-level
groups:

.. code-block:: text

   /metadata
   /metadata/numerics
   /input
   /options
   /results
   /diagnostics
   /events

``metadata``
   Program, package version, scientific module, method, creation information,
   and shared schema version.

``metadata/numerics``
   Fixed ``float64``/``complex128`` working and storage policy.

``input``
   Source, optional raw content, and normalized input mapping.

``options``
   Resolved workflow options.

``results``
   Module-specific payload.

``diagnostics``
   Warnings, optional deterministic report text, and module diagnostics.

``events``
   Persistent non-progress event records.

Module-specific scientific payloads remain in the module package. The shared
I/O layer handles recursive mappings, common metadata, options, diagnostics,
and events.

Writer structure
----------------

A module writer normally:

#. validates that the envelope belongs to the expected module;
#. opens the destination with an explicit overwrite policy handled by the
   caller;
#. writes shared metadata and precision information;
#. writes normalized input and resolved options;
#. writes the module payload with units and descriptions;
#. writes diagnostics and persistent events;
#. closes the file before returning the final :class:`pathlib.Path`.

Use shared helpers from ``quantas.io.hdf5`` for common values. Keep payload
layout decisions inside the module.

Dataset requirements
--------------------

A numerical dataset should define, where applicable:

* stable path and name;
* dtype;
* shape and axis order;
* physical unit;
* human-readable description;
* normalization basis;
* masks or status values aligned with the data;
* compression/chunking chosen for access pattern rather than appearance.

Do not store a formatted table as the only representation of scientific values.

Recursive values and ``None``
-----------------------------

The shared mapping writer supports dataclasses, mappings, sequences, enums,
paths, timestamps, scalars, and arrays. ``None`` is represented explicitly
through an internal kind marker rather than being confused with zero or an
empty string.

For module payloads with scientifically important optional fields, an explicit
``_is_none`` attribute or equivalent schema marker may be preferable. The
reader and format documentation must agree.

Reader structure
----------------

A current-schema reader:

#. opens the file;
#. validates program, shared schema, and expected module;
#. reads shared metadata, input, options, warnings, and events;
#. validates module payload version when present;
#. validates required datasets, units, shapes, and enum values;
#. reconstructs the typed result;
#. returns a complete ``ResultData`` envelope.

Readers inspect metadata. They do not infer file type from the extension or
filename.

Package version versus schema version
-------------------------------------

The Quantas package version describes Python code and public API. The shared
result schema and module payload schemas describe persisted data. They change
for different reasons and must not be coupled mechanically.

Increment a schema when the stored structure or scientific interpretation
changes incompatibly. Adding an optional field that old readers may ignore can
be backward compatible; changing units, normalization, axis order, or meaning
is not.

Migration policy
----------------

New writers emit only the current schema. Readers may support older layouts
when migration is unambiguous and tested.

A migration should:

* identify the historical layout explicitly;
* document every assumption;
* convert values and units without display rounding;
* add a warning or migration metadata;
* preserve source information;
* have a frozen historical fixture and expected current result.

Do not guess a missing normalization or unit when more than one interpretation
is plausible. Reject the file with an actionable error.

Module payload versions
-----------------------

A module may maintain its own payload version under the shared envelope. Use
this when the module evolves independently but the top-level envelope remains
compatible.

Document:

* current writer version;
* older versions accepted by the reader;
* unsupported versions;
* migration rules;
* fields added or reinterpreted.

EOS archive exception
---------------------

EOS uses a persistent archive containing datasets, immutable fit records,
accepted/candidate state, and history. It is not a conventional single-shot
``ResultData`` payload. Its reader, writer, and session lifecycle remain
module-specific, while metadata-based dispatch still identifies it as a
Quantas EOS archive.

Round-trip testing
------------------

For each module test:

.. code-block:: text

   construct/run result
       -> write HDF5
       -> read HDF5
       -> compare typed result

Compare:

* metadata needed for dispatch;
* arrays and dtypes;
* masks;
* units and axis semantics;
* optional fields;
* warnings and persistent events;
* covariance and diagnostics;
* source/option provenance.

Also test malformed files:

* wrong module;
* unsupported schema;
* missing required dataset;
* incorrect shape;
* wrong unit;
* invalid enum/status value.

User inspection and derived files
---------------------------------

Users may inspect native files with :mod:`quantas.api.registry` or ``h5py``.
They should not modify a native result in place. A custom analysis should read
the file and write a new derived artifact with its own provenance.
