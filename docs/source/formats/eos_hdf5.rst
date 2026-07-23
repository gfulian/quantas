EOS native HDF5 archive
=======================

EOS persists datasets, every executed fit, user or workflow decisions, and the
currently accepted scientific results in a native HDF5 archive. The archive is
not a terminal log and does not depend on the batch, interactive, CLI, or GUI
frontend that created it.

The current EOS archive schema is ``1.1``. Quantas continues to
read schema ``1.0`` archives from earlier phases; missing selection nodes are
interpreted as an all-observation default mask.  The version increase prevents
older readers from silently discarding selection semantics.  Native
floating-point arrays are stored as IEEE-754 ``float64`` values, and the
numerical precision policy is written under ``/metadata/numerics``.

Scientific state model
----------------------

Three layers are deliberately separate.

Input datasets
   Raw values, normalized values, units, uncertainties, comments, provenance,
   source text, source hash, and coordinate classification are stored once.

Immutable fit records
   Every attempted fit has a monotonic global identifier. Successful, failed,
   and invalid-input results are retained. A fit record contains its request,
   initial/fixed/implied parameter declarations, solver settings, result,
   covariance, residuals, predictions, warnings, and provenance.

Accepted state
   Acceptance is not inferred from record order. A later rejected attempt does
   not replace an earlier accepted fit. The accepted state for each
   ``domain/target`` slot is a reference to a record and a compact materialized
   result copy.

An available property that has not been fitted is represented by an explicit
``not_processed`` slot. A property that is absent from the archived datasets is
not registered and is therefore reported as unavailable rather than unfitted. EOS does not create fake numerical fit results for such
properties. A deliberate decision to leave a property untreated can be stored
as a ``no_operation`` state event.

Schema overview
---------------

.. code-block:: text

   /metadata
       archive_kind = "quantas-eos-archive"
       program, version, module, method
       schema_version
       created_at, created_by
       /numerics

   /input
       /datasets/000001
           /source_text
           /raw_data
           /normalized_data
           /uncertainties
           /selection
               /default_mask
               /groups             # optional
           /units
           /raw_units
           /provenance
           /metadata
           /classification
           /column_order

   /session
       next_dataset_id, next_record_id, next_event_id
       /records/000001
           /request
           /result
           /provenance
       /state_events/000001
           /metadata
       /current/pv/volume
       /current/pv/a
       /current/vt/volume
       ...

   /results
       /accepted/pv/volume
           /request
           /result
       ...

   /events

Record identifiers are archive-global and zero-padded only in HDF5 group names.
The numerical identifier remains the integer ``1``, ``2``, and so on.

Input preservation
------------------

When the source file is available at archive creation, its UTF-8 text and a
SHA-256 digest are stored. This is in addition to the parsed raw and normalized
columns. Therefore reopening the archive does not require the original input
file to remain at the same filesystem path.

Both raw and normalized columns contain all reported uncertainty columns. The
``/uncertainties`` group provides a convenient dedicated view of the normalized
``sigma_*`` arrays without making it the authoritative copy.

``/selection/default_mask`` stores the input-level non-destructive selection
derived from ``USE`` and trailing ``*`` markers. ``/selection/groups`` is
present when the source ``FORMAT`` contains ``GROUP``. Every observation remains
in ``raw_data`` and ``normalized_data`` regardless of selection. A fit record
stores its resolved final mask in ``/request``; therefore two jobs can use
different subsets while sharing one immutable input dataset.

Fit-record immutability
-----------------------

A record is never modified after creation. Decisions are appended as state
events:

``record_created``
   A numerical attempt was persisted.

``record_accepted``
   The record became the accepted result for its slot. The event records the
   previously accepted record, when one existed.

``record_unaccepted``
   The current acceptance was explicitly revoked. The immutable fit and its
   earlier acceptance event remain in the history, while the materialized
   accepted-result view is removed.

``record_rejected`` and ``record_candidate``
   Scientific evaluation without deletion or replacement of the record.

``record_used_as_initial_guess``
   A later request reused this result. The child fit also records the source as
   ``parent_record_id``.

``note_added``
   A later annotation that does not mutate the creation note.

``no_operation``
   A registered property was deliberately left unfitted.

The mutable ``/session/current`` groups are compact indices. The append-only
records and state events remain the provenance authority.

A currently accepted record cannot be rejected directly. Acceptance must first
be revoked, then a separate rejection may be appended. This prevents one record
from being both the current accepted result and the latest explicit rejection.

Accepted-result materialization
-------------------------------

``/results/accepted/<domain>/<target>`` contains a compact copy of the selected
request and result plus the source ``record_id``. This controlled duplication
is intentional: readers and future calculators can access the current EOS
without replaying the complete session history. The immutable source record
remains authoritative.

Python API
----------

Create an archive containing input but no fits:

.. code-block:: python

   from quantas.api import eos

   dataset = eos.read_input("PV_topaz.dat")
   with eos.Archive.create("topaz.hdf5", dataset=dataset) as archive:
       print(archive.slot_state("pv/c").status)
       # EOSSlotStatus.NOT_PROCESSED

Persist and accept one script-driven fit:

.. code-block:: python

   from quantas.api.eos import WLSOptions
   from quantas.api import eos

   dataset = eos.read_input("PV_quartz.dat")
   request = eos.FitRequest(
       model="BM3",
       options=eos.FitOptions(solver_options=WLSOptions()),
   )
   result = eos.fit(dataset, request)

   with eos.Archive.create("quartz.hdf5", dataset=dataset) as archive:
       record = archive.store_fit(1, request, result)
       print(record.record_id)

Record a second attempt without replacing the accepted result:

.. code-block:: python

   second = archive.append_fit(
       1,
       second_request,
       second_result,
       parent_record_id=record.record_id,
   )
   archive.reject(second.record_id, note="model not retained")

   accepted = archive.accepted("pv/volume")
   assert accepted.record_id == record.record_id

The explicit ``append_fit``/``accept`` separation is the primitive used by a
session controller. ``store_fit`` is a convenience operation for scripts and
ordinary batch jobs that accept each successful result immediately.

Resumable session and inspection API
------------------------------------

The archive can be opened directly for update, but ordinary guided workflows
should use :class:`quantas.api.eos.Session`:

.. code-block:: python

   from quantas.api import eos

   with eos.Session.resume("quartz.hdf5") as session:
       inspection = session.inspect()
       print(inspection.record_count)

``EOSArchive.inspect()``, ``inspect_slot()``, and ``inspect_record()`` return
passive views that derive current record dispositions from state events. The
derived status is not written back to HDF5 and therefore cannot drift from the
provenance authority.

``EOSArchive.size_info()`` and the session controller expose the physical file
size. The session default warning threshold is 100 MiB. A warning is advisory:
the archive has no hard size ceiling, and the HDF5 schema remains valid above
the threshold. Large histories are expected only when many attempts contain
large prediction or residual arrays.

The calculator, diagnostics, CSV exporters, and static plotter derive their
outputs from immutable records without introducing additional persisted sources
of truth. A separate future Quantas GUI package will consume the same session,
inspection, calculator, and neutral plotting contracts; a prompt-driven
Click/Rich EOS session is not part of the current roadmap.


Batch manifest
--------------

A file produced by :command:`quantas eos run` stores the normalized declarative
plan under ``/session/batch_plan``. The manifest records job order, domain,
target, model, typed solver options, parameter constraints, acceptance policy,
failure policy, and CLI provenance. It is written before numerical execution.
Consequently an interrupted or partially failed batch still documents what was
requested, while immutable fit records document what was actually attempted.
The accepted-result indices remain independent of job order; replacement of a
previously accepted result must be declared explicitly by the plan.
