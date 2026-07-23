Events, observers, and progress
===============================

Quantas events decouple workflow communication from presentation. Calculators
emit structured events; observers decide whether to display, collect, forward,
or ignore them.

Event contracts
---------------

An operational :class:`quantas.core.events.Event` contains:

* message;
* :class:`~quantas.core.events.EventLevel`;
* optional progress in the interval ``[0, 1]``;
* lightweight structured data;
* timestamp.

Supported levels are:

``DEBUG``
   Detailed diagnostic information.

``INFO``
   Normal workflow information.

``WARNING``
   A non-fatal condition requiring attention.

``ERROR``
   Failure information emitted before an exception is propagated.

``PROGRESS``
   Operational progress for live observers.

``RESULT``
   Availability or completion of a structured result.

Observers
---------

An observer is any callable satisfying:

.. code-block:: python

   def observer(event: Event) -> None:
       ...

Shared implementations include:

``NullObserver``
   Ignore events. This is the default for silent library use.

``ListObserver``
   Collect operational events for tests or custom processing.

``CallbackObserver``
   Forward events to an application callback.

Frontends may implement additional observers without changing calculators.

Calculator emission
-------------------

:class:`quantas.models.BasicCalculator` provides ``emit`` and ``add_warning``.
A calculator should emit meaningful stage boundaries and structured results,
not every internal scalar operation.

.. code-block:: python

   self.emit("Fitting static equation of state")
   self.emit(
       "Static EOS fit completed",
       level=EventLevel.RESULT,
       data={"model": model_name, "r_squared": r_squared},
   )

Event data should remain lightweight. Large arrays are summarized before
persistence, because the authoritative numerical arrays already belong in the
result payload.

Numerical progress
------------------

Core numerical functions do not emit Quantas events. They accept a callback
such as ``callback(current, total)``. The calculator converts that callback to:

.. code-block:: python

   self.emit(
       "Evaluating states",
       level=EventLevel.PROGRESS,
       progress=current / total,
   )

This keeps core routines independent from observers and lets the same algorithm
run inside tests, notebooks, the CLI, or a GUI.

Persistent and operational events
---------------------------------

``PROGRESS`` events are observer-only and are intentionally excluded from
``ResultData.events`` and native HDF5. Persisting thousands of loop updates
would inflate files without adding scientific provenance.

Meaningful messages, warnings, errors, and result events are converted to
serializable :class:`quantas.core.events.EventRecord` objects. Array payloads
are represented by type/shape/dtype summaries.

Warnings
--------

Use ``add_warning`` when a calculation completes but interpretation is limited,
for example:

* extrapolation;
* a poorly conditioned fit;
* an unresolved polarization branch;
* incomplete optional metadata;
* a fallback explicitly allowed by the method.

The warning is both emitted and stored. Do not emit a warning for a state that
must instead raise an exception.

CLI observer behaviour
----------------------

The CLI translates events into Rich presentation:

* warnings in yellow;
* errors in red;
* transient progress bars by default;
* no live progress with ``--no-progress``;
* ordinary text using the terminal theme.

Rich writes directly to standard output/error. Redirected or unsupported
terminals must receive plain text without ANSI sequences.

Testing events
--------------

Test:

* event ordering for meaningful stages;
* progress normalization and final completion;
* warning persistence;
* exclusion of progress from persistent histories;
* serializable summarization of event data;
* observer independence of numerical results;
* CLI rendering separately from calculator behaviour.

Do not freeze incidental debug wording more tightly than necessary. Freeze
levels, structured data, and scientifically meaningful messages.
