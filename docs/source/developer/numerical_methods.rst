Changing numerical and scientific code
======================================

Scientific reliability has priority over cleanup, abstraction, and speed. A
numerical refactor is successful only when it preserves the intended physical
quantity, units, shapes, conventions, tolerances, and supported domain.

Start from current behaviour
----------------------------

Before changing a formula or algorithm:

#. identify the current implementation and every caller;
#. capture current behaviour with focused characterization tests;
#. record units, shapes, axis order, dtype, and edge-case behaviour;
#. identify an analytical, legacy, published, experimental, or independent
   computational reference;
#. make the smallest possible change;
#. compare old and new results over ordinary and difficult cases;
#. optimize only after the validated implementation is stable.

A characterization test is not automatically a scientific validation. It
protects against accidental behavioural drift while the correctness of that
behaviour is assessed separately.

Units, conventions, and shapes
------------------------------

Every public numerical function should make the following unambiguous:

* units of every input and output;
* whether values are per cell, formula unit, atom, mole, or unit volume;
* tensor convention and coordinate frame;
* axis order and array shape;
* whether masks and uncertainties align with the leading dimensions;
* expected dtype.

For example, ``(nT, nP, 6, 6)`` and ``(nP, nT, 6, 6)`` are not interchangeable
merely because they contain the same number of values. Axis order is part of
the scientific contract.

Precision policy
----------------

Use the central precision helpers in ``quantas.core.numerics`` when generic
serialization or normalization needs to preserve the Quantas policy.
Scientific arrays use ``float64`` and complex arrays use ``complex128``.

Do not:

* introduce a user-selectable single-precision mode;
* downcast HDF5 values to save space;
* round arrays before storing them;
* use report precision as numerical precision;
* mix dtype changes with an unrelated refactor.

Tolerances are scientific decisions
-----------------------------------

A tolerance should have:

* a name;
* a unit or a clear dimensionless meaning;
* a documented default;
* a test at and around the threshold;
* a rationale connected to numerical noise or physical resolution;
* an explicit consequence when the threshold is exceeded.

Do not use a larger tolerance to hide an unstable tensor, an ill-conditioned
fit, a discontinuous mode assignment, or an invalid external-code output.

Validation, extrapolation, and failure
--------------------------------------

Distinguish these states:

``valid``
   The requested quantity is defined and inside the supported numerical domain.

``extrapolated``
   A value was calculated outside the sampled support according to a documented
   policy. Preserve a mask and warning.

``not applicable``
   The quantity has no physical definition for the selected method or state.

``unresolved``
   The base quantity may be valid, but a branch, derivative, or decomposition is
   ambiguous, for example a shear polarization at a degeneracy.

``invalid``
   The requested result cannot be interpreted safely.

A finite number is not necessarily a valid result. Do not replace invalid
states with zero, clip unstable eigenvalues, or silently copy a fallback value
unless the method explicitly defines and records that fallback.

Progress callbacks
------------------

Core and numerical objects do not emit Quantas events. A long loop may accept a
simple callback:

.. code-block:: python

   def solve_states(
       states: list[object],
       progress_callback: Callable[[int, int], None] | None = None,
   ) -> list[object]:
       results = []
       for current, state in enumerate(states, start=1):
           results.append(solve_one(state))
           if progress_callback is not None:
               progress_callback(current, len(states))
       return results

The calculator translates this into ``EventLevel.PROGRESS`` with a normalized
fraction. This preserves reuse in tests, notebooks, and future frontends.

Uncertainties and covariance
----------------------------

When a result includes uncertainty:

* state whether it is observational, fitted, propagated, or heuristic;
* preserve covariance when derived quantities share parameters;
* use Jacobians or resampling consistently with the documented method;
* distinguish statistical uncertainty from model discrepancy;
* preserve masks when an uncertainty is unavailable;
* test deterministic seeds for stochastic propagation.

Do not create a zero uncertainty for a missing covariance. ``Zero`` means a
known exact value; ``unavailable`` means something different.

Optimization and performance
----------------------------

Use NumPy and SciPy first. Before adding an accelerator:

#. profile a realistic workflow;
#. identify the dominant operation;
#. establish a numerical baseline;
#. improve array layout, repeated allocation, or algorithmic complexity;
#. measure wall time and memory again;
#. verify identical scientific results.

Numba or another accelerator must not be introduced without measured benefit
and a maintenance plan. A faster algorithm that changes branch selection,
tolerances, or reproducibility is a scientific change, not a transparent
optimization.

Batch sizes and sampling densities are different
------------------------------------------------

A batch size generally controls memory and throughput. A sampling density
controls scientific resolution. Document them separately.

Changing a batch size should leave results unchanged. Changing an angular,
pressure, temperature, or interpolation grid may alter extrema or derivatives
and requires convergence testing.

Recommended scientific test ladder
----------------------------------

For a new numerical method, build evidence in this order:

#. exact identities and input validation;
#. analytical or manufactured cases;
#. invariance tests, such as rotation or unit conversion;
#. comparison with the legacy implementation;
#. comparison with an independent package or published reference;
#. real-data regression with frozen provenance;
#. CLI/API and HDF5 round-trip equivalence;
#. performance measurement after correctness is established.

Document tolerances separately for numerical regression and scientific
acceptance. A test tolerance chosen to accommodate platform floating-point
noise is not automatically an acceptable disagreement with experiment.
