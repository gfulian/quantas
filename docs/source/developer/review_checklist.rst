Change review checklist
=======================

Use the relevant parts of this checklist before creating a checkpoint, opening
a pull request, or handing work to another contributor.

Scientific contract
-------------------

* [ ] The physical quantity and approximation are stated.
* [ ] Units, normalization, tensor convention, axis order, and dtype are explicit.
* [ ] Current behaviour was characterized before a numerical refactor.
* [ ] Analytical, legacy, published, experimental, or independent references
      were identified where required.
* [ ] Tolerances have rationale and boundary tests.
* [ ] Invalid, unresolved, not-applicable, and extrapolated states are distinct.
* [ ] Uncertainty and covariance semantics are documented.
* [ ] No scientific data are rounded for storage or internal use.

Architecture
------------

* [ ] Core code has no frontend or workflow imports.
* [ ] Scientific modules do not import one another.
* [ ] Active objects are classes; passive contracts are dataclasses.
* [ ] Long numerical loops use simple callbacks rather than Quantas events.
* [ ] Module plot builders do not import Matplotlib.
* [ ] The CLI calls :mod:`quantas.api`.
* [ ] Optional plotting dependencies do not leak into base calculations.

Inputs and interfaces
---------------------

* [ ] File recognition and completion checks are separate.
* [ ] External units and conventions are converted explicitly.
* [ ] Source provenance and relevant external-code settings are preserved.
* [ ] Complete, incomplete, malformed, and variant fixtures are tested.
* [ ] Normalized inputs have documented shapes and units.

Results and persistence
-----------------------

* [ ] The typed result contains full-precision values and aligned masks.
* [ ] New HDF5 fields have unit and description metadata.
* [ ] Writer and reader were updated together.
* [ ] Current and supported historical schemas are explicit.
* [ ] Write/read round-trip tests pass.
* [ ] Missing values are not confused with zero.
* [ ] Metadata-based dispatch still works.

Events and errors
-----------------

* [ ] Meaningful stages, warnings, and failures emit appropriate levels.
* [ ] Progress events are not persisted.
* [ ] Event data are lightweight and serializable.
* [ ] Errors are actionable and not hidden by broad exception handling.
* [ ] A warning is not used where the result should be rejected.

Reports, plots, and exports
---------------------------

* [ ] Report tables contain raw values, units, and deterministic ordering.
* [ ] Plain reports contain no ANSI or live-progress artifacts.
* [ ] Plot specifications are frontend-neutral.
* [ ] Plot resolution and renderer quality controls are distinguished from
      scientific convergence controls.
* [ ] Exports preserve documented precision, units, and provenance.

Public API and CLI
------------------

* [ ] Public names are reviewed and added explicitly to ``__all__``.
* [ ] Every public symbol has a technical docstring and one API-reference entry.
* [ ] Registry capabilities and types are current.
* [ ] CLI help explains each new option.
* [ ] API and CLI execute the same numerical path.
* [ ] Output paths, overwrite behaviour, and exit status are tested.

Tests
-----

* [ ] Focused unit/module tests pass.
* [ ] Architecture and public-surface tests pass.
* [ ] HDF5 and export tests pass.
* [ ] CLI/API equivalence is covered.
* [ ] Plot-specification and renderer tests are separated.
* [ ] A real-data or scientific regression exists when appropriate.
* [ ] Ruff, Mypy, and ``compileall`` pass.
* [ ] The complete staged test runner passes.

Documentation and citations
---------------------------

* [ ] Theory, workflow, tutorial, format, CLI, and API pages were inspected.
* [ ] Scientific references are in the canonical registry.
* [ ] Bibliography and tutorial assets were regenerated when affected.
* [ ] RST links, downloads, images, and title underlines are valid.
* [ ] Sphinx builds with warnings as errors.
* [ ] Commands and Python examples were executed against the snapshot.

Packaging
---------

* [ ] New files are included in the intended distribution artifact.
* [ ] Wheel and source distribution build successfully.
* [ ] Twine validation passes.
* [ ] Clean-install distribution checks pass.
* [ ] ``py.typed`` and the console entry point remain available.
* [ ] Platform-specific paths and subprocesses are portable.

A checkbox is not evidence by itself. Record the command, test, dataset, or
reference that demonstrates each important scientific and release claim.
