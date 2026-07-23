Common change recipes
=====================

This page is a practical map for changes that cross several Quantas layers.
Use it as a checklist, then consult the focused architecture pages for details.

Add a scientific option
-----------------------

An option that changes a formula, approximation, domain, tolerance, sampling
resolution, or failure policy is a scientific option.

#. Add a typed field to the passive module ``Options`` contract.
#. Validate range, unit, allowed values, and cross-option consistency.
#. Add focused tests for default, valid alternatives, and invalid values.
#. Thread the option into the analysis layer; do not read Click state there.
#. Preserve the resolved option in ``ResultData.options`` and native HDF5.
#. Include it in reports when it affects interpretation.
#. Expose it through the public API contract.
#. Add the Click option and a meaningful help string.
#. Update Implementation and Workflows, CLI reference, API reference, and the
   relevant tutorial when users need guidance.
#. Compare CLI and API results with the option enabled.

Presentation-only controls such as DPI, colormap, and report precision belong
to renderers or CLI rendering options, not the scientific ``Options`` dataclass.

Add a result field
------------------

#. Define the field on the typed ``Result`` contract with shape and unit in its
   docstring.
#. Calculate it in the analysis layer using full precision.
#. Define invalid/unavailable semantics and any aligned mask.
#. Include uncertainty or covariance when scientifically supported.
#. Add the HDF5 dataset with ``unit`` and ``description`` attributes.
#. Read and validate it in the payload reader.
#. Add round-trip and malformed-shape/unit tests.
#. Decide whether it belongs in reports, exports, plots, and interoperability.
#. Document it in the format specification and API reference.
#. Add a regression test for at least one meaningful state.

Do not reuse an existing field name for a quantity with a different physical
normalization.

Add a CLI option
----------------

#. Confirm that a public API option or operation already exists.
#. Add the Click declaration in the relevant adapter.
#. Group it with related options.
#. Let Click perform basic type conversion and choice validation.
#. Construct the public options/request object.
#. Add or update the centralized descriptive help text when used by the command
   reference generator.
#. Test ``--help``, valid use, invalid use, and output-path behaviour.

A CLI-only scientific option is an architectural error because notebooks and a
GUI could not request the same calculation.

Add a report quantity
---------------------

#. Read the value from the typed result; do not recalculate it in the renderer.
#. Keep the cell value numerical.
#. Add unit and numeric-format metadata to the ``ReportTable``.
#. Update plain-text and Rich snapshots only when their semantic content
   changes.
#. Check redirected output for absence of ANSI sequences.

Add a plot
----------

#. Decide the scientific question answered by the figure.
#. Build required arrays and masks from the typed result.
#. Represent them with existing neutral plot contracts.
#. Add the module plot builder and plot-selection API.
#. Add renderer support only if a genuinely new neutral primitive is required.
#. Add CLI selection and renderer options.
#. Test the plot specification separately from Matplotlib output.
#. Add one renderer smoke or image-structure test.
#. Document scientific interpretation and convergence controls.

Do not place Matplotlib calls in a module plot package.

Add an HDF5 field or change a schema
------------------------------------

A backward-compatible optional dataset may not require a shared envelope
version change, but it still needs a module payload version or explicit reader
policy when interpretation changes.

#. Define whether the change is optional, required, or incompatible.
#. Update the writer to emit the current layout only.
#. Update the reader to validate and migrate supported older layouts.
#. Never infer a missing unit or module solely from the filename.
#. Preserve ``None`` and unavailable states explicitly.
#. Add old-file and new-file tests.
#. Update the format documentation.
#. Update registry-based opening if metadata or dispatch changes.
#. Decide whether a schema version increment is required.

See :doc:`persistence`.

Add a public API symbol
-----------------------

#. Add a stable wrapper or alias under ``quantas.api``.
#. Add it to the namespace ``__all__``.
#. Add a complete public docstring.
#. Add it once to the curated API reference.
#. Add surface and usage tests.
#. Add a registry capability only when frontend discovery needs it.
#. Ensure the CLI uses the public operation when applicable.

Do not promote a calculator or reader simply to avoid writing a wrapper.

Add an external-code parser
---------------------------

#. Place code-specific syntax under ``quantas.interfaces.<code>``.
#. Define file-recognition and completion checks.
#. Parse into explicit physical fields with units.
#. Preserve source provenance and relevant code settings.
#. Convert to a normalized module input through a separate adapter.
#. Test complete, incomplete, malformed, and version-variant fixtures.
#. Keep parser errors contextual and actionable.

See :doc:`interfaces`.

Add a canonical citation
------------------------

#. Add one immutable ``Citation`` to ``quantas.references.registry``.
#. Use a stable lowercase key; store the DOI without a URL prefix.
#. Add the record to ``CITATIONS``.
#. Associate the key with module/method sets when appropriate.
#. Cite the key in the scientific page at the relevant statement.
#. Regenerate bibliography fragments.
#. Update report citation sets if the implemented method changed.
#. Run citation and documentation tests.

See :doc:`citation_registry`.

Change a numerical default
--------------------------

A default change can alter scientific results for every user. Treat it as a
scientific change:

#. characterize the old default;
#. explain the failure or limitation motivating the new value;
#. run convergence or sensitivity studies;
#. update Options, CLI, workflow documentation, tutorial commands, reports,
   and frozen examples;
#. preserve the old value as an explicit selectable option when scientifically
   meaningful;
#. record the change in release notes.

Do not change a default merely because a larger grid or higher polynomial order
sounds more accurate.
