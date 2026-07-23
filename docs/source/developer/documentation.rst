Documentation workflow
======================

Documentation is part of the tested package. A change that modifies public
options, result fields, units, schemas, scientific defaults, or module
contracts must update the relevant documentation in the same change.

Build locally
-------------

Install the documentation environment:

.. code-block:: console

   python -m pip install -e ".[docs]"

Build HTML on Linux or macOS:

.. code-block:: console

   make -C docs clean html SPHINXOPTS="-W --keep-going"

On Windows:

.. code-block:: doscon

   .\docs\make.bat clean
   .\docs\make.bat html

For release validation, warnings should fail the build. ``--keep-going`` is
useful because it reports multiple independent problems in one pass.

Manual structure
----------------

The public documentation separates:

``Scientific Background``
   Theory, equations, assumptions, and references.

``Implementation and Workflows``
   Quantas-specific numerical choices, defaults, limitations, diagnostics, and
   performance controls.

``Tutorials``
   Reproducible CLI and API analyses using curated data.

``Input and Output Formats``
   Data contracts, shapes, units, schemas, and independent HDF5 inspection.

``Command Reference``
   Syntax generated from Click plus contextual guidance.

``API Reference``
   Curated public ``quantas.api`` symbols.

``Scientific Validation``
   External and regression evidence.

``Development Guide``
   Architecture and contributor procedures.

Do not copy the same long explanation into every section. Link to the
authoritative discussion and repeat only the context required for readability.

Docstrings
----------

Every public module, class, dataclass, function, and method needs a technical
NumPy/Sphinx-compatible docstring. Include:

* concise description;
* parameters and units;
* returns and shapes;
* raises;
* notes on conventions or normalization;
* examples only when they add real value.

A public function that returns ``ndarray`` should state the array shape and
axis order. A tolerance should state its unit and effect.

CLI reference
-------------

Command syntax is generated with ``sphinx-click`` from the actual Click tree.
Option help is therefore part of both terminal UX and the manual.

Describe what an option controls and when it should change. Do not write help
such as “Set the tolerance” when the user needs to know which comparison the
tolerance affects.

API reference
-------------

API pages are curated from explicit ``__all__`` inventories. Do not use
``:imported-members:`` or broad ``automodule`` blocks that expose accidental
implementation names.

When a public symbol is added, document it exactly once and update the coverage
tests.

Scientific citations
--------------------

Theory pages use labelled auto-numbered footnotes backed by the canonical
registry. Regenerate fragments with:

.. code-block:: console

   python docs/tools/generate_theory_bibliographies.py

See :doc:`citation_registry`.

Tutorial assets
---------------

Checked-in tutorial figures and downloadable tables are reproducible from
scripts under ``docs/tools``. Generated HDF5 archives should normally remain
temporary; commit selected figures, reports, CSV files, specifications, and
public-API scripts required by the tutorial.

Asset generators must:

* use public APIs for user-facing workflows;
* be deterministic;
* clean temporary files;
* create parent directories;
* preserve source provenance;
* be covered by manifest or existence tests.

The EOS asset generator is intentionally explicit because full MGD fitting is
not appropriate for every Sphinx startup.

RST quality
-----------

Check:

* title underlines match or exceed title length;
* all ``:doc:`` targets exist;
* every RST page is in one toctree or marked ``:orphan:``;
* downloads, literal includes, and images exist;
* no duplicate toctree targets;
* no unresolved references;
* no trailing whitespace;
* code blocks use the correct language.

Infrastructure tests cover many of these conditions, but a strict Sphinx build
is still required.

Writing examples
----------------

Examples should:

* import from :mod:`quantas.api` only;
* use paths relative to a documented working directory;
* show output creation explicitly;
* avoid hidden state;
* use real or scientifically meaningful data;
* include reproducibility checkpoints;
* distinguish expected values from illustrative placeholders.

Do not publish an example that has not been executed against the snapshot.

Updating documentation for a code change
----------------------------------------

Use this map:

.. list-table::
   :header-rows: 1

   * - Change
     - Documentation to inspect
   * - Formula or approximation
     - Scientific Background, workflow, validation, citation
   * - Scientific option/default
     - Workflow, CLI, API, tutorial
   * - Input field
     - Format, API, tutorial/input generator
   * - Result/HDF5 field
     - Format, API, report/export/plot, validation
   * - CLI option
     - Click help, command context, workflow
   * - Public symbol
     - API reference, examples, stability notes
   * - Plot type
     - Workflow, tutorial, renderer/API reference
   * - Parser/interface
     - Format, interface developer page, tutorial provenance

Documentation review
--------------------

Before merging:

#. run documentation infrastructure tests;
#. regenerate affected assets and bibliographies;
#. build Sphinx with warnings as errors;
#. inspect the changed pages in a browser at desktop and narrow widths;
#. test all copied commands in a clean working directory;
#. verify downloads from the built HTML;
#. check external DOI links during release preparation.
