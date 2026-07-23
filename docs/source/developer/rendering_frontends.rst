Rendering and frontend integration
==================================

Quantas separates scientific preparation from concrete presentation. Reports
and plots are first represented by passive frontend-neutral contracts, then
rendered by Rich, plain text, Matplotlib, or a future GUI.

Report pipeline
---------------

.. code-block:: text

   typed Result
       |
       v
   module report builder
       |
       v
   list[ReportTable]
       |
       +--> plain-text renderer
       +--> Rich renderer
       +--> CSV/custom frontend

A :class:`quantas.models.ReportTable` stores:

* title;
* ordered columns;
* rows containing raw values;
* metadata for units, alignments, formats, notes, and provenance.

The module decides which quantities are scientifically relevant. The renderer
decides spacing, terminal borders, colors, and textual formatting.

Numeric formatting
------------------

Formatting profiles live in the table renderer. Stored and in-memory values
remain full precision. A report profile may display elastic constants with a
fixed number of decimals or energies in scientific notation without changing
HDF5 values.

Do not convert all numerical cells to strings inside the report builder.

Plain-text reports
------------------

Persistent reports must be deterministic and contain:

* no ANSI escapes;
* no box-drawing dependence;
* no transient progress artifacts;
* stable ordering;
* documented units;
* values formatted for reading, not for re-import as the authoritative data.

The native HDF5 result remains the source for subsequent numerical analysis.

Rich terminal rendering
-----------------------

Rich owns interactive terminal presentation, including:

* tables;
* warning/error colors;
* transient progress;
* terminal-aware wrapping.

Rich writes directly to ``sys.stdout`` or ``sys.stderr`` rather than Click
stream wrappers. Redirected output must fall back to plain text without ANSI
sequences.

Plot pipeline
-------------

.. code-block:: text

   typed Result
       |
       v
   module plot builder
       |
       v
   PlotCollection
       |
       +--> Matplotlib renderer
       +--> future GUI renderer

Neutral plot contracts describe axes, series, bands, contours, spherical maps,
3D surfaces, vector overlays, scalar backgrounds, and panel layouts. They
contain prepared scientific data and portable style hints, not concrete artist
objects.

Module plot builders
--------------------

A module builder owns:

* result-field selection;
* scientific unit conversion for the displayed quantity;
* masks and invalid-state handling;
* branch selection;
* physical versus normalized geometry;
* labels and legends;
* uncertainty bands;
* plot metadata required for interpretation.

It must not import Matplotlib.

Matplotlib backend
------------------

The backend owns:

* figure and axes creation;
* artist construction;
* colormap resolution;
* file format and DPI;
* final layout;
* renderer-specific fallbacks.

Matplotlib is an optional dependency and is confined to
``quantas.renderers.plots.matplotlib``.

Adding a plot primitive
-----------------------

Before adding a new neutral class, check whether existing series, bands,
contours, colored paths, surfaces, or vector layers can represent the data.

When a new primitive is necessary:

#. define a passive dataclass in ``quantas.models.plot``;
#. validate aligned shapes and values in ``__post_init__``;
#. add it to the appropriate union/collection contract;
#. add backend dispatch and rendering;
#. test the model independently from Matplotlib;
#. test one rendered figure;
#. update public API/rendering documentation when exposed.

Future GUI requirements
-----------------------

A GUI should:

* call :mod:`quantas.api`;
* use public observer callbacks for progress;
* consume typed results, ``ReportTable``, and ``PlotCollection``;
* inspect module capabilities through the registry;
* keep GUI state and widgets outside Quantas core and modules;
* never reproduce a scientific formula merely to update a display.

The GUI may render a neutral plot with Plotly or another backend, but it must
respect the same masks, units, and branch semantics prepared by the module.

Testing rendering boundaries
----------------------------

Architecture tests verify that:

* Matplotlib imports are confined to the backend;
* module plot builders do not import concrete renderers;
* plot specification classes have one authoritative definition;
* ``ReportTable`` has one authoritative definition;
* public API rendering exposes bridges, not concrete renderer internals.

Renderer tests should not be used as substitutes for scientific result tests.
Test the arrays and plot specification before testing pixels or artists.
