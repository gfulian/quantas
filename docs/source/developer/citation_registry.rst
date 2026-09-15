Citation registry
=================

Scientific references are centralized under :mod:`quantas.references`. This
prevents inconsistent author lists, titles, journal details, and DOI strings
across reports and documentation.

Data model
----------

A canonical :class:`quantas.references.Citation` stores:

* stable key;
* ordered authors;
* title and year;
* article, book, chapter, preprint, report, or software kind;
* journal, volume, pages, publisher, or parent-book metadata as appropriate;
* a report number for technical reports when applicable;
* DOI without a URL prefix;
* external URL only when a DOI is unavailable.

Citation objects are immutable dataclasses.

Canonical bibliographic style
-----------------------------

The registry stores scientific metadata independently from the surface used to
render it. Canonical records follow these rules:

* authors and editors use ``initials + surname`` in publication order, for
  example ``R. J. Angel`` or ``F.-X. Coudert``;
* scientific spelling and Unicode typography are retained in canonical data,
  for example ``G. Valdrè`` and ``Zeitschrift für Kristallographie``;
* DOI values are stored without ``doi:`` or URL prefixes;
* page ranges use a simple ASCII hyphen in machine data;
* book chapters and technical reports use their own record kinds rather than
  being represented as journal articles;
* a URL is stored only when no DOI is available.

The plain-text renderer transliterates canonical Unicode metadata to portable
ASCII for terminal reports and HDF5-embedded text. The RST renderer preserves
the canonical scientific typography. Do not degrade the registry itself merely
to satisfy a terminal encoding constraint.

Key conventions
---------------

Use a lowercase stable identifier such as:

.. code-block:: text

   stixrude_lithgow_bertelloni_2005
   eosfit7_angel_gonzalez_platas_alvaro_2014

The key is an internal scientific identifier used by tests, report sets, and
RST footnotes. Do not change it merely to adjust typography.

Register a reference
--------------------

#. Add one ``Citation`` constant in ``quantas.references.registry``.
#. Add it to the canonical ``CITATIONS`` mapping.
#. Store the DOI as, for example, ``10.1111/...``, not
   ``https://doi.org/10.1111/...``.
#. Add or update registry tests.
#. Associate it with module or method citation sets when appropriate.

Module and method sets
----------------------

``MODULE_CITATION_KEYS`` identifies the standard references associated with a
workflow. ``METHOD_CITATION_KEYS`` identifies references for a particular
scientific method such as quasi-harmonic thermodynamics, Christoffel acoustics,
or cold finite strain.

A module report may combine its general citation set with method-specific
references selected by the actual options. Do not cite a method that was not
used merely because the module supports it.

Report rendering
----------------

The citation renderer produces deterministic plain text and DOI URLs. This
supports report footers and HDF5-embedded report text without duplicating
bibliographic formatting logic.

Scientific-background pages
---------------------------

Theory pages use labelled auto-numbered footnotes:

.. code-block:: rst

   ... Eulerian finite strain [#stixrude_lithgow_bertelloni_2005]_.

The displayed number is local to the page and follows first appearance. The
stable label is the canonical registry key.

Bibliography fragments are generated with:

.. code-block:: console

   python docs/tools/generate_theory_bibliographies.py

The generated entry links the DOI to ``https://doi.org/<DOI>``. Do not write a
second free-form copy of the same reference in the page source.

Adding a citation to documentation
----------------------------------

#. Register the citation.
#. Add its key to the page list used by the bibliography generator.
#. Cite the key at the statement it supports.
#. Regenerate the fragment.
#. Check first-appearance ordering.
#. Run citation and documentation-tree tests.

Changing a reference
--------------------

Correct a canonical record only after checking the publication. A correction
changes every generated occurrence. Keep the key stable unless the original
record represented the wrong publication.

Validation
----------

Tests should verify:

* unique keys and DOI values;
* canonical ``initials + surname`` author/editor formatting;
* valid year and record kind;
* record-specific metadata for chapters and reports;
* DOI stored without URL prefix and no redundant URL when a DOI exists;
* referenced keys exist;
* module/method sets preserve intended order;
* documentation footnotes match the registry;
* generated fragments are current;
* DOI links use the canonical resolver.
