Developer setup
===============

Quantas is a scientific library first. Numerical behaviour, units, shapes, and
native HDF5 contracts must remain independent from the command line, Quantas GUI,
and other graphical frontends.

Development environment
-----------------------

Use Python 3.10 through 3.13 in a fresh virtual environment. Python 3.14
is not yet supported because the complete scientific dependency stack and
validation suite have not been certified on that runtime. From the
repository root:

.. code-block:: console

   python -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install -e ".[dev]"

On Windows PowerShell, activate with:

.. code-block:: powershell

   .venv\Scripts\Activate.ps1

The ``dev`` extra installs the scientific runtime, tests, static-analysis tools,
packaging tools, Sphinx, and the plotting backend. The smaller extras are useful
when reproducing an end-user installation:

``plot``
   Static Matplotlib rendering.

``test``
   Pytest and test-only dependencies.

``docs``
   Sphinx, numpydoc, ``sphinx-click``, and the Read the Docs theme.

``typecheck``
   Mypy and the supported NumPy range used for type checking.

GUI frameworks are intentionally not Quantas dependencies. A separate GUI
application should install Quantas and consume :mod:`quantas.api`.

Verify the checkout
-------------------

Before editing code, establish a clean baseline:

.. code-block:: console

   ruff check src tests tools
   mypy
   python -m compileall -q src tests tools
   python tools/run_tests.py all -- -q

The complete test runner uses isolated subprocesses, disables unrelated pytest
plugin autoloading, gives each stage its own Matplotlib configuration directory,
and limits numerical library thread counts. A failure in a clean baseline must
be understood before it is mixed with a new change.

Use focused loops while developing
----------------------------------

Run the smallest relevant suite during implementation:

.. code-block:: console

   pytest --suite ha -q
   pytest --suite elasticity -q
   pytest --suite cli -q
   pytest tests/infrastructure/test_architecture.py -q

Then run the complete staged matrix before producing a checkpoint:

.. code-block:: console

   python tools/run_tests.py all -- -q

Do not repeatedly run every plotting suite while changing a scalar numerical
function. Conversely, do not consider a numerical change complete merely
because one unit test passes.

Repository orientation
----------------------

``src/quantas/core``
   Reusable physics, mathematics, chemistry, geometry, numerical policies, and
   frontend-neutral events.

``src/quantas/models``
   Passive contracts and shared active base classes.

``src/quantas/interfaces``
   Parsers and adapters for external scientific codes.

``src/quantas/modules``
   Complete scientific workflows, reports, plot builders, persistence adapters,
   and exports.

``src/quantas/api``
   Supported application-facing namespaces.

``src/quantas/renderers``
   Concrete table and plot presentation.

``src/quantas/cli``
   Click option parsing, validation, dispatch, Rich presentation, and output
   naming.

``tests``
   Scientific, architectural, frontend, persistence, example, and regression
   tests.

``docs``
   Sphinx sources and reproducible documentation asset generators.

``examples``
   Curated scientific inputs used by tutorials and smoke tests.

``tools``
   Test orchestration, distribution validation, and repository checks.

A safe first contribution
-------------------------

For a first code contribution, choose a change that does not alter formulas or
persisted schemas. Good examples are:

* improve one public docstring and its API-reference wording;
* add a missing validation test for an existing option;
* add a report row using a value that is already present in the result payload;
* improve a CLI help string without changing the option contract.

A useful workflow is:

#. identify the authoritative contract and its existing tests;
#. write or extend a focused test;
#. make the smallest change;
#. run lint, type checking, and the focused test;
#. inspect the rendered output or report;
#. run the relevant staged suite.

Avoid beginning with a new Click option or a new HDF5 field. Both are
cross-cutting changes and are better attempted after reading
:doc:`change_recipes`.

Working rules
-------------

* Use ``from __future__ import annotations`` in public modules.
* Use modern Python 3.10 type hints and :class:`pathlib.Path`.
* Use dataclasses for passive data, not for active workflow controllers.
* Keep real calculations and native floating HDF5 values in ``float64``;
  complex values use ``complex128``.
* Do not add a runtime precision selector.
* Do not print from scientific code.
* Do not import Click, Rich, Matplotlib, or GUI objects into numerical code.
* Add technical NumPy/Sphinx-compatible docstrings with parameters, returns,
  and raises sections to every public object.
* Update tests and documentation in the same change as the code.
