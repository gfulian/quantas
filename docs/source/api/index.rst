API overview
============

The supported Python contract begins at :mod:`quantas.api`.  Applications,
notebooks, services, the command-line interface, Quantas GUI, and other
graphical frontends should import from these namespaces rather than from implementation packages.  Those packages may change as the
code evolves and are not application contracts.

The API is organized by scientific domain rather than as one flat collection
of functions:

.. code-block:: python

   from quantas.api import (
       common,
       elasticity,
       eos,
       ha,
       interop,
       profiles,
       plotting,
       qha,
       registry,
       rendering,
       seismic,
       thermoelasticity,
   )

Only names exported through the ``__all__`` attribute of these modules are
part of the supported public surface.  A class may be implemented internally
in another package, but users should import and reference it through its
``quantas.api`` alias.

The reference pages use explicit autodoc directives for those aliases rather
than expanding every imported module member.  Automated tests require every
exported name to be documented exactly once and reject undocumented or
accidentally promoted symbols.

API lifecycle patterns
----------------------

Most scientific modules follow a single-shot lifecycle:

.. code-block:: text

   input or path
   -> normalize_input
   -> run
   -> ResultData
   -> typed module payload
   -> report / plots / export / HDF5

HA, QHA, Elasticity, SEISMIC, and Thermoelasticity expose passive ``Input``,
``Options``, and ``Result`` contracts around this pattern.  EOS is deliberately
different: one dataset can produce many immutable fit records, so its public
surface is organized around fit requests, batch plans, archives, sessions,
diagnostics, and post-fit calculations.

Result envelopes and typed payloads
-----------------------------------

Single-shot workflows return :class:`quantas.api.common.ResultData`.  The
envelope contains common metadata, normalized input, resolved options,
warnings, persistent events, and one module-specific payload.  The
module-specific ``get_result`` function validates the envelope before
returning the typed payload:

.. code-block:: python

   from quantas.api import qha

   result_data = qha.run("phonons.yaml", options=qha.Options())
   result = qha.get_result(result_data)
   print(result.equilibrium_volume.shape)

Do not extract a payload by assuming a string key when a typed accessor is
available.

Frontend-neutral reporting and plotting
---------------------------------------

Calculations build :class:`quantas.api.common.ReportTable` and
:class:`quantas.api.plotting.PlotCollection` objects.  Concrete plot types and
primitives are public through :mod:`quantas.api.plotting`; plain-text and
Matplotlib output is produced separately by :mod:`quantas.api.rendering`.
This keeps numerical workflows independent from terminals, notebooks, web
applications, and graphical interfaces.

Observers and progress
----------------------

Long-running public operations accept an :class:`quantas.api.common.Observer`.
Events contain messages, levels, structured data, and optional normalized
progress.  The calculator does not create a Rich progress bar or GUI widget;
the frontend decides how events are presented.

Runtime dependencies and optional features
------------------------------------------

Importing :mod:`quantas` remains lightweight.  The base scientific runtime
includes ``odrpack`` for orthogonal-distance EOS regression and ``spglib`` for
crystallographic symmetry and structure-normalization operations.  These
packages are required dependencies even though their functionality is imported
only when the corresponding workflow is used.

Static plotting remains optional and requires the ``plot`` extra.  A missing
required runtime backend indicates an incomplete or broken installation rather
than an intentionally unavailable Quantas capability.

Public namespaces
-----------------

.. list-table::
   :header-rows: 1
   :widths: 22 48 30

   * - Namespace
     - Purpose
     - Continue with
   * - :mod:`quantas.api.common`
     - Shared envelopes, neutral report/plot contracts, and events.
     - :doc:`common`
   * - :mod:`quantas.api.ha`
     - Harmonic thermodynamic calculations.
     - :doc:`ha`
   * - :mod:`quantas.api.qha`
     - Quasi-harmonic inspection, calculation, validation, and comparison.
     - :doc:`qha`
   * - :mod:`quantas.api.elasticity`
     - Second-order elastic analysis and directional surfaces.
     - :doc:`elasticity`
   * - :mod:`quantas.api.seismic`
     - Christoffel phase/group analysis, polarization, and enhancement.
     - :doc:`seismic`
   * - :mod:`quantas.api.eos`
     - EOS datasets, fitting, archives, diagnostics, and calculations.
     - :doc:`eos`
   * - :mod:`quantas.api.thermoelasticity`
     - QSA calibration, P--T reconstruction, and geological profiles.
     - :doc:`thermoelasticity`
   * - :mod:`quantas.api.profiles`
     - Frontend-neutral terrestrial profile specifications.
     - :doc:`profiles`
   * - :mod:`quantas.api.plotting`
     - Frontend-neutral plot specifications and reusable plotting primitives.
     - :doc:`plotting`
   * - :mod:`quantas.api.rendering`
     - Supported rendering bridge for neutral tables and plots.
     - :doc:`rendering`
   * - :mod:`quantas.api.registry`
     - Capability-based discovery for frontends and services.
     - :doc:`registry`
   * - :mod:`quantas.api.interop`
     - Explicit transformations between scientific workflows.
     - :doc:`interoperability`

Related documentation
---------------------

- :doc:`../getting_started/python_api` introduces the common programming
  patterns.
- :doc:`../workflows/index` explains scientific implementation choices and
  numerical limitations.
- :doc:`../tutorials/index` provides complete CLI/API examples.
- :doc:`../formats/index` specifies persisted inputs and results.
