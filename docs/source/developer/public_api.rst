Maintaining the public API
==========================

The supported application contract begins at :mod:`quantas.api`. Internal
packages may be imported by Quantas itself, but notebooks, the CLI, and future
frontends should not depend on them.

Supported namespaces
--------------------

The public package exposes domain namespaces lazily:

* :mod:`quantas.api.common`;
* :mod:`quantas.api.elasticity`;
* :mod:`quantas.api.seismic`;
* :mod:`quantas.api.ha`;
* :mod:`quantas.api.qha`;
* :mod:`quantas.api.eos`;
* :mod:`quantas.api.thermoelasticity`;
* :mod:`quantas.api.profiles`;
* :mod:`quantas.api.rendering`;
* :mod:`quantas.api.registry`;
* :mod:`quantas.api.interop`.

``import quantas`` remains lightweight. ``quantas.api`` loads one namespace
only when it is first accessed.

What belongs in the public API
------------------------------

Good public objects are stable scientific contracts or operations:

* passive ``Input``, ``Options``, and ``Result`` types;
* enums and immutable requests needed to configure supported operations;
* normalized readers and constructors;
* workflow operations such as ``run``, ``fit``, or ``analyze_grid``;
* result accessors;
* HDF5 read/write operations;
* report, plot, and export builders;
* registry and interoperability operations intended for applications.

The following normally remain internal:

* concrete calculators, readers, exporters, and solvers;
* implementation controllers;
* CLI and renderer classes;
* low-level HDF5 layout helpers;
* convenience functions that expose unstable internal structures;
* broad re-exports from ``quantas.modules``.

Facade pattern
--------------

A public facade imports internal objects under private aliases and exposes a
small wrapper with a technical docstring:

.. code-block:: python

   from quantas.modules.example.api import run_example as _run
   from quantas.modules.example.models import (
       ExampleInput as Input,
       ExampleOptions as Options,
       ExampleResult as Result,
   )

   def run(
       input_data: Input | str | Path,
       options: Options | None = None,
       observer: Observer | None = None,
   ) -> ResultData:
       """Run the supported Example workflow."""
       return _run(input_data, options=options, observer=observer)

   __all__ = ["Input", "Options", "Result", "run"]

The wrapper is not pointless duplication. It fixes names, type hints,
docstrings, and the supported boundary while internal organization remains
free to change.

``__all__`` is the contract
---------------------------

Every public namespace defines an explicit ``__all__``. The curated API
reference is generated from that inventory, and tests require each public name
to be documented exactly once.

Do not rely on ``from module import *`` behaviour or on names that happen to be
present in a module global namespace.

Adding a public symbol
----------------------

Before promotion, answer:

#. Is the name scientifically meaningful and stable?
#. Can it be represented without exposing a concrete calculator or archive
   implementation?
#. Are units, shapes, exceptions, and lifecycle documented?
#. Does it work without Click, Rich, or Matplotlib?
#. Does it require a capability entry in the registry?
#. Is it needed by notebooks, a frontend, or another supported public
   operation?

Then update, in the same change:

#. the public facade and ``__all__``;
#. the API reference page;
#. public-surface tests;
#. representative usage tests;
#. registry descriptors when applications need discovery;
#. the CLI if the operation is command-line accessible;
#. migration or deprecation documentation when replacing an existing name.

Closed public annotations
-------------------------

A public dataclass is not a closed contract when one of its fields requires an
enum, protocol, or passive structure that can only be imported from
``quantas.core``, ``quantas.models``, or ``quantas.modules``.  Re-export the
minimum stable type through the relevant public namespace or
:mod:`quantas.api.common`; do not re-export the calculator or implementation
package that owns it.

Tests resolve type hints for public ``Input``, ``Options``, request, and context
contracts and require every Quantas-defined referenced type to have a public
identity.

Capability registry
-------------------

:mod:`quantas.api.registry` describes supported operations without importing
all modules eagerly. A :class:`quantas.api.registry.ModuleDescriptor` records:

* stable module name and title;
* public namespace import path;
* result payload key, when applicable;
* supported capabilities;
* canonical public operation names;
* named operation descriptors for capabilities with several supported actions;
* public input, options, and result type aliases.

Add a capability only when the corresponding operation is part of the
supported API. Modules are not required to advertise identical lifecycles.
When one capability has several operations, add stable
:class:`quantas.api.registry.OperationDescriptor` entries instead of inventing
one generic callable with an untyped context mapping.

Result accessors
----------------

Single-shot facades expose ``get_result`` to validate the generic envelope and
return the typed payload. This is preferable to asking users to know a private
``results`` dictionary key.

A new payload type or result key must be updated consistently in:

* the module contract;
* the public accessor;
* the capability descriptor;
* HDF5 read/write code;
* registry-based file opening;
* API and round-trip tests.

Versioning policy
-----------------

The Python package version governs the public API. Persisted text and HDF5
formats have independent schema versions because files may outlive the import
surface.

Before the first stable 2.0 release, reviewed breaking corrections may be made
without maintaining deprecated aliases. After a stable public release,
incompatible removal requires a documented migration and deprecation period.

Application and CLI equivalence
-------------------------------

The CLI must call public API operations. It may perform frontend tasks before
and after the call, but it must not use a different calculator or reconstruct a
scientific result independently. Representative tests should compare API and
CLI results for exact or tolerance-controlled equality.

Contract verification
---------------------

The tests freeze:

* public namespace names;
* exact ``__all__`` inventories;
* lazy imports;
* absence of calculators/readers/exporters/renderers from the supported surface;
* registry capability mappings;
* metadata-based HDF5 dispatch;
* API-reference coverage;
* representative workflows through public facades.

Treat every addition to ``__all__`` as a reviewed API decision.
