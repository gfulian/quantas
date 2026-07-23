Packaging and platform validation
=================================

Quantas is distributed as a pure-Python wheel and source distribution. The
runtime supports Python 3.10 and newer on Linux, Windows, and macOS. Quantas
does not ship its own compiled extension, but NumPy, SciPy, h5py, spglib, and
ODRPACK may install platform-specific wheels.

Build artifacts
---------------

From a clean checkout:

.. code-block:: console

   rm -rf build dist *.egg-info
   python -m build
   python -m twine check dist/*

On Windows, remove the same directories with the native shell or Python rather
than assuming Unix commands.

The source distribution should contain:

* license and README;
* source package;
* tests and developer tools required by project policy;
* examples and provenance;
* documentation sources and selected generated assets.

The wheel contains runtime code, package metadata, and the ``py.typed`` marker.

Clean-install validation
------------------------

Editable installs can hide packaging mistakes. Validate built artifacts in
fresh virtual environments:

.. code-block:: console

   python tools/check_distribution.py dist

The distribution checker verifies package metadata, command entry point,
organized public API namespaces, capability registry, ``py.typed``, and absence
of unintended GUI dependencies.

Manual smoke tests should also cover:

.. code-block:: console

   quantas --version
   quantas --help
   python -c "import quantas; from quantas.api import registry; print(registry.list_modules())"

Optional dependencies
---------------------

``plot``
   Matplotlib static rendering.

``docs``
   Sphinx documentation toolchain.

``test``
   Pytest and test-only dependencies.

``typecheck``
   Mypy environment.

``dev``
   Complete contributor environment.

ODRPACK is a base scientific dependency because ODR is part of the supported
EOS fitting surface. spglib is also a base scientific dependency because
crystallographic symmetry analysis and structure normalization are supported
runtime capabilities. GUI frameworks remain outside Quantas and belong to a
separate application distribution.

Package data and examples
-------------------------

When adding a non-Python file, decide deliberately whether it belongs in:

* the runtime wheel;
* source distribution only;
* documentation build only;
* examples only;
* tests only.

Update ``MANIFEST.in``, setuptools package data, example manifests, and
distribution tests as required. Do not assume a file present in the checkout is
present after installation.

Versioning
----------

The package version is defined once and exposed through package metadata. Text
and HDF5 schemas use independent versions. Do not increment a persisted schema
merely because the package version changed.

Pre-release corrections may remove unsupported historical aliases according to
project policy. Stable releases require explicit deprecation and migration for
public API changes.

Cross-platform concerns
-----------------------

Test:

* :class:`pathlib.Path` handling;
* UTF-8 text and newline differences;
* shell scripts and Windows ``make.bat``;
* subprocess invocation without shell-specific quoting;
* temporary file closing before reopen/delete;
* Matplotlib backend isolation;
* console entry-point installation;
* compiled dependency availability on supported Python versions.

Continuous integration
----------------------

The intended matrix covers:

* lint and complete type checking;
* supported Python versions;
* staged scientific and frontend tests;
* documentation build with warnings as errors;
* wheel and source distribution checks;
* clean-install tests;
* Windows and macOS platform jobs.

Release checkpoint
------------------

Before a release candidate:

#. start from a clean checkout;
#. run the complete staged test runner;
#. run Ruff and Mypy;
#. build documentation strictly;
#. regenerate and verify examples/assets;
#. build wheel and sdist;
#. validate both distributions in clean environments;
#. inspect package contents;
#. verify version and changelog;
#. record checksums and the scientific validation matrix.
